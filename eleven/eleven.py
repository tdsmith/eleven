#!/usr/bin/env python
# coding=utf-8

from math import isnan
from numpy import mean, std, power, asarray, log
from scipy.stats.mstats import gmean
from warnings import warn
from types import *
from itertools import repeat
import pandas as pd

log2 = lambda x: log(x)/log(2)

def average_cq(seq, efficiency=1.0):
    """Given a set of Cq values, return the Cq value that represents the
    average expression level of the input.

    The intent is to average the expression levels of the samples,
    since the average of Cq values is not biologically meaningful.

    :param iterable seq: A sequence (e.g. list, array, or Series) of Cq values.
    :param float efficiency: The fractional efficiency of the PCR reaction; i.e.
        1.0 is 100% efficiency, producing 2 copies per amplicon per cycle.
    :return: Cq value representing average expression level
    :rtype: float
    """
    denominator = sum( [pow(2.0*efficiency, -Ci) for Ci in seq] )
    return log(len(seq)/denominator)/log(2.0*efficiency)

def validate_sample_frame(sample_frame):
    """Makes sure that `sample_frame` has the columns we expect.

    :param DataFrame sample_frame: A sample data frame.
    :return: True (or raises an exception)
    :raises TypeError: if sample_frame is not a pandas DataFrame
    :raises ValueError: if columns are missing or the wrong type
    """
    if not isinstance(sample_frame, pd.core.frame.DataFrame):
        raise TypeError("Expected a pandas DataFrame, received {}".format(type(sample_frame)))
    for col in ['Sample', 'Target', 'Cq']:
        if col not in sample_frame:
            raise ValueError("Missing column {} in sample frame".format(col))
    if sample_frame['Cq'].dtype.kind != 'f':
        raise ValueError("Expected Cq column to have float type; has type {} instead".format(str(sample_frame['Cq'].dtype)))
    return True

def censor_background(sample_frame, ntc_samples=['NTC'], margin=log2(10)):
    """Selects rows from the sample data frame that fall `margin` or greater
    cycles earlier than the NTC for that target. NTC wells are recognized by
    string matching against the Sample column.

    :param DataFrame sample_frame: A sample data frame.
    :param iterable ntc_samples: A sequence of strings giving the sample names of your NTC wells, i.e. ['NTC']
    :param float margin: The number of cycles earlier than the NTC for a "good" sample, i.e. log2(10)
    :return: a view of the sample data frame containing only non-background rows
    :rtype: DataFrame
    """
    ntcs = sample_frame.loc[ sample_frame['Sample'].apply(lambda x: x in ntc_samples), ]
    if ntcs.empty:
        return sample_frame
    g = ntcs.groupby('Target')
    min_ntcs = g['Cq'].min()
    # if a target has no NTC, min_ntcs.loc[sample] is NaN
    # we should retain all values from targets with no NTC
    # all comparisons with NaN are false
    # so we test for the "wrong" condition and invert the result
    censored = sample_frame.loc[ ~(sample_frame['Cq'] > (min_ntcs.loc[sample_frame['Target']] - margin)) ]
    return censored

def expression_ddcq(sample_frame, ref_target, ref_sample):
    """Calculates expression of samples in a sample data frame relative to a
    single reference gene and reference sample using the ∆∆Cq method.

    For best results, the ref_sample should be defined for all targets and the
    ref_target should be defined for all samples, or else the series you get
    back will have lots of NaNs.
    
    :param DataFrame sample_frame: A sample data frame.
    :param string ref_target: A string matching an entry of the Target column;
        the target to use as the reference target (e.g. 'Gapdh')
    :param string ref_sample: A string matching an entry of the Sample column.
    :return: a Series of expression values for each row of the sample data
        frame.
    :rtype: Series
    """
    # It might be more correct to replace asarray calls (to discard indexes)
    # with proper joins.
  
    ref_target_df = sample_frame.ix[sample_frame['Target'] == ref_target, ['Sample', 'Cq']]
    ref_target_grouped = ref_target_df.groupby('Sample')
    ref_target_mean_by_sample = ref_target_grouped['Cq'].aggregate(average_cq)
    ref_target_mean_list = ref_target_mean_by_sample.ix[sample_frame['Sample']]
    ref_target_delta = asarray(ref_target_mean_list - ref_target_mean_by_sample[ref_sample])

    ref_sample_df = sample_frame.ix[sample_frame['Sample'] == ref_sample, ['Target', 'Cq']]
    ref_sample_grouped = ref_sample_df.groupby('Target')
    ref_sample_mean_by_target = ref_sample_grouped['Cq'].aggregate(average_cq)
    ref_sample_delta = asarray(sample_frame['Cq'] - asarray(ref_sample_mean_by_target.ix[sample_frame['Target']]))

    rel_exp = pd.Series(
            power(2, ref_target_delta - ref_sample_delta),
            index = sample_frame.index)

    return rel_exp
    
def expression_nf(sample_frame, nf_n, ref_sample):
    """Calculates expression of samples in a sample data frame relative to
    pre-computed normalization factors.

    ref_sample should be defined for all targets or the result will contain
    many NaNs.

    :param DataFrame sample_frame: A sample data frame.
    :param Series nf_n: A Series of normalization factors indexed by sample.
        You probably got this from `compute_nf`.
    :param string ref_sample: The name of the sample to normalize against,
        which should match a value in the sample_frame Sample column.
    :return: a Series of expression values for each row in the sample data
        frame.
    :rtype: Series
    """
    ref_sample_df = sample_frame.ix[sample_frame['Sample'] == ref_sample, ['Target', 'Cq']]
    ref_sample_cq = ref_sample_df.groupby('Target')['Cq'].aggregate(average_cq)

    delta = -sample_frame['Cq'] + asarray(ref_sample_cq.ix[sample_frame['Target']])
    rel = power(2, delta) / asarray(nf_n.ix[sample_frame['Sample']])
    return rel

def collect_expression(sample_frame, ref_targets, ref_sample):
    """Calculates the expression of all rows in the sample_frame relative to
    each of the ref_targets. Used in rank_targets.

    :param DataFrame sample_frame: A sample data frame.
    :param iterable ref_targets: A sequence of targets from the Target column of
        the sample frame.
    :param string ref_sample: The name of the sample to which expression should
        be referenced.
    :return: a DataFrame of relative expression; rows represent rows of the
        sample_frame and columns represent each of the ref_targets.
    :rtype: DataFrame
    """
    by_gene = {'Sample': sample_frame['Sample'], 'Target': sample_frame['Target']}
    for target in ref_targets:
        by_gene[target] = expression_ddcq(sample_frame, target, ref_sample)
    return pd.DataFrame(by_gene)

def rank_targets(sample_frame, ref_targets, ref_sample):
    """Uses the geNorm algorithm to determine the most stably expressed
    genes from amongst ref_targets in your sample.

    See Vandesompele et al.'s 2002 Genome Biology paper for information about
    the algorithm: http://dx.doi.org/10.1186/gb-2002-3-7-research0034

    :param DataFrame sample_frame: A sample data frame.
    :param iterable ref_targets: A sequence of targets from the Target column
        of sample_frame to consider for ranking.
    :param string ref_sample: The name of a sample from the Sample
        column of sample_frame. It doesn't really matter what it is but it
        should exist for every target.
    :return: a sorted DataFrame with two columns, 'Target' and 'M' (the
        relative stability; lower means more stable).
    :rtype: DataFrame
    """
    table = collect_expression(sample_frame, ref_targets, ref_sample)
    all_samples = sample_frame['Sample'].unique()
    t = table.groupby(['Sample', 'Target']).mean()
    logt = log2(t)
    ref_targets = set(ref_targets)

    worst = []
    worst_m = []
    while len(ref_targets) - len(worst) > 1:
        M = []
        for test_target in ref_targets:
            if test_target in worst: continue
            Vs = []
            for ref_target in ref_targets:
                if ref_target == test_target or ref_target in worst: continue
                A = logt.ix[zip(all_samples, repeat(test_target)), ref_target]
                Vs.append(A.std())
            M.append( (sum(Vs)/(len(ref_targets)-len(worst)-1), test_target) )
        worst.append(max(M)[1])
        worst_m.append(max(M)[0])
    best = ref_targets - set(worst)
    worst.reverse()
    worst_m.reverse()
    worst_m = [worst_m[0]] + worst_m
    return pd.DataFrame({'Target': list(best) + worst, 'M': worst_m}, columns=['Target', 'M'])

def calculate_all_nfs(sample_frame, ranked_targets, ref_sample):
    """For a set of n ranked_genes, calculates normalization factors NF_1,
    NF_2, ..., NF_n. NF_i represents the normalization factor generated by
    considering the first i targets in ranked_targets.
    
    calculate_nf (which returns only NF_n) is probably more
    useful for routine analysis.

    :param DataFrame sample_frame: A sample data frame.
    :param iterable ranked_targets: A list or Series of target names, in order
        of descending stability (ascending M).
    :param string ref_sample: The name of the sample to normalize against.
    :return: a DataFrame with columns 1, 2, ..., n containing normalization
        factors NF_1, ..., NF_n for each sample, indexed by sample name.
    :rtype: DataFrame
    """

    # Returns a DataFrame, where rows represent samples and columns represent a number of reference genes.
    grouped = sample_frame.groupby(['Target', 'Sample'])['Cq'].aggregate(average_cq)
    samples = sample_frame['Sample'].unique()
    nfs = {}
    for i in xrange(1, len(ranked_targets)+1):
        nfs[i] = gmean([pow(2, -grouped.ix[zip(repeat(ref_gene), samples)] + grouped.ix[ref_gene, ref_sample]) for ref_gene in ranked_targets[:i]])
    return pd.DataFrame(nfs, index=samples)

def calculate_nf(sample_frame, ref_targets, ref_sample):
    """Calculates a normalization factor from the geometric mean of the
    expression of all ref_targets, normalized to a reference sample.

    :param DataFrame sample_frame: A sample data frame.
    :param iterable ref_targets: A list or Series of target names.
    :param string ref_sample: The name of the sample to normalize against.
    :return: a Series indexed by sample name containing normalization factors
        for each sample.
    """
    grouped = sample_frame.groupby(['Target', 'Sample'])['Cq'].aggregate(average_cq)
    samples = sample_frame['Sample'].unique()
    nfs = gmean([pow(2, -grouped.ix[zip(repeat(ref_gene), samples)] + grouped.ix[ref_gene, ref_sample]) for ref_gene in ref_targets])
    return pd.Series(nfs, index=samples)

def calculate_v(nfs):
    """Calculates V(n+1/n) values. Useful for establishing the quality of
    your normalization regime. See Vandesompele 2002 for advice on
    interpretation.

    :param DataFrame nfs: A matrix of all normalization factors, produced by
       `calculate_all_nfs`.
    :return: a Series of values [V(2/1), V(3/2), V(4/3), ...].
    """
    v = []
    if (nfs.columns != range(1, nfs.columns[-1]+1)).any():
        raise ValueError("Column names invalid in nf_v_frame")
    for i in nfs.columns[:-1]:
        v.append(std(log2(nfs[i]/nfs[i+1]), ddof=1))
    return pd.Series(v, index=nfs.columns[:-1])

"""
This function may return, someday. But not yet.

def recommend_refset(sample_list, ref_genes, ref_sample):
    ranked_genes = rank_genes(sample_list, ref_genes, ref_sample)
    nfs = calculate_all_nfs(sample_list, ref_genes, ref_sample)
    vs = nf_v(nfs)
    rec = [(ranked_genes[0], 0)]
    for v in sorted(vs.index):
        if v > 3 and vs[v-1] < 0.15: break
        rec.append((ranked_genes[v-1], vs[v]))
    return rec
"""
