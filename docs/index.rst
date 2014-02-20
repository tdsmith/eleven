eleven: RT-qPCR normalization in Python
=======================================
Tim D. Smith, `@biotimylated <https://twitter.com/biotimylated>`_.

Eleven is a Python library for performing multi-gene RT-qPCR gene expression normalization. It is a free, open-source implementation of the `GeNorm algorithm <http://dx.doi.org/10.1186/gb-2002-3-7-research0034>`_ described by Vandesompele et al. in 2002.

Eleven requires Python 2.7. Earlier versions will not be supported. Python 3.x support is on the roadmap. You will need a Scientific Python stack, including pandas and scipy. If you don't have these, you can install the free version of the `Anaconda environment <https://store.continuum.io/cshop/anaconda/>`_, which has everything you need.

A sample analysis session looks like this::

    # Read PCR data into a pandas DataFrame. You want a data file where each
    # row corresponds to a separate well, with columns for the sample name,
    # target name, and Cq value. NTC wells should have the sample name set to
    # a value like 'NTC'.
    >> df = pd.read_csv('my_data.csv')

    # If your Sample, Target, and Cq columns are called other things, they
    # should be renamed to Sample, Target, and Cq.
    >> df = df.rename(columns={'Gene': 'Target', 'Ct': 'Cq'})

    # Drop the wells that are too close to the NTC for that target.
    >> censored = eleven.censor_background(df)

    # Rank your candidate reference genes.
    >> ranked = eleven.rank_targets(censored, ['Gapdh', 'Rn18s', 'Hprt',
        'Ubc', 'Actb'], 'Control')

    # Normalize your data by your most stable genes and compute normalization
    # factors (NFs).
    >> nf = eleven.calculate_nf(censored, ranked.ix['Target', 0:3], 'Control')

    # Now, normalize all of your expression data.
    >> censored['RelExp'] = eleven.expression_nf(censored, nf, 'Control')

Wasn't that easy? This adds the relative expression of each well as a column of the data frame. Now you can use regular pandas tools for handling the data, so ::

    censored.groupby(['Sample', 'Target'])['RelExp'].aggregate(['mean', 'std'])

gives you a nice table of means and standard deviations for each target in each sample.

.. toctree::
   :maxdepth: 2

   eleven


