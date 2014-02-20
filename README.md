# eleven
Tim D. Smith, [@biotimylated](https://twitter.com/biotimylated)

Eleven is a Python library for performing multi-gene RT-qPCR gene expression normalization. It is a free, open-source implementation of the [GeNorm algorithm](http://dx.doi.org/10.1186/gb-2002-3-7-research0034) described by Vandesompele et al. in 2002.

[Documentation](http://eleven.readthedocs.org) is hosted at Read the Docs.

## How do I use eleven?

Eleven requires Python 2.7. Earlier versions will not be supported. Python 3.x support is on the roadmap. You will need a Scientific Python stack, including pandas and scipy. If you don't have these, you can install the free version of the [Anaconda environment](https://store.continuum.io/cshop/anaconda/), which has everything you need.

A sample analysis session looks like this:

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

Wasn't that easy? This adds the relative expression of each well as a column of the data frame. Now you can use regular pandas tools for handling the data, so `censored.groupby(['Sample', 'Target'])['RelExp'].aggregate(['mean', 'std'])` gives you a nice table of means and standard deviations for each target in each sample.

## Isn't Gapdh/Actb/Rn18s good enough?

If you're expecting 40-fold changes in your experiments, normalizing against a single "usual suspect" reference gene will probably do it.

But if you're interested in reliably measuring smaller changes, remember that the quality of your results cannot be better than the quality of your normalization. Without at least assessing the stability of your favorite reference gene under your experimental conditions against a panel of other genes that are likely to be more or less stably expressed, the systematic error of your comparison is totally uncontrolled. __Unless you show your reference gene is quantitatively stable, you have no evidence you are running a quantitative experiment.__

## Why GeNorm?

Several algorithms have been proposed and are in use for selecting an ensemble of stably expressed targets from a panel of candidate reference genes. GeNorm is one of the older and more popular algorithms. A [2009 review by Vandesompele, Kubista, and Pfaffl](http://www.gene-quantification.de/Vandesompele-Kubista-Pfaffl-real-time-PCR-chapter-4.pdf) explains the mathematical basis behind several normalization algorithms and concludes while "every scientist should at least validate their reference genes, the actual method used \[to normalize genes\] is less critical" since they give "highly similar rankings."

Adding other algorithms isn't a priority for me but I'll gladly accept pull requests supported by regression tests.

## Why should I use eleven?

Eleven has a simple, clean interface and uses familiar data structures. Also, I think we're the only game in town for PCR analysis in Python.

There are other options in R; [SLqPCR](http://www.bioconductor.org/packages/devel/bioc/html/SLqPCR.html) is probably the most kindred to eleven. [qpcR](http://www.dr-spiess.de/qpcR.html) does a number of very sophisticated things but I found it correspondingly mysterious. [But I don't like R](http://tim-smith.us/arrgh/).

## Why is it named eleven?

PCR is amplification based. [Our amplifier goes to 11.](https://en.wikipedia.org/wiki/Up_to_eleven)

