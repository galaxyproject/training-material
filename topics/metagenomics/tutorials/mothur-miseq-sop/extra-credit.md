### Can't get enough?

Below are some more exercises for those who wish to go into more detail about statistical significance testing and population-level analysis.

Or, click [here](#wrap-up) to jump to the end and skip this extra credit section.


# Extra Credit

## Determine statistical significance of clusterings

We can perform a test to determine whether the clustering within the tree is statistically significant or not
using by choosing from the `parsimony`, `unifrac.unweighted`, or `unifrac.weighted` commands. To run these we
will first need to create a design file that indicates which treatment each sample belongs to.

> <hands-on-title>Obtain design file</hands-on-title>
>
> - Find the file `mouse.time.design` in your history (you imported this file at the start of this tutorial)
> - Make sure the datatype is set to `mothur.design`.
>
> > <tip-title>Changing datatype of a datasets</tip-title>
> >  - Click on the **pencil icon** of the dataset
> >  - Click on the **Datatypes** tab
> >  - Select the new datatype from dropdown menu
> >  - Click **Save**
> {: .tip}
{: .hands_on}


The design file look something like this:

```
group    time
F3D0     Early
F3D1     Early
F3D141   Late
F3D142   Late
F3D143   Late
F3D144   Late
F3D145   Late
F3D146   Late
F3D147   Late
F3D148   Late
F3D149   Late
F3D150   Late
F3D2     Early
F3D3     Early
F3D5     Early
F3D6     Early
F3D7     Early
F3D8     Early
F3D9     Early
```

Using the `parsimony` command let's look at the pairwise comparisons. Specifically, let's focus on the
early vs. late comparisons for each mouse:

> <hands-on-title>Compare Early-vs-Late</hands-on-title>
> - **Parsimony** {% icon tool %} with the following parameters
>   - "tree" to the `tre` output from Tree.Shared (collection)
>   - "group" to the design file described above
>   - "output logfile?" to `yes`
{: .hands_on}

In the logfile for `thetayc.0.03.lt.ave` we see

```
Tree#   Groups      ParsScore   ParsSig
1       Early-Late  1           <0.001
```

There was clearly a significant difference between the clustering of the early and late time points.
Recall that this method ignores the branch length.

The two distance matrices that we generated earlier (i.e. `jclass.0.03.lt.ave.dist` and
    `thetayc.0.03.lt.ave.dist`) can then be visualized using the pcoa or nmds plots.

Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional
data in as few dimensions as possible. Our data is highly dimensional (~9 dimensions).

> <hands-on-title>PCoA</hands-on-title>
>
> - **Pcoa** {% icon tool %} with the following parameters
>   - "phylip" to dist files from Dist.shared (collection)
{: .hands_on}

The loadings files will tell you what fraction of the total variance in the data are represented
by each of the axes. For instance the loading file for `thetayc.0.03.lt.ave` looks something like:

```
axis  loading
1     45.354207
2     13.526582
3     11.791424
4     4.493544
5     4.012474
...
```

In this case the first and second axis represent about 45 and 14% of the variation (59% of the total)
for the thetaYC distances. The output to the logfile:

```
Processing...
Rsq 1 axis: 0.736369
Rsq 2 axis: 0.882025
Rsq 3 axis: 0.978093
```

indicates that the R-squared between the original distance matrix and the distance between the points in 2D
PCoA space was 0.88, but that if you add a third dimension the R-squared value increases to 0.98. All in all,
not bad.

Alternatively, non-metric multidimensional scaling (NMDS) tries to preserve the distance between samples using
a user defined number of dimensions. We can run our data through NMDS with 2 dimensions with the following
tool:

> <hands-on-title>Nmds</hands-on-title>
>
> - **Nmds** {% icon tool %} with the following parameters
>   - "phylip" to dist files from Dist.shared (collection)
>   - "output logfile?" to `yes`
>
> Opening the `stress` file for `thetayc.0.03.lt.ave` we can inspect the stress and R^2 values, which describe
> the quality of the ordination. Each line in this file represents a different iteration and the configuration
> obtained in the iteration with the lowest stress is reported in the `axes` file. In the logfile:
>
> ```
> Number of dimensions:           2
> Lowest stress :                 0.113657
> R-squared for configuration:    0.947622
> ```
>
> We find that the lowest stress value was 0.11 with an R-squared value of 0.95; that stress level is
> actually pretty good. You can test what happens with three dimensions in the following way:
>
> - **Nmds** {% icon tool %} with the following parameters
>   - "phylip" to dist files collection from Dist.shared
>   - "mindim" to `3`
>   - "maxdim" to `3`
>   - "output logfile?" to `yes`
>
> > <question-title></question-title>
> >
> > What are stress and R-squared values when using 3 dimensions?
> >
> > > <solution-title></solution-title>
> > > The stress value drops to 0.05 and the R2 value goes up to 0.99 (see logfile). Not bad.
> > {: .solution }
> {: .question}
{: .hands_on}



In general, we would like a stress value below 0.20 and a value below 0.10 is even better. Thus, we can conclude that,
NMDS is better than PCoA. We can plot the three dimensions of the NMDS data by plotting the contents of the `axes`
file. <!-- TODO: tool for 3D plots in Galaxy? -->

Again, it is clear that the early and late samples cluster separately from each other. Ultimately, ordination
is a data visualization tool. We might ask if the spatial separation that we see between the early and late
plots in the NMDS plot is statistically significant. To do this we have two statistical tools at our disposal.
The first analysis of molecular variance (AMOVA), tests whether the centers of the clouds representing a group
are more separated than the variation among samples of the same treatment. This is done using the distance
matrices we created earlier and does not actually use ordination.

> <hands-on-title>Amova</hands-on-title>
>
> - **Amova** {% icon tool %} with the following parameters
>   - "phylip" to dist files from Dist.shared (collection)
>   - "design" to mouse.time.design file from your history
>   - "output logfile?" to `yes`
{: .hands_on}

in logfile for thetaYC we find:

```
Early-Late    Among       Within     Total
SS            0.628379    0.552221   1.1806
df            1           17         18
MS    0.628379    0.0324836

Fs:    19.3445
p-value: <0.001*
```

Here we see from the AMOVA that the "cloud" early and late time points has a significantly different centroid
for this mouse. Thus, the observed separation in early and late samples is statistically significant. We can
also see whether the variation in the early samples is significantly different from the variation in the late
samples using the `Homova` command:

> <hands-on-title>Homova</hands-on-title>
>
> - **Homova** {% icon tool %} with the following parameters
>   - "phylip" to dist files from Dist.shared (collection)
>   - "design" to mouse.time.design file from your history
>   - "output logfile?" to `yes`
{: .hands_on}

```
HOMOVA        BValue     P-value    SSwithin/(Ni-1)_values
Early-Late    7.51408    <0.001*    0.0603208    0.00773943
```

We see that there is a significant difference in the variation with the early samples having a larger amount
of variation (0.061) than the late samples (0.008). This was what we found in the original study - the early
samples were less stable than the late samples.

Next, we might ask which OTUs are responsible for shifting the samples along the two axes. We can determine
this by measuring the correlation of the relative abundance of each OTU with the two axes in the NMDS dataset.
We do this with the `corr.axes` tool:

> <hands-on-title>Correlation</hands-on-title>
>
> - **Corr.axes** {% icon tool %} with the following parameters
>   - "axes" to axes output from Nmds in 3 dimension (collection)
>   - "shared" to shared output from collapse collection on Sub.sample
>   - "method" to `Spearman`
>   - "numaxes" to `3`
{: .hands_on}

Examining the axes output, we see the data for the first five OTUs look something like this..

```
OTU         axis1       p-value      axis2       p-value     axis3       p-value     length
Otu0001     0.285213    0.226258    -0.742431    0.000272    0.676613    0.001466    1.044201
Otu0002     0.283582    0.228923    -0.636524    0.003387    0.873574    0.000001    1.117458
Otu0003     0.461270    0.046828    -0.586271    0.008337    0.767610    0.000125    1.070378
Otu0004    -0.131579    0.576679    -0.240351    0.307860    0.408772    0.082266    0.492114
Otu0005    -0.315327    0.180955     0.046553    0.843432    0.097497    0.679135    0.333323
...
```

What these results show is that OTUs 1 and 2 are responsible for moving points in a negative direction along
axis 2. Recalling that we classified each OTU earlier (see taxonomy output from `Classify.otu`), we can see
that these first five OTUs are mainly members of the Porphyromonadaceae:

```
OTU        Size   Taxonomy
Otu0001    12329   Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
Otu0002    8912    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
Otu0003    7857    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
Otu0004    7483    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);Barnesiella(100);
Otu0005    7479    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
...
```

This helps to illustrate the power of OTUs over phylotypes since each of these OTUs is behaving differently.
These data can be plotted in what's known as a biplot where lines radiating from the origin (axis1=0, axis2=0,
axis3=0) to the correlation values with each axis are mapped on top of the PCoA or NMDS plots.
<!-- TODO: make this plot? -->

Later, using the metastats command, we will see another method for describing which populations are
responsible for differences seen between specific treatments.

An alternative approach to building a biplot would be to provide data indicating metadata about each sample.
For example, we may know the weight, height, blood pressure, etc. of the subjects in these samples. For
discussion purposes the file `mouse.dpw.metadata` is provided and looks something like this:

```
group    dpw
F3D0     0
F3D1     1
F3D141   141
F3D142   142
F3D143   143
F3D144   144
F3D145   145
F3D146   146
F3D147   147
F3D148   148
F3D149   149
F3D150   150
F3D2     2
F3D3     3
F3D5     5
F3D6     6
F3D7     7
F3D8     8
F3D9     9
```

> <hands-on-title>Hands-on</hands-on-title>
>
> - **Corr.axes** {% icon tool %} with the following parameters
>   - "axes" to axes output from Nmds in 3 dimension
>   - "Generate Collector Curvers for" to Metadata table
>   - "metadata table" to `mouse.dpw.metadata`
>   - "method" to `Spearman`
>   - "numaxes" to `3`
>
> This will output a file like the following:
>
> ```
> Feature    axis1       p-value      axis2       p-value     axis3       p-value     length
> dpw        0.205263    0.383832    -0.292982    0.213861    0.821053    0.000016    0.895600
> ```
>
> Indicating that as the dpw increases, the communities shift to in the positive direction along axis 3.
>
> Another tool we can use is `get.communitytype` to see whether our data can be partitioned in to separate
> community types
>
> <!-- TODO: add this tool to mothur suite -->
> - **Get.communitytype** {% icon tool %} with the following parameters
>   - "shared" to Subsample.shared file
>   - "output logfile?" to `yes`
>
{: .hands_on}

In logfile we find the following output:

```
K    NLE        logDet    BIC         AIC         Laplace
1    9612.15    522.97    10070.01    9923.15     9587.84
2    9688.76    464.05    10605.95    10311.76    9348.28
3    10329.39   329.18    11705.91    11264.39    9634.77
4    11026.12   97.78     12861.98    12273.12    9929.10
5    11662.52  -250.61    13957.71    13221.52    10104.59
```

We see that the minimum Laplace value is for a K value of 2 (9348.28). This indicates that our samples
belonged to two community types. Opening the `design` output we see that all of the late samples and the Day 0
sample belonged to Partition_1 and the other early samples belonged to Partition_2. We can look at the
`summary` output to see which OTUs were most responsible for separating the communities:

```
OTU        P0.mean  P1.mean  P1.lci  P1.uci  P2.mean  P2.lci  P2.uci  Difference   CumFraction
Otu0006    3.36     10.48    9.17    11.97   0.46     0.28    0.78    10.01        0.15
Otu0014    6.17     8.45     7.35    9.72    3.76     2.98    4.73    4.70         0.22
Otu0002    5.63     7.14     6.17    8.25    3.83     3.05    4.81    3.31         0.27
Otu0008    4.01     2.92     2.41    3.54    5.85     4.80    7.12    2.92         0.31
Otu0019    2.07     3.48     2.90    4.18    0.94     0.63    1.40    2.54         0.35
...
```

Again we can cross reference these OTU labels with the consensus classifications in the taxonomy file to get
the names of these organisms.

> <question-title></question-title>
>
> What organisms were the top 5 contributing OTUs classified as?
>
> > <solution-title></solution-title>
> > Note down the names of the top 5 OTUs as output by thesummary output of get.communitytype.
> > Then look at the taxonomy file output by Classify.otu.
> >
> > In our example these top 5 OTUs were classified
> > as belonging to Porphyromonadaceae (top 3 OTUs), Alistipes and Lactobacillus.
> {: .solution }
{: .question}

## Population-level Analysis

In addition to the use of `corr.axes` and `get.communitytype` we have several tools to differentiate between
different groupings of samples. The first we'll demonstrate is `metastats`, which is a non-parametric T-test
that determines whether there are any OTUs that are differentially represented between the samples from early and late in this study.

> <hands-on-title>T-test</hands-on-title>
>
> - **Metastats** {% icon tool %} with the following parameters
>   - "shared" to Subsample.shared
>   - "design" to `mouse.time.design`
{: .hands_on}

Looking at the first 5 OTUs from `Late-Early` output file we see the following:

```
OTU        mean(group1)  variance(group1)  stderr(group1)  mean(group2)  variance(group2)  stderr(group2)  p-value
Otu0001    0.026104      0.000079          0.002807        0.011304      0.000031          0.001856        0.000999
Otu0002    0.072869      0.000101          0.003176        0.041946      0.000208          0.004805        0.000999
Otu0003    0.015261      0.000023          0.001531        0.002182      0.000003          0.000539        0.000999
Otu0004    0.029451      0.000064          0.002536        0.020427      0.000140          0.003947        0.074925
Otu0005    0.068139      0.000087          0.002957        0.070058      0.000163          0.004254        0.729271
```

These data tell us that OTUs 1, 2, and 3 was significantly different between the early and late samples.

> <question-title></question-title>
>
>  Which of the top 10 OTUs in your output were significantly different between early and late samples?
>
> > <solution-title></solution-title>
> > Looking at the p-value cut-off and using your favorite cutoff threshold (say 0.01).
> > Answer to the question is all OTUs with a value lower than this threshold. Note that these OTU labels may
> > be different for you and may very between one repetition of this tutorial to the next, and therefore may
> > vary between you and your neighbour as well.
> {: .solution }
{: .question}

Another non-parametric tool we can use as an alternative to metastats is lefse:

> <hands-on-title>Lefse</hands-on-title>
>
> - **Lefse** {% icon tool %} with the following parameters
>   - "shared" to Subsample.shared
>   - "design" to `mouse.time.design`
{: .hands_on}

Looking at the top of the lefse summary file we see:

```
OTU        LogMaxMean  Class   LDA         pValue
Otu0001    4.41671     Late    3.91585    0.000601825
Otu0002    4.86254     Late    4.20329    0.000695271
Otu0003    4.18358     Late    3.82749    0.00022674
Otu0004    4.4691      -
Otu0005    4.84546     -
```

Again, OTUs 1, 2, and 3 are significantly different between the two groups and are significantly elevated in the
late samples

Finally, Mothur has an implementation of the random forest algorithm build into her as classify.rf. This will tell
us which features (i.e. OTUs) are useful in discriminating between the two groups of samples:

> <hands-on-title>Classify.rf</hands-on-title>
>
> - **Classify.rf** {% icon tool %} with the following parameters
>   - "shared" to Subsample.shared
>   - "design" to `mouse.time.design`
{: .hands_on}

in the logfile we see:

```
Creating 100 (th) Decision tree
numCorrect = 19
forrestErrorRate = 0
confusion matrix:
        Early    Late    time
Early   9        0       0
Late    0        10      0
time    0        0       0
```

We can ignore the time row and column and see that our samples were all correctly assigned to the proper groups.
Looking at `summary` output, we see the top 10 OTUs that resulted in the greatest mean decrease in activity were:

```
OTU        Mean decrease accuracy
Otu0038    0.21
Otu0003    0.15
Otu0091    0.14
Otu0096    0.13
Otu0024    0.12
Otu0006    0.1
Otu0011    0.1
Otu0015    0.09
Otu0082    0.08
Otu0042    0.07
```

# Wrap-up
{:.no_toc}

