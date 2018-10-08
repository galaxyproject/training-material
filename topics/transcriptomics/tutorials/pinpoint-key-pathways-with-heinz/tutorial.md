---
layout: tutorial_hands_on

title: "Metatranscriptomics: Pinpoint key pathways using Heinz"
zenodo_link: "https://doi.org/10.5281/zenodo.1344105"
questions:
    - "Which pathways are potentially contributing to dental caries?"
objectives:
    - "Analyze metatranscriptomics data using Heinz in Galaxy to pinpoint the optimal scoring subnetwork."
time_estimation: "1-2h"
key_points:
    - "Analyzing differential expression for (meta)transcriptomics data."
    - "Validating the p-value distribution of the differential expression analysis."
    - "Finding the most differentially expressed network in a gene functional network using Heinz."
    - "Analyzing the pinpointed network."
contributors:
    - cicozhang
    - sanneabeln
---

# Overview
{:.no_toc}

The human microbiome plays a key role in health and disease. Thanks to comparative metatranscriptomics,
the cellular functions that are deregulated by the microbiome in disease can now be computationally
explored. Unlike gene-centric approaches, pathway-based methods provide a systemic view of such
functions; however, they typically consider each pathway in isolation and in its entirety.
They can therefore overlook the key differences that (i) span multiple pathways, (ii) contain
bidirectionally deregulated components, (iii) are confined to a pathway region. To capture these
properties, computational methods that reach beyond the scope of predefined pathways are needed.

In this tutorial, we will perform a network analysis using [Heinz](https://github.com/ls-cwi/heinz) in Galaxy. The data comes from the study [May et al.](https://academic.oup.com/bioinformatics/article/32/11/1678/2240171), and we will reproduce some of the computational steps from this study with simplified data and parameters to speed up the analysis for the purposes of this tutorial.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


> ### {% icon comment %} Note
> Your results may be slightly different from the ones presented in this tutorial due to the tool versions or
> reference data versions or stochastic processes in the algorithms.
{: .comment}

In this tutorial, we will create a Heinz workflow step by step, as the picture below shows.

![Heinz workflow](../../images/heinz-workflow.png)

# Obtaining and preparing data  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1344105.svg)](https://doi.org/10.5281/zenodo.1344105)


The study [May et al.](https://academic.oup.com/bioinformatics/article/32/11/1678/2240171) includes the computation
steps starting from the raw RNAseq datasets. The operations that processed raw data into the interpreted data are
beyond the scope of this tutorial. To learn that, please refer to the relevant topics in the Galaxy training
material. In this tutorial, we start with the interpreted data, which are KO
([KEGG Orthology](https://www.genome.jp/kegg/ko.html)) count data. All the data needed for this tutorial are
available from Zenodo.

## Understanding our input data

According to the study [May et al.](https://academic.oup.com/bioinformatics/article/32/11/1678/2240171),
dental caries (DC) dataset in this tutorial came from an experiment that comprised of supragingival plaque
samples collected from all dental surfaces of 36 individuals who had either a caries-positive (disease) or
a caries-negative (health) oral profile. Each of the 36 samples was sequenced, pre-processed and transformed
into KO counts. We will use these count data as the starting point to perform the network analysis.

> ### {% icon comment %} Dataset details
> The count data of the 36 samples are separated into 36 files, organized into two groups:
> 'CP' (caries-positive) and 'CN' (caries-negative).
{: .comment}


## Importing the data into Galaxy

After knowing what our input data are like, let's get them into Galaxy history:

> ### {% icon hands_on %} Hands-on: Obtaining our data
>
> 1. Make sure you have an empty Galaxy history. Give it a sensible name.
>
>    {% include snippets/history_create_new.md %}
>
> 2. **Upload Disease Dataset**
>    - Open the file upload menu
>    - Click on  **Collection** tab
>    - Click on the **Paste/Fetch data** button
>    - Copy the Zenodo links for the Disease Datasets
>      > ### {% icon details %} View list of Zenodo URLs for Dental Caries Dataset (Disease, CP)
>      > ```
>      > https://zenodo.org/record/1344105/files/2241_CP_DZ_PairTo_2242.txt
>      > https://zenodo.org/record/1344105/files/2126_CP_MZ_PairTo_2125.txt
>      > https://zenodo.org/record/1344105/files/2991_CP_DZ_PairTo_2992.txt
>      > https://zenodo.org/record/1344105/files/2931_CP_DZ_PairTo_2930.txt
>      > https://zenodo.org/record/1344105/files/2284_CP_DZ_PairTo_2283.txt
>      > https://zenodo.org/record/1344105/files/2125_CP_MZ_PairTo_2126.txt
>      > https://zenodo.org/record/1344105/files/4131_CP_DZ_PairTo_4132.txt
>      > https://zenodo.org/record/1344105/files/2954_CP_DZ_PairTo_2955.txt
>      > https://zenodo.org/record/1344105/files/2170_CP_MZ_PairTo_2169.txt
>      > https://zenodo.org/record/1344105/files/2955_CP_DZ_PairTo_2954.txt
>      > https://zenodo.org/record/1344105/files/2011_CP_DZ_PairTo_2012.txt
>      > https://zenodo.org/record/1344105/files/2012_CP_DZ_PairTo_2011.txt
>      > https://zenodo.org/record/1344105/files/2269_CP_DZ_PairTo_2270.txt
>      > https://zenodo.org/record/1344105/files/3215_CP_MZ_PairTo_3214.txt
>      > https://zenodo.org/record/1344105/files/2354_CP_DZ_PairTo_2355.txt
>      > https://zenodo.org/record/1344105/files/3306_CP_DZ_PairTo_3307.txt
>      > https://zenodo.org/record/1344105/files/2061_CP_DZ_PairTo_2062.txt
>      > https://zenodo.org/record/1344105/files/2355_CP_DZ_PairTo_2354.txt
>      > https://zenodo.org/record/1344105/files/2242_CP_DZ_PairTo_2241.txt
>      > ```
>      {: .details}
>    - Click on **Start** to Upload the files
>    - Click **Build** once upload has completed
>    - Enter **Name**: `CN`
>    - Click **Create list** to make the collection
>
> 3. **Upload control (healthy) datasets**
>    - Repeat the previous steps with the samples from healthy individuals:
>      > ### {% icon details %} View list of Zenodo URLs for Dental Caries Dataset (Healthy, CN)
>      > ```
>      > https://zenodo.org/record/1344105/files/2310_CN_DZ_PairTo_2309.txt
>      > https://zenodo.org/record/1344105/files/2062_CN_DZ_PairTo_2061.txt
>      > https://zenodo.org/record/1344105/files/2191_CN_MZ_PairTo_2192.txt
>      > https://zenodo.org/record/1344105/files/2052_CN_MZ_PairTo_2051.txt
>      > https://zenodo.org/record/1344105/files/2051_CN_MZ_PairTo_2052.txt
>      > https://zenodo.org/record/1344105/files/2192_CN_MZ_PairTo_2191.txt
>      > https://zenodo.org/record/1344105/files/2234_CN_DZ_PairTo_2233.txt
>      > https://zenodo.org/record/1344105/files/2233_CN_DZ_PairTo_2234.txt
>      > https://zenodo.org/record/1344105/files/2270_CN_DZ_PairTo_2269.txt
>      > https://zenodo.org/record/1344105/files/2225_CN_MZ_PairTo_2226.txt
>      > https://zenodo.org/record/1344105/files/4132_CN_DZ_PairTo_4131.txt
>      > https://zenodo.org/record/1344105/files/2309_CN_DZ_PairTo_2310.txt
>      > https://zenodo.org/record/1344105/files/2992_CN_DZ_PairTo_2991.txt
>      > https://zenodo.org/record/1344105/files/3214_CN_MZ_PairTo_3215.txt
>      > https://zenodo.org/record/1344105/files/2169_CN_MZ_PairTo_2170.txt
>      > https://zenodo.org/record/1344105/files/2930_CN_DZ_PairTo_2931.txt
>      > https://zenodo.org/record/1344105/files/3307_CN_DZ_PairTo_3306.txt
>      > ```
>      {: .details}
>    - Name this collection: `CN`
>
>    > ### {% icon question %} Question
>    >
>    > 1. How many samples do you have in your disease collection (CP)? How many healthy samples (CN)?
>    > 2. How many columns in each file? What are these columns?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. You should have 19 samples in the disease collection (CP), and 17 in the negative collection (CN).
>    > > 2. There are two columns, one is the KO IDs, the other is the count.
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

> ### {% icon tip %} Tip: Creating a collection from files already in your history
>
> *Dataset collections* enables us to easily run tools on multiple datasets at once, you probably
> have done that in the file upload menu using the `Collection` tab. If not, you can still create
> *dataset collections* manually.
>
> > ### {% icon hands_on %} Manually organizing our data into a collection
> >
> > If you have all the datasets in the history, but they are not organized into a collection yet,
> > you can follow these steps to create a collection:
> >
> > 1. Click on the **checkmark icon** at top of your history.
> >   ![Checkmark icon in history menu](../../../../shared/images/history_menu_buttons2.png)
> >
> > 2. Select all the files whose name contains `CP`, then click on **for all selected..** and select
> >   **Build Dataset List** from the dropdown menu.
> >
> > 3. In the next dialog window, you need to give a name, here we just set it to `CP`, then click **Create list**.
> >   ![List of suggested paired datasets](../../images/create_collection.png)
> >
> > 4. Hidden these selected files by clicking on **for all selected..** and selecting **Hidden datasets**.
> >   **Note:** This step is optional, we do it here to keep Galaxy history clean.
> >
> > 5. Redo the Step 2, 3, 4 for CN, set the name of the data list as 'CN'.
> {: .details}
{: .tip}



# Differential Expression Analysis (DEA) by DESeq2

## What is differential expression analysis?

> ### {% icon comment %} A defintion
>
> The definition of differential expression analysis given by
> [EBI](https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene)
> means taking the normalised read count data and performing statistical analysis to discover quantitative
> changes in expression levels between experimental groups. For example, we use statistical testing to decide
> whether, for a given gene, an observed difference in read counts is significant, that is, whether it is greater
> than what would be expected just due to natural random variation.
>
{: .comment}

In principle, DEA is a causal analysis; but in reality, it is hampered by the complexity of the experimental
situation and measurement. Back to our datasets, CP and CN, they are from two experimental groups. By DEA,
we hope to pinpoint the candidate genes relevant to dental caries first, then we will use Heinz to infer the
related pathways.

## Which tools are available for DEA?

There are a few canned tools commonly used for DEA, like [Limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). If you are interested, you
may look up the pros and cons of each tool. Here we use DESeq2.


## Conduct differential expression analysis

After learning about differential expression analysis, let's get some hands-on experience.

> ### {% icon hands_on %} Hands-on: DEA via DESeq2
>
> - **DESeq2** {% icon tool %} with the following parameters
>   - "Specify a factor name" to `dental_caries` (in "1: Factor")
>   - "Specify a factor level" to `CP`, "Counts file(s)" to `CP` by clicking the "Dataset collections" icon (in "1: Factor level")
>   - "Specify a factor level" to `CN`, "Counts file(s)" to `CN` by clicking the "Dataset collections" icon (in "2: Factor level")
>   - "Files have header" to `No`
>   - "Visualising the analysis results" to `No`
>   - Leave all other parameters to the default settings <br><br>
>
{: .hands_on}

It usually takes 10 - 15 minutes to finish DESseq2. After the analysis is finished, have a look at the file, it should look like something this:

```
GeneID  Base mean   log2(FC)	        StdErr	            Wald-Stats	        P-value	                P-adj
K15661  55.0733128361562    2.49892393773464	0.508475451930939	4.91454194739384	8.89902710599613e-07	0.00491938218419466
K15666  48.4532561772761    2.22805029428311	0.493391718412049	4.51578372951607	6.30830206230236e-06	0.00820667379798327
K10213  25.1966742274619    -2.45624670858868	0.550895251183889	-4.45864563782342	8.24791475064012e-06	0.00820667379798327
K03732  563.63634258472 -1.41961721984814	0.316992240900184	-4.47839737596338	7.52055094378821e-06	0.00820667379798327
K01792  112.923146882195    -1.26925892659617	0.285732234964578	-4.44212717810268	8.90738834802815e-06	0.00820667379798327
```


# Fit a BUM model (a mixture model)

From a statistical point of view, p-values are uniformly distributed under null hypothesis; in other words,
under alternative hypothesis, the noise component (which holds under null hypothesis) will be adequately
modeled by a uniform distribution. With this knowledge, we can fit these p-values to a mixture model
(BUM model), as the figure below shows.

![p-values are fitted to a mixture model](../../images/bum.jpeg){:width="50%"}

Before fitting to BUM model in Galaxy, we need to prepare the input data for the tool
**Fit a BUM model**, that's a file that only contains p-values.

> ### {% icon hands_on %} Hands-on: extract p-values from DESeq2 output
>
> - **cut** {% icon tool %} with the following parameters
>   - "Cut columns" to `c6` (in "1: Factor")
>   - "Delimited by" to `TAB`
>   - "From" to the output of **DESeq2** <br><br>
>
{: .hands_on}

Then we can **Fit a BUM model** now.

> ### {% icon hands_on %} Hands-on: fit the BUM model
>
> - **Fit a BUM model** {% icon tool %} with the following parameters
>   - "Input file" to the output of **cut** <br><br>
>
{: .hands_on}

# Pinpoint the key pathways with Heinz

After getting the parameters of the BUM model from the last step, we will use Heinz to pinpoint
the key pathways. Before we continue, let's figure out what Heinz is actually doing.

Heinz is an algorithm in searching an optimal subnetwork from a bigger network. You may wonder what
the networks are here. Through the previous steps, we have got a list of identities, that is a list
of gene IDs with p-values, which form the nodes  of 'the bigger network', the relations between the
nodes, that is the edges, need to be obtained from a background network, which represents a pathway
relation databases, such as [Reactome](https://reactome.org/) and [STRING](https://string-db.org/).
In this tutorial, we only use a small sample background network for demonstration purposes. The
background network is represented as edges in a txt file where each line denotes an edge as follows:

```
ACTR1B	ACVR2B
ZSWIM9	FOXP3
LGALS4	PRKX
NPTX1	CIAO1
```

Upload this edge file (hereafter we call it edge file) into the Galaxy instance.

> ### {% icon details %} View the Zenodo URL for edge file
> ```
> https://zenodo.org/record/1344105/files/edge.txt
> ```
{: .details}


## Calculate Heinz scores

As the first step, we need to calculate a Heinz score for each node, using the BUM model parameters
we obtained; meanwhile, we also need to specify an FDR value as input.

> ### {% icon comment %} What is FDR value?
>
> FDR is short for false discovery rate, which is a method of conceptualizing the rate of type I errors
> in null hypothesis testing when conducting multiple comparisons, if you are interested, view the detail
> in [Wikipedia](https://en.wikipedia.org/wiki/False_discovery_rate).
>
{: .comment}

In our case, the higher an FDR value is, the more positive nodes (regarding the Heinz scores) we get,
which means it may include a lot of false positive nodes. For different datasets and problems, we
probably need different FDR values. Here we set FDR to 0.11.

Similar to **Fit a BUM model**, we also need to prepare the input data for the tool
**Calculate a Heinz score**.

> ### {% icon question %} Question
>
> What is the requirement of the input data format for **Calculate a Heinz score**?
>
> > ### {% icon solution %} Solution
> >
> > In the user interface of the tool "Calculate a Heinz score", we see that "A node file with p-values" is needed.
> > It should contain two columns delimited by tab, one is KO ID; the other, p-value.
> >
> {: .solution}
{: .question}

> ### {% icon hands_on %} Hands-on: extract geneID and p-values from DESeq2 output
>
> - **cut** {% icon tool %} with the following parameters
>   - "Cut columns" to `c1,c6` (in "1: Factor")
>   - "Delimited by" to `TAB`
>   - "From" to the output of **DESeq2** <br><br>
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on: calculate Heinz scores
>
> - **Calculate a Heinz score** {% icon tool %} with the following parameters
>   - "A node file with p-values" to the output of **cut** (a different **cut**, see below)
>   - "FDR value" to 0.11
>   - "Choose your input type for BUM parameters" to "The output file of BUM model"
>   - "Output file of BUM model as input: lambda on the first line and alpha, the second" to the output of **Fit a BUM model** <br><br>
{: .hands_on}

## Run Heinz: pinpoint the optimal subnetwork

After getting Heinz scores, let's run Heinz program to find the optimal subnetwork from the background network which we mentioned earlier.

> ### {% icon hands_on %} Hands-on: pinpoint the optimal subnetwork
>
> - **Identify optimal scoring subnetwork** {% icon tool %} with the following parameters
>   - "File containing Heinz scores" to the output of **Calculate a Heinz score**
>   - "Edge file" to the edge file uploaded <br><br>
>
{: .hands_on}

It usually takes a few minutes to get the result, but mind you, for some tasks, it might take a few hours to get a result in practice.

> ### {% icon tip %} Tip: the running time of the program
>
> * Graph problem is way more complicated than we thought.
> * It might take a much longer time for some complicated datasets.
> * We can use multiple CPUs to accelerate the computation (for now, this function is not available in Galaxy yet),
>   but, to use that, you can install Heinz directly via [Bioconda](https://anaconda.org/bioconda/heinz) in a Linux environment.
{: .tip}

## Visualize the output: visualize the optimal subnetwork

The result we got from the last step is not very human readable. Therefore we need to visualize the output by making it
into graphs. Except the tool we will use in Galaxy, you may consider using eXamine plugin in Cytoscape for a richer
visualization.

> ### {% icon hands_on %} Hands-on: visualize the optimal subnetwork
>
> - **Visualize** {% icon tool %} with the following parameters
>   - "Heinz output file" to the output of **Identify optimal scoring subnetwork** <br><br>
>
{: .hands_on}

# Save the history into a workflow

Congrats! You have finished all the tools in Heinz workflow!

At the end of the tutorial, as a self-practice, you may save all of your correct operations into a workflow,
which you can reuse for different datasets next time.
