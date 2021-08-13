---
layout: tutorial_hands_on
enable: false

title: 'Mass spectrometry: LC-MS data processing'
zenodo_link: 'https://zenodo.org/record/5179809'
level: Introductory
questions:
- Why should one consider performing data processing steps when dealing with Metabolomics data?
- How to conduct LCMS-based untargeted metabolomic data processing using Galaxy?
objectives:
- To comprehend the key role of data processing
- To comprehend the diversity of steps necessary to perform untargeted LC-MS metabolomic data processing.
- To be familiar with the "identify then perform" approach necessary to deploy relevant data processing strategies.
time_estimation: 1H
key_points:
- Data processing is a key step in untargeted Metabolomics analyses. The question of data filtering and correction must be addressed in all projects, even thought in some cases it may lead to the decision of no action on data. In particular, blank filtering, pool variation study and signal drift correction are common aspects to consider when dealing with LC-MS.
- Although some main steps are standard, various ways to combine tools exist. Remember that depending on your context (type of samples, protocol specificities...) specific filters/normalisations may be needed, independently of standards ones.
- Tools are available in Galaxy, but do not forget that you need appropriate knowledge to decide what to use depending on your data.
contributors:
- melpetera
- workflow4metabolomics

---


# Introduction
{:.no_toc}

Metabolomics is a *-omic* science known for being one of the most closely related to phenotypes.
It involves the study of different types of matrices, such as blood, urine, tissues, in various organisms including plants.
It focuses on studying the very small molecules which are called *metabolites*, to better understand matters linked to the metabolism.

Metabolomics analyses can be quite complex to conduct, especially when dealing with untargeted approaches. 
**Liquid-Chromatography Mass Spectrometry** (LC-MS) is one of the three main technologies used to perform this kind of approach. 
Data analysis for this technology requires a large variety of steps, ranging from extracting information from the raw data, to statistical analysis and annotation. 
One of these steps is called "data processing". It takes place after the pre-processing step (extraction of the peak list from raw data) and before any statistical analysis.
You can get an overview of a complete LC-MS untargeted metabolomic workflow by following the dedicated 
training material [here](https://training.galaxyproject.org/training-material/topics/metabolomics/tutorials/lcms/tutorial.html).

After the pre-processing step, what you have at your disposal is a list of ions (variableMetadata file) and corresponding intensities (dataMatrix file).
What you may want now is to get some relevant information from your tables. However, your data may not be suitable yet for statistical analysis.
What should you do to ensure the quality of your tables? This tutorial will show you what the usual quality steps are.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Overview and data upload

Data processing covers a large range of actions. They can generally be described as transformation (*e.g.* normalisation) or filter (*e.g.* removal of unwanted ions). 
In this tutorial we will focus on three main types of processing:
- Removing "trash" signals
- Correcting intensities
- Removing signals of insufficient quality

## Galaxy tools' workflow

The different Galaxy tools that will be used in this tutorial are given in the following workflow picture.

![The picture is composed of 3 main boxes linked by arrows from left to right. In each of the first and last boxes, labelled '3 tabulars', are found 3 files named sampleMetadata, variableMetadata and dataMatrix. In the middle box labelled 'Quality Control' is found the workflow used in this tutorial, represented as tiny boxes with tool names, linked with arrows following the tutorial order.](../../images/tutorial-lcms-proc-wf.png "The full tutorial workflow")

All these modules are part of the [Wokflow4Metabolomics](http://workflow4metabolomics.org/) tool suit ({% cite Giacomoni2014 %}, {% cite Guitton2017 %}).
They are compatible with the whole data analysis solution maintained by the W4M team. 

> ### {% icon comment %} Workflow4Metabolomics public history
>
> This training material can be followed running it on any Galaxy instance holding the Galaxy modules needed.
> Nonetheless, if you happen to be a W4M user and do not want to run the hands-on yourself, please note that
> you can find the entire history in the 'published histories' section:
> [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****)
>
{: .comment}

## Dataset description

To illustrate the steps in this tutorial, a dataset has been built purposely.
It is composed of 3 files. They are text files of tables with tabulation as separator. 

The *dataMatrix* file is a table containing the intensities of measured variables (ions, in lines) for every samples (in column).
The first column is for ions' identifiers while the first line is for samples' identifiers.

The *variableMetadata* file is a table containing information about the ions. 
The first column is for ions's identifiers while the other columns gather information about m/z and retention time (rt). 

The *sampleMetadata* file is a table containing information about the samples. 
The first column is for samples's identifiers while the other columns gather analytical and biological information
such as the order of injection in the analytical sequence and the biological groups of interest for the supposed study. 

The simulated design is composed of 30 biological samples (tagged "sample" in the *sampleType* column),
completed with 8 Quality-control pooled samples (tagged "pool" in the *sampleType* column)
and 6 extraction solvent samples (tagged "blank" in the *sampleType* column). 
The samples have been supposedly injected in two distinct sequences (tagged "B1" and "B2" in the *batch* column),
the injection order being given in the *injectionOrder* column. 
Two sample characteristics are given:
- The *Group* column represents two groups "A" and "B", supposedly two biological groups (*e.g.* phenotypes, treatment groups...).
- the *Osmo* column represents a measurement of supposed osmolarity, imagining that the samples may be urine samples. 

## Data upload

To perform the different exercices of this tutorial, you need to create a new history and upload the dedicated dataset. 

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the 3 starting files in your history. Two possibilities: 
>    - Option 1: from a shared data library (ask your instructor)
>    - Option 2: from [Zenodo](https://zenodo.org/record/5179809) using the URLs given below:
>
>    ```
>   https://zenodo.org/record/5179809/files/Dataprocessing_dataMatrix.txt
>   https://zenodo.org/record/5179809/files/Dataprocessing_sampleMetadata.txt
>   https://zenodo.org/record/5179809/files/Dataprocessing_variableMetadata.txt 
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets for their names to be shorter: "dataMatrix", "sampleMetadata", "variableMetadata".
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Check that the datatype is "tabular". If not, you may change it. 
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> > ### {% icon tip %} Comment to W4M users
> >
> > If you happen to be a W4M user, please note that you can find at the following link a ready-to-start history:
> > [GTN_LCMSprocessing_start](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****).
> >
> > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets ****TOUPDATE****.
> {: .tip}
>
{: .hands_on}

# Removing "trash" signals

Data are often affected by various sources of unwanted variability. This can be found in your tables in different ways, 
two common ones being the presence of unwanted ions, and the effect of biological or analytical variables on intensity measures. 
This unwanted information can limit the effectiveness of statistical methods, leading sometimes to difficulties in revealing investigated effects. 
Thus, identifying such variability can help analysing your data at its full potential. 
Yet, getting rid of such information may not be a trivial problem. It requires different elements to be completed successfully. 

In this section we will adress the question of "trash" signal filtering. 
By this we mean ions that are present in the extracted peak list we have, but that do not correspond to relevant compounds to analyse.
This can be for example noise, or ions from compounds that are not present in original biological samples. 

To make it clearer, we will illustrate this with two examples: a "simple" one first and a more advanced one in a second time. 

## Filtering signals at given retention times

When using a chromatography column for MS analysis, you may want to exclude some time ranges where you know the information found there is not of interest.
For example, you may want to exclude the dead volume, a calibration zone at the begining or the end, or to exclude a column flush.

In this tutorial, let's suppose data are from some LC-QTOF analysis with a dead volume between 0 and 0.4 minutes and a column flush from 16 minutes.
Then we may want to exclude ions that may be found at theses specific retention time (rt) ranges. 
A quick check at the variableMetadata reveals that a retention time column is available ("rt") with values in minutes.
We can then use this column to filter the dataset. 

> ### {% icon hands_on %} Hands-on: Using **Generic_filter** to filter ions found at specific retention times
>
> 1. {% tool [Generic_Filter](toolshed.g2.bx.psu.edu/repos/melpetera/generic_filter/generic_filter/2020.01) %} with the following parameters:
>    - {% icon param-file %} *"Data Matrix file"*: the `dataMatrix` file
>    - {% icon param-file %} *"Sample metadata file"*: the `sampleMetadata` file
>    - {% icon param-file %} *"Variable metadata file"*: the `variableMetadata` file
>    - *"Deleting samples and/or variables according to Numerical values"*: `yes`
>        - In *"Identify the parameter to filter "*:
>            - *"1: Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `rt`
>                - *"Interval of values to remove"*: `extremity`
>                    - *"Remove all values lower than"*: `0.4`
>                    - *"And upper than"*: `16.0`
>    - *"Deleting samples and/or variables according to Qualitative values"*: `no`
>
>
>    > ### {% icon tip %} Comment to W4M users
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets number ****TOUPDATE****.
>    {: .tip}
>
{: .hands_on}

The **Generic Filter** {% icon tool %} tool generates 3 tables. They correspond to the 3 original tables, 
except the content has been filtered according to the specified parameters. 
By "filtering", it means removing from the dataset some variables (ions) and/or samples according to the defined filters. 

> ### {% icon question %} Question
>
> What have changed between the 3 input tables and the 3 output ones?
>
> > ### {% icon solution %} Solution
> >
> >  The sampleMetadata output is identical to the input one since no filter on samples has been applied.
> >  The variableMetadata output contains less lines than the input one: some ions has been removed from the table.
The same happened with the dataMatrix table, corresponding ions being removed from the file.
Thus, these two files went from 201 lines (header included) to 172 lines, meaning only 171 ions remained in the dataset.
> >
> {: .solution}
>
{: .question}

This was a relatively easy task to do. The key point to perform the filter was to know the ranges of rt that needed to be filter,
since the information about rt values was already explicitly found in the dataset. 
However, sometimes filtering requires more steps to be able to perform the wanted processing.
This will be illustrated in the next example.

## Using blanks to filter noise signals

As mentioned before, measured signal using mass spectrometry may not always be of interest.
It can be noise, or compounds not characteristic of the analysed biological samples.
There are several ways to reduce the impact on gathered data.
But one key point is always, as a starting point, to identify the issue. 

In the previous example we knew there were retention time values where signals were not relevant.
But they may also be some signals that represent noise, found at no specific retention times. 
So a question can be "how do we identify these signals?". 

One possible alternative is the use of blanks to estimate the noise, as a reference. 
The idea is to compare blanks’ intensities with other samples’ intensities.
If there is no subtantial difference, we can assume that the concerned signal is noise.

Of course, to be able to do so, you need to inject reference blanks along with your biological samples to get these intensities in your dataMatrix.
Thus, you need to anticipate it when you define your injection sequence. 
Ideally, the blanks to use are extraction blanks, but you can also use injection solvent depending on your protocols. 

When blanks are available in your dataset, the last thing to consider is "how do I formalise the information to be able to filter?".
One common way to compare may be to set a minimum difference between means or medians, or to test for significant difference with a statistical test.
In this tutorial, we will choose to calculate a mean fold change between blanks and non-blanks samples, and to set a threshold value for filtering.  

The mean fold change ("fold") for each ion can be calculated using the **Intensity Check** {% icon tool %} tool. 

> ### {% icon hands_on %} Hands-on: Using **Intensity Check** to generate the information needed to filter
>
> 1. {% tool [Intensity Check](toolshed.g2.bx.psu.edu/repos/melpetera/intensity_checks/intens_check/1.2.8) %} with the following parameters:
>    - {% icon param-file %} *"Data Matrix file"*: `GF_dataMatrix` (output of the previous **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file"*: `GF_sampleMetadata` (output of the previous **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file"*: `GF_variableMetadata` (output of the previous **Generic Filter** {% icon tool %} job)
>    - *"Computation method"*: `Between one class and all the remaining samples`
>        - *"Class column"*: `c5: sampleType`
>        - *"Selected class"*: `blank`
>        - *"Calculate the mean fold change"*: `Yes`
>            - *"Where should the class be placed for the mean fold change calculation?"*: `Denominator (Bottom)`
>
>
>    > ### {% icon tip %} Comment to W4M users
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets number ****TOUPDATE****.
>    {: .tip}
>
{: .hands_on}

This module generates two outputs: a pdf file for plots (that is not of interest in our example)
and a table corresponding to the variableMetadata file used as input, completed with new columns depending on the selected parameters. 
In our case, it generated a column named *fold_Other_VS_blank* that we will use for filtering. 

> ### {% icon question %} Question
>
> What does a value of "4" mean in the *fold_Other_VS_blank* column?
Remember that we used as parameter *"Selected class"*=`blank` and
*"Where should the class be placed for the mean fold change calculation?"*=`Denominator (Bottom)`.
>
> > ### {% icon solution %} Solution
> >
> > Since *"Selected class"*=`blank`, the samples are devided in two groups: the blank samples on one hand
and all the other samples on the other hand. 
To calculate a mean fold change (*i.e.* a ratio of means) between the two classes, we need to define which mean will be used as numerator 
and which one will be used as denominator. Since we defined that the selected class should be the denominator, the values we will get are 
(mean of non-blank samples)/(mean of blank samples). 
Thus, a value of "4" for a given ion means that the ion has a mean 4 times higher in non-blank samples compared to blank samples. 
> >
> {: .solution}
>
{: .question}

The **Intensity Check** {% icon tool %} tool generates additional information about your dataset, but do not perform any filter.
We can now define a threshold for filtering, and use the **Generic Filter** {% icon tool %} tool again to remove noise signal. 
Here we will use a threshold value of "4". What we want is to remove ions having a fold value lower than 4,
meaning that the mean values of biological samples' intensities for theses ions are not sufficiently high compared to blanks to be considered
to be resulting from relevant compounds. 

> ### {% icon hands_on %} Hands-on: Using **Generic_filter** to filter ions with insuffisant mean contrast with blank samples
>
> 1. {% tool [Generic_Filter](toolshed.g2.bx.psu.edu/repos/melpetera/generic_filter/generic_filter/2020.01) %} with the following parameters:
>    - {% icon param-file %} *"Data Matrix file"*: `GF_dataMatrix` (output of the previous **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file"*: `GF_sampleMetadata` (output of the previous **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file"*: `IC_variableMetadata` (output of the previous **Intensity Check** {% icon tool %} job)
>    - *"Deleting samples and/or variables according to Numerical values"*: `yes`
>        - In *"Identify the parameter to filter "*:
>            - *"1: Insert Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `fold_Other_VS_blank`
>                - *"Interval of values to remove"*: `lower`
>                    - *"Remove all values lower than"*: `4.0`
>    - *"Deleting samples and/or variables according to Qualitative values"*: `no`
>
>
>    > ### {% icon tip %} Comment to W4M users
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets number ****TOUPDATE****.
>    {: .tip}
>
{: .hands_on}

As for the previous use of **Generic Filter**, we now have a dataset filtered from noise signals. 

> ### {% icon question %} Questions
>
> What is the current stage of our trining dataset?
>
> > ### {% icon solution %} Solution
> >
> > From a dataset containing originally 200 ions, the two successive filters lead to a dataset containing only 131 ions.
The number of samples did not change, so we still have 30 biological samples, 8 QC pools and 6 blanks. 
> >
> {: .solution}
>
{: .question}

At that point we will not use the blank samples anymore. However, before removing them from the dataset, it is always interesting to 
get an overview of the dataset including blanks. 

## Filtered dataset overview

The idea here is to have a glance at what the dataset looks like at a macro scale.
For this tutorial, we consider that at this step what we have in our dataset is ions only resulting from compounds originally present in the biological samples.
Thus, there are some assumptions we can begin to make, that we can try to check using graphical tools. 

Here, we will use a tool that is called **Quality Metrics** to have an overview of our dataset through the generation of a pdf file containing some plots. 

> ### {% icon hands_on %} Hands-on: Using **Quality Metrics** to get an overview of the dataset
>
> 1. {% tool [Quality Metrics](toolshed.g2.bx.psu.edu/repos/ethevenot/qualitymetrics/quality_metrics/2.2.8) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix file"*: `GF_dataMatrix` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file"*: `GF_sampleMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file"*: `GF_variableMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - *"Coefficient of Variation"*: `no`
>    - *"Advanced parameters"*: `Use default`
>
>
>    > ### {% icon tip %} Comment to W4M users
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets number ****TOUPDATE****.
>    {: .tip}
>
{: .hands_on}

This tool generates several outputs (some of them are going to be used in a later use of this module), however for now only the pdf output is of interest.

![Screenshot of the PDF output from Quality Metrics. It is a composition of several distinct plots. In particular there are a PCA plot and a total intensity per sample plot.](../../images/lcmsproc_QM.png "PDF output from Quality Metrics")

On the top left of the picture, we can see the two first components of a Principal Component Analysis (PCA) with the projection of samples on it, colored by sample type.
We can see that the main variability in the dataset distinguishes the blank samples (on the left in blank) from the other ones (red for pools and green for samples).
This is awaited since blanks are supposed to have very low intensities compared to biological samples. 
This observation is consistant with the top middle plot which represents the sum of intensities for each sample (plotting according to the injection order),
where blank samples have very low values. 

This kind of plots is a good way to detect samples that may have abnormally low profiles: 
if a sample is positioned at the same area as the blanks, it is suspicious and a special attention should be given to the concerned sample.

In this tutorial dataset, no atypical sample is observed, so no special attention needs to be paid on specific samples.

Since the blank samples are of no use anymore, we can remove them from the dataset. 
Again, this can be done running the **Generic Filter** {% icon tool %} tool, using the *sampletype* column of the sampleMetadata table. 

> ### {% icon hands_on %} Hands-on: Using **Generic_filter** to remove blank samples from the dataset
>
> 1. {% tool [Generic_Filter](toolshed.g2.bx.psu.edu/repos/melpetera/generic_filter/generic_filter/2020.01) %} with the following parameters:
>    - {% icon param-file %} *"Data Matrix file"*: `GF_dataMatrix` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file"*: `GF_sampleMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file"*: `GF_variableMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - *"Deleting samples and/or variables according to Numerical values"*: `no`
>    - *"Deleting samples and/or variables according to Qualitative values"*: `yes`
>        - In *"Removing a level in factor"*:
>            - *"1: Insert Removing a level in factor"*
>                - *"Name of the column to filter"*: `sampleType`
>                - *"Remove factor when"*: `blank`
>
>    > ### {% icon comment %} Comment
>    >
>    > The tabular outputs of the previous **Quality Metrics** job are not needed here.
This is why what we use as input here are the outputs from the last **Generic Filter** job.
>    {: .comment}
>
>    > ### {% icon tip %} Comment to W4M users
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets number ****TOUPDATE****.
>    {: .tip}
>
{: .hands_on}

This step leads to a dataset containing 131 ions and 38 samples (from which 8 pools). 

You may have noticed that from module to module, output names tend to become longer and longer.
To prevent very long and not-so-informative names due to successive use of modules, we highly recommand to regularly rename the outputs.
Here, we completed a first row of data processing with the filter of trash signals. Now is a good time for a little renaming.

> ### {% icon hands_on %} Hands-on: Rename the last 3 tables
>
> 1. Rename the tables for their names to be shorter: "RT_blank_Filter_dataMatrix", "RT_blank_Filter_sampleMetadata", "RT_blank_Filter_variableMetadata".
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
>
{: .hands_on}

Now we are ready to investigate another type of data processing: correcting intensities

# Signal drift and batch effect correction

Here we will illustrate an example of data processing that is not about filtering, but about correcting. 
Indeed, when it comes to perform statistics comparing biological samples, it is crucial for the intensity values used to be reflecting
relevant variability between samples. 

In untargeted Metabolomics studies, we manipulate measures that are relative abundancies.
Although we have no unit attached to the intensities, we at least assume that, for a given ion, an intensity value for one sample being higher
than the one from another sample means that the compound from which the ion is generated is original found in higher abundance in the biological sample from
the first sample compared to the other. 
This assumption may seem trivial, but truth is it is not when dealing with LC-MS data.

It is known that when injecting successively a large number of samples, the LC-MS system tends to get dirty. This may cause a measure drift. 
To prevent inability to catch signal anymore, in case of large injection series, the sequence is generally divided into several batches and the source is cleaned between batches.
Unfortunately, these signal drift and batch design can add significant variability in the data. 
It makes sample comparison complicated, since the assumption stated previously about abundance in biological samples may not be true anymore. 
In case data is impacted by these effects, we need a way to normalise the data to obtain something reliable for statistical analysis. 

## [Optional step] Checking batch effect on data

You may have noticed in the **Quality Metrics** output PDF that appart from blanks, samples seemed to be seperated in two groups. 
Truth is it is indeed. However, the best would be to confirm it, in particular by finding out what these groups are linked to.
Here dices are already rolled, and given the fact that in this section we are adressing signal drift and batch effects, 
one could suppose that the groups may have something to do with it. 

Let's confirm it by performing a PCA, with the specificity to colour the sample projections according to the supposed effect: the batch information. 
This can be done using the **Multivariate** {% icon tool %} tool. 

> ### {% icon hands_on %} Hands-on: Using **Multivariate** to get a coloured score plot from a PCA
>
> 1. {% tool [Multivariate](toolshed.g2.bx.psu.edu/repos/ethevenot/multivariate/Multivariate/2.3.10) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix file"*: `RT_blank_Filter_dataMatrix` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file"*: `RT_blank_Filter_sampleMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file"*: `RT_blank_Filter_variableMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - *"Number of predictive components"*: `2`
>    - *"Advanced graphical parameters"*: `Full parameter list`
>        - *"Sample colors"*: `batch`
>        - *"Amount by which plotting text should be magnified relative to the default"*: `0.4`
>    - *"Advanced computational parameters"*: `Use default`
>
>    > ### {% icon comment %} Comment
>    >
>    > Changing the `0.8` value to `0.4` makes the text size of labels on the score plot smaller.
Since the plot box size by default is tiny, this enhances the readability of the plot. 
>    {: .comment}
>
>    > ### {% icon tip %} Comment to W4M users
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets number ****TOUPDATE****.
>    {: .tip}
>
{: .hands_on}

This output of **Multivariate** we will be interested in here is only the PDF file (Multivariate_figure.pdf). 
Among the four plots displayed in the file, we can see at the bottom left corner that the first component clearly reveals two distinct groups.
These groups perfectly match the batch information we coloured the sample projections with. 

This confirms that we indeed have at least a batch effect in the data. 
To note, this comment is relevant thanks to the fact that supposedly we randomised the injection sequence. 
This is essential to be able to efficiently separate any analytical effect from known biological variables when searching for effects. 
If you look at the sampleMetadata file we originally imported, you can notice that the biological group (*A* and *B*) are equally found in the two batches,
with alternation all through the injection order. 

## Performing the correction process

We previously saw that in data processing issues, it is sometimes crucial to anticipate what needs to be done before defining the injection sequence.
This is even more the case when dealing with batch effects and signal drifts.

There is a strategy, described initially by Van Der Kloet in 2009 ({% cite VdK2009 %}), that has made its way to nowadays procedures ({% cite Dunn2011 %}).
The idea is to model the signal drift and the batch effect level by using a reference that can be representative for every ions in the dataset.
Indeed, signal drift and batch effects can have very different impacts accross ions in a same dataset, so it is crucial to have a reference that is reliable
for each ion independantly. 

Here the "universal" reference is obtained by using samples that are made by pooling together an extract of every samples in study. 
Thus, these Quality-control pooled samples ('pools') contain all the compounds that are originally found in the samples, being a reference for a very wide range of ions.
By injecting these pools all through the injection sequence of the study samples, we obtain a reference for which the main variability observed is composed of
the analytical effects we want to correct, since biologically they are supposed to be identical.

> ### {% icon details %} More details about the theory
>
> The procedure that is used to correct data in this tutorial is the following.
>
> For each ion independently, the normalisation process works as described in the following picture:
>
> ![An example plot with the normalisation formula](../../images/lcms_BC_theo.png "How this works")
>
> You can see a plot representing 8 sample measures (blue points) in a batch for a given extracted ion. 
The yellow line represents a model for the signal drift that can be used to normalise the data.
This line is determined using the pools only (red squares). 
With the given formula, we can correct the signal drift. 
>
> This work has to be done for each batch. 
Thus, if your sequence is divided into several batches, the idea is to obtain something similar to the following picture:
>
> ![A before/after plot showing an example of intensities before correction, with clear signal drift and batch effects, and after correction where the effects have been erased thanks to the correction process*](../../images/lcms_BC_theo2.png "Before/after picture")
>
{: .details}

In this tutorial, we have all the information we need to perform the correction:
- Pools representing the molecule diversity have been injected regularly all through the sequences.
- They are numerous enough in each batch for the regression to be reliable (well, at least we have enough of them to perform a linear modeling).
- The dataset contains the mandatory information needed in the sampleMetadata file: the injection order, the batches of analysis and the sample type (pool or sample).

We can then use the **Batch correction** {% icon tool %} tool to perform the correction. 

> ### {% icon hands_on %} Hands-on: Using **Batch correction** to correct the data from signal drift and batch effet
>
> 1. {% tool [Batch_correction](toolshed.g2.bx.psu.edu/repos/melpetera/batchcorrection/Batch_correction/3.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Data Matrix file "*: `RT_blank_Filter_dataMatrix` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file "*: `RT_blank_Filter_sampleMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file "*: `RT_blank_Filter_variableMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - *"Type of regression model "*: `linear`
>        - *"Factor of interest "*: `Group`
>        - *"Level of details for plots "*: `complete`
>
>
>    > ### {% icon tip %} Comment to W4M users
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/****TOUPDATE****) history,
this step corresponds to the datasets number ****TOUPDATE****.
>    {: .tip}
>
{: .hands_on}

The tools generates 3 outputs: ...

Here I stopped for now
**TODO** to complete this section about BC of course, and to completed the ones after (not forgetting the "renaming hands:on" as the one previously)
**NOTTOFORGET** Links and numbers from the reference history; requesting the GTN category for the Zenodo repo




# Filtering signals of insufficient quality






> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Quality Metrics](toolshed.g2.bx.psu.edu/repos/ethevenot/qualitymetrics/quality_metrics/2.2.8) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix file"*: `BC_linear_dataMatrix` (output of the **Batch_correction** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file"*: `RT_blank_Filter_sampleMetadata` (output of the last **Generic Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file"*: `BC_linear_variableMetadata` (output of the **Batch_correction** {% icon tool %} job)
>    - *"Coefficient of Variation"*: `no`
>    - *"Advanced parameters"*: `Use default`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}





> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Generic_Filter](toolshed.g2.bx.psu.edu/repos/melpetera/generic_filter/generic_filter/2020.01) %} with the following parameters:
>    - {% icon param-file %} *"Data Matrix file"*: `dataMatrix_out` (output of **Batch_correction** {% icon tool %})
>    - {% icon param-file %} *"Sample metadata file"*: `sampleMetadata_out` (output of **Quality Metrics** {% icon tool %})
>    - {% icon param-file %} *"Variable metadata file"*: `variableMetadata_out` (output of **Quality Metrics** {% icon tool %})
>    - *"Deleting samples and/or variables according to Numerical values"*: `yes`
>        - In *"Identify the parameter to filter "*:
>            - {% icon param-repeat %} *"Insert Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `poolCV_over_sampleCV`
>                - *"Interval of values to remove"*: `upper`
>                    - *"Remove all values upper than"*: `1.0`
>            - {% icon param-repeat %} *"Insert Identify the parameter to filter "*
>                - *"On file"*: `Variable metadata`
>                - *"Name of the column to filter"*: `pool_CV`
>                - *"Interval of values to remove"*: `upper`
>                    - *"Remove all values upper than"*: `0.3`
>    - *"Deleting samples and/or variables according to Qualitative values"*: `yes`
>        - In *"Removing a level in factor"*:
>            - {% icon param-repeat %} *"Insert Removing a level in factor"*
>                - *"Name of the column to filter"*: `sampleType`
>                - *"Remove factor when"*: `pool`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}







# Study-specific data processing: example normalising data according to a non-analytical effect






## Sub-step with **Multivariate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Multivariate](toolshed.g2.bx.psu.edu/repos/ethevenot/multivariate/Multivariate/2.3.10) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix file"*: `dataMatrix_out` (output of **Generic_Filter** {% icon tool %})
>    - {% icon param-file %} *"Sample metadata file"*: `sampleMetadata_out` (output of **Generic_Filter** {% icon tool %})
>    - {% icon param-file %} *"Variable metadata file"*: `variableMetadata_out` (output of **Generic_Filter** {% icon tool %})
>    - *"Number of predictive components"*: `2`
>    - *"Advanced graphical parameters"*: `Full parameter list`
>        - *"Sample colors"*: `Osmo`
>        - *"Amount by which plotting text should be magnified relative to the default"*: `0.4`
>    - *"Advanced computational parameters"*: `Use default`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Normalization**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Normalization](toolshed.g2.bx.psu.edu/repos/marie-tremblay-metatoul/normalization/normalization/1.0.7) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix of preprocessed data"*: `dataMatrix_out` (output of **Generic_Filter** {% icon tool %})
>    - *"Normalization method"*: `Quantitative variable`
>        - {% icon param-file %} *"Sample metadata matrix"*: `sampleMetadata_out` (output of **Generic_Filter** {% icon tool %})
>        - *"Name of the column of the numerical variable for normalization (weight, osmolality, ...)"*: `Osmo`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Multivariate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Multivariate](toolshed.g2.bx.psu.edu/repos/ethevenot/multivariate/Multivariate/2.3.10) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix file"*: `dataMatrixOut` (output of **Normalization** {% icon tool %})
>    - {% icon param-file %} *"Sample metadata file"*: `sampleMetadata_out` (output of **Generic_Filter** {% icon tool %})
>    - {% icon param-file %} *"Variable metadata file"*: `variableMetadata_out` (output of **Generic_Filter** {% icon tool %})
>    - *"Number of predictive components"*: `2`
>    - *"Advanced graphical parameters"*: `Full parameter list`
>        - *"Sample colors"*: `Osmo`
>        - *"Amount by which plotting text should be magnified relative to the default"*: `0.4`
>    - *"Advanced computational parameters"*: `Use default`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Multivariate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Multivariate](toolshed.g2.bx.psu.edu/repos/ethevenot/multivariate/Multivariate/2.3.10) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix file"*: `dataMatrixOut` (output of **Normalization** {% icon tool %})
>    - {% icon param-file %} *"Sample metadata file"*: `sampleMetadata_out` (output of **Generic_Filter** {% icon tool %})
>    - {% icon param-file %} *"Variable metadata file"*: `variableMetadata_out` (output of **Generic_Filter** {% icon tool %})
>    - *"Number of predictive components"*: `2`
>    - *"Advanced graphical parameters"*: `Full parameter list`
>        - *"Sample colors"*: `Group`
>        - *"Amount by which plotting text should be magnified relative to the default"*: `0.4`
>    - *"Advanced computational parameters"*: `Use default`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}













# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.