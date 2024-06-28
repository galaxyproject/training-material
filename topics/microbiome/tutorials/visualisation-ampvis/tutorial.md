---
layout: tutorial_hands_on


title: Divers and Adaptable Visualisations of Metabarcoding Data Using ampvis2
level: Intermediate
zenodo_link: https://zenodo.org/records/11281381
questions:
- How can we adapt the plots to our research data?
- How can we filter the data to show only significant information?
- How can we compare multiple visualisation methods?
- How can we use numerical and categorical metadata for amplicon visualization?
objectives:
- Use heatmap workflow to analyse and visualise amplicon data
- Use ungrouped or grouped data or grouped data with facets
- Use ordination plot, or boxplot, or rarefaction curve, or timeseries
time_estimation: 2H
key_points:
- Using various visualisation methods can present data from different perspectives
- With sufficient metadata, the data can be visualised in relation to different groups and facets 
contributors:
- lenaarenot
- paulzierep

---




Microbiome analysis using amplicon sequencing is central to many ecological studies.
The produced amplicon sequencing data are converted into OTU tables and represent the input 
for the ampvis2 tool, where they can be visualised in various ways {% cite Andersen2018 %}.
If you already have amplicon data produced and ready to visualise it, you can start with this tutorial. 
You can use your own data or download the data we used and follow us step-by-step though the tutorial. 

First of all you can put your data into a 
rarefaction curve to explore species richness. Then you can input your
data into subsets and finally create a heatmap, a boxplot, an ordination plot
or even a timeseries plot. Most of these visualisation methods are described in 
[Introduction to ampvis2](https://kasperskytte.github.io/ampvis2/articles/ampvis2.html#heatmap).
![overview of visualisation methods](./images/overview.png 
"Overview of posible visualisation methods (taken from: Introduction to ampvis2 by Kasper Skytte Andersen)")
Your data need to be in an acceptable format for the ampvis_load tool. The tool 
requires an OTU table and accepts the following formats: _phyloseq_, _biom_, 
_dada2_sequencetable_ or _tabular_. The OTU table is the only mandatory input for 
ampvis_load, but you can also input _sample_metadata_ (in _tabular_ or _tsv_ formats), 
_taxonomy_table_ (in _tabular_ format), _fasta_file_ (in _fasta_ format)
and _phylogenetic_tree_ (in _newick_ format), as well as various combinations thereof.

For this tutorial we chose to demonstrate all visualisation tools using a combination 
of 3 inputs:
OTU table, sample metadata and taxonomy table, all in _tabular_ format.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Rarefaction Curve

As first exploration of your data, you can start with rarefaction curve.
Rarefaction curves are a refined version of accumulation curves. They help visualize 
the number of species (species richness) in your samples by showing the average number 
of species observed as more samples are added. This method pools all the samples 
together and calculates the mean species richness for different sample sizes, resulting 
in a smooth curve. Rarefaction curves are particularly useful for comparing different datasets, 
as they provide a standard way to assess species richness regardless of sample size differences {% cite Gotelli2001 %}.
> <comment-title> </comment-title>
> - for this part, we need 'raw' data, it should not be normalised
> - for this section, we used a different dataset than for the rest of the tutorial
{: .comment}

## Get Data
If you don't use your own dataset, you can also go to **Zenodo** and find a dataset to download.
We looked for a dataset marked "open" and used the following:
[Environmental DNA metabarcoding data from Marina di Camerota coast (Italy) based on citizen science sampling](https://zenodo.org/records/10362755).
There are 2 datasets: **"V4-18S"** and **"COI"**, of which we used the one named **"COI"**.

### Sub-step: Generate **Uploadable Datasets from Downloaded Excel Sheet**
All data (OTU, metadata and tax table) we need to upload to Galaxy separately are combined in one excel sheet and looks like this...
![data to separate](./images/all_data.png 
"OTU, metadata and tax table is combined in one sheet and needs to be separated")

> <hands-on-title> Generate separated datasets for Galaxy upload </hands-on-title>
>
> 1. for OTU table:
>    - keep the sample names in the first row
>    - keep the asv+number in the first column
>    - the first cell (A1) needs to reed ASV or OTU 
>
>    > <comment-title> Attention! </comment-title>
>    >
>    > make sure sample names have no blank spaces
>    {: .comment}
>    
>    > <details-title> How it will look like </details-title>
>    >
>    > After separeting OTUs from the main data sheet it will look like this...
>    > 
>    >![separated OTU table](./images/otu.png "Separated OTU table")
>    {: .details}
>
> 2. for metadata table:
>    - copy the metadata (marked in blue) and the sample names to a new sheet
>    
>    > <details-title> How it will look like </details-title>
>    >
>    > The copied set of metadata will look like this...
>    > 
>    >![copied metadata](./images/meta_row.png "Copied set of metadata")
>    {: .details}
>
>    - transpose the dataset and copy to a new sheet 
>    - remove blank spaces from sample names 
>    
>    > <details-title> How it will look like </details-title>
>    >
>    > The transposed set of metadata will look like this...
>    > 
>    >![transposed metadata](./images/meta.png "Transposed metadata")
>    {: .details}
>
> 3. for tax table:
>    - keep the asv+number in the first column
>    - keep the last column "lineage"
>    - split the "lineage"-column by the delimeter __semicolon__
>    - give all columns a name 
>    
>    > <details-title> How it will look like </details-title>
>    >
>    > After separeting the taxa into diferent columns and renaming, it will look like this...
>    > 
>    >![separated and renamed tax table](./images/taxa.png "Separated and renamed tax table")
>    {: .details}
>
> 4. save all 3 data sheets separately
>
> 5. convert to tsv format
>    - e.g. use any free available tool online to convert xlsx to tsv
>    
{: .hands_on}


> <question-title></question-title>
> 
> 1. Why do you need to ensure that sample names have no blank spaces in the OTU table?
> 2. Why don't you use the metadata in its original untransposed form?
>
> > <solution-title></solution-title>
> >
> > 1. Many bioinformatics tools assume that sample names are continuous strings, and spaces can be 
interpreted as delimiters or end-of-string characters, leading to incorrect data parsing or analysis failures.
> > 2. You need the sample names as a column. Bioinformatics tools expect the metadata to be organised such that 
each metadata attribute has its own column.
> >
> {: .solution}
> 
{: .question}

## create a rarefaction curve
You can find the workflow "ampvis2 rarefaction v1.0 " on Galaxy and use it for the tutorial.

We have pre-selected the "step size", "colour curves by" and set __"free__ __scale"__ for "scales of the facets".

> <comment-title></comment-title>
> - first you need to upload the freshly generated dataset to Galaxy
> - the steps from heatmap hands-on box (next section) might be helpful
>
> > <details-title> forgot how to upload data? </details-title>
> >
> > If you forgot how to upload data to Galaxy, here is a nice tutorial:
> >
> > [Data Manipulation Olympics](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/data-manipulation-olympics/tutorial.html#upload-data)
> >
> {: .details}
>
{: .comment}

> <details-title> How it will look like </details-title>
>
> Result of the rarefaction curve.
> 
>![Result of the rarefaction curve](./images/rarefaction.png "Result of the rarefaction curve")
{: .details}

> <question-title></question-title>
> 
> 1. If you run the workflow "ampvis2 rarefaction v1.1 (with subset)" on Galaxy and use the following metadata for the subset: 
	metadata variable = sample_id and metadata values = (select all possible samples).
	What is the difference from the rarefaction curve you generated with the workflow without subsets? 
> 2. If you consider the output of the rarefaction curve before, one sample has a high curve and the rest are close to each other.
	Can you make the rest more "visible"?
>
> > <solution-title></solution-title>
> >
> > 1. There is no difference;, it's the same output 
> > 2. Yes, if you run the workflow (with subsets) and select all samples except of __"COI-B2b"__
> >
> > > <details-title> How it will look like </details-title>
> > > 
> > >  Result of this rarefaction curve.
> > > 
> > >![Result of this rarefaction curve](./images/rarefaction_without.png "Result of the rarefaction curve without "COI-B2b" ")
> > > 
> > {: .details}
> > 
> {: .solution}
>
{: .question}

# Use Case 1: Heatmap, Ordination Plot or Boxploot

To create a heatmap, ordination plot, or boxplot you can continue with your dataset or use the same as we do for the next sections.

> <comment-title> </comment-title>
> - we now use normalised data and a different dataset than for the rarefaction curve (as it has more metadata)
{: .comment}

## Get Data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/records/11281381/files/Galaxy11-[MiDAS_otushort_table.tsv].mothur.axes/content
>    https://zenodo.org/api/records/11281381/files/Galaxy1-[MiDAS_metadata.tsv].tabular/content
>    https://zenodo.org/api/records/11281381/files/Galaxy3-[MiDAS_taxtable.tsv].tabular/content
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype is in the right format
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 6. When you have your history created and ready on the right side, click on "workflow" above.
>	 Choose the needed workflow by clicking on the blue play button in its box(if you hover over it, 
>	 it says: Run workflow).
>
> 7. The workflow will ask you to input mandatory parameters. After doing so, click the blue button above
>	 "Run Workflow".
>    See on the next picture how for how it looks.
>
{: .hands_on}

![Running the workflow](./images/choose_parameters.png "Running the workflow, choose the right datasets and mandatory parameters")

## Heatmaps
Now, we can use our data, put them into subsets, and create ungrouped or grouped outputs, including those with facets. 
The subsets are based on variables we define and are available in the metadata {% cite Andersen2018 %}.
> <comment-title> </comment-title>
> - in the following sections, we provide prepared workflows on Galaxy and a set of parameters for running the indicated workflow
> - some parameters are pre-selected for you, such as the taxonomic level to aggregate the OTUs
{: .comment}

### Heatmap (ungrouped)
You can find the workflow "ampvis2 heatmap v3.0 (no group)" on Galaxy and use it for the tutorial.

We used the following metadata for this subset: metadata variable = Plant and metadata values = Aalborg East & Aalborg West.

Choose the metadata variable as Plant and metadata values as Aalborg East & Aalborg West.
In the next box you can see the resulting heatmap.

> <details-title> How it will look like </details-title>
>
> Result of the heatmap created with ungrouped data.
> 
>![Result of the heatmap](./images/heatmap_no_group.png "Result of the heatmap")
{: .details}

> <details-title> Error with metadata values while running the workflow </details-title>
>
> If a set in history shows red, indicating an error, click on the set to expand it.
> If it says: "parameter 'vals': an invalid option", then click on the "Run
> job again" button. The 'values' box will be highlighted in blue; simply select the values again from dropdown
> menue and click "Run Tool" button above.
>
> Additionally, note that the heatmap generation is paused. Once all sets in the history except the heapmap turns
> green, expand the heapmap set and click "Run job again". You may need to choose the metadata
> list again, make sure you choose the recently generated one. Click "Run Tool" again.
> 
{: .details}

### Heatmap (grouped)
You can find the workflow "ampvis2 heatmap v2.0 (only group)" on Galaxy and use it for the tutorial.

We used 2 different metadata subsets:
- 1) Metadata used for this subset: metadata variable = Plant, metadata values = Aalborg East & Aalborg West, grouped by = Plant
- 2) Metadata used for this subset: metadata variable = Period, metadata values = Winter & Summer, grouped by = Year

> <hands-on-title> Run a workflow </hands-on-title>
>
> 1. Create a new history (you can use the previous, but if you run multiple workflows, it might
>  	 become difficult to find your heatmap later)
>  
> 2. Use the same data and rename if you wish following the hands-on sections above
>
> 3. Choose the needed workflow 
>
> 4. Choose the mandarory parameters and click the "Run Workflow" button.
>    See on the next pictures how it looks like.
>
{: .hands_on}

> <details-title> How it will look like </details-title>
>
> Result of the first metadata subset heatmap created with grouped by Plant data.
> 
>![Result of the heatmap](./images/heatmap_gr_by_plant.png "Result of the heatmap created with grouped by _Plant_ ")
> 
> Result of the second metadata subset heatmap created with grouped by Year data.
> 
>![Result of the heatmap](./images/heatmap_gr_by_year.png "Result of the heatmap created with grouped by _Year_")
{: .details}

> <question-title></question-title>
>
> 1. Can you create a heatmap that shows only the first and the last year of data collection?
> 2. Can you create a heatmap using the following settings: metadata variable = Year and 
	metadata value = Date plus grouped by = Year?
>
> > <solution-title></solution-title>
> >
> > 1. Yes, with the following settings: metadata variable = Year, 
   metadata values = 2006 & 2015, grouped by = Year
> > 2. No, the metadata values must correspond to the metadata variable options
> >
> {: .solution}
>
{: .question}

### Heatmap (grouped with facets)
You can find the workflow "ampvis2 heatmap v1.0 (group+facet)" on Galaxy and use it for the tutorial.

We used 2 different metadata subsets:
- 1) Metadata used for this subset: metadata variable = Plant, metadata values = Aalborg East & Aalborg West, 
	grouped by = Plant, facet by = Period 
- 2) Metadata used for this subset: metadata variable = Period, metadata values = Winter & Summer, 
   grouped by = Year, facet by = Period 
   

> <hands-on-title> Run a workflow </hands-on-title>
>
> 1. Create a new history (if you wish)
>  
> 2. Use the same data and rename if you wish following the hands-on sections above
>
> 3. Choose the needed workflow 
>
> 4. Choose the mandarory parameters and click the "Run Workflow" button.
>    See on the next pictures how it looks like.
>
{: .hands_on}

> <details-title> How it will look like </details-title>
>
> Result of the first metadata subset heatmap created with grouped by Plant data and facet by Period.
> 
>![Result of the heatmap](./images/heatmap_plant_period.png "Result of the heatmap created with grouped by _Plant_ and facet by _Period_ ")
> 
> Result of the second metadata subset heatmap created with grouped by Year data.
> 
>![Result of the heatmap](./images/heatmap_year_period.png "Result of the heatmap created with grouped by _Year_ and facet by _Period_ ")
> 
{: .details}

## Ordination Plots
We can now use our data, generate subsets, and create different plots by applying various ordination methods. 
As with heatmaps, the subsets are based on variables we define and are available in the metadata {% cite Andersen2018 %}.

> <comment-title> </comment-title>
> - in the following sections, we provide prepared workflows on Galaxy along with the set of parameters to select
for running each workflow
> - some parameters are pre-selected for you, such as ordination method, 
transformation (if used) and options to colour and label the points or frames
{: .comment}

### Ordination Method: PCA
You can find the workflow "ampvis2 ordination plot v1.0 (pca)" on Galaxy and use it for the tutorial.

Metadata used for this subset: metadata variable = Plant and metadata values = Aalborg East & Aalborg West.

> <comment-title></comment-title>
> - use the same data set as for heatmaps
> - the steps from heatmap hands-on boxes remain the same
{: .comment}

> <details-title> How it will look like </details-title>
>
> Result of the ordination plot created with the PCA method.
> 
>![Result of the ordiantion plot](./images/ordination_pca.png "Result of the ordination plot created with the PCA method")
{: .details}

### Ordination Method: PCA plus Trajectory: _date_
You can find the workflow "ampvis2 ordination plot v1.1 (pca+trajectory_date)" on Galaxy and use it for the tutorial.

Metadata used for this subset: metadata variable = Plant and metadata values = Aalborg East & Aalborg West.

> <comment-title></comment-title>
> - use the same data set as for heatmaps
> - the steps from heatmap hands-on boxes remain the same
{: .comment}

> <details-title> How it will look like </details-title>
>
> Result of the ordination plot created with the PCA method plus using the trajectory date.
> 
>![Result of the ordiantion plot](./images/ordination_pca_date.png "Result of the ordination plot created with the PCA method plus using the trajectory date")
{: .details}

> <question-title></question-title>
>
>  Does it make sense to run the following settings: metadata variable = Period and metadata values = Winter & Summer?
>
> > <solution-title></solution-title>
> >
> >  No, you will get a very messy bundle of colours.
> >
> {: .solution}
>
{: .question}

### Ordination Method: CCA
You can find the workflow "ampvis2 ordination plot v1.2 (cca transform_hellinger)" on Galaxy and use it for the tutorial.

Metadata used for this subset: metadata variable = Plant and metadata values = Aalborg East & Aalborg West.

> <comment-title></comment-title>
> - use the same data set as for heatmaps
> - the steps from heatmap hands-on boxes remain the same
{: .comment}

> <details-title> How it will look like </details-title>
>
> Result of the ordination plot created with the CCA method and the Hellinger transformation.
> 
>![Result of the ordiantion plot](./images/ordination_cca_hellinger.png "Result of the ordination plot created with the CCA method and the Hellinger transformation")
{: .details}

> <question-title></question-title>
>
>  If you use the CCA ordination method with the following settings: metadata variable = Period and metadata values = Winter & Summer.
	What do you need to remove from pre-selected parameters to keep the ordination plot stays readable?
>
> > <solution-title></solution-title>
> >
> >  When you expand the ordination plot set in your history, you see options for colour, shape, frame, and label by. Select colour by _Period_ 
and frame by _Period_ and deselect the other mentioned options above, so they read _"Nothing_ _selected"_ .
> >
> {: .solution}
>
{: .question}

## Boxplot
We can now use our data, put them into subsets and create a boxplot. 
As with heatmaps, the subsets are based on variables we define and are available in the metadata {% cite Andersen2018 %}.
> <comment-title> </comment-title>
> - in the prepared workflow on Galaxy provided in this tutorial, some parameters are pre-selected for you, 
such as the number of taxa to show
> - The samples are grouped by _Period_
{: .comment}

Metadata used for this subset: metadata variable = Plant and metadata values = Aalborg East & Aalborg West.

> <comment-title></comment-title>
> - use the same data set as for heatmaps
> - the steps from heatmap hands-on boxes remain the same
{: .comment}

> <details-title> How it will look like </details-title>
>
> Result of the boxplot grouped by Period.
> 
>![Result of the boxplot](./images/boxplot_period.png "Result of the boxplot grouped by Period")
{: .details}

> <tip-title>Create a different boxplot</tip-title>
>
> * use the same data set
> * set metadata variable = Period and metadata values = Summer & Winter
>
> > <details-title> How it will look like </details-title>
> >
> > Result of this boxplot.
> >
> > ![Result of the boxplot](./images/boxplot_other.png "Result of this boxplot")
> >
> {: .details}
>
{: .tip}

> <question-title></question-title>
>
> 1. Can you create an output where only odd years are considered?
> 2. Do you need to change any pre-selected parameter for question 1?
>
> > <solution-title></solution-title>
> >
> > 1. Yes, if you set metadata variable = Year and metadata values = 2007, 2009, 2011, 2013, 2015 
> > 2. Yes, set "group the sample" to _Year_
> >
> {: .solution}
>
{: .question}

# Use Case 2: Time Series Plot

Time series analysis is primarily known for forecasting. A time series can be seen as 
an example of a random or stochastic process, which we can use to visualise seasonal 
differences {% cite DeGooijer2006 %}. 

In our dataset, and with the settings listed below, we can observe the 
temporal evolution of the 3 most common microorganisms in the plants Aalborg East and Aalborg West 
over the entire period data was collected.

## Create a Time Series Plot
You can find the workflow "ampvis2 timeseries v1.0" on Galaxy and use it for the tutorial.

Metadata used for this subset: metadata variable = Plant and metadata values = Aalborg East & Aalborg West.

The time variable is mandatory here, and as **_Date_** is the only valid time variable we can select, it has already been pre-selected. 
Since displaying a large number of taxa can make the plot messy (for this data set at least), we have chosen to show only the top 3 taxa.

> <comment-title></comment-title>
> - use the same data set as for heatmaps
> - the steps from heatmap hands-on boxes remain the same
{: .comment}

> <details-title> How it will look like </details-title>
>
> Result of the time series plot.
> 
>![Result of the time series plot](./images/timeseries.png "Result of the time series plot")
{: .details}

> <question-title></question-title>
>
> 1. If you run the following settings: metadata variable = Period and metadata values = Winter & Summer. 
	Do you need to change the time variable in the pre-selected parameters?
> 2. Can you adjust the settings to separate the 3 curves into distinct curves for each period, corresponding to the shown Phylum?
>
> > <solution-title></solution-title>
> >
> > 1. No, _Date_ is still the only possible option 
> > 2. Yes, you can separate the curves by expanding the time series set in the history and running it again with "group the sample by" _Period_
> >
> {: .solution}
>
{: .question}


# Conclusion
This tutorial has equipped you with essential skills for visualising and interpreting microbiome data using ampvis2. 
By exploring diverse visualisation methods and leveraging metadata effectively, you've learned how to gain valuable insights 
from complex datasets. Don't hesitate to explore further in Galaxy's toolbox for additional data manipulation tools and 
continue refining your analysis techniques to enhance your research capabilities. 

Happy exploring!