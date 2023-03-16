---
layout: tutorial_hands_on
subtopic: deconvo
priority: 3

title: Creating the bulk RNA-seq dataset for deconvolution
zenodo_link: 'https://zenodo.org/record/7319173'
questions:
- Where can I find good quality RNA-seq datasets?
- How can I reformat and manipulate these downloads to create the right format for MuSiC?
objectives:
- You will retrieve raw data from the EMBL-EBI Expression Atlas.
- You will manipulate the metadata and matrix files.
- You will combine the metadata and matrix files into an ESet object for MuSiC deconvolution.
- You will create multiple ESet objects - both combined and separated out by disease phenotype for your bulk dataset.
time_estimation: 1H
key_points:
- The EMBL-EBI Expression Atlas contains high quality datasets.
- Metadata manipulation is key for generating the correctly formatted resource.
contributions:
  authorship:
    - nomadscientist
    - mtekman
  testing:
    - MarisaJL

tags:
  - single-cell
  - human
  - deconvolution
  - bulk
  - transcriptomics

requirements:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
      - bulk-music
      - bulk-music-2-preparescref

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - bulk-music-4-compare

---


# Introduction

After completing the [MuSiC deconvolution tutorial](https://xuranw.github.io/MuSiC/articles/MuSiC.html) ({% cite wang2019bulk %}), you are hopefully excited to apply this analysis to data of your choice. Annoyingly, getting data in the right format is often what prevents us from being able to successfully apply analyses. This tutorial is all about reformatting a raw bulk RNA-seq dataset pulled from a public resource (the EMBL-EBI Expression atlas ({% cite Moreno2021 %}).  Let's get started!

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Metadata Manipulation

Just as in our scRNA-dataset preparation tutorial, we will tackle the metadata first. We are roughly following the same concept as in the previous bulk deconvolution tutorial, by comparing human pancreas data across a disease variable (type II diabetes vs healthy), but using public datasets to do it.

## Find the data
We explored the [expression atlas](https://www.ebi.ac.uk/gxa/experiments), browsing experiments in order to find the bulk RNA-seq pancreas dataset ({% cite segerstolpe2016single %}). You can [explore this dataset here](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5060/Results) using their browser. These cells come from 7 healthy individuals and 4 individuals with Type II diabetes, so we will create reference Expression Set objects for the total as well as separating out by phenotype, as you may have reason to do this in your analysis (or you may not!). This dataset is from the same lab that we built our scRNA-seq reference from, so we should get quite accurate results given the same lab made both datasets!

> <hands-on-title> Data upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    {{ page.zenodo_link }}/files/E-MTAB-5060-experiment-design.tsv
>    ```
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets as needed
> 
> 4. Check that the datatype is tabular
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 5. Add to `experiment-design` the following tags `#metadata #bulk #ebi`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

As before, the metadata object annoyingly has a bunch of unnecessary columns. You can examine this with the {% icon galaxy-eye %} in the Galaxy history. Let's remove them!

![Columns in a table where some contain run info or Sample Characteristic[age] while others are empty.](../../images/bulk-music/bulk-metadata-annoying.png "Ridiculous metadata columns and labels")

{% snippet faqs/galaxy/tutorial_mode.md %}

> <hands-on-title> Remove unnecessary columns </hands-on-title>
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `output` (Input dataset)
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"List of Fields"*: `3 5 7 8 9 10 11 12 13 15 16 17 18`
>
>    > <comment-title></comment-title>
>    >
>    > Advanced cut works slightly differently in a workflow versus running the tool independently. Independently, there is a list and you can click through the list to note your columns, while in a workflow it appears as a text option and you put each column on a different line. The point is, each number above represents a column, so remove them!
>    {: .comment}
>
{: .hands_on}

Now let's take care of the excessively wordy header titles - and note that oftentimes various programmes struggle with titles or cells that have any spaces ` ` in them, so removing those now often saves hassle later.

> <comment-title></comment-title>
> You might also remember in the MuSiC tutorial that we can analyse numeric parameters in the metadata (in that case, hbac1c content). Reformatting to ensure numerical values in these columns (i.e. taking the ` years` out of an age cell) is helpful then too.
{: .comment}

> <hands-on-title> Fixing titles </hands-on-title>
>
> 1. {% tool [Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.2) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Advanced Cut** {% icon tool %})
>    - In *"Check"*:
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[age\]`
>            - *"Replacement"*: `Age`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `year`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[body mass index\]`
>            - *"Replacement"*: `BMI`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[disease\]`
>            - *"Replacement"*: `Disease`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[individual\]`
>            - *"Replacement"*: `Individual`
>        - {% icon param-repeat %} *"Insert Check"*
>            - *"Find Regex"*: `Sample Characteristic\[sex\]`
>            - *"Replacement"*: `Sex`
> 
> 2. Change the datatype to tabular 
>    
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
{: .hands_on}

Now examine {% icon galaxy-eye %} your resultant metadata file in the Galaxy history. Better, right?

![5 columns with numerical or string information on Run, Age, BMI, Disease and Sex](../../images/bulk-music/bulk-metadata-pretty.png "Look at the pretty metadata")

This is ready to go, so now we'll reformat the matrix!

# Manipulate the expression matrix

Let's upload the dataset.

> <hands-on-title> Data upload </hands-on-title>
>
> 1. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    {{ page.zenodo_link }}/files/E-MTAB-5060-raw-counts.tsv
>    ```
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. Rename the dataset as needed
> 3. Check that the datatype is tabular
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 4. Add to `raw-counts` the following tags `#matrix #bulk #ebi`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

Now examine {% icon galaxy-eye %} your raw counts file in the Galaxy history.

> <question-title></question-title>
>
> 1. Are samples in the rows or columns?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > ![Column 1 contains Gene ID followed by many lines of ENSG####. Column 2 contains the gene names. The following columns contain numerous iterations of ERR#####](../../images/bulk-music/raw-matrix.png "Gene info")
> >
> > 1. By examining the matrix, you can find that genes are the rows while samples are the columns.
> >
> {: .solution}
{: .question}

While it's awesome that there's a gene name column, unfortunately the gene names will be duplicated - different ENS IDs can refer to the same Gene Name. This going to be a problem later. So we need to get this in a format to collapse the ENS IDs, just as we did previously in the scRNA-seq data reference preparation. Sadly, we'll start by removing the column of gene names to prepare for the ENS ID collapse.

> <hands-on-title> Remove gene names column </hands-on-title>
>
> 1. {% tool [Remove columns](toolshed.g2.bx.psu.edu/repos/iuc/column_remove_by_header/column_remove_by_header/1.0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular file"*: `raw-counts` (Input dataset)
>    - In *"Select Columns"*:
>        - {% icon param-repeat %} *"Insert Select Columns"*
>            - *"Header name"*: `Gene Name`
>
{: .hands_on}

Now that your data is in a format of having a rows of ENS IDs and samples as columns, you can apply the handy ENS ID collapsing workflow as we did in the scRNA-seq reference. If you have already imported this workflow during the first tutorial, then you can use it again now. 

> <hands-on-title> Convert from Ensembl to GeneSymbol using workflow </hands-on-title>
>
> 1. Import this [workflow](https://usegalaxy.eu/u/wendi.bacon.training/w/convert-from-ensembl-to-genesymbol-summing-duplicate-genes).
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Run the workflow on your sample with the following parameters:
>    - *"Organism"*: `Human`
>    - {% icon param-file %} *"Expression Matrix (Gene Rows)"*: `output` (output of **Remove columns** {% icon tool %})
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
{: .hands_on}

The output will likely be called **Text transformation** and will look like this:

![Alphabetised gene symbols appear in column one with integers in the following columns corresponding to samples](../../images/bulk-music/ensidcollapse_bulk.png "Output of the ENS ID collapsing workflow for bulk dataset")

Success! You've now prepared your metadata and your matrix. It's time to put it together to create the Expression Set objects needed for MuSiC deconvolution.

# Construct Expression Set Objects

We have three more tasks to do: first, we need to create the expression set object with all the phenotypes combined. Then, we will create the two objects we actually need - one for healthy and one for diseased.

> <hands-on-title> Creating the combined object </hands-on-title>
>
> 1. {% tool [Construct Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_construct_eset/music_construct_eset/0.1.1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Assay Data"*: `out_file` #matrix (output of **Text transformation** {% icon tool %})
>    - {% icon param-file %} *"Phenotype Data"*: `out_file1` #metadata (output of **Regex Find And Replace** {% icon tool %})
>
> 2. Remove the `#metadata #matrix` tags from the output **RData ESet Object**
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many genes are in your object?
> 2. How many samples?
> 3. What metadata categories are there?
>
> > <div id="solution-1" class="box-title"><button type="button" aria-controls="solution-1-contents" aria-expanded="true" aria-label="Toggle solution box: "><i class="far fa-eye" aria-hidden="true"></i><span class="visually-hidden"></span> Solution<span role="button" class="fold-unfold fa fa-minus-square"></span></button></div>
> >
> > The trick with all of these questions is to examine {% icon galaxy-eye %} the `General info` output {% icon param-file %} of the **Construct Expression Set Object** tool.
> > ![Lines showing ExpressionSet'; assayData: 34997 features, 7 samples; protocolData: none; phenoData; sampleNames: ERR### (7 total); varLabels: Age BMI Disease Sex; varMetadata: labelDescription; and 3 more useless lines ](../../images/bulk-music/generalinfo.png "General info output")
> >
> > 1. There are `34997` features, which are the genes.
> > 2. There are `7` samples.
> > 3. The metadata categories are the same you prepared earlier, shown here in a category of phenoData: `Age BMI Disease Sex`
> >
> {: .solution}
{: .question}

> <hands-on-title> Creating the disease-only object </hands-on-title>
>
> 1. {% tool [Manipulate Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_manipulate_eset/music_manipulate_eset/0.1.1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Expression Set Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Concatenate other Expression Set objects?"*: `No`
>    - *"Subset the dataset?"*: `Yes`
>        - *"By"*: `Filter Samples and Genes by Phenotype Values`
>            - In *"Filter Samples by Condition"*:
>                - {% icon param-repeat %} *"Insert Filter Samples by Condition"*
>                    - *"Name of phenotype column"*: `Disease`
>                    - *"List of values in this column to filter for, comma-delimited"*: `type II diabetes mellitus`
>
> 2. Add the tag `#T2D` to the output **RData ESet Object**
>
{: .hands_on}

You can either re-run this tool or set it up again to create the healthy-only object.

> <hands-on-title> Creating the healthy-only object </hands-on-title>
>
> 1. {% tool [Manipulate Expression Set Object](toolshed.g2.bx.psu.edu/repos/bgruening/music_manipulate_eset/music_manipulate_eset/0.1.1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"Expression Set Dataset"*: `out_rds` (output of **Construct Expression Set Object** {% icon tool %})
>    - *"Concatenate other Expression Set objects?"*: `No`
>    - *"Subset the dataset?"*: `Yes`
>        - *"By"*: `Filter Samples and Genes by Phenotype Values`
>            - In *"Filter Samples by Condition"*:
>                - {% icon param-repeat %} *"Insert Filter Samples by Condition"*
>                    - *"Name of phenotype column"*: `Disease`
>                    - *"List of values in this column to filter for, comma-delimited"*: `normal`
>
> 2. Add the tag `#healthy` to the output **RData ESet Object**
>
{: .hands_on}

# Conclusion

{% icon congratulations %} Congrats! You have successfully reformatted the RNA-seq samples into two ESet objects consisting of disease-only or healthy-only samples. You're ready to take all this hard work and start comparing cell compositions in the next tutorial.

You can find the [workflow for generating the ESet object](https://usegalaxy.eu/u/wendi.bacon.training/w/music-deconvolution-data-generation--bulk--eset) and the [answer key history](https://usegalaxy.eu/u/wendi.bacon.training/h/music-deconvolution-data-generation--bulk--eset).

![7 boxes in the workflow editor and a subworkflow box for converting Ensembl to GeneSymbol](../../images/bulk-music/workflow-bulk.png "Workflow: Generating the bulk ESet Objects")

This tutorial is part of the [https://singlecell.usegalaxy.eu](https://singlecell.usegalaxy.eu) portal ({% cite tekman2020single %}).
