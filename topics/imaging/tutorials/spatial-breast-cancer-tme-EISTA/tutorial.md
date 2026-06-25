---
layout: tutorial_hands_on

title: Spatial transcriptomics of the breast cancer tumour microenvironment using
  EISTA Galaxy
zenodo_link: https://zenodo.org/records/15129356
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

---


# Introduction

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/api/records/15129356/files/BreastCancer1.zip/content
>    ```
>    The data on Zenodo are provided as a compressed archive. When you paste the link
>    into Galaxy's upload tool, Galaxy will show an **Export** option. Select it and
>    click **Next** to import the archive.
>
>    The archive contains the following files:
>    - `V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5`
>    - `V1_Breast_Cancer_Block_A_Section_1_image.tif`
>    - `scalefactors_json.json`
>    - `tissue_hires_image.png`
>    - `tissue_lowres_image.png`
>    - `tissue_positions_list.csv`
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Once the data are imported, check the datatype of each dataset and set it according
>    to the file extension:
>    - `.csv` → `csv`
>    - `.tiff` → `tiff`
>    - `.h5` → `h5`
>    - `.png` → `png`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 4. Rename the datasets
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 6. Create a text file named `Resolution.txt` with the following contents and upload it to Galaxy:
>    ```
>    0.2
>    0.5
>    0.8
>    1.0
>    1.2
>    ```
>
>    Use the Galaxy upload tool or the paste/fetch data option to create the file, and
>    follow the standard upload steps shown in the Galaxy FAQ snippets.
>
>    {% snippet faqs/galaxy/datasets_upload_paste_fetched_data.md %}
>
> 7. Run {% tool [SpatialData IO](toolshed.g2.bx.psu.edu/repos/iuc/spatialdata_io/spatialdata_io/0.7.2+galaxy1) %} with the imported files to create the SpatialData object:
>    - *"Spatial technology"*: `10x Genomics Visium`
>    - {% icon param-file %} *"Feature BC matrix (Counts file)"*: `V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5`
>    - *"Is the matrix input raw?"*: `No`
>    - {% icon param-file %} *"Scale factors file"*: `scalefactors_json.json`
>    - {% icon param-file %} *"Full resolution image"*: `V1_Breast_Cancer_Block_A_Section_1_image.tif`
>    - {% icon param-file %} *"Tissue high resolution image"*: `tissue_hires_image.png`
>    - {% icon param-file %} *"Tissue low resolution image"*: `tissue_lowres_image.png`
>    - {% icon param-file %} *"Tissue positions file"*: `tissue_positions_list.csv`
>
>    This step creates the SpatialData object that will be used in the rest of the tutorial.
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **SpatialData Operations**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SpatialData Operations](toolshed.g2.bx.psu.edu/repos/iuc/spatialdata_operation/spatialdata_operation/0.7.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"SpatialData object"*: `output` (Input dataset)
>    - *"Operation"*: `Export the table of a SpatialData object to anndata`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Split file**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Split file](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_to_collection/split_file_to_collection/0.5.2) %} with the following parameters:
>    - *"Select the file type to split"*: `Text files`
>        - {% icon param-file %} *"Text file to split"*: `output` (Input dataset)
>        - *"Specify number of output files or number of records per file?"*: `Number of records per file ('chunk mode')`
>        - *"Method to allocate records to new files"*: `Maintain record order`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Replace Text**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output` (Input dataset)
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `^`
>            - *"Replace with:"*: `leiden_res_`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Map parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Map parameter value](toolshed.g2.bx.psu.edu/repos/iuc/map_param_value/map_param_value/0.2.0) %} with the following parameters:
>    - *"Select type of input parameter to match"*: `Integer`
>        - *"Value to map"*: `{'id': 12, 'output_name': 'output'}`
>    - *"Select how to handle unmapped values"*: `Use unmodified input parameter value if input parameter value not found in mappings`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Map parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Map parameter value](toolshed.g2.bx.psu.edu/repos/iuc/map_param_value/map_param_value/0.2.0) %} with the following parameters:
>    - *"Select type of input parameter to match"*: `Boolean`
>        - *"Value to map"*: `Yes`
>        - In *"Add value mapping"*:
>            - {% icon param-repeat %} *"Insert Add value mapping"*
>                - *"to this value"*: `True`
>    - *"Select how to handle unmapped values"*: `Use unmodified input parameter value if input parameter value not found in mappings`
>    - *"Select type of parameter to output"*: `Boolean`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `result_adata` (output of **SpatialData Operations** {% icon tool %})
>    - *"Method used for inspecting"*: `Calculate quality control metrics, using 'pp.calculate_qc_metrics'`
>        - *"Proportions of top genes to cover"*: `{'id': 1, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Parse parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Parse parameter value](param_value_from_file) %} with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: `list_output_txt` (output of **Split file** {% icon tool %})
>    - *"Select type of parameter to parse"*: `Float`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Transpose**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Transpose](toolshed.g2.bx.psu.edu/repos/iuc/datamash_transpose/datamash_transpose/1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `outfile` (output of **Replace Text** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `total_counts`
>            - *"y coordinate"*: `n_genes_by_counts`
>            - *"Color by"*: `volume`
>            - *"Use the layers attribute?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `n_genes_by_counts`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `total_counts`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `volume`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for filtering"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - *"Filter"*: `Minimum number of genes expressed`
>            - *"Minimum number of genes expressed required for a cell to pass filtering"*: `{'id': 2, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Compose text parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Compose text parameter value](toolshed.g2.bx.psu.edu/repos/iuc/compose_text_param/compose_text_param/0.1.1) %} with the following parameters:
>    - In *"components"*:
>        - {% icon param-repeat %} *"Insert components"*
>            - *"Choose the type of parameter for this field"*: `Text Parameter`
>                - *"Enter text that should be part of the computed value"*: `leiden_res_`
>        - {% icon param-repeat %} *"Insert components"*
>            - *"Choose the type of parameter for this field"*: `Float Parameter`
>                - *"Enter float that should be part of the computed value"*: `{'id': 24, 'output_name': 'float_param'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Replace Text**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file` (output of **Transpose** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `\t`
>            - *"Replace with:"*: `,`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - *"Filter"*: `Minimum number of counts`
>            - *"Minimum number of counts required for a cell to pass filtering"*: `{'id': 4, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Parse parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Parse parameter value](param_value_from_file) %} with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: `outfile` (output of **Replace Text** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter genes based on number of cells or counts, using 'pp.filter_genes'`
>        - *"Filter"*: `Minimum number of cells expressed`
>            - *"Minimum number of cells expressed required for a gene to pass filtering"*: `{'id': 5, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter genes based on number of cells or counts, using 'pp.filter_genes'`
>        - *"Filter"*: `Minimum number of counts`
>            - *"Minimum number of counts required for a gene to pass filtering"*: `{'id': 6, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - *"Filter"*: `Maximum number of counts`
>            - *"Maximum number of counts required for a cell to pass filtering"*: `{'id': 8, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter cell outliers based on counts and numbers of genes expressed, using 'pp.filter_cells'`
>        - *"Filter"*: `Maximum number of genes expressed`
>            - *"Maximum number of genes expressed required for a cell to pass filtering"*: `{'id': 10, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for inspecting"*: `Calculate quality control metrics, using 'pp.calculate_qc_metrics'`
>        - *"Proportions of top genes to cover"*: `{'id': 1, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `total_counts`
>            - *"y coordinate"*: `n_genes_by_counts`
>            - *"Color by"*: `volume`
>            - *"Use the layers attribute?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `n_genes_by_counts`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `total_counts`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `volume`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for filtering"*: `Filter on any column of observations or variables`
>        - *"What to filter?"*: `Observations (obs)`
>        - *"Type of filtering?"*: `By key (column) values`
>            - *"Key to filter"*: `volume`
>            - *"Type of value to filter"*: `Number`
>                - *"Filter"*: `greater than`
>                - *"Value"*: `{'id': 11, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for filtering"*: `Filter on any column of observations or variables`
>        - *"What to filter?"*: `Observations (obs)`
>        - *"Type of filtering?"*: `By key (column) values`
>            - *"Key to filter"*: `volume`
>            - *"Type of value to filter"*: `Number`
>                - *"Filter"*: `less than`
>                - *"Value"*: `{'id': 13, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>    - *"Method used for inspecting"*: `Calculate quality control metrics, using 'pp.calculate_qc_metrics'`
>        - *"Proportions of top genes to cover"*: `{'id': 1, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Scatter plot along observations or variables axes, using 'pl.scatter'`
>        - *"Plotting tool that computed coordinates"*: `Using coordinates`
>            - *"x coordinate"*: `total_counts`
>            - *"y coordinate"*: `n_genes_by_counts`
>            - *"Color by"*: `volume`
>            - *"Use the layers attribute?"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `n_genes_by_counts`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `total_counts`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Generic: Violin plot, using 'pl.violin'`
>        - *"Keys for accessing variables"*: `Subset of variables in 'adata.var_names' or fields of '.obs'`
>            - *"Keys for accessing variables"*: `volume`
>        - In *"Violin plot attributes"*:
>            - *"Add a stripplot on top of the violin plot"*: `Yes`
>                - *"Add a jitter to the stripplot"*: `Yes`
>                    - *"Size of the jitter points"*: `0.4`
>            - *"Display keys in multiple panels"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Manipulate AnnData**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.11.4+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Function to manipulate the object"*: `Copy layers from a different anndata object`
>        - {% icon param-file %} *"Source anndata object"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>        - In *"Layers to copy"*:
>            - {% icon param-repeat %} *"Insert Layers to copy"*
>                - *"Layer to be copied from the source anndata"*: `X`
>                - *"Target layer name"*: `counts`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy normalize**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy normalize](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_normalize/scanpy_normalize/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used for normalization"*: `Normalize counts per cell, using 'pp.normalize_total'`
>        - *"Target sum"*: `10000.0`
>        - *"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell"*: `No`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy normalize** {% icon tool %})
>    - *"Method used for inspecting"*: `Logarithmize the data matrix, using 'pp.log1p'`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for filtering"*: `Annotate (and filter) highly variable genes, using 'pp.highly_variable_genes'`
>        - *"Choose the flavor for identifying highly variable genes"*: `Cell Ranger`
>            - *"Number of highly-variable genes to keep"*: `4000`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for filtering"*: `Annotate (and filter) highly variable genes, using 'pp.highly_variable_genes'`
>        - *"Choose the flavor for identifying highly variable genes"*: `Seurat`
>            - *"Maximal mean cutoff"*: `1000000000.0`
>            - *"Maximal normalized dispersion cutoff"*: `50.0`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Pick parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Pick parameter value](toolshed.g2.bx.psu.edu/repos/iuc/pick_value/pick_value/0.2.0) %} with the following parameters:
>    - *"Picking behavior"*: `Pick first value`
>        - *"Select type of parameter to select from"*: `Dataset`
>            - In *"Pick from"*:
>                - {% icon param-repeat %} *"Insert Pick from"*
>                    - {% icon param-file %} *"Value"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy remove confounders**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy remove confounders](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_remove_confounders/scanpy_remove_confounders/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `data_param` (output of **Pick parameter value** {% icon tool %})
>    - *"Method used for plotting"*: `Regress out unwanted sources of variation, using 'pp.regress_out'`
>        - *"Keys for observation annotation on which to regress on"*: `total_counts,volume`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `data_param` (output of **Pick parameter value** {% icon tool %})
>    - *"Method used for plotting"*: `Preprocessing: Plot dispersions versus means for genes, using 'pl.highly_variable_genes'`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Pick parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Pick parameter value](toolshed.g2.bx.psu.edu/repos/iuc/pick_value/pick_value/0.2.0) %} with the following parameters:
>    - *"Picking behavior"*: `Pick first value`
>        - *"Select type of parameter to select from"*: `Dataset`
>            - In *"Pick from"*:
>                - {% icon param-repeat %} *"Insert Pick from"*
>                    - {% icon param-file %} *"Value"*: `anndata_out` (output of **Scanpy remove confounders** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `data_param` (output of **Pick parameter value** {% icon tool %})
>    - *"Method used for inspecting"*: `Scale data to unit variance and zero mean, using 'pp.scale'`
>        - *"Maximum value"*: `10.0`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Pick parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Pick parameter value](toolshed.g2.bx.psu.edu/repos/iuc/pick_value/pick_value/0.2.0) %} with the following parameters:
>    - *"Picking behavior"*: `Pick first value`
>        - *"Select type of parameter to select from"*: `Dataset`
>            - In *"Pick from"*:
>                - {% icon param-repeat %} *"Insert Pick from"*
>                    - {% icon param-file %} *"Value"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy cluster, embed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `data_param` (output of **Pick parameter value** {% icon tool %})
>    - *"Method used"*: `Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using 'pp.pca'`
>        - *"Type of PCA?"*: `Full PCA`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used for plotting"*: `PCA: Scatter plot in PCA coordinates, using 'pl.pca'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `log1p_total_counts,log1p_n_genes_by_counts,volume,total_counts`
>    - In *"Advanced Options"*:
>        - *"Output annotated data matrix?"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy plot** {% icon tool %})
>    - *"Method used for inspecting"*: `Compute a neighborhood graph of observations, using 'pp.neighbors'`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy cluster, embed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used"*: `Embed the neighborhood graph using UMAP, using 'tl.umap'`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `log1p_total_counts,log1p_n_genes_by_counts,volume`
>        - *"Show edges?"*: `No`
>    - In *"Advanced Options"*:
>        - *"Output annotated data matrix?"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy cluster, embed**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy cluster, embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy plot** {% icon tool %})
>    - *"Method used"*: `Cluster cells into subgroups, using 'tl.leiden'`
>        - *"Coarseness of the clusterin"*: `{'id': 24, 'output_name': 'float_param'}`
>        - *"Key under which to add the cluster labels"*: `{'id': 31, 'output_name': 'out1'}`
>        - *"How many iterations of the Leiden clustering algorithm to perform."*: `2`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Inspect AnnData**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.11.4+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy Inspect and manipulate**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy cluster, embed** {% icon tool %})
>    - *"Method used for inspecting"*: `Rank genes for characterizing groups, using 'tl.rank_genes_groups'`
>        - *"Get ranked genes as a Tabular file?"*: `True`
>        - *"The key of the observations grouping to consider"*: `{'id': 31, 'output_name': 'out1'}`
>        - *"Comparison"*: `Compare each group to the union of the rest of the group`
>        - *"Method"*: `Wilcoxon-Rank-Sum`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Text reformatting**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Text reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `obs` (output of **Inspect AnnData** {% icon tool %})
>    - *"AWK Program"*: `awk 'BEGIN {FS=\t} {print $1, $NF}'`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy Inspect and manipulate** {% icon tool %})
>    - *"Method used for plotting"*: `Marker genes: Plot ranking of genes using dotplot plot, using 'pl.rank_genes_groups'`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Column join**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Column join](toolshed.g2.bx.psu.edu/repos/iuc/collection_column_join/collection_column_join/0.0.3) %} with the following parameters:
>    - {% icon param-file %} *"Tabular files"*: `outfile` (output of **Text reformatting** {% icon tool %})
>    - *"Number of header lines in each input file"*: `1`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Advanced Cut**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Advanced Cut](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to cut"*: `tabular_output` (output of **Column join** {% icon tool %})
>    - *"Operation"*: `Discard`
>    - *"Cut by"*: `fields`
>        - *"Is there a header for the data's columns ?"*: `Yes`
>            - *"List of Fields"*: `c['1']`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Replace**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output` (output of **Advanced Cut** {% icon tool %})
>    - In *"Find and Replace"*:
>        - {% icon param-repeat %} *"Insert Find and Replace"*
>            - *"Find pattern"*: `split_file_\d+.txt_`
>            - *"Find-Pattern is a regular expression"*: `Yes`
>            - *"Replace all occurences of the pattern"*: `Yes`
>            - *"Find and Replace text in"*: `entire line`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Replace**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `outfile` (output of **Replace** {% icon tool %})
>    - In *"Find and Replace"*:
>        - {% icon param-repeat %} *"Insert Find and Replace"*
>            - *"Find pattern"*: `(?<!\.)\b(\d+)\b`
>            - *"Replace with"*: `c_\1`
>            - *"Find-Pattern is a regular expression"*: `Yes`
>            - *"Replace all occurences of the pattern"*: `Yes`
>            - *"Find and Replace text in"*: `entire line`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Table Compute**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy2) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `outfile` (output of **Replace** {% icon tool %})
>        - *"Input data has"*: ``
>        - *"Type of table operation"*: `Compute expression across rows or columns`
>            - *"Calculate"*: `Number of Unique Observations`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Manipulate AnnData**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.11.4+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy plot** {% icon tool %})
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>        - *"What to annotate?"*: `Observations (obs)`
>        - {% icon param-file %} *"Table with new annotations"*: `outfile` (output of **Replace** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Add column**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Add column](addValue) %} with the following parameters:
>    - *"Add this value"*: `{'id': 21, 'output_name': 'output_param_text'}`
>    - {% icon param-file %} *"to Dataset"*: `table` (output of **Table Compute** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Manipulate AnnData**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.11.4+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Function to manipulate the object"*: `Transform string annotations to categoricals`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Compute**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `out_file1` (output of **Add column** {% icon tool %})
>    - *"Input has a header line with column names?"*: `Yes`
>        - In *"Expressions"*:
>            - {% icon param-repeat %} *"Insert Expressions"*
>                - *"Add expression"*: `abs(c2-c3)`
>                - *"The new column name"*: `difference`
>    - In *"Error handling"*:
>        - *"If an expression cannot be computed for a row"*: `Fail the entire tool run`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy plot**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used for plotting"*: `Embeddings: Scatter plot in UMAP basis, using 'pl.umap'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `{'id': 34, 'output_name': 'text_param'}`
>        - *"Show edges?"*: `No`
>    - In *"Advanced Options"*:
>        - *"Output annotated data matrix?"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Sort**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Sort](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sort_header_tool/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Sort Query"*: `out_file1` (output of **Compute** {% icon tool %})
>    - *"Number of header lines"*: `1`
>    - In *"Column selections"*:
>        - {% icon param-repeat %} *"Insert Column selections"*
>            - *"Sort on column"*: `c4`
>        - {% icon param-repeat %} *"Insert Column selections"*
>            - *"Sort on column"*: `c1`
>            - *"in"*: `Descending order`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Scanpy filter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Scanpy filter](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_filter/scanpy_filter/1.11.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Scanpy plot** {% icon tool %})
>    - *"Method used for filtering"*: `Sample observations or variables with or without replacement, using 'pp.sample'`
>        - *"Type of sampling"*: `By fraction`
>            - *"Sample to this 'fraction' of the number of observations"*: `{'id': 15, 'output_name': 'output'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Table Compute**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy2) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `outfile` (output of **Sort** {% icon tool %})
>        - *"Type of table operation"*: `Drop, keep or duplicate rows and columns`
>            - *"List of rows to select"*: `2`
>    - *"Output formatting options"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Pick parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Pick parameter value](toolshed.g2.bx.psu.edu/repos/iuc/pick_value/pick_value/0.2.0) %} with the following parameters:
>    - *"Picking behavior"*: `Pick first value`
>        - *"Select type of parameter to select from"*: `Dataset`
>            - In *"Pick from"*:
>                - {% icon param-repeat %} *"Insert Pick from"*
>                    - {% icon param-file %} *"Value"*: `anndata_out` (output of **Scanpy filter** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Cut**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `table` (output of **Table Compute** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Analyze and visualize spatial multi-omics data**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Analyze and visualize spatial multi-omics data](toolshed.g2.bx.psu.edu/repos/goeckslab/squidpy/squidpy_spatial/1.5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input anndata"*: `data_param` (output of **Pick parameter value** {% icon tool %})
>    - *"Select an analysis"*: `Spatial neighbors -- Create a graph from spatial coordinates`
>        - In *"Advanced Graph Options"*:
>            - *"Type of coordinate system"*: `generic`
>            - *"Whether to compute the graph from Delaunay triangulation"*: `Yes`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Parse parameter value**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Parse parameter value](param_value_from_file) %} with the following parameters:
>    - {% icon param-file %} *"Input file containing parameter to parse out of"*: `out_file1` (output of **Cut** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Analyze and visualize spatial multi-omics data**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Analyze and visualize spatial multi-omics data](toolshed.g2.bx.psu.edu/repos/goeckslab/squidpy/squidpy_spatial/1.5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input anndata"*: `output` (output of **Analyze and visualize spatial multi-omics data** {% icon tool %})
>    - *"Select an analysis"*: `centrality_scores -- Compute centrality scores per cluster or cell type`
>        - *"Key in anndata.AnnData.obs where clustering is stored"*: `{'id': 88, 'output_name': 'text_param'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Analyze and visualize spatial multi-omics data**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Analyze and visualize spatial multi-omics data](toolshed.g2.bx.psu.edu/repos/goeckslab/squidpy/squidpy_spatial/1.5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input anndata"*: `output` (output of **Analyze and visualize spatial multi-omics data** {% icon tool %})
>    - *"Select an analysis"*: `centrality_scores -- Compute centrality scores per cluster or cell type`
>        - *"Key in anndata.AnnData.obs where clustering is stored"*: `{'id': 88, 'output_name': 'text_param'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Analyze and visualize spatial multi-omics data**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Analyze and visualize spatial multi-omics data](toolshed.g2.bx.psu.edu/repos/goeckslab/squidpy/squidpy_spatial/1.5.0+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select the input anndata"*: `output` (output of **Analyze and visualize spatial multi-omics data** {% icon tool %})
>    - *"Select an analysis"*: `spatial_autocorr -- Calculate Global Autocorrelation Statistic (Moran’s I or Geary’s C)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **CellTypist**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [CellTypist](toolshed.g2.bx.psu.edu/repos/iuc/celltypist/celltypist/1.7.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input AnnData file"*: `output` (output of **Analyze and visualize spatial multi-omics data** {% icon tool %})
>    - *"Select model from"*: `Cached`
>        - *"Choose CellTypist model"*: ``
>    - *"Refine the predicted labels by running the majority voting classifier after over-clustering"*: `Yes`
>    - *"Generate a dotplot of the predicted cell types"*: `Yes`
>        - *"Reference column in AnnData.obs for dotplot"*: `{'id': 88, 'output_name': 'text_param'}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Liana methods**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Liana methods](toolshed.g2.bx.psu.edu/repos/iuc/liana_methods/liana_methods/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **CellTypist** {% icon tool %})
>    - *"Method for ligand-receptor inference"*: `Aggregate ligand-receptor scores from multiple methods (rank_aggregate)`
>        - *"Group By"*: `{'id': 88, 'output_name': 'text_param'}`
>        - *"Interaction source"*: `Use built-in database`
>            - *"Resource source"*: `Download from LIANA API`
>        - *"Subset cell type pairs"*: `Use all possible combinations`
>        - *"Consensus options"*: `Default (Specificity and Magnitude)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **SpatialData Operations**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [SpatialData Operations](toolshed.g2.bx.psu.edu/repos/iuc/spatialdata_operation/spatialdata_operation/0.7.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"SpatialData object"*: `output` (Input dataset)
>    - *"Operation"*: `Import anndata table to a SpatialData object`
>        - {% icon param-file %} *"annotated data object to add"*: `anndata_out` (output of **Liana methods** {% icon tool %})
>        - *"Table name"*: `table_processed`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.