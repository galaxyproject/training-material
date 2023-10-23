---
layout: tutorial_hands_on

title: Converting NCBI data to the AnnData Format
subtopic: single-cell-CS-code
priority: 2
zenodo_link: 'https://zenodo.org/record/7053673'

questions:
- How do i understand NCBI data?
- How can i convert raw gene data to the AnnData format?
- How do i manually and automatically add metadata to my AnnData object?
objectives:
- add objectives here
time_estimation: 1H
key_points:
- Single cell data is huge, and must have its many (# genes) dimensions reduced for analysis
- Analysis is more subjective than we think, and biological understanding of the samples as well as many iterations of analysis are important to give us our best change of attaining real biological insights

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin
        - scrna-case_alevin-combine-datasets
tags:
- 10x
- paper-replication

contributions:
  authorship:
    - hexhowells

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_JUPYTER-trajectories
        - scrna-case_monocle3-trajectories

---

> <hands-on-title>Convert raw data to AnnData</hands-on-title>
>
> 1. {% tool [Import AnnData and loom](toolshed.g2.bx.psu.edu/repos/iuc/anndata_import/anndata_import/0.7.5+galaxy1) %} with the following parameters:
>    - *"hd5 format to be created"*: `Anndata file`
>         - *"Format for the annotated data matrix?"*: `Tabular, CSV, TSV`
>         - {% icon param-file %} *"Annotated data matrix"*: `Select all imported files`
>         - *"Does the first column store the row names?"*: `Yes`
>
{: .hands_on}

> <hands-on-title>Transpose AnnData objects</hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Select all AnnData files`
>    - *"Function to manipulate the object"*: `Transpose the data matrix, leaving observations and variables interchanged`
>
{: .hands_on}

> <hands-on-title>Combine AnnData objects</hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Select first Manipulate AnnData (transpose) output`
>    - *"Function to manipulate the object"*: `Concatenate along the observations axis`
>         - *"Annotated data matrix to add"*: `Select all other Manipulate AnnData (transpose) outputs`
>         - *"Join method"*: `Intersection of variables`
>         - *"Key to add the batch annotation to obs"*: `batch`
>         - *"Separator to join the existing index names with the batch category"*: `-`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Combined Object`
>
{: .hands_on}

> <hands-on-title></hands-on-title>
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined Object`
>    - *"What to inspect?"*: `Key-indexed observations annotation (obs)`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Combined Object`
>
{: .hands_on}

> <hands-on-title>Create replicate metadata</hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Observation data`
>    - *"1: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `2|4|8`
>         - *"Replace with"*: `poolA`
>    - **+ Insert Replacement**
>    - *"2: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `3|5|9`
>         - *"Replace with"*: `poolB`
>    - **+ Insert Replacement**
>    - *"3: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `6`
>         - *"Replace with"*: `poolC`
>    - **+ Insert Replacement**
>    - *"4: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `0|1|7`
>         - *"Replace with"*: `NA`
>    - **+ Insert Replacement**
>    - *"5: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `replicate`
>
> 2. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c2`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: `output of Replace Text`
>
> 3. **Rename** {% icon galaxy-pencil %} output `Replicate Metadata`
>
{: .hands_on}

> <hands-on-title>Create patient data</hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Observation data`
>    - *"1: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `(0$)|(1$)`
>         - *"Replace with"*: `patient1`
>    - **+ Insert Replacement**
>    - *"2: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `(2$)|(3$)|(4$)|(5$)|(6$)`
>         - *"Replace with"*: `patient2`
>    - **+ Insert Replacement**
>    - *"3: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `(7$)|(8$)|(9$)`
>         - *"Replace with"*: `patient3`
>    - **+ Insert Replacement**
>    - *"5: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `patient`
>
> 2. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c2`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: `output of Replace Text`
>
> 3. **Rename** {% icon galaxy-pencil %} output `Patient Metadata`
>
{: .hands_on}

> <hands-on-title>Create patient metadata</hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Observation data`
>    - *"1: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `0`
>         - *"Replace with"*: `AUG_PB1A`
>    - **+ Insert Replacement**
>    - *"2: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `1$`
>         - *"Replace with"*: `AUG_PB1B`
>    - **+ Insert Replacement**
>    - *"3: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `(2$)|3`
>         - *"Replace with"*: `MAY_PB1A`
>    - **+ Insert Replacement**
>    - *"4: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `4|5|6`
>         - *"Replace with"*: `MAY_PB1B`
>    - **+ Insert Replacement**
>    - *"5: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `7`
>         - *"Replace with"*: `MAY_PB2A`
>    - **+ Insert Replacement**
>    - *"6: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `8|9`
>         - *"Replace with"*: `MAY_PB2B`
>    - **+ Insert Replacement**
>    - *"7: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `specimenID`
>
> 2. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c2`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: `output of Replace Text`
>
> 3. **Rename** {% icon galaxy-pencil %} output `Specimen Metadata`
>
{: .hands_on}

> <hands-on-title>Create tumor metadata</hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Observation data`
>    - *"1: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `0`
>         - *"Replace with"*: `left-mid`
>    - **+ Insert Replacement**
>    - *"2: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `1|2|3|8|9`
>         - *"Replace with"*: `right-mid`
>    - **+ Insert Replacement**
>    - *"3: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `4|5|6`
>         - *"Replace with"*: `right-apex`
>    - **+ Insert Replacement**
>    - *"4: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `7`
>         - *"Replace with"*: `right-anterior`
>    - **+ Insert Replacement**
>    - *"5: Replacement"*
>         - *"in column"*: `Column: 2`
>         - *"Find pattern"*: `batch`
>         - *"Replace with"*: `tumorSpecimen`
>
> 2. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c2`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: `output of Replace Text`
>
> 3. **Rename** {% icon galaxy-pencil %} output `Tumor Metadata`
>
{: .hands_on}

> <hands-on-title>Combine metadata</hands-on-title>
>
> 1. {% tool [Paste](Paste1) %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: `Replicate Metadata`
>    - {% icon param-file %} *"and"*: `Patient Metadata`
>    - *"Delimit by"*: `Tab`
>
> 2. {% tool [Paste](Paste1) %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: `Output of previous Paste`
>    - {% icon param-file %} *"and"*: `Specimen Metadata`
>    - *"Delimit by"*: `Tab`
>
> 3. {% tool [Paste](Paste1) %} with the following parameters:
>    - {% icon param-file %} *"Paste"*: `Output of previous Paste`
>    - {% icon param-file %} *"and"*: `Tumor Metadata`
>    - *"Delimit by"*: `Tab`
>
> 4. **Rename** {% icon galaxy-pencil %} output `Cell Metadata`
>
{: .hands_on}

> <hands-on-title>Add metadata to AnnData object</hands-on-title>
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `Combined Object`
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations of variables`
>         - *"What to annotate?"*: `Observations (obs)`
>         - {% icon param-file %} *"Table with new annotations"*: `Cell Metadata`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Annotated Object`
>
{: .hands_on}

> <hands-on-title>Add initial metadata</hands-on-title>
>
> 1. {% tool [Scanpy FilterCells](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_filter_cells/scanpy_filter_cells/1.8.1+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `Annotated Object`
>    - *"Format of input object"*: `AnnData format hdf5`
>    - *"Format of output object"*: `AnnData format`
>    - *"Name of the column in `anndata.var` that contains gene name"*: `_index`
>    - **+ Insert Parameters to select cells to keep**
>    - *"1: Parameters to select cells to keep"*
>         - *"Name of parameter to filter on"*: `n_genes`
>         - *"Min value"*: `0.0`
>         - *"Max value"*: `1000000000.0`
>    - *"Force recalculation of QC vars"*: `No`
>
{: .hands_on}

> <hands-on-title>Add final metadata</hands-on-title>
>
> 1. {% tool [AnnData Operations](toolshed.g2.bx.psu.edu/repos/ebi-gxa/anndata_ops/anndata_ops/1.8.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in hdf5 AnnData format"*: `Output of Scanpy FilterCells`
>    - *"Format of output object"*: `AnnData format`
>    - *"Copy AnnData to .raw"*: `No`
>    - *"Gene symbols field in AnnData"*: `index`
>    - **+ Insert Flag genes that start with these names**
>    - *"1: Parameters to select cells to keep"*
>         - *"starts withn"*: `MT-`
>         - *"Var name"*: `mito`
>    - *"Number of top genes"*: `50`
>
> 2. **Rename** {% icon galaxy-pencil %} output `Final Object`
>
{: .hands_on}