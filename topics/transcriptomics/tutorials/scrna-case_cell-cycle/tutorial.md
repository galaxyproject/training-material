---
layout: tutorial_hands_on

title: Removing the effects of cell cycle genes
zenodo_link: https://zenodo.org/record/7311628/
questions:
- How can I reduce the effects of the cell cycle on my scRNA-seq data?
objectives:
- Identify and regress out the effects of cell cycle genes
- Create PCA plots to understand the impact of the regression
requirements:
- 
time_estimation: 1H
key_points:
- Cell cycle genes can conceal what is happening in your data 
- Identifying the cell cycle genes and regressing out their effects can reveal underlying patterns in the data
contributors:
- MarisaJL

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

Single-cell RNA sequencing can be sensitive to both biological and technical
variation, which is why preparing your data carefully is an important part of
the analysis. You want the results to reflect real differences in expression
between cells that relates to their type or state. Other sources of variation
can conceal or confound this variation, making it harder for you to see what
is really going on. 

One common biological confounder is the cell cycle. Cells express different
genes during different parts of the cell cycle, depending on whether they are
in their growing phase or duplicating their DNA in preparation for cell
division. If these cell cycle genes are having a big impact on your data,
then you could end up with separate clusters that actually represent cells of
the same type that are just at different stages of the cycle. 

In this tutorial, we will identify the genes whose expression varies during
the cell cycle so that we can regress out (or remove) their effects on the
data. 

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

> ### Agenda
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

The data used in this tutorial is from a mouse dataset of fetal growth restriction {% raw %} `{% cite Bacon2018 %}`{% endraw %}. 

If you've been working through the Filter, Plot and Explore Single-cell 
RNA-seq Data tutorial then you can use your dataset here. Cell cycle
regression should be performed after the data has been filtered, normalised,
and scaled. At the end of this tutorial, you can return to the main tutorial
to plot and explore your data with reduced effects from the cell cycle genes.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the dataset Processed_AnnData
> 4. Check that the datatype is h5ad
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

In addition to the scRNA-seq dataset, we will also need lists of the genes
that are expressed at different points in the cell cycle. The lists used in
this tutorial were prepared by.......... and can be downloaded from Zenodo
below. 

The first list includes genes that are expressed during S Phase. The second
list contains genes that are expressed during the G2/M Phases. We don't need
a list of genes that are expressed in the G1 Phase because ........

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 2. Check that the datatype for both is tabular
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}



# Cell Cycle Scoring

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Inspect and manipulate](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_inspect/scanpy_inspect/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"Method used for inspecting"*: `Score cell cycle genes, using 'tl.score_genes_cell_cycle'`
>        - *"Format for the list of genes associated with S phase"*: `File`
>            - {% icon param-file %} *"File with the list of genes associated with S phase"*: `output` (Input dataset)
>        - *"Format for the list of genes associated with G2M phase"*: `File`
>            - {% icon param-file %} *"File with the list of genes associated with G2M phase"*: `output` (Input dataset)
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

# Cell Cycle Regression

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Scanpy RegressOut](toolshed.g2.bx.psu.edu/repos/ebi-gxa/scanpy_regress_variable/scanpy_regress_variable/1.8.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input object in AnnData/Loom format"*: `anndata_out` (output of **Inspect and manipulate** {% icon tool %})
>    - *"Variables to regress out"*: `phase`
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

# Plotting the Effects of Cell Cycle Regression

Your data is now ready for further analysis, so you can return to the Filter,
Plot and Explore Single-cell RNA-seq Data tutorial and move on to the
Preparing coordinates step there. However, if you want to understand how cell
cycle regression has affected your data then you might want to work through
the following steps first. 

In order to do this, we need to identify the cell cycle genes in our AnnData
dataset so that we can select them for plotting. We need to create a new
column to annotate the AnnData - you might find it easier to do this using a
spreadsheet and then upload the column as a tabular dataset, but it is
possible to complete all the steps on Galaxy. 

## Prepare a table of cell cycle genes
>If we're going to mark all the cell cycle genes, we'll need a list of all 97
>genes instead of the two separate lists for S Phase and G2/M Phase. We'll 
>join the two lists together and then add another column that just reads TRUE,
>which we'll use later to mark these as cell cycle genes in the main dataset.
>
> ### {% icon hands_on %} Hands-on: Create a list of all cell cycle genes
>
> 1. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `output` (Input dataset)
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `output` (Input dataset)
>
> 2. {% tool [Add column](toolshed.g2.bx.psu.edu/repos/devteam/add_value/addValue/1.0.0) %} with the following parameters:
>    - *"Add this value"*: `TRUE`
>    - {% icon param-file %} *"to Dataset"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>
> 3. Rename the dataset CC_Genes
{: .hands_on}

## Create an ordered list of gene names
> Next, we'll need a list of all the genes in our dataset, so that we can mark
> the ones that are in our cell cycle list. We'll also add a column of 
> numbers as this will help us keep the gene names in order.
> 
> ### {% icon hands_on %} Hands-on: Get the gene names from your dataset
>
> 1. {% tool [Inspect AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_inspect/anndata_inspect/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output` (Input dataset)
>    - *"What to inspect?"*: `Key-indexed annotation of variables/features (var)`
>
> 2. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `var` (output of **Inspect AnnData** {% icon tool %})
>        - *"Type of table operation"*: `Drop, keep or duplicate rows and columns`
>            - *"List of columns to select"*: `1`
>            - *"List of rows to select"*: `2:15396`
>    - *"Output formatting options"*: ``
>
>
>    > ### {% icon comment %} Comment
>    >
>    > Since we don't want to keep the header, we select the rows from 2 to the end
>    >  of the dataset. If you were using a dataset of a different size, you would
>    >  need to change this parameter to include all the genes. 
>    {: .comment}
>
> 3. {% tool [Add column](toolshed.g2.bx.psu.edu/repos/devteam/add_value/addValue/1.0.0) %} with the following parameters:
>    - {% icon param-file %} *"to Dataset"*: `table` (output of **Table Compute** {% icon tool %})
>    - *"Iterate?"*: `YES`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > Adding these numbers will enable us to keep the genes in their original
>    > order. This is essential for adding the cell cycle gene annotation back
>    > into the AnnData. 
>    {: .comment}
>
{: .hands_on}


## Mark the cell cycle genes
> We can now combine our table of cell cycle genes with the table of gene
> names. 
>  
> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Join](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_easyjoin_tool/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"1st file"*: `out_file1` (output of **Add column** {% icon tool %})
>    - *"Column to use from 1st file"*: `c1`
>    - {% icon param-file %} *"2nd File"*: `out_file1` (output of **Add column** {% icon tool %})
>    - *"Column to use from 2nd file"*: `c1`
>    - *"Output lines appearing in"*: `All lines [-a 1 -a 2]`
>    - *"Value to put in unpaired (empty) fields"*: `FALSE`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > When we do this, we'll ask for any empty fields to be filled in with
>    > FALSE. The cell cycle gene table has an extra column where they are
>    >  all marked as TRUE - they will retain this marking when we join the
>    >   tables but since there are no entries for the rest of the genes,
>    >    their rows will be filled in as FALSE. This will enable us to pick
>    >     out the cell cycle genes later. 
>    {: .comment}
>
> 2. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `output` (output of **Join** {% icon tool %})
>    - *"on column"*: `c2`
>    - *"everything in"*: `Ascending order`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > Sorting the genes using the numbers we added earlier will put them back in their original order - make sure to sort them in ascending order, otherwise they'll end up the opposite way around. 
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. What would happen if any of the cell cycle genes were not present in the dataset?
> 2. How would we remove these genes from the table?
>
> > ### {% icon solution %} Solution
> >
> > 1. Any cell cycle genes that weren't in the dataset would have an empty
> >  field in numbered column, which would be filled in with FALSE when we
> >   created the table. These rows would appear at the top of the table
> >    after it was sorted. 
> > 2. We should check the first rows of the table for any unnumbered genes
> >  and then cut these rows out in the next step. 
> >
> {: .solution}
>
{: .question}

## Add an annotation to the AnnData
We now have a table with all the gene names in the same order as the main
dataset and a column indicating which ones are cell cycle genes. If we cut
this column out of the table then we can add it as a new annotation to the
main dataset. We'll also need to add a column header, which will be used as
the key for this annotation in the AnnData. 

> ### {% icon hands_on %} Hands-on: Add the cell cycle annotation 
>
> 1. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} with the following parameters:
>    - *"Input Single or Multiple Tables"*: `Single Table`
>        - {% icon param-file %} *"Table"*: `out_file1` (output of **Sort** {% icon tool %})
>        - *"Type of table operation"*: `Drop, keep or duplicate rows and columns`
>            - *"List of columns to select"*: `3`
>            - *"List of rows to select"*: `1:15395`
>    - *"Output formatting options"*: ``
>
>
>    > ### {% icon comment %} Comment
>    >
>    > If there were any cell cycle genes that weren't present in the main
>    > dataset, we could remove them at this stage by excluding them from the
>    >  List of rows to select. 
>    {: .comment}
> 
> 2. Create a header for the new column 
>
>
> 3. {% tool [Concatenate datasets](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `output` (Input dataset)
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `table` (output of **Table Compute** {% icon tool %})
>
>
{: .hands_on}


## Sub-step with **Manipulate AnnData**
We will need to add the annotation to both the original dataset and to the one that we created by regressing out the cell cycle genes. This will allow us to plot the cell cycle genes before and after regression. 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Inspect and manipulate** {% icon tool %})
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>        - {% icon param-file %} *"Table with new annotations"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>
>
>
> 2. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `output_h5ad` (output of **Scanpy RegressOut** {% icon tool %})
>    - *"Function to manipulate the object"*: `Add new annotation(s) for observations or variables`
>        - {% icon param-file %} *"Table with new annotations"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>
## Filter the cell cycle genes
For both the nely annotated datasets
> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Function to manipulate the object"*: `Filter observations or variables`
>        - *"Type of filtering?"*: `By key (column) values`
>            - *"Key to filter"*: `CC_genes`
>            - *"Type of value to filter"*: `Boolean`
>
> 2. {% tool [Manipulate AnnData](toolshed.g2.bx.psu.edu/repos/iuc/anndata_manipulate/anndata_manipulate/0.7.5+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Function to manipulate the object"*: `Filter observations or variables`
>        - *"Type of filtering?"*: `By key (column) values`
>            - *"Key to filter"*: `CC_genes`
>            - *"Type of value to filter"*: `Boolean`
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

## Plot the cell cycle genes before regression

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Cluster, infer trajectories and embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used"*: `Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using 'tl.pca'`
>        - *"Type of PCA?"*: `Full PCA`
>
> 2. {% tool [Cluster, infer trajectories and embed](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_cluster_reduce_dimension/scanpy_cluster_reduce_dimension/1.7.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata` (output of **Manipulate AnnData** {% icon tool %})
>    - *"Method used"*: `Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using 'tl.pca'`
>        - *"Type of PCA?"*: `Full PCA`
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

## Plot the cell cycle genes after regression

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.7.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Cluster, infer trajectories and embed** {% icon tool %})
>    - *"Method used for plotting"*: `PCA: Plot PCA results, using 'pl.pca_overview'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `phase`
>        - In *"Plot attributes"*:
>            - *"Colors to use for plotting categorical annotation groups"*: `rainbow (Miscellaneous)`
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
> 2. {% tool [Plot](toolshed.g2.bx.psu.edu/repos/iuc/scanpy_plot/scanpy_plot/1.7.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Annotated data matrix"*: `anndata_out` (output of **Cluster, infer trajectories and embed** {% icon tool %})
>    - *"Method used for plotting"*: `PCA: Plot PCA results, using 'pl.pca_overview'`
>        - *"Keys for annotations of observations/cells or variables/genes"*: `phase`
>        - In *"Plot attributes"*:
>            - *"Colors to use for plotting categorical annotation groups"*: `rainbow (Miscellaneous)`
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
