---
layout: tutorial_hands_on

title: 'Trajectory Analysis: Monocle3 in RStudio'
subtopic: single-cell-CS
priority: 6
zenodo_link: 'https://zenodo.org/record/7455590'

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

time_estimation: 1H

key_points:
- The take-home messages
- They will appear at the end of the tutorial

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin
        - scrna-case_alevin-combine-datasets
        - scrna-case_basic-pipeline
        - scrna-case_JUPYTER-trajectories
        - scrna-case_monocle3-trajectories
        
tags:
- single-cell
- trajectory-analysis
- paper-replication

contributions:
  authorship:
    - wee-snufkin

notebook:
  language: r
  snippet: topics/single-cell/tutorials/scrna-case_monocle3-rstudio/preamble.md 
---

# Working in RStudio

## Uploading files

If you are using RStudio locally, then you don’t have to bother about uploading the files – just skip to ```file.choose()``` and navigate to your files using the pop-up window. 

If you are working in RStudio Cloud, you have to download the generated files from your history first. To do so, just click on the {% icon galaxy-save %} save icon for `Cell metadata (obs)`, `Gene metadata (var)` and `Expression matrix`. Then, return to the RStudio and click on ‘Upload’ button in the right bottom window toolbar and choose already downloaded files to upload. You should now see all three filed in this window. You might want to rename the files to make their names shorter.

![Screenshot of Files tab in RStudio, highlighting 'Upload' and 'Rename' buttons and listing three uploaded and renamed files: 'cell_metadata', 'gene_metadata', 'expression_matrix'.](../../images/scrna-casestudy-monocle/r_files_tab.png "The view of the Files tab with uploaded files and highlighted relevant buttons.")

If you are using RStudio Galaxy tool, you can get data directly from your history by running:
```r
file.copy(gx_get(2), "cell_metadata")
file.copy(gx_get(3), "gene_metadata")
file.copy(gx_get(4), "expression_matrix")
```
The number in the brackets corresponds to the dataset number in your history, so make sure you put the right number for the corresponding file. We can specify the name of the fetched files in the quotation marks, as shown above. All three files should appear in the Files tab window.

> <question-title></question-title>
>
> What is the datatype of the uploaded files?
>
> > <solution-title></solution-title>
> >
> > If you first downloaded the files from Galaxy and then uploaded them into RStudio, you should be able to see the extension `.tabular`. You might not see it if you fetched data directly from your history, but you can easily check the type of data in Galaxy, and - what's more - change it there. But today we're focusing on coding!
> >
> {: .solution}
>
{: .question}

Once we have our data loaded, let's specify paths to access the files. There is a strightforward way of getting the correct path without the concern of making typos or getting the path wrong. Just run `file.choose()` and choose the corresponding file which you want to get path to:
```r
cells_path <- file.choose() 
genes_path <- file.choose() 
expression_path <- file.choose() 
```
You should now see the new variables in the Environment tab window. 

As mentioned above, the datatype of our files is tabular, so we will use ```read.delim()``` function to read them in. The first argument is the file path and the second one, `row.names=1` takes the column number of the data file from which to take the row names. 
```r
cell_metadata <- read.delim(cells_path, row.names=1)
gene_metadata <- read.delim(genes_path, row.names=1)
expression_matrix <- read.delim(expression_path, row.names=1)
```

We have now three dataframes that we will use to generate cell_data_set object.

> <question-title></question-title>
>
> Why should we set `row.names=1`?
>
> > <solution-title></solution-title>
> >
> > This allows us to ensure that the expression value matrix has the same number of columns as the `cell_metadata` has rows and the same number of rows as the `gene_metadata` has rows. Importantly, row names of the `cell_metadata` object should match the column names of the expression matrix and row names of the `gene_metadata` object should match row names of the expression matrix.
> >
> {: .solution}
>
{: .question}

## View and modify files

According to [Monocle3 documentation](https://cole-trapnell-lab.github.io/monocle3/docs/starting/), `expression_matrix` should have genes as rows and cells as columns. Let's check if that's the case here.
```r
View(expression_matrix)
```
`View()` opens a new tab with a preview of the content of the file. We can see that in our matrix rows are cells and genes are columns, so we have to transpose the matrix simply using function `t()`. But before doing so, we will change its type from dataframe to matrix - this is Monocle's requirement to generate cell_data_set afterwards.
```r
expression_matrix <- as.matrix(expression_matrix)
expression_matrix <- t(expression_matrix)
```
Another condition we have to satisfy if that one of the columns of the `gene_metadata` should be named "gene_short_name", which represents the gene symbol for each gene. Some functions won't work without that. Do we have such a column? Let's check.
```r
View(gene_metadata)
```

The second column indeed contains gene symbols, but is called "Symbol" instead of "gene_short_name". That can be easily changed by a simple assignment, as long as we know the number of the column that we want to modify. We can access the column names by `colnames()`
```r
colnames(gene_metadata)[2] <- 'gene_short_name'
```

You can now switch to the `gene_metadata` tab and check if the name has changed. 

It looks like our data fulfils the requirements to generate the cell_data_set file, but before that…


## Loading Monocle3

First things first, we need to load Monocle3! Generally, it is a good practice to load all the packages at the very beginning of the script, but for the purpose of this tutorial, we will load the needed packages on the way, to make you aware of which library is needed and when. 
Monocle 3 runs in the R statistical computing environment. You will need R version 4.1.0 or higher, Bioconductor version 3.14, and monocle3 1.2.7 or higher to have access to the latest features.
```r
# Install Bioconductor and some of required dependencies

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

# Install monocle3 through the cole-trapnell-lab GitHub:
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
```

> <warning-title>Installation errors</warning-title>
> Sometimes there might be some changes in repositories, packages or dependencies and you might also use different versions of them. Therefore, it might happen that you would need to carefully read the error messages and follow the suggestions. If you are facing any difficulties with installation process, it is recommended that you consult your problem with additional online resources. It is more likely that RStudio Cloud or Galaxy tool would fail rather than local RStudio.
>
{: .warning}

## Generating CDS object
Now let’s store our files in one object – the cell_data_set. This is the main class used by Monocle to hold single cell expression data. The class is derived from the Bioconductor SingleCellExperiment class. It's similar to Python's AnnData storing a data matrix together with annotations of observations and variables. There are three ways of creating CDS object in monocle:
-	Using ```new_cell_data_set() ``` function with three data frames as arguments (not their paths!): expression matrix (can also be a sparseMatrix), cell metadata and gene metadata
-	Using ```load_cellranger_data() ``` function and providing the path to the folder containing 10X Genomics Cell Ranger output files. This function takes an argument `umi_cutoff` that determines how many reads a cell must have to be included
-	Using ```load_mm_data() ``` function providing the paths to matrix file and two metadata files (features and cell information)
In this tutorial we will use the first option:
```r
cds <- new_cell_data_set(expression_matrix, cells_metadata, genes_metadata)
```
We are now ready to process our data!

> <details-title>Format conversion</details-title>
>
>  Since Monocle’s CDS object is analogous to Python's AnnData, why don’t we use some kind of conversion between those two formats? There is indeed a package called `sceasy` that helps easy conversion of different single-cell data formats to each other. However, when we tested this conversion on our dataset and then used Monocle to plot the expression of genes, the plots were not correct – the expression was shown to be ideantical throughout the sample. For comparison, Seurat did well when plotting gene expression of the same converted object! Although conversion functions are very handy, you have to be aware that their output might be interpreted differently by certain packages. Therefore, to make sure that the analysis is reliable, we decided to generate CDS object directly using Monocle’s function. 
> ![Comparison between plots of gene expression generated by Monocle and Seurat using the CDS object that was converted from AnnData using SCEasy tool. Gene expression in Monocle plots is identical throughout the sample, while the expression of the same genes plotted by Seurat is noticeable only in specific clusters.](../../images/scrna-casestudy-monocle/monocle_seurat.png "Comparison between plots of gene expression generated by Monocle and Seurat using the CDS object that was converted from AnnData using SCEasy tool.")
>
{: .details}

# Additional step: adding genes symbols based on their IDs
> <warning-title>Additional step</warning-title>
>  This step is not necessary for the dataset we are working on but some users might find it helpful when analysing their own data.
>
{: .warning}

If you remember the very first tutorial, we were starting with gene IDs and later on adding gene symbols based on the Ensembl GTF file.  
But what if we didn’t have the genes symbols in our CDS object and wanted to add them now? Of course - it's possible! We will also base this annotation on Ensembl - the genome database – with the use of the library BioMart. 
```r
cds_extra <- cds		# assign our CDS to a new object for the demonstration purpose 
library("biomaRt")		# load the BioMart library
ensembl.ids <- rownames(fData(cds_extra))		# fData() allows to access cds rowData table and the rownames are stored in ensembl.ids
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL") # connect to a specified BioMart database and dataset hosted by Ensembl
ensembl_m = useMart("ensembl", dataset="mmusculus_gene_ensembl", 
                    host='https://nov2020.archive.ensembl.org')
#ensembl_m = useMart("ensembl", dataset="mmusculus_gene_ensembl") # connect to a specified BioMart database and dataset within this database; in our case we choose the mus musculus database
genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
               filters = 'ensembl_gene_id', 
               values = ensembl.ids, 
               mart = ensembl_m) # retrieve the specified attributes from the connected BioMart database; 'ensembl_gene_id' are genes IDs, 'external_gene_name' are the genes symbols that we want to get for our values stored in ‘ensembl.ids’

listDatasets(mart)

gene_names <- rownames(fData(cds_a)) #genes IDs
#replace IDs for gene names
count = 1
for (geneID in gene_names)
{
  print(geneID)
  index <- which(genes==geneID)
  print(index)
  gene_names[count] <- genes$external_gene_name[index]
  count = count + 1
}

fData(cds_extra)$gene_short_name2 <- gene_names
```



# Monocle workflow
Do you remember the Monocle workflow introduced in the previous tutorial? Here is a recap:
![Monocle workflow: scRNA-seq dataset, pre-process data (normalise, remove batch effects), non-linear dimensionality reduction (t-SNE, UMAP), cluster cells, compare clusters (identify top markers, targeted contrasts), trajectory analysis](../../images/scrna-casestudy-monocle/monocle3_new_workflow.png "Workflow provided by Monocle3 documentation")

Therefore, let’s start with normalisation and pre-processing that can be performed using the function `preprocess_cds()`. The argument `num_dim` is the number of principal components that will be computed when using PCA during normalisation. The overall ‘shape’ of the plot changes slightly as the value of `num_dim` changes, but for this demonstration we will use the value of 210, which, compared to the results from the previous tutorial, makes the most sense for our dataset. Then you can check that you're using enough PCs to capture most of the variation in gene expression across all the cells in the data set. 

```r
cds <- preprocess_cds(cds, num_dim = 210)
plot_pc_variance_explained(cds)
```

![Six plots showing only the shaded shape of how the cells are clustered depending on the num_dim argument. The general trend is maintained though](../../images/scrna-casestudy-monocle/num_dim.jpg "The ‘shape’ of the plot showing how the cells are clustered depending on the num_dim argument.")



# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

><hands-on-title> Data upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from
>
{: .hands_on}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
