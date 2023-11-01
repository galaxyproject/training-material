---
layout: tutorial_hands_on

title: 'Cite-Seq Tool Data Processing into RStudio Visualization (Cite-Seq, Seurat, R)'
subtopic: single-cell-CS-code
priority: 2

questions:
- How can I use Seurat's Cite-Seq capabilities?
- How can I visualize and interpret multimodal data in Seurat?

objectives:
- Learn to use Galaxy's Seurat tool with Cite-seq capabilities to create a Seurat Object 
- Understand the parameters of the Seurat tool 
- Move between Galaxy and RStudio to holistically explore Cite-Seq data

time_estimation: 3H

key_points:
- Being able to switch between Galaxy and RStudio when analyzing datasets can be useful when looking to adjust default parameters within Seurat's functions and workflow.
- Seurat in RStudio gives more flexibility and specificity of analyses, but Galaxy offers great reproducibility and ease of analysis.
- Beginning to learn the syntax and use of R will expand your access to bioinformatic analyses. 

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        
tags:
- Cite-Seq
- RStudio

contributions:
  authorship:
    - Camila-goclowski

notebook:
  language: r
  snippet: topics/single-cell/tutorials/scCiteSeq-Tool-to-RStudio/preamble.md
---

{% snippet topics/single-cell/faqs/notebook_warning.md %}

Before we can do any real biological investigation, we need to understand what each of the outputs from our Seurat tool are. Maybe you've already begun to dissect what's what, but just in case, let's run through each of the datasets here, together. 

><comment-title>gx_get()</comment-title>
> RStudio in galaxy comes with a gx_get() function. This function is critical to understand and be able to use in order to move datasets from your history and into RStudio. The function will output the file path with which you can access the data via RStudio.
> To use it, simply use the numbered location of the dataset you are looking to import. For example: 
> If we want to find the first dataset we imported, simply run the following command: 
```{r}
gx_get(1)
```
>The result of this command will be a file path to the first dataset in your galaxy history. Use that file path for importing purposes. 
{: .comment}

To take a look at the pre analysis RNA-seq matrix, use the following commands: 
```{r}
gx_get(1)
RNA<-read.csv('/import/1')
```
Note that the dataset we are using also contains ~5% of mouse cells, which we can use as negative controls for the cell surface protein measurements. As such, the RNA expression matrix has "HUMAN_" or "MOUSE_" appended to each gene. 

Now let's take a look at what's in here. 
```{r}
view(RNA)
```
If you're familiar with scRNA-seq matrices, this may look familiar to you. That's because it is exactly that--an RNA-seq matrix! In these matrices we have genes as row names and cell barcodes as column names. The values within the matrix denote the number of transcripts from a given gene within a given cell.

You may have noticed there are TONS of zero values in this matrix. You may also be thinking, "Won't that create noise in the dataset??" The answer is yes, and these zeros are one of the first things that the Seurat preprocessing tool will accomplish. This matrix that we've labelled as RNA is *not* what we will be analyzing further into this tutorial. We are simply taking a look to ground ourselves in what the data looked like *before* preprocessing. 

We can do the same thing with the pre-analysis protein matrix. We'll call it the ADT matrix for now, since that is how Seurat recognizes it! 
```{r}
gx_get(2)
ADT<-read.csv('/import/2')
```
Again, let's take a look at what's in here:
```{r}
view(ADT)
```
Looks shockingly similar, doesn't it?!

In the ADT matrix, we have cell surface proteins (instead of gene names) as row names and the same cell barcodes as column names. 

If you ran the same parameters as I did, the next output (number 3 in our history) will be Seurat's run log. This is unfortunately not super easy to import into RStudio since it comes as an html format. It contains all of the run information from the background coding done by the tool. Any warnings, errors, or progress bars will be present in here and are often useful for troubleshooting in case something goes awry. Because of the html formatting, we will not look at this output together, but feel free to explore it on your own using the view (eye) icon in your history. 

The next output in my galaxy history are protein markers! Let's take a look: 
```{r}
gx_get(4)
protein_markers<-read.table('/import/4', header = T)
view(protein_markers)
```
There are tons of markers in this list and if you look closely, you'll see that some are not statistically significant. Let's take care of that and filter out any marker that has an adjusted p-value above 0.045:
```{r}
protein_markers<-subset(protein_markers, p_val_adj < 0.045)
```

Now we have a statistically signficant list of protein markers per cluster! There are a number of statistics that are included here, if you're interested in better understanding them, take a look at [Seurat's documentation of FindAllMarkers] (https://satijalab.org/seurat/reference/findallmarkers) for more details. 

The next dataset in our history should be RNA markers. Let's import them, remove the statistically insignifcant ones, and take a look: 
```{r}
gx_get(5)
rna_markers<-read.table('/import/5', header = T)
rna_markers<-subset(rna_markers, p_val_adj < 0.045)
view(rna_markers)
```

Just like the RNA and ADT matrices looked quite similar, the protein and rna markers will as well. This is because Seurat is interpretting and analyzing the RNA and ADT assays in the same manner, with the same tools. So once again, if you're interested in what some of the statistic on the rna_markers file mean, take a look at the [Seurat documentation of FindAllMarkers] (https://satijalab.org/seurat/reference/findallmarkers). 

