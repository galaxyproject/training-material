---
layout: tutorial_hands_on

title: Filter, plot and explore single-cell RNA-seq data (Scanpy, Python)
subtopic: single-cell-CS-code
priority: 2
zenodo_link: 'https://zenodo.org/record/7053673'

questions:
- Is my single cell dataset a quality dataset?
- How do I generate and annotate cell clusters?
- How do I pick thresholds and parameters in my analysis? What's a "reasonable" number, and will the world collapse if I pick the wrong one?
objectives:
- Interpret quality control plots to direct parameter decisions
- Repeat analysis from matrix to clustering
- Identify decision-making points
- Appraise data outputs and decisions
- Explain why single cell analysis is an iterative (i.e. the first plots you generate are not final, but rather you go back and re-analyse your data repeatedly) process
time_estimation: 3H
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
- Python

contributions:
  authorship:
    - hexhowells
    - nomadscientist

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_JUPYTER-trajectories
        - scrna-case_monocle3-trajectories

notebook:
  language: python
  snippet: topics/single-cell/tutorials/scrna-case-jupyter_basic-pipeline/preamble.md

---

{% snippet topics/single-cell/faqs/notebook_warning.md %}

# Install libraries

This tutorial requies some libraries to be installed which is done below (igraph and louvain are not used directly and are just required for plotting). The ```-q``` parameter hides most of the outputs of the installation in order to make the notebook a bit cleaner. If there are any issues with the installation, then removing this parameter may give you more information about the issue.

```python
pip install scanpy -q
```

```python
pip install igraph -q
```

```python
pip install louvain -q
```

```python
pip install pandas -q
```

---

We can now import the two libraries that we will be using, **scanpy** is the primary library that we will use and will handle all the plotting and data processing. Meanwhile, **pandas** is used briefly for some manual data manipulation.

```python
import scanpy as sc
import pandas as pd
```

# Load Data

You can import files from your Galaxy history directly using the following code. This will depend on what number in your history the final annotated object is. If your object is dataset #4 in your history, then you import it with the following:

```python
mito_counted_anndata = get(4)
```

You now need to read it in as a h5ad object.

```python
adata = sc.read_h5ad(mito_counted_anndata)
```

# Filtering

You have generated an annotated AnnData object from your raw scRNA-seq fastq files. However, you have only completed a ‘rough’ filter of your dataset - there will still be a number of ‘cells’ that are actually just background from empty droplets or simply low-quality. There will also be genes that could be sequencing artifacts or that appear with such low frequency that statistical tools will fail to analyse them. This background garbage of both cells and genes not only makes it harder to distinguish real biological information from the noise, but also makes it computationally heavy to analyse. These spurious reads take a lot of computational power to analyse! First on our agenda is to filter this matrix to give us cleaner data to extract meaningful insight from, and to allow faster analysis.

> <question-title></question-title>
>
> 1. What information is stored in your AnnData object? The last tool to generate this object counted the mitochondrial associated genes in your matrix. Where is that data stored?
> 2. While you are figuring that out, how many genes and cells are in your object?
>
>   > <tip-title>Hint</tip-title>
>   > Inspect the Anndata object by printing it with:
>   > ```python
>   > print(adata)
>   >
>   > print(adata.obs)
>   >
>   > print(adata.var)
>   > ```
>   {: .tip}
> > <solution-title></solution-title>
> >
> > 1. If you examine your AnnData object, you'll find a number of different quality control metrics for both cells (**obs**) and genes (**var**).
> >   - For instance, you can see a `n_cells` under **var**, which counts the number of cells that gene appears in.
> >   - In the **obs**, you have both discrete and log-based metrics for `n_genes`, how many genes are counted in a cell, and `n_counts`, how many UMIs are counted per cell. So, for instance, you might count multiple GAPDHs in a cell. Your `n_counts` should thus be higher than `n_genes`.
> >   - But what about the mitochondria?? Within the cells information **obs**, the `total_counts_mito`,  `log1p_total_counts_mito`, and `pct_counts_mito` has been calculated for each cell.
> > 2. You can see by printing the object that the matrix is `31178 x 35734`. This is `obs x vars`, or rather, `cells x genes`, so there are `31178 cells` and `35734 genes` in the matrix.
> >
> {: .solution}
>
{: .question}

## Generate QC Plots

We want to filter our cells, but first we need to know what our data looks like. There are a number of subjective choices to make within scRNA-seq analysis, for instance we now need to make our best informed decisions about where to set our thresholds (more on that soon!). We’re going to plot our data a few different ways. Different bioinformaticians might prefer to see the data in different ways, and here we are only generating some of the myriad of plots you can use. Ultimately you need to go with what makes the most sense to you.

## Creating the Plots

```python
# Violin - genotype - log
sc.pl.violin(
  adata,
  keys=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mito'],
  groupby='genotype',
  save='-genotype-log.png'
)
```


```python
# Violin - sex - log
sc.pl.violin(
  adata,
  keys=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mito'],
  groupby='sex',
  save='-sex-log.png'
)
```

```python
# Violin - batch - log
sc.pl.violin(
  adata,
  keys=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mito'],
  groupby='batch',
  save='-batch-log.png'
)
```

```python
# Scatter - mito x UMIs
sc.pl.scatter(
  adata,
  x='log1p_total_counts',
  y='pct_counts_mito',
  save='-mitoxUMIs.png'
)
```

```python
# Scatter - mito x genes
sc.pl.scatter(
  adata,
  x='log1p_n_genes_by_counts',
  y='pct_counts_mito',
  save='-mitoxgenes.png'
)
```

```python
# Scatter - genes x UMIs
sc.pl.scatter(
  adata,
  x='log1p_total_counts',
  y='log1p_n_genes_by_counts',
  color='pct_counts_mito',
  save='-genesxUMIs.png'
)
```

## Analysing the plots

That’s a lot of information! Let’s attack this in sections and see what questions these plots can help us answer.

> <question-title>Batch Variation</question-title>
>
> Are there differences in sequencing depth across the samples?
> 1. Which plot(s) addresses this?
> 2. How do you interpret it?
>
> > <solution-title></solution-title>
> >
> > 1. The plot `violin - batch - log` will have what you're looking for!
> >     ![Violin - batch - log](../../images/scrna-casestudy-jupyter/violin-batch-log.png "Violin - batch - log (Raw)")
> > 2. Keeping in mind that this is a log scale - which means that small differences can mean large differences - the violin plots probably look pretty similar.
> >    - `N703` and `N707` might be a bit lower on genes and counts (or UMIs), but the differences aren't catastrophic.
> >    - The `pct_counts_mito` looks pretty similar across the batches, so this also looks good.
> >    - Nothing here would cause us to eliminate a sample from our analysis, but if you see a sample looking completely different from the rest, you would need to question why that is and consider eliminating it from your experiment!
> >
> {: .solution}
>
{: .question}

> <question-title>Biological Variables</question-title>
>
> Are there differences in sequencing depth across sex? Genotype?
> 1. Which plot(s) addresses this?
> 2. How do you interpret the `sex` differences?
> 3. How do you interpret the `genotype` differences?
>
> > <solution-title></solution-title>
> >
> > 1. Similar to above, the plots `violin - sex - log` and `violin - genotype - log` will have what you're looking for!
> >      ![Violin - sex - log](../../images/scrna-casestudy-jupyter/violin-sex-log.png "Violin - sex - log (Raw)")
> >      ![Violin - genotype - log](../../images/scrna-casestudy-jupyter/violin-genotype-log.png "Violin - genotype - log (Raw)")
> > 2. There isn’t a major difference in sequencing depth across sex, I would say - though you are welcome to disagree!
> >    - It is clear there are far fewer female cells, which makes sense given that only one sample was female. *Note - that was an unfortunate discovery made long after generating libraries. It's quite hard to identify the sex of a neonate in the lab! In practice, try hard to not let such a confounding factor into your data! You could consider re-running all the following analysis without that female sample, if you wish.*
> > 3. In `Violin - genotype - log`, however, we can see there is a difference. The `knockout` samples clearly have fewer genes and counts. From an experimental point of view, we can consider, does this make sense?
> >    - Would we biologically expect that those cells would be smaller or having fewer transcripts? Possibly, in this case, given that these cells were generated by growth restricted neonatal mice, and in which case we don’t need to worry about our good data, but rather keep this in mind when generating clusters, as we don’t want depth to define clusters, we want biology to!
> >    - On the other hand, it may be that those cells didn’t survive dissociation as well as the healthy ones (in which case we’d expect higher mitochondrial-associated genes, which we don’t see, so we can rule that out!).
> >    - Maybe we unluckily poorly prepared libraries for specifically those knockout samples. There are only three, so maybe those samples are under-sequenced.
> >    - So what do we do about all of this?
> >        - Ideally, we consider re-sequencing all the samples but with a higher concentration of the knockout samples in the library. Any bioinformatician will tell you that the best way to get clean data is in the lab, not the computer! Sadly, absolute best practice isn’t necessarily always a realistic option in the lab - for instance, that mouse line was long gone! - so sometimes, we have to make the best of it. There are options to try and address such discrepancy in sequencing depth. Thus, we’re going to take these samples forward and see if we can find biological insight despite the technical differences.
> >
> {: .solution}
>
{: .question}

Now that we’ve assessed the differences in our samples, we will look at the libraries overall to identify appropriate thresholds for our analysis.

> <question-title>Filter Thresholds: genes</question-title>
>
> What threshold should you set for `log1p_n_genes_by_counts`?
> 1. Which plot(s) addresses this?
> 2. What number would you pick?
>
> > <solution-title></solution-title>
> >
> > 1. Any plot with `log1p_n_genes_by_counts` would do here, actually! Some people prefer scatterplots to violins.
> >    ![Scatter-genesxmito](../../images/scrna-casestudy-jupyter/scatter-mito-genes.png "Scatter - mito x genes (Raw)")
> >
> > 2. In `Scatter - mito x genes` you can see how cells with `log1p_n_genes_by_counts` up to around, perhaps, `5.7` (around 300 genes) often have high `pct_counts_mito`.
> >   - You can plot this as just `n_counts` and see this same trend at around 300 genes, but with this data the log format is clearer so that's how we're presenting it.
> >   - You could also use the violin plots to come up with the threshold, and thus also take batch into account. It's good to look at the violins as well, because you don't want to accidentally cut out an entire sample (i.e. N703 and N707).
> >   - Some bioinformaticians would recommend filtering each sample individually, but this is difficult in larger scale and in this case (you're welcome to give it a go! You'd have to filter separately and then concatenate), it won't make a notable difference in the final interpretation.
> >
> {: .solution}
>
{: .question}

> <question-title>Filter Thresholds: UMIs</question-title>
>
> What threshold should you set for `log1p_total_counts`?
> 1. Which plot(s) addresses this?
> 2. What number would you pick?
>
> > <solution-title></solution-title>
> >
> > 1. As before, any plot with `log1p_total_counts` will do! Again, we'll use a scatterplot here, but you can use a violin plot if you wish!
> >    ![Scatter-countsxmito](../../images/scrna-casestudy-jupyter/scatter-mito-umis.png "Scatterplot - mito x UMIs (Raw)")
> >
> > 2. We can see that we will need to set a higher threshold (which makes sense, as you'd expect more UMI's per cell rather than unique genes!). Again, perhaps being a bit aggressive in our threshold, we might choose `6.3`, for instance (which amounts to around 500 counts/cell).
> >   - In an ideal world, you'll see a clear population of real cells separated from a clear population of debris. Many samples, like this one, are under-sequenced, and such separation would likely be seen after deeper sequencing!
> >
> {: .solution}
>
{: .question}

> <question-title>Filter Thresholds: mito</question-title>
>
> What threshold should you set for `pct_counts_mito`?
> 1. Which plot(s) addresses this?
> 2. What number would you pick?
>
> > <solution-title></solution-title>
> >
> > 1. Any plot with `pct_counts_mito` would do here, however the scatterplots are likely the easiest to interpret. We'll use the same as last time.
> >    ![Scatter-countsxmito](../../images/scrna-casestudy-jupyter/scatter-mito-umis.png "Scatterplot - mito x UMIs (Raw)")
> >
> > 2. We can see a clear trend wherein cells that have around 5% mito counts or higher also have far fewer total counts. These cells are low quality, will muddy our data, and are likely stressed or ruptured prior to encapsulation in a droplet. While 5% is quite a common cut-off, this is quite messy data, so just for kicks we'll go more aggressive with a `4.5%`.
> >    - In general, you must adapt all cut-offs to your data - metabolically active cells might have higher mitochondrial RNA in general, and you don't want to lose a cell population because of a cut-off.
> >
> {: .solution}
>
{: .question}

## Applying the Thresholds

It’s now time to apply these thresholds to our data! First, a reminder of how many cells and genes are in your object: ```31178 cells``` and ```35734 genes```. Let’s see how that changes each time!

> <details-title>Working in a group? Decision-time!</details-title>
> If you are working in a group, you can now divide up a decision here with one *control* and the rest varied numbers so that you can compare results throughout the tutorials.
> - Control
>      - **log1p_n_genes_by_counts** > `5.7`
>      - **log1p_total_counts** > `6.3`
>      - **pct_counts_mito** < `4.5%`
> - Everyone else: Choose your own thresholds and compare results!
{: .details}

We will plot the raw data before applying any filters so that we can more clearly see the changes we will make.

```python
# Raw
sc.pl.violin(
  adata,
  keys=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mito'],
  groupby='genotype',
  save='-raw.png'
)
```

```python
genes_filtered_obj = adata[adata.obs['log1p_n_genes_by_counts'] >= 5.7]
genes_filtered_obj = genes_filtered_obj[genes_filtered_obj.obs['log1p_n_genes_by_counts'] <= 20.0]

# Violin - Filterbygenes
sc.pl.violin(
  genes_filtered_obj,
  keys=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mito'],
  groupby='genotype',
  save='-Filterbygenes.png'
)
```

```python
print(genes_filtered_obj)
```

> <question-title></question-title>
>
> 1. Interpret the violin plot
> 2. How many genes & cells do you have in your object now?
>
> > <solution-title></solution-title>
> >
> > ![Violinplot-filteronce](../../images/scrna-casestudy-jupyter/violin-raw-filteredgenes.png "Raw vs 1st filter - genes/cell")
> > 1. The only part that seems to change is the `log1p_n_genes_by_counts`.  You can see a flatter bottom to the violin plot - this is the lower threshold set. Ideally, this would create a beautiful violin plot because there would be a clear population of low-gene number cells. Sadly not the case here, but still a reasonable filter.
> > 2. In the printed AnnData information, you can see you now have `17,040 cells x 35,734 genes`.
> >
> {: .solution}
>
{: .question}


```python
counts_filtered_obj = genes_filtered_obj[genes_filtered_obj.obs['log1p_total_counts'] >=  6.3]
counts_filtered_obj = counts_filtered_obj[counts_filtered_obj.obs['log1p_total_counts'] <= 20.0]

# Violin - Filterbycounts
sc.pl.violin(
  counts_filtered_obj,
  keys=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mito'],
  groupby='genotype',
  save='-Filterbycounts.png'
)
```

```python
print(counts_filtered_obj)
```

> <question-title></question-title>
>
> 1. Interpret the violin plot
> 2. How many genes & cells do you have in your object now?
>
> > <solution-title></solution-title>
> >
> > ![Violinplot-filtertwice](../../images/scrna-casestudy-jupyter/violin-filteredgenesxfilteredcounts.png "1st filter vs 2nd filter - counts/cell")
> > 1. We will focus on the `log1p_total_counts` as that shows the biggest change. Similar to above, the bottom of the violin shape has flattered due to the threshold.
> > 2. In the printed AnnData information, you can see you now have `8,678 cells x 35,734 genes`.
> >
> {: .solution}
>
{: .question}

```python
mito_filtered_obj = counts_filtered_obj[counts_filtered_obj.obs['pct_counts_mito'] >=  0]
mito_filtered_obj = mito_filtered_obj[mito_filtered_obj.obs['pct_counts_mito'] <= 4.5]

# Violin - Filterbymito
sc.pl.violin(
  mito_filtered_obj,
  keys=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mito'],
  groupby='genotype',
  save='-Filterbymito.png'
)
```

```python
print(mito_filtered_obj)
```

> <question-title></question-title>
>
> 1. Interpret the violin plot
> 2. How many genes & cells do you have in your object now?
>
> > <solution-title></solution-title>
> >
> > ![Violinplot-filtermito](../../images/scrna-casestudy-jupyter/violin-mitofilter.png "Violin plots after filtering genes, counts, and mito content/cell")
> > 1. If we carefully check the axes, we can see that the `pct_counts_mito` has shrunk.
> > 2. In the printed AnnData information, you can see you now have `8,605 cells x 35,734 genes`.
> >
> {: .solution}
>
{: .question}

Here’s a quick overall summary for easy visualisation if you fancy it.

![Filtering summary](../../images/scrna-casestudy-jupyter/violins2gether.png "Filtering summary")

Fantastic work! However, you’ve now removed a whole heap of cells, and since the captured genes are sporadic (i.e. a small percentage of the overall transcriptome per cell) this means there are a number of genes in your matrix that are currently not in any of the remaining cells. Genes that do not appear in any cell, or even in only 1 or 2 cells, will make some analytical tools break and overall will not be biologically informative. So let’s remove them! Note that ```3``` is not necessarily the best number, rather it is a fairly conservative threshold. You could go as high as 10 or more.

> <details-title>Working in a group? Decision-time!</details-title>
> If you are working in a group, you can now divide up a decision here with one *control* and the rest varied numbers so that you can compare results throughout the tutorials.
> - **min_cells** = `3`
> - Everyone else: Choose your own thresholds and compare results! Note if you go less than 3 (or even remove this step entirely), future tools are likely to fail due to empty gene data.
{: .details}

```python
filtered_obj = mito_filtered_obj.copy()

sc.pp.filter_genes(filtered_obj, min_cells=3)
sc.pp.filter_genes(filtered_obj, max_cells=1000000000)

print(filtered_obj)
```

In practice, you’ll likely choose your thresholds then set up all these filters to run without checking plots in between each one. But it’s nice to see how they work!

Using the final ```filtered_object```, we can summarise the results of our filtering:

|                    | Description | Genes |
| ------------------ | ----------- | ----- |
| Raw                | 31178       | 35734 |
| Filter genes/cell  | 17040       | 35734 |
| Filter counts/cell | 8678        | 35734 |
| Filter mito/cell   | 8605        | 35734 |
| Filter cells/gene  | 8605        | 15395 |

{% icon congratulations %} Congratulations! You have filtered your object! Now it should be a lot easier to analyse.

# Processing

So currently, you have a matrix that is 8605 cells by 15395 genes. This is still quite big data. We have two issues here - firstly, you already know there are differences in how many transcripts and genes have been counted per cell. This technical variable can obscure biological differences. Secondly, we like to plot things on x/y plots, so for instance Gapdh could be on one axis, and Actin can be on another, and you plot cells on that 2-dimensional axis based on how many of each transcript they possess. While that would be fine, adding in a 3rd dimension (or, indeed, in this case, 15393 more dimensions), is a bit trickier! So our next steps are to transform our big data object into something that is easy to analyse and easy to visualise.

```python
output_h5ad = filtered_obj.copy()
sc.pp.normalize_total(output_h5ad)
```

Normalisation helps reduce the differences between gene and UMI counts by fitting total counts to 10,000 per cell. The inherent log-transform (by log(count+1)) aligns the gene expression level better with a normal distribution. This is fairly standard to prepare for any future dimensionality reduction.

Now we need to look at reducing our gene dimensions. We have loads of genes, but not all of them are different from cell to cell. For instance, housekeeping genes are defined as not changing much from cell to cell, so we could remove these from our data to simplify the dataset. We will flag genes that vary across the cells for future analysis.

```python
output_h5ad = sc.pp.log1p(output_h5ad, copy=True)  # below function requires log scaled data
sc.pp.highly_variable_genes(output_h5ad)
```

Next up, we’re going to scale our data so that all genes have the same variance and a zero mean. This is important to set up our data for further dimensionality reduction. It also helps negate sequencing depth differences between samples, since the gene levels across the cells become comparable. Note, that the differences from scaling etc. are not the values you have at the end - i.e. if your cell has average GAPDH levels, it will not appear as a ‘0’ when you calculate gene differences between clusters.

```python
scaled_data = sc.pp.scale(output_h5ad, max_value=10.0, copy=True)
```

{% icon congratulations %} Congratulations! You have processed your object!

> <comment-title></comment-title>
> At this point, we might want to remove or regress out the effects of unwanted variation on our data. A common example of this is the cell cycle, which can affect which genes are expressed and how much material is present in our cells. If you’re interested in learning how to do this, then you can move over to the [Removing the Effects of the Cell Cycle]({% link topics/single-cell/tutorials/scrna-case_cell-cycle/tutorial.md %}) tutorial now – then return here to complete your analysis.
{: .comment}

# Preparing coordinates

We still have too many dimensions. Transcript changes are not usually singular - which is to say, genes were in pathways and in groups. It would be easier to analyse our data if we could more easily group these changes.

## Principal components

Principal components are calculated from highly dimensional data to find the most spread in the dataset. So in our, ```1982``` highly variable gene dimensions, there will be one line (axis) that yields the most spread and variation across the cells. That will be our first principal component. We can calculate the first ```x``` principal components in our data to drastically reduce the number of dimensions.

> <comment-title>1982???</comment-title>
> Where did the `1982` come from?
>
> The quickest way to figure out how many highly variable genes you have, in my opinion, is to re-run ```sc.pp.highly_variable_genes``` function with the added parameter ```subset=True```, therefore: ```sc.pp.highly_variable_genes(output_h5ad, subset=True)```. This subsetting removes any nonvariable genes.
>
> Then you can ```print(output_h5ad)``` and you'll see only 1982 genes. The following processing steps will use only the highly variable genes for their calculations, but depend on keeping all genes in the object. Thus, please use the original output of your ```sc.pp.highly_variable_genes``` function with far more than 1982 genes!, currently stored as ```scaled_data```.
{: .comment}

> <warning-title>Check your AnnData object!</warning-title>
> Run ```print(scaled_data)```
> Your AnnData object should have far more than 1982 genes in it (if you followed our settings and tool versions, you'd have a matrix 8605 × 15395 (cells x genes). Make sure to use that AnnData object output from FindVariableGenes, rather than the 1982 from your testing in the section above labelled '1982'.
{: .warning}

```python
pca_components = sc.tl.pca(scaled_data, n_comps=50, copy=True)
```

Why 50 principal components you ask? Well, we’re pretty confident 50 is an over-estimate. Let's visualise the variance of each principal component.

```python
sc.pl.pca_variance_ratio(pca_components, n_pcs=50, save='-variance-ratio.png')
```

![Variance ratio](../../images/scrna-casestudy-jupyter/pca-variance.png "Variance ratio")

We can see that there is really not much variation explained past component 19. So we might save ourselves a great deal of time and muddied data by focusing on the top ```20``` PCs.

## Neighborhood graph

We’re still looking at around 20 dimensions at this point. We need to identify how similar a cell is to another cell, across every cell across these dimensions. For this, we will use the k-nearest neighbor (kNN) graph, to identify which cells are close together and which are not. The kNN graph plots connections between cells if their distance (when plotted in this 20 dimensional space!) is amonst the k-th smallest distances from that cell to other cells. This will be crucial for identifying clusters, and is necessary for plotting a UMAP. From <a href="https://github.com/lmcinnes/umap">UMAP developers</a>: “Larger neighbor values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default”.

> <details-title>Working in a group? Decision-time!</details-title>
> If you are working in a group, you can now divide up a decision here with one *control* and the rest varied numbers so that you can compare results throughout the tutorials.
> - Control
>      - **Number of PCs to use** = `20`
>      - **Maximum number of neighbours used** = `15`
> - Everyone else: Use the PC variance plot to pick your own PC number, and choose your own neighbour maximum as well!
{: .details}

```python
neighbours = sc.pp.neighbors(pca_components, n_neighbors=15, use_rep='X_pca', n_pcs=20, copy=True)
```

## Dimensionality reduction for visualisation

Two major visualisations for this data are tSNE and UMAP. We must calculate the coordinates for both prior to visualisation. For tSNE, the parameter [perplexity](https://www.nature.com/articles/s41467-019-13056-x) can be changed to best represent the data, while for UMAP the main change would be to change the kNN graph above itself, by changing the <b>neighbours.</b>

> <details-title>Working in a group? Decision-time!</details-title>
> If you are working in a group, you can now divide up a decision here with one *control* and the rest varied numbers so that you can compare results throughout the tutorials.
> - Control
>      - **Perplexity** = `30`
> - Everyone else: Choose your own perplexity, between 5 and 50!
{: .details}

```python
tsne_components = sc.tl.tsne(neighbours, use_rep='X_pca', perplexity=30, copy=True)
```

```python
umap_components = sc.tl.umap(tsne_components, copy=True)
```

{% icon congratulations %} Congratulations! You have prepared your object and created neighborhood coordinates. We can now use those to call some clusters!

# Cell clusters & gene markers

> <question-title></question-title>
>
> Let's take a step back here. What is it, exactly, that you are trying to get from your data? What do you want to visualise, and what information do you need from your data to gain insight?
>
> > <solution-title></solution-title>
> >
> > Really we need two things - firstly, we need to make sure our experiment was set up well. This is to say, our biological replicates should overlap and our variables should, ideally, show some difference. Secondly, we want insight - we want to know which cell types are in our data, which genes drive those cell types, and in this case, how they might be affected by our biological variable of growth restriction. How does this affect the developing cells, and what genes drive this? So let's add in information about cell clusters and gene markers!
> >
> {: .solution}
>
{: .question}

Finally, let’s identify clusters! Unfortunately, it’s not as majestic as biologists often think - the maths doesn’t necessarily identify true cell clusters. Every algorithm for identifying cell clusters falls short of a biologist knowing their data, knowing what cells should be there, and proving it in the lab. Sigh. So, we’re going to make the best of it as a starting point and see what happens! We will define clusters from the kNN graph, based on how many connections cells have with one another. Roughly, this will depend on a resolution parameter for how granular you want to be.

> <details-title>Working in a group? Decision-time!</details-title>
> Oh yes, yet another decision! Single cell analysis is sadly not straight forward.
> - Control
>      - **Resolution, high value for more and smaller clusters** = `0.6`
> - Everyone else: Pick your own number. If it helps, this sample should have a lot of very similar cells in it. It contains developing T-cells, so you aren't expecting massive differences between cells, like you would in, say, an entire embryo, with all sorts of unrelated cell types.
{: .details}

```python
# Find Clusters
clusters = sc.tl.louvain(umap_components, resolution=0.6, copy=True)
```

Nearly plotting time! But one final piece is to add in SOME gene information. Let’s focus on genes driving the clusters.

# FindMarkers

```python
markers_cluster = sc.tl.rank_genes_groups(clusters, groupby="louvain", method='t-test_overestim_var', n_genes=50, copy=True)
```

But we are also interested in differences across genotype, so let’s also check that (note that in this case, it’s turning it almost into bulk RNA-seq, because you’re comparing all cells of a certain genotype against all cells of the other)

```python
markers_genotype = sc.tl.rank_genes_groups(markers_cluster, groupby="genotype", method='t-test_overestim_var', n_genes=50, copy=True)
```

**Note:** The function ```rank_genes_groups``` does not return a DataFrame that we can use but instead metadata about the marker table, so first we need to construct the marker table using this generated metadata. This is done using the following function, however it's not too important to understand what this code does!

```python
def generate_marker_table(adata):
    # extract marker table metadata
    res = adata.uns['rank_genes_groups']

    # generate DataFrame from metadata
    res_df = pd.DataFrame({
                "genes": pd.DataFrame(res["names"]).stack(),
                "scores": pd.DataFrame(res["scores"]).stack(),
                "logfoldchanges": pd.DataFrame(res["logfoldchanges"]).stack(),
                "pvals": pd.DataFrame(res["pvals"]).stack(),
                "pvals_adj": pd.DataFrame(res["pvals_adj"]).stack(),
            })

    # convert row names to columns
    res_df.index.name = 'newhead'
    res_df.reset_index(inplace=True)

    # rename generic column names
    res_df = res_df.rename(columns={'level_0': 'rank', 'level_1':'cluster'})

    # reorder columns
    res_df = res_df.reindex(columns=['cluster', 'rank', 'genes', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'])

    # insert ref column
    res_df.insert(2, 'ref', 'rest')

    return res_df
```

Now we can generate our marker tables!

```python
# Generate marker tables
cluster_marker_table = generate_marker_table(markers_cluster)
genotype_marker_table = generate_marker_table(markers_genotype)

display(cluster_marker_table.head(4))
display(genotype_marker_table.head(4))
```

Now, there’s a small problem here, which is that if you inspect the output marker tables (above tables), you won’t see gene names, you’ll see Ensembl IDs. While this is a more bioinformatically accurate way of doing this (not every ID has a gene name!), we might want to look at more well-recognised gene names, so let’s pop some of that information in!

```python
# Join two datasets

cluster_joined = pd.merge(cluster_marker_table, markers_cluster.var, left_on='genes', right_on='ID')
genotype_joined = pd.merge(genotype_marker_table, markers_genotype.var, left_on='genes', right_on='ID')

display(cluster_joined.head(5))
display(genotype_joined.head(5))
```

```python
# Cut columns from tables

cluster_markers_named = cluster_joined[['cluster', 'ref', 'rank', 'genes', 'Symbol', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']]
genotype_markers_named = genotype_joined[['cluster', 'ref', 'rank', 'genes', 'Symbol', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']]

display(cluster_markers_named.head(5))
display(genotype_markers_named.head(5))
```

Well done! It’s time for the best bit, the plotting!

# Plotting!

It’s time! Let’s plot it all! But first, let’s pick some marker genes from the ```markers_cluster``` list that you made as well. I’ll be honest, in practice, you’d now be spending a lot of time looking up what each gene does (thank you google!). There are burgeoning automated-annotation tools, however, so long as you have a good reference (a well annotated dataset that you’ll use as the ideal). In the mean time, let’s do this the old-fashioned way, and just copy a bunch of the markers in the original paper.

```python
# PCA
sc.pl.embedding(
    markers_cluster,
    basis='pca',  
    color=['louvain','sex','batch','genotype','Il2ra','Cd8b1','Cd8a','Cd4','Itm2a','Aif1','log1p_total_counts'],
    gene_symbols='Symbol',
    use_raw=False,
    save='.png'
)
```

```python
# TSNE
sc.pl.embedding(
    markers_cluster,
    basis='tsne',  
    color=['louvain','sex','batch','genotype','Il2ra','Cd8b1','Cd8a','Cd4','Itm2a','Aif1','log1p_total_counts'],
    gene_symbols='Symbol',
    use_raw=False,
    save='.png'
)
```

```python
# UMAP
sc.pl.embedding(
    markers_cluster,
    basis='umap',  
    color=['louvain','sex','batch','genotype','Il2ra','Cd8b1','Cd8a','Cd4','Itm2a','Aif1','log1p_total_counts'],
    gene_symbols='Symbol',
    use_raw=False,
    save='.png'
)
```

{% icon congratulations %} Congratulations! You now have plots galore!

# Insights into the beyond

Now it’s the fun bit! We can see where genes are expressed, and start considering and interpreting the biology of it. At this point, it’s really about what information you want to get from your data - the following is only the tip of the iceberg. However, a brief exploration is good, because it may help give you ideas going forward with for your own data. Let us start interrogating our data!

## Biological Interpretation

> <question-title>Appearance is everything</question-title>
>
> Which visualisation is the most useful for getting an overview of our data, *pca*, *tsne*, or *umap*?
> ![PCA-tSNE-UMAP](../../images/scrna-casestudy-jupyter/louvain_clustering.png "Louvain clustering by dimension reduction")
>
> > <solution-title></solution-title>
> >
> > You can see why a PCA is generally not enough to see clusters in samples - keep in mind, you're only seeing components 1 and 2! - and therefore why the tSNE and UMAP visualisation dimensionality reductions are so useful. But there is not necessarily a clear winner between tSNE and UMAP, but I think UMAP is slightly clearer with its clusters, so we'll stick with that for the rest of the analysis.
> >
> {: .solution}
>
{: .question}

Note that the cluster numbering is based on size alone - clusters 0 and 1 are not necessarily related, they are just the clusters containing the most cells. It would be nice to know what exactly these cells are. This analysis (googling all of the marker genes, both checking where the ones you know are as well as going through the marker tables you generated!) is a fun task for any individual experiment, so we’re going to speed past that and nab the assessment from the original paper!

| Clusters   | Marker                  | Cell type                            |
| ---------- | ----------------------- | ------------------------------------ |
| 4          | Il2ra                   | Double negative (early T-cell)       |
| 0,1,2,6    | Cd8b1, Cd8a, Cd4        | Double positive (middle T-cell)      |
| 5          | Cd8b1, Cd8a, Cd4 - high | Double positive (late middle T-cell) |
| 3          | Itm2a                   | Mature T-cell                        |
| 7          | Aif1                    | Macrophages                          |


The authors weren’t interested in further annotation of the DP cells, so neither are we. Sometimes that just happens. The maths tries to call similar (ish) sized clusters, whether it is biologically relevant or not. Or, the question being asked doesn’t really require such granularity of clusters.

![Known marker gene locations](../../images/scrna-casestudy-jupyter/known_marker_gene_locations.png "Known marker genes locations")

> <details-title>Working in a group? Important!</details-title>
> If you have deviated from any of the original parameters in this tutorial, you will likely have a different number of clusters. You will, therefore, need to change the ```categories``` parameter in ```rename_categories``` accordingly. Best of luck!
>
{: .details}

## Annotating Clusters

```python
# Add meaningful names to each category
markers_cluster.rename_categories(key='louvain', categories=['DP-M4','DP-M3','T-mat','DN','DP-M2','DP-M1','DP-L','Macrophages'])

# Copy AnnData object
markers_cluster_copy = markers_cluster.copy()

# Rename 'louvain' column
markers_cluster_copy.obs = markers_cluster_copy.obs.rename(columns={'louvain': 'cell_type'})

# Scanpy - plot updated object
sc.pl.embedding(
    markers_cluster_copy,
    basis='umap',  
    color=['cell_type','sex','batch','genotype','Il2ra','Cd8b1','Cd8a','Cd4','Itm2a','Aif1','Hba-a1','log1p_total_counts'],
    gene_symbols='Symbol',
    use_raw=False,
    save='-annotated.png'
)
```

![Annotated Embeddings](../../images/scrna-casestudy-jupyter/annotated_umap.png "Our Annotated UMAP")

Now that we know what we’re dealing with, let’s examine the effect of our variable, proper science!

> <question-title>Genotype</question-title>
>
> Are there any differences in genotype? Or in biological terms, is there an impact of growth restriction on T-cell development in the thymus?
>
> ![Genotype Images](../../images/scrna-casestudy-jupyter/genotype_differences.png "Genotype differences")
>
> > <solution-title></solution-title>
> >
> > We can see that DP-L, which seems to be extending away from the DP-M bunch, as well as the mature T-cells (or particularly the top half) are missing some knockout cells. Perhaps there is some sort of inhibition here? INTERESTING! What next? We might look further at the transcripts present in both those populations, and perhaps also look at the genotype marker table... So much to investigate! But before we set you off to explore to your heart's delight, let's also look at this a bit more technically.
> >
> {: .solution}
>
{: .question}

## Technical Assessment

Is our analysis real? Is it right? Well, we can assess that a little bit.

> <question-title>Batch effect</question-title>
>
> Is there a batch effect?
>
> ![Batch effect](../../images/scrna-casestudy-jupyter/batch_effect.png "Batch effect?")
>
> > <solution-title></solution-title>
> >
> > While some shifts are expected and nothing to be concerned about, DP-L looks to be mainly comprised of N705. There might be a bit of batch effect, so you could consider using batch correction on this dataset. However, if we focus our attention on the other cluster - mature T-cells -  where there is batch mixing, we can still assess this biologically even without batch correction.
> > Additionally, we will also look at the confounding effect of sex.
> >
> > ![Sex effect](../../images/scrna-casestudy-jupyter/sex_differences.png "Sex differences")
> >
> > We note that the one female sample - unfortunately one of the mere three knockout samples - seems to be distributed in the same areas as the knockout samples at large, so luckily, this doesn't seem to be a confounding factor and we can still learn from our data. Ideally, this experiment would be re-run with either more female samples all around or swapping out this female from the male sample.
> >
> {: .solution}
>
{: .question}

> <question-title>Depth effect</question-title>
>
> Are there any clusters or differences being driven by sequencing depth, a technical and random factor?
>
> ![Sequencing depth](../../images/scrna-casestudy-jupyter/counts_across_clusters.png "Counts across clusters")
>
> > <solution-title></solution-title>
> >
> > Eureka! This explains the odd DP shift between wildtype and knockout cells - the left side of the DP cells simply have a higher sequencing depth (UMIs/cell) than the ones on the right side. Well, that explains some of the sub-cluster that we're seeing in that splurge. Importantly, we don't see that the DP-L or (mostly) the mature T-cell clusters are similarly affected. So, whilst again, this variable of sequencing depth might be something to regress out somehow, it doesn't seem to be impacting our dataset. The less you can regress/modify your data, in general, the better - you want to stay as true as you can to the raw data, and only use maths to correct your data when you really need to (and not to create insights where there are none!).
> >
> {: .solution}
>
{: .question}

> <question-title>Sample purity</question-title>
>
> Do you think we processed these samples well enough?
>
> ![Sequencing depth](../../images/scrna-casestudy-jupyter/hemoglobin.png "Hemoglobin across clusters")
>
> > <solution-title></solution-title>
> >
> > We have seen in the previous images that these clusters are not very tight or distinct, so we could consider stronger filtering. Additionally, hemoglobin - a red blood cell marker that should NOT be found in T-cells - appears throughout the entire sample in low numbers. This suggests some background in the media the cells were in, and we might consider in the wet lab trying to get a purer, happier sample, or in the dry lab, techniques such as SoupX or others to remove this background. Playing with filtering settings (increasing minimum counts/cell, etc.) is often the place to start in these scenarios.
> >
> {: .solution}
>
{: .question}

> <question-title>Clustering resolution</question-title>
>
> Do you think the clustering is appropriate? i.e. are there single clusters that you think should be separate, and multiple clusters that could be combined?
>
> ![Itm2a Expression](../../images/scrna-casestudy-jupyter/itm2a_across_clusters.png "Itm2a across clusters")
>
> > <solution-title></solution-title>
> >
> > Important to note, lest all bioinformaticians combine forces to attack the biologists: just because a cluster doesn't look like a cluster by eye is NOT enough to say it's not a cluster! But looking at the biology here, we struggled to find marker genes to distinguish the DP population, which we know is also affected by depth of sequencing. That's a reasonable argument that DP-M1, DP-M2, and DP-M3 might not be all that different. Maybe we need more depth of sequencing across all the DP cells, or to compare these explicitly to each other (consider variations on FindMarkers!). However, DP-L is both seemingly leaving the DP cluster and also has fewer knockout cells, so we might go and look at what DP-L is expressing in the marker genes. If we look at T-mat further, we can see that its marker gene - Itm2a - is only expressed in half of the cluster. You might consider sub-clustering this to investigate further, either through changing the resolution or through analysing this cluster alone.
> > If we look at the differences between genotypes alone (so the pseudo-bulk), we can see that most of the genes in that list are actually ribosomal. This might be a housekeeping background, this might be cell cycle related, this might be biological, or all three. You might consider investigating the cycling status of the cells, or even regressing this out (which is what the authors did).
> {: .solution}
>
{: .question}

Ultimately, there are quite a lot ways to analyse the data, both within the confines of this tutorial (the many parameters that could be changed throughout) and outside of it (batch correction, sub-clustering, cell-cycle scoring, inferred trajectories, etc.) Most analyses will still yield the same general output, though: there are fewer knockout cells in the mature T-cell population.

{% icon congratulations %} Congratulations! You have interpreted your plots in several important ways!

# Export your data, and notebook, and figures

It’s now time to export your data! First, we need to get Jupyter to see it as a file.

```adata.write('markers_cluster_copy')```

Now you can export it.

```put("MarkersCluster.h5ad")```

To export your notebook to your Galaxy history, you can use the following. Change the text to be your notebook name. Do not use spaces!

```put("yourtitlehere.ipynb")```    

Want to export some plots? Choose any (or all) of the plots you saved as files in the folder at the left and put their titles in the following. You can run multiple exports at the same time.

```python
put("figures/plotname.png")
put("figures/plotname.png")
put("figures/plotname.png")
```

# Conclusion

{% icon congratulations %} Congratulations! You’ve made it to the end!

In this tutorial, you moved from technical processing to biological exploration. By analysing real data - both the exciting and the messy! - you have, hopefully, experienced what it’s like to analyse and question a dataset, potentially without clear cut-offs or clear answers. If you were working in a group, you each analysed the data in different ways, and most likely found similar insights. One of the biggest problems in analysing scRNA-seq is the lack of a clearly defined pathway or parameters. You have to make the best call you can as you move through your analysis, and ultimately, when in doubt, try it multiple ways and see what happens!
