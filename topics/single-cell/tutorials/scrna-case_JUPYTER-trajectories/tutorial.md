---
layout: tutorial_hands_on

title: Inferring single cell trajectories (Scanpy, Python)
subtopic: single-cell-CS-code
priority: 3
zenodo_link: 'https://zenodo.org/record/7075718'
redirect_from:
- /topics/transcriptomics/tutorials/scrna-JUPYTER-trajectories/tutorial
- /topics/transcriptomics/tutorials/scrna-case_JUPYTER-trajectories/tutorial
questions:
- How can I infer lineage relationships between single cells based on their RNA, without a time series?
objectives:
- Execute multiple plotting methods designed to maintain lineage relationships between cells
- Interpret these plots
time_estimation: 2H
key_points:
- Trajectory analysis is less robust than pure plotting methods, as such 'inferred relationships' are a bigger mathematical leap
- As always with single-cell analysis, you must know enough biology to deduce if your analysis is reasonable, before exploring or deducing novel insight
requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin
        - scrna-case_alevin-combine-datasets
        - scrna-case_basic-pipeline
-
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
        - galaxy-intro-jupyter
tags:
- 10x
- paper-replication
- Python

contributions:
  authorship:
    - nomadscientist
    - wee-snufkin
    - mtekman
  editing:
    - hexylena

notebook:
  language: python
  snippet: topics/single-cell/tutorials/scrna-case_JUPYTER-trajectories/preamble.md

---

# Run the tutorial!

From now on, you can view this tutorial in the Jupyter notebook, which will allow you to read the material and simultaneously execute the code cells! You may have to change certain numbers in the code blocks, so do read carefully. The tutorial is adapted from the [Scanpy Trajectory inference tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html).

{% snippet topics/single-cell/faqs/notebook_warning.md %}

## Install modules & activate them

```python
pip install scanpy
```
```python
pip install fa2
```
```python
pip install igraph
```
```python
pip install louvain
```
```python

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
```

## Import dataset

You can now import files from your Galaxy history directly using the following code. This will depend on what number in your history the final annotated object is. If your object is dataset #2 in your history, then you import it as following:

```python
thymusobject = get(2)
```

You now you need to read it in as an h5ad object.

```python
adata = sc.read_h5ad(thymusobject)
```

## Draw force-directed graph

First, we will calculate a [force-directed graph](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.draw_graph.html), as an alternate to tSNE, which will likely work better for trajectory analysis.

```python
sc.tl.draw_graph(adata)
```

And now time to plot it!
*Note: We're saving to **png**, but you can also choose pdf*


```python
sc.pl.draw_graph(adata, color='cell_type', legend_loc='on data', save = 'Plot1.png')
```

![Plot1-Force-Directed Graph](../../images/scrna-casestudy/draw_graph_faPlot1.png "Plot1-Force-Directed Graph")

Well now this is exciting! Our DP-late is more clearly separating, and we might also suppose that DP-M1, DP-M2, and DP-M3 are actually earlier on in the differentiation towards mature T-cells. And we're only just getting started!


## Diffusion maps

We'll now perform an *optional step*, that basically takes the place of the PCA. Instead of using PCs, we can use [diffusion maps](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.diffmap.html).


```python
sc.tl.diffmap(adata)
```

Now that we have our diffusion map, we need to re-calculate neighbors using the diffusion map instead of the PCs. Then we re-draw and plot a new force directed graph using the new neighbors.

```python
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_diffmap')
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='cell_type', legend_loc='on data', save = 'Plot2.png')
```
![Diffusion Map](../../images/scrna-casestudy/draw_graph_faPlot2.png "Diffusion Map")


Oh dear! This doesn't look great. Maybe the DP-M4 cells are a whole other trajectory? That doesn't seem right. Saying that, this spreads out our T-mature cells, which makes a lot more sense when it comes to T-cell biology (we expect T-cells to differentiate into two types of T-cells, Cd8+Cd4- and Cd4+Cd8-). If you wanted to, you could also re-cluster your cells (since you've changed the neighborhood graph on which the clusterisation depends). You could use this:
`sc.tl.louvain(adata, resolution=0.6)`
However, we tried that, and it called far too many clusters given the depth of sequencing in this dataset. Let's stick with our known cell types and move from there.


## Working in a group? Decision-time!
If you are working in a group, you can now divide up a decision here with one *control* and the rest can vary numbers so that you can compare results throughout the tutorials.
- Control
   - Go straight to the PAGA section
- Everyone else:
   - you could re-call clusters `sc.tl.louvain(adata, resolution=0.6` or use other resolutions! (Tip, go low!)
        - Please note that in this case, you will want to change the PAGA step `sc.pl.paga` to group by `louvain` rather than `cell_type`. You can certainly still plot both, we only didn't because with using our old Louvain calls, the cell_type and louvain categories are identical.
   - you could undo the diffusion map step by running the following
        `sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')`
        `sc.tl.draw_graph(adata)`
   - you could also change the number of neighbors used in the pp.neighbors step (this is the same as the Galaxy tool **Scanpy ComputeGraph**

- Everyone else: You will want to compare FREQUENTLY with your control team member.

## PAGA

[PAGA](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.paga.html) is used to generalise relationships between groups, or likely clusters, in this case.


```python
sc.tl.paga(adata, groups='cell_type')
```

Now we want to plot our PAGA, but we might also be interested in colouring our plot by genes as well. In this case, remembering that we are dutifully counting our genes by their EnsemblIDs rather than Symbols (which do not exist for all EnsemblIDs), we have to look up our gene of interest (CD4, CD8a) and plot the corresponding IDs.


```python
sc.pl.paga(adata, color=['cell_type', 'ENSMUSG00000023274', 'ENSMUSG00000053977'], title=['Cell type', 'CD4', 'Cd8a'], save = 'Plot4.png')
```

![PAGA](../../images/scrna-casestudy/pagaPlot4.png "PAGA")


Well now that is interesting! This analysis would find that DP-M1 and DP-M4 are both driving towards differentiation, which is not something we had necessarily been able to specify before by just looking at our cluster graphs or applying our biological knowledge.

## Re-draw force-directed graph

Force directed graphs can be initialised randomly, or we can prod it in the right direction. We'll prod it with our PAGA calculations. Note that you could also try prodding it with tSNE or UMAP. A lot of these tools can be used on top of each other or with each other in different ways, this tutorial is just one example. Similarly, you could be using any **obs** information for grouping, so could do this for *louvain* or *cell_type* for instance.

```python
sc.tl.draw_graph(adata, init_pos='paga')

sc.pl.draw_graph(adata, color=['cell_type'], title=['Cluster'], legend_loc='on data', save = 'Plot5.png')
sc.pl.draw_graph(adata, color=['genotype'], title=['Genotype'], save = 'Plot6.png')
sc.pl.draw_graph(adata, color=['ENSMUSG00000023274', 'ENSMUSG00000053977'], title=['CD4', 'Cd8a'], save = 'Plot7.png')
```

![Force-Directed + PAGA - Cell type](../../images/scrna-casestudy/draw_graph_faPlot5.png "Force-Directed + PAGA - Cell type")

![Force-Directed + PAGA - Genotype](../../images/scrna-casestudy/draw_graph_faPlot6.png "Force-Directed + PAGA - Genotype")

![Force-Directed + PAGA - Markers](../../images/scrna-casestudy/draw_graph_faPlot7.png "Force-Directed + PAGA - Markers")

**Note** - we are aware that something about these graphs has gotten a bit odd in the recent Scanpy updates. Watch this space for a fix!

Well aren't those charts interesting! Using the diffusion map to drive the force-directed graph, we see correct ordering of our cells (from DN to DP to T-mature, which was lost with the diffusion map alone) as well as two apparent branches leaving the mature T-cell population, which is what we'd biologically expect. In terms of our experiment, we're seeing a clear trajectory issue whereby the knockout cells are not found along the trajectory into T-mature (which, well, we kind of already figured out with just the cluster analysis, but we can feel even more confident about our results!) More importantly, we can see the T-mature population dividing itself, which we did not see in the clustering via UMAP/tSNE alone, and we can verify that as the leftmost branch has CD4 but the rightmost branch does not. This is suggesting our branchpoint from to CD4+ and CD8+ single positive cells. Exciting! However, it is important to note, that the branches there are quite small and sparsely populated, which can indicate artifact branches (i.e. trajectory analysis does its best to find branches, particularly diffusion map, so you can pretty easily force branches to appear even if they are not biologically real!). However, to be frank, we were surprised not to find this clearer in the main cluster map, as we know that the T-cells should diverge at that point, so if anything this is a relief that our data is believable!

And now, just for fun, we can compare the scatter graph with our PAGA side by side.


```python
sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)
```

![PAGA Compare](../../images/scrna-casestudy/paga_compare.png "PAGA Compare")

**Note** - we are aware that something about these graphs has gotten a bit odd in the recent Scanpy updates. Watch this space for a fix!

## Diffusion pseudotime

We know that our cells are initialising at DN. We can use feed that information into our algorithms to then calculate a trajectory.

First, let's name our 'root'.

```python
adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type']  == 'DN')[0]
```

## Working in a group? Decision-time!
If you called new clusters using the louvain algorithm, you might want to choose one of those clusters to be your root cell instead, so change the `cell_type` above for `louvain` and then name the cluster number. Use the plots you created to help you pick the number!


Onto the [diffusion pseudotime](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dpt.html), where we are infer multiple time points within the same piece of data!

```python
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color=['cell_type', 'dpt_pseudotime'], legend_loc='on data', save = 'Plot8.png')
```

![Force-Directed + Pseudotime](../../images/scrna-casestudy/draw_graph_faPlot8.png "Force-Directed + Pseudotime")

This is nice, as it supports our conclusions thus far on the trajectory of the T-cell differentiation. With single-cell, the more ways you can prove to yourself what you're seeing is real, the better! If we did not find consistent results, we would need to delve in further to see if the algorithm (not all algorithms fit all data!) or the biology.

Where might we go from here? We might consider playing with our louvain resolutions, to get the two branches to be called as different clusters, and then comparing them to each other for gene differences or genotype differences. We might also use different objects (for instance, what if we regressed out cell cycle genes?) and see if that changes the results. Perhaps we would eliminate the DN double-branch input. Or perhaps that's real, and we should investigate that. What would you do?



## Working in a group? The finale!
Look at each others images! How do yours differ, what decisions were made? Previously, when calling clusters in the 'Filter, Plot and Explore Single-cell RNA-seq Data', the interpretation at the end is largely consistent, no matter what decisions are made throughout (mostly!). Is this the case with your trajectory analyses? You may find that it is not, which is why pseudotime analysis even more crucially depends on your understanding of the underlying biology (we have to choose the root cells, for instance, or recognise that DN cells should not be found in the middle of the DPs) as well as choosing the right analysis. That's why it is a huge field! With analysing scRNA-seq data, it's almost like you need to know about 75% of your data and make sure your analysis shows that, for you to then identify the 25% new information.

# Export your data, figures, and notebook

It's now time to export your data! First, we need to get it Jupyter to see it as a file.


```python
adata.write('Trajectorythymus.h5ad')
```

Now you can export it, as well as all your lovely plots! If you go into the *figures* folder at the left, you'll see your lovely plots and can choose which ones to export. The following code will push them into your galaxy history. You can also directly download them onto your computer from the file window at the left.

```python
put("Trajectorythymus.h5ad")
```
```python
put("figures/draw_graph_faPlot1.png")
put("figures/draw_graph_faPlot2.png")
put("figures/draw_graph_faPlot5.png")
put("figures/draw_graph_faPlot6.png")
put("figures/draw_graph_faPlot7.png")
put("figures/draw_graph_faPlot8.png")
put("figures/paga_compare.pdf")
put("figures/pagaPlot4.png")
```
The cell below will only work if you haven't changed the name of the notebook. If you renamed it, simply type its new name in the parenthesis.
```python
put("single-cell-scrna-case_JUPYTER-trajectories.ipynb")       
```

This may take a moment, so go check your Galaxy history to make sure your images, anndata object, and notebook (.ipynb) have all made it back into your Galaxy history. Once they are all there, you can exit this browser and return to the Galaxy tutorial!

If things have gone wrong, you can also download this [answer key tutorial]({{ page.zenodo_link }}/files/Trajectories_Answer_Key.ipynb).

# Citation

Please note, this is largely based on the trajectories tutorial found on the Scanpy site itself [https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html](https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html).



# After Jupyter

Congratulations! You've made it through Jupyter!

> <hands-on-title>Closing JupyterLab</hands-on-title>
> 1. Click **User**: **Active Interactive Tools**
> 2. Tick the box of your Jupyter Interactive Tool, and click **Stop**
{: .hands_on}

If you want to run this notebook again, or share it with others, it now exists in your history. You can use this 'finished' version just the same way as you downloaded the directions file and uploaded it into the Jupyter environment.


# Conclusion

Congratulations! You've made it to the end! You might be interested in the [Answer Key History](https://usegalaxy.eu/u/wendi.bacon.training/h/cs4inferring-trajectories-using-python-in-galaxyanswer-key) or the [Answer Key Jupyter Notebook](https://zenodo.org/record/7054806/files/Trajectories_Answer_Key.ipynb?download=1).

In this tutorial, you moved from called clusters to inferred relationships and trajectories using pseudotime analysis. You found an alternative to PCA (diffusion map), an alternative to tSNE (force-directed graph), a means of identifying cluster relationships (PAGA), and a metric for pseudotime (diffusion pseudotime) to identify early and late cells. If you were working in a group, you found that such analysis is slightly more sensitive to your decisions than the simpler filtering/plotting/clustering is. We are inferring and assuming relationships and time, so that makes sense!

To discuss with like-minded scientists, join our [Gitter](https://gitter.im/Galaxy-Training-Network/galaxy-single-cell?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) channel for all things Galaxy-single cell!
