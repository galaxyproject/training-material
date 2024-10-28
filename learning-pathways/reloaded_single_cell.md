---
layout: learning-pathway
tags: [advanced]
cover-image: assets/images/wab-annotatedcells-2.png
cover-image-alt: "Image of cells in different coloured clusters"
type: use

editorial_board:
- nomadscientist
- pavanvidem

title: Applying single-cell RNA-seq analysis in Coding Environments
description: |
  Gone is the pre-annotated, high quality tutorial data - now you have real, messy data to deal with. You have decisions to make and parameters to decide. This learning pathway challenges you to replicate a published analysis as if this were your own dataset. You will perform this analysis in coding environments hosted on Galaxy, instead of Galaxy's button-based tool interface.

  The data is messy. The decisions are tough. The interpretation is meaningful. Come here to advance your single cell skills! Note that you get two options: performing the analysis predominantly in R or in Python.

  For support throughout these tutorials, join our Galaxy [single cell chat group on Matrix](https://matrix.to/#/#Galaxy-Training-Network_galaxy-single-cell:gitter.im) to ask questions!



pathway:
  - section: "Module 1: Coding environments in Galaxy"
    description: |
      Let's start with the basics of running a coding environments in Galaxy.
    tutorials:
      - name: jupyterlab
        topic: galaxy-interface
      - name: galaxy-intro-jupyter
        topic: galaxy-interface
      - name: rstudio
        topic: galaxy-interface

  - section: "Module 2: Preparing the dataset"
    description: |
      These tutorials take you from raw scRNA sequencing reads to a matrix ready for downstream analysis. Galaxy coding environments don't have the same level of computational power as the easy-to-use Galaxy tools, so in practice, dataset preparation is usually performed in the Galaxy user interface to process the dataset into something smaller, which can then be analysed in the coding environment. Nevertheless, the whole process can be performed in a coding environment.
    tutorials:
      - name: alevin-commandline
        topic: single-cell

  - section: "Module 2: Generating cluster plots"
    description: |
              These tutorials take you from the pre-processed matrix to cluster plots and gene expression values. You can pick whether to follow the Python (Scanpy) or R (Seurat) tutorial.
    tutorials:
      - name: scrna-case_jupyter_basic-pipeline
        topic: single-cell
      - name: scrna-case_FilterPlotandExplore_RStudio
        topic: single-cell

  - section: "Module 3: Inferring trajectories"
    description: |
              This isn't strictly necessary, but if you want to infer trajectories - pseudotime relationships between cells - you can try out these tutorials with the same dataset. Again, you can choose whether to follow the Python (Scanpy) or R (Monocle) tutorial.
    tutorials:
     - name: scrna-case_JUPYTER-trajectories
       topic: single-cell
     - name: scrna-case_monocle3-rstudio
       topic: single-cell
 
  - section: "The End!"
    description: |
      And now you're done! You will find more features, tips and tricks in our general [Galaxy Single-cell Training page](/training-material/topics/single-cell/index.html).
---

Want to try scRNA-seq analysis in a coding environment? Follow this learning path!
