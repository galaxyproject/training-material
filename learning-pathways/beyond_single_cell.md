---
layout: learning-pathway
tags: [intermediate]
cover-image: assets/images/wab-annotatedcells-2.png
cover-image-alt: "Image of cells in different coloured clusters"
type: use

editorial_board:
- nomadscientist
- pavanvidem

title: Applying single-cell RNA-seq analysis
description: |
  Gone is the pre-annotated, high quality tutorial data - now you have real, messy data to deal with. You have decisions to make and parameters to decide. This learning pathway challenges you to replicate a published analysis as if this were your own dataset. You will be introduced to a few more tools available for scRNA-seq in Galaxy. Finally, if our tool offerings are not enough for you, you will be directed towards how to use coding notebooks within Galaxy, setting you up to analyse scRNA-seq in R or python notebooks.

  The data is messy. The decisions are tough. The interpretation is meaningful. Come here to advance your single cell skills! Note that you get two options for inferring trajectories.

  For support throughout these tutorials, join our Galaxy [single cell chat group on Matrix](https://matrix.to/#/#Galaxy-Training-Network_galaxy-single-cell:gitter.im) to ask questions!

pathway:
  - section: "Module 1: Case study"
    description: |
      These tutorials take you from raw scRNA sequencing reads to inferred trajectories to replicate a published analysis. Note that you get two options for inferring trajectories.
    tutorials:
      - name: scrna-case_alevin
        topic: single-cell
      - name: name-tags
        topic: single-cell

  - section: "Module 2: Theory of Single-Cell RNA-seq"
    description: |
      When analysing sequencing data, you should always start with a quality control step to clean your data and make sure your data is good enough to answer your research question. After this step, you will often proceed with a mapping (alignment) or genome assembly step, depending on whether you have a reference genome to work with.
    tutorials:
      - name: scrna-intro
        topic: single-cell

  - section: "Module 3: Time to analyse data!"
    description: |
      It's time to apply your skills! You'll now analyse some clean data from the 10X Chromium platform.
    tutorials:
      - name: scrna-preprocessing-tenx
        topic: single-cell
      - name: scrna-scanpy-pbmc3k
        topic: single-cell

  - section: "The End!"
    description: |
      And now you're done! There are still loads of resources to take you from basic analysis to more difficult decision-making, deconvolution, multiomics, or ingesting from different data sources. See the [Galaxy Single Cell Training page](/training-material/topics/single-cell/index.html) for more!
---

New to Galaxy and/or the field of scRNA-seq? Follow this learning path to get familiar with the basics!
