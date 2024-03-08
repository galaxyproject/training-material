---
layout: learning-pathway
tags: [intermediate]
cover-image: assets/images/wab-annotatedcells-2.png
cover-image-alt: "Image of cells in different coloured clusters"
type: use

editorial_board:
- nomadscientist
- pavanvidem
- pcm32

title: Applying single-cell RNA-seq analysis
description: |
  Gone is the pre-annotated, high quality tutorial data - now you have real, messy data to deal with. You have decisions to make and parameters to decide. This learning pathway challenges you to replicate a published analysis as if this were your own dataset. You will be introduced to a few more tools available for scRNA-seq in Galaxy. Finally, if our tool offerings are not enough for you, you will be directed towards how to use coding notebooks within Galaxy, setting you up to analyse scRNA-seq in R or python notebooks.

  The data is messy. The decisions are tough. The interpretation is meaningful. Come here to advance your single cell skills! Note that you get two options for inferring trajectories.

  For support throughout these tutorials, join our Galaxy [single cell chat group on Matrix](https://matrix.to/#/#Galaxy-Training-Network_galaxy-single-cell:gitter.im) to ask questions!

pathway:
  - section: "Module 1: Case study"
    description: |
      These tutorials take you from raw scRNA sequencing reads to cell cluster plots to replicate a published analysis.
    tutorials:
      - name: scrna-case_alevin
        topic: single-cell
      - name: scrna-case_alevin-combine-datasets
        topic: single-cell
      - name: scrna-case_basic-pipeline
        topic: single-cell

  - section: "Module 2: Inferring trajectories"
    description: |
      This isn't strictly necessary, but if you want to infer trajectories - pseudotime relationships between cells - you can try out these tutorials with the same dataset.  Note that you get two options for inferring trajectories, you can choose either.
    tutorials:
      - name: scrna-case_trajectories
        topic: single-cell
      - name: scrna-case_monocle3-trajectories
        topic: single-cell

  - section: "Module 3: Moving into coding environments"
    description: |
      Did you know Galaxy can host coding environments? They don't have the same level of computational power as the easy-to-use Galaxy tools, but you can unlock the full freedom in your data analysis. You can install your favourite single-cell tool suite that is not available on Galaxy, export your data into these coding environments and run your analysis there. If you want your favourite tool suite as a Galaxy tool, you can always request [here](https://docs.google.com/spreadsheets/d/15hqgqA-RMDhXR-ylKhRF-Dab9Ij2arYSKiEVoPl2df4/edit?usp=sharing). Let's start with the basics of running these environments in Galaxy.
    tutorials:
      - name: jupyterlab
        topic: galaxy-interface
      - name: galaxy-intro-jupyter
        topic: galaxy-interface
      - name: rstudio
        topic: galaxy-interface

  - section: "The End!"
    description: |
      And now you're done! If you are interested in trying out the case study analyses in a coding environment, try out our ["Case study: Reloaded" series](/training-material/topics/single-cell#st-single-cell-cs-code) next! Otherwise, you will find more features, tips and tricks in our general [Galaxy Single-cell Training page](/training-material/topics/single-cell/index.html).
---

New to Galaxy and/or the field of scRNA-seq? Follow this learning path to get familiar with the basics!
