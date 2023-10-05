---
layout: learning-pathway
tags: [beginner]
cover-image: assets/images/wab-annotatedcells-2.png
cover-image-alt: "Image of cells in different coloured clusters"
type: use

editorial_board:
- nomadscientist
- pavanvidem

title: Introduction to Galaxy and Single Cell RNA Sequence analysis
description: |
  This learning path aims to teach you the basics of Galaxy and analysis of Single Cell RNA-seq data.
  You will learn how to use Galaxy for analysis, and an important Galaxy feature for iterative single cell analysis. You'll tbe guided through the general theory of single analysis and then perform a basic analysis of 10X chromium data. For support throughout these tutorials, join our Galaxy [single cell chat group on Matrix](https://matrix.to/#/#Galaxy-Training-Network_galaxy-single-cell:gitter.im) to ask questions!

priority: 9

pathway:
  - section: "Module 1: Introduction to Galaxy"
    description: |
      Get a first look at the Galaxy platform for data analysis. We start with a
      short introduction (video slides & practical) to familiarize you with the Galaxy
      interface, and then proceed with a short tutorial of how to tag - and organise! - your history.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
      - name: name-tags
        topic: galaxy-interface

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
