---
layout: learning-pathway
tags: [beginner, genomic-intervals, epigenetics]
cover-image: assets/images/chipseq-viz.png
cover-image-alt: "Plot of ChIP-Seq signal from multiple samples across a genome region"
type: use

editorial_board:
- wm75
- pavanvidem

title: Epigenetics data analysis with Galaxy
description: |
  Familiarize yourself with the Galaxy user interface and its most important features,
  then move on to analysis of genomic interval and deep-sequencing data
  with tools and tool suites like MACS2, bedtools and deepTools.

priority: 9

pathway:
  - section: "Module 1: Introduction to Galaxy"
    description: |
      Get a first look at the Galaxy platform for data analysis. We start with a
      short introduction (video slides & practical) to familiarize you with the Galaxy
      interface, then proceed with slightly longer, but still introductory tutorials that
      demonstrate the handling of genomic intervals data in Galaxy.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
      - name: galaxy-intro-101
        topic: introduction
      - name: galaxy-intro-peaks2genes
        topic: introduction

  - section: "Module 2: Analysis of ChIP-Seq, CUT&RUN and ATAC-Seq data"
    description: |
      Learn how to perform mapping, peak calling with **MACS-2**,
      read depth analysis and correlation with the **deepTools** suite,
      and results vizualization with **IGV** and other tools
      for various chromatin selection assays.

    tutorials:
      - name: formation_of_super-structures_on_xi
        topic: epigenetics
      - name: cut_and_run
        topic: epigenetics
      - name: atac-seq
        topic: epigenetics

  - section: "Module 3: Analysis of DNA modification"
    description: |
      Learn how to detect DNA methylation through sequencing and to analyze its distribution across genome regions.
    tutorials:
      - name: introduction-dna-methylation
        topic: epigenetics
      - name: methylation-seq
        topic: epigenetics

---

