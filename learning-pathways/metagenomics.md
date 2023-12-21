---
layout: learning-pathway  # (uncomment this line to activate it)
type: use
cover-image: assets/images/microgalaxy-logo.png
cover-image-alt: "microgalaxy logo"

editorial_board:
- bebatut
funding:
- gallantries

title: Metagenomics data processing and analysis for microbiome
description: |
  This learning path aims to teach you the basics of Galaxy and analysis of metagenomics data.
  You will learn how to use Galaxy for analysis, and will be guided through the common
  steps of microbiome data analysis: quality control, taxonomic profiling, taxonomic binning, assembly, functional profiling, and also some applications
tags: [microbiome]

pathway:
  - section: "Module 1: Introduction to metagenomics"
    description: Why study the microbiome? What are the different approaches for metagenomics? This module will give you a short introduction to metagenomics.
    tutorials:
      - name: introduction
        topic: metagenomics

  - section: "Module 2: Introduction to Galaxy"
    description: |
      Get a first look at the Galaxy platform for data analysis. We start with a
      short introduction (video slides & practical) to familiarize you with the Galaxy
      interface, and then proceed with a slightly longer introduction tutorials where
      you perform a first, very simple, analysis.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
      - name: galaxy-intro-101
        topic: introduction

  - section: "Module 3: Quality control"
    description: When analysing sequencing data, you should always start with a quality control step to clean your data and make sure your data is good enough to answer your research question.
    tutorials:
      - name: quality-control
        topic: sequence-analysis

  - section: "Module 4: Assembly"
    description: |
        This module covers the following questions:
        - Why metagenomic data should be assembled?
        - What is the difference between co-assembly and individual assembly?
        - What is the difference between reads, contigs and scaffolds?
        - How tools based on De Bruijn graph work?
        - How to assess the quality of metagenomic data assembly?
    tutorials:
      - name: metagenomics-assembly
        topic: metagenomics

  - section: "Module 5: Community taxonomic profiling"
    description: |
        This module covers the following questions:
        - Which species (or genera, families, ...) are present in my sample?
        - What are the different approaches and tools to get the community profile of my sample?
        - How can we visualize and compare community profiles?

        This module will cover taxonomic profiling in theory and also with an example tutorial.
    tutorials:
      - name: taxonomic-profiling
        topic: metagenomics
        #description: More theoretical tutorial
      - name: beer-data-analysis
        topic: metagenomics
        #description: Application tutorial

  #- section: "Module 6: Taxonomic binning"
  #  description:
  #
  #  tutorials:
  #    - name: learner_participation_engagement
  #      topic: metagenomics


  #- section: "Module 7: Community Biodiversity"
  #  description:
  #
  #  tutorials:
  #    - name: learner_participation_engagement
  #      topic: metagenomics

  - section: "Module 8: Applying concepts to metatranscriptomics data"
    description: |
        This module covers the following questions:
        - How to analyze metatranscriptomics data?
        - What information can be extracted of metatranscriptomics data?
        - How to assign taxa and function to the identified sequences?
    tutorials:
      - name: metatranscriptomics
        topic: metagenomics

  - section: "Recommended follow-up tutorials"
    tutorials:
      - name: pathogen-detection-from-nanopore-foodborne-data
        topic: metagenomics
---


