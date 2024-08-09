---
layout: learning-pathway
tags: [beginner]
type: use


title: Clinical metaproteomics workflows within Galaxy
description: |
  This learning path aims to teach you the basics of how to perform metaproteomics analysis of the clinical data within the Galaxy platform. You will learn how to use Galaxy for analysis, and will be guided through the most common first steps of any metaproteomics database generation to searching the database, verifying the proteins/peptides, and data analysis.

cover-image: shared/images/proteomics.png
cover-image-alt: image of a 3D protein folding structure


editorial_board:
- subinamehta

pathway:
  - section: "Module 1: Database generation"
    description: |
      Get a first look at the Galaxy platform for data analysis. We start with a
      short introduction to familiarize you with the Galaxy
      interface, and then proceed with understanding how to generate a customized database for clinical metaproteomics
    tutorials:
      - name: clinical-mp-1-database-generation
        topic: proteomics

  - section: "Module 2: Discovery"
    description: |
      This section helps to guide the users through MSMS dataset search against the compact database generated in the first module.
      The identified peptides and proteins from various softwares will be combined later to perform verification.
    tutorials:
      - name: clinical-mp-2-discovery
        topic: proteomics


  - section: "Module 3: Verification"
    description: |
      Here we use the PepQuery tool to verify the presence of the peptides as well as validate that the peptides/proteins
      identified are indeed of microbial origin.
    tutorials:
      - name: clinical-mp-3-verification
        topic: proteomics

  - section: "Module 4: Quantitation"
    description: |
      In this module, we perform quantitative analysis of our data using MaxQuant. Quantitative analysis will help us identify
      differertially abundant proteins present in the sample and their abundance in various conditions.
    tutorials:
      - name: clinical-mp-4-quantitation
        topic: proteomics

  - section: "Module 5: Data Interpretation"
    description: |
      We perform statistical analysis of the quantified peptides using MS stats and also used Unipept to perform taxonomic classification.
    tutorials:
      - name: clinical-mp-5-data-interpretation
        topic: proteomics
---

New to Galaxy and/or the field of metaproteomics? Follow this learning path to get familiar with the basics!

