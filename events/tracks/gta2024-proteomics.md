---
layout: event-track

title: Proteomics
description: Learn all about proteomics and how to analyze data in Galaxy

program:
  - section: "Introduction"  # section title is optional
    description: |
      Get started
    tutorials:
      - name: introduction
        topic: proteomics
      - name: maxquant-msstats-dda-lfq
        topic: proteomics

  - section: "Advanced analysis"
    description: |
      Some more advanced analysis
    tutorials:
      - name: encyclopedia
        topic: proteomics
      - name: proteogenomics-dbcreation
        topic: proteomics
      - name: proteogenomics-novel-peptide-analysis
        topic: proteomics
      - name: metaproteomics
        topic: proteomics

  - section: "Clinical Metaproteomics Workflow"
    description: |
        This learning path aims to teach you the basics of how to perform metaproteomics analysis of the clinical data within the Galaxy platform. You will learn how to use Galaxy for analysis and will be guided through the most common first steps of any metaproteomics database generation to searching the database, verifying the proteins/peptides, and data analysis.
 
 
 - section: "Module 1: Database generation"
    description: |
      Get a first look at the Galaxy platform for data analysis. We start with a short introduction to familiarize you with the Galaxy
      interface, and then proceed with understanding how to generate a customized database for clinical metaproteomics
    tutorials:
      - name: clinical-mp-1-database-generation
        topic: proteomics

  - section: "Module 2: Discovery"
    description: |
      This section helps to guide the users through the MSMS dataset search against the compact database generated in the first module.
      The identified peptides and proteins from various software will be combined later to perform verification.
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
      differentially abundant proteins present in the sample and their abundance in various conditions.
    tutorials:
      - name: clinical-mp-4-quantitation
        topic: proteomics  

  - section: "Module 5: Data Interpretation"
    description: |
      We performed statistical analysis of the quantified peptides using MS stats and Unipept to perform taxonomic classification.
    tutorials:
      - name: clinical-mp-5-data-interpretation
        topic: proteomics

---
