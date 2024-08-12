---
layout: learning-pathway
tags: [beginner]
type: use


title: Proteogenomics
description: |
  This learning path aims to teach you the basics of how to perform proteogenomics analysis of the Mass spectrometry data within the Galaxy platform. You will learn how to use Galaxy for analysis and will be guided through the most common first steps of any proteogenomics database generation to searching the database, followed by novel peptide data analysis.

cover-image: shared/images/proteomics.png
cover-image-alt: image of a 3D protein folding structure

editorial_board:
- subinamehta

pathway:
  - section: "Module 1: Database generation"
    description: |
      Get a first look at the Galaxy platform for data analysis. We start with a short introduction to familiarize you with the Galaxy interface, and then proceed with understanding how to generate a customized database for proteogenomics.
    tutorials:
      - name: proteogenomics-dbcreation
        topic: proteomics

  - section: "Module 2: Database searching"
    description: |
      This section helps to guide the users through an MSMS dataset search against the customized database generated in the first module. The identified peptides and proteins will be then analyzed later in the novel peptide analysis.
    tutorials:
      - name: proteogenomics-dbsearch
        topic: proteomics

  - section: "Module 3: Novel Peptide Analysis"
    description: |
      The last module in the proteogenomics tutorial is to identify "novel peptides" using BlastP and to localize the peptides to their genomic coordinates. Both inputs from modules 1 and 2 are required to run this tutorial.
    tutorials:
      - name: proteogenomics-novel-peptide-analysis
        topic: proteomics

---

New to Galaxy and/or the field of metaproteomics? Follow this learning path to get familiar with the basics!

