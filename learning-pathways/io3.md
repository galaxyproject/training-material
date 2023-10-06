---
layout: learning-pathway
tags: [beginner]
type: use

editorial_board:
- abretaud
- hexylena
- shiltemann
funding:
- erasmusplus

title: Gallantries Grant - Intellectual Output 3 - Data stewardship, federation, standardisation, and collaboration

priority: 1
draft: true

pathway:
  - section: "Year 1: Introduction to genomics and genome annotation"
    description: |
       This will give students a good basic knowledge in the application domain of this IO and give them their first taste of data management [SC3.1,SC3.3,SC3.5]
    tutorials:
      - name: introduction
        topic: genome-annotation

  - section: "Year 1: Prokaryotic annotation"
    description: |
       This module will cover the background relevant to annotating prokaryotic genomes in Galaxy (one of the two main classes of genomes), and collaborative curation with Apollo, as well as further exploration of annotation from code. [SC1.5, SC3.1-4]
    tutorials:
      - name: annotation-with-prokka
        topic: genome-annotation
      - name: apollo
        topic: genome-annotation

  - section: "Year 2: FAIR Data"
    description: |
       This submodule will focus specifically on how learners can make their data more FAIR (findable, accessible, interoperable, and reusable) [SC3.5]
    tutorials:
    - topic: fair
      name: fair-intro
    - topic: fair
      name: fair-gtn
    - topic: fair
      name: data-management
    - topic: fair
      name: bioimage-metadata
    - name: ro-crate-intro
      topic: fair
    - name: ro-crate-in-galaxy
      topic: fair
    - name: ro-crate-in-python
      topic: fair
    - name: ro-crate-galaxy-best-practices
      topic: fair
    - name: ro-crate-workflow-run-ro-crate
      topic: fair

  - section: "Year 2: Automatic Annotation"
    description: |
       Building on the modules developed in the previous years, this will be further automated giving students the tools required to scale genome annotation regardless of the size of their organism. [SC1.1, SC1.6, SC2.1, SC3.1, SC3.3]
    tutorials:
      - name: funannotate
        topic: genome-annotation

  - section: "Year 3: Eukaryotic annotation"
    description: |
       This module will cover the background relevant to annotating eukaryotic genomes in Galaxy (the second of the two main genome classes), and collaborative curation with Apollo. Additionally students will learn about automating this annotation process using Galaxy and code. [SC1.5, SC2.1, SC3.1-4]
    tutorials:
      - name: repeatmasker
        topic: genome-annotation
      - name: lncrna
        topic: genome-annotation
      - name: apollo-euk
        topic: genome-annotation

  - section: "Year 3: Official Gene Set"
    description: |
       ⚠️  **TODO** ⚠️One of the key tasks in annotation is producing an official gene set (OGS), and ensuring integrity and validation of all of the curated annotations. This will also further familiarise students with public databases and the process for submitting datasets. [SC3.1, SC3.5]


---
