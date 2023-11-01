---
layout: learning-pathway
tags: [beginner]
type: use

editorial_board:
- abretaud
- hexylena
- shiltemann
funding:
- gallantries

title: Gallantries Grant - Intellectual Output 3 - Data stewardship, federation, standardisation, and collaboration

description: |
  This Learning Pathway collects the results of Intellectual Output 3 in the Gallantries Project

cover-image: ./shared/images/Gallantries_logo.png
cover-image-alt: "Gallantries logo with the carpentries wrench in galaxy 2 stripes 1 strip colour scheme."

priority: 5
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
       One of the key tasks in annotation is producing an official gene set (OGS), and ensuring integrity and validation of all of the curated annotations. This will also further familiarise students with public databases and the process for submitting datasets. [SC3.1, SC3.5]
    tutorials:
      - name: official-gene-set
        topic: genome-annotation


---


Success Criteria:

- SC3.1) Data stewardship. This competency will provide learners with the necessary skills to evaluate data owned by their organisation, identify key metadata and content requirements, and establish controls and assurance metrics to ensure new data produced by their organisation is of sufficient quality.
- SC3.2) Data federation. While hosting data internally is good, sharing data with teams across the Union and around the world is even better. Students need to understand how to achieve data federation and interchange with other organisations.
- SC3.3) Data standardisation. A key component of stewardship and federation, standardisation of data allows it to be reused internally by common pipelines, but more importantly to be submitted to external databases. This is a key requirement of many data related projects.
- SC3.4) Data collaboration. Many projects now scale beyond the ability of an individual to work on alone. Students need to learn how to work together with large datasets.
- SC3.5) FAIR Data. For high quality, reproducible science, datasets should be FAIR; findable, accessible, interoperable, and reusable. This competency aids learners in ensuring their data is FAIR.
