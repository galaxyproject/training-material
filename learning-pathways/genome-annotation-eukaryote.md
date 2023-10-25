---
layout: learning-pathway
title: Genome annotation for eukaryotes
description: |
  Learn how to annotate an eukaryotic genome sequence: identify repeated regions, find the position and function of genes, and even set up a manual curation environment with Apollo.
type: use
tags: [genome annotation, eukaryote]

cover-image: assets/images/gga.png
cover-image-alt: "Galaxy Genome Annotation logo"
editorial_board:
- abretaud
funding:
- gallantries



pathway:
  - section: "Module 1: Introduction"
    description: |
      General introduction to genome annotation
    tutorials:
      - name: introduction
        topic: genome-annotation

  - section: "Module 2: Repeat Masking"
    description: |
      Learn how to identifying, and "mask", repeated regions, a first step before annotating genes
    tutorials:
      - name: repeatmasker
        topic: genome-annotation

  - section: "Module 3: Gene calling and functional annotation"
    description: |
      In this module you will learn to run the Funannotate tool suite to find the position of genes and to functionally annotate them. Optionally you can also identify long non-coding RNAs.
    tutorials:
      - name: funannotate
        topic: genome-annotation
      - name: lncrna
        topic: genome-annotation

  - section: "Module 4: Manual curation"
    description: |
      Automatic annotation is rarely perfect, in this module you will learn how to start a collaborative manual curation project using Galaxy and Apollo
    tutorials:
      - name: apollo-euk
        topic: genome-annotation

---
