---
layout: learning-pathway
title: Genome annotation for prokaryotes
description: |
  Learn how to annotate a prokaryotic genome sequence: find the position and function of genes, and even set up a manual curation environment with Apollo.
type: use
tags: [genome annotation, prokaryote]

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

  - section: "Module 2: Gene calling and functional annotation"
    description: |
      In this module you will learn to run the Funannotate tool suite to find the position of genes and to functionally annotate them
    tutorials:
      - name: annotation-with-prokka
        topic: genome-annotation

  - section: "Module 3: Manual curation"
    description: |
      Automatic annotation is rarely perfect, in this module you will learn how to start a collaborative manual curation project using Galaxy and Apollo
    tutorials:
      - name: apollo
        topic: genome-annotation

---
