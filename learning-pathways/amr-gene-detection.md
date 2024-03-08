---
layout: learning-pathway
title: Detection of AMR genes in bacterial genomes
description: |
  This learning path aims to teach you the basic steps to detect and check Antimicrobial resistance (AMR) genes in bacterial genomes using Galaxy.
type: use
tags: [amr, bacteria, microgalaxy, one-health]

editorial_board:
- bebatut
funding:
- abromics

pathway:
#  - section: "Module 1: Introduction"
#    description: |
#      General introduction to AMR detection
#    #tutorials:
#    #  - name: introduction
#    #    topic: genome-annotation

  - section: "Module: Species and contamination checking"
    description: |
      Taxonomic assignation is useful in AMR detection to check contamination and confirm species
    tutorials:
      - name: bacterial-isolate-species-contamination-checking
        topic: ecology

  - section: "Module: Assembly"
    description: |
      Assembly is a major step in the process of detecting AMR genes as it combines sequenced reads into contigs, longer sequences where it will be easier to identify genes and in particular AMR genes
    tutorials:
      - name: mrsa-illumina
        topic: assembly
      - name: mrsa-nanopore
        topic: assembly

  - section: "Module: Genome annotation"
    description: |
      The generated contigs can be annotated to detect genes, potential plasmids, etc. This will help the AMR gene detection process, especially the verification and visualization
    tutorials:
      - name: bacterial-genome-annotation
        topic: genome-annotation

  - section: "Module: AMR gene detection"
    description: |
      AMR gene content can be assessed from the contigs to detect known resistance mechanisms and potentially identify novel mechanisms.
    tutorials:
      - name: amr-gene-detection
        topic: genome-annotation

  - section: "Recommended follow-up tutorials"
    tutorials:
      - name: pathogen-detection-from-nanopore-foodborne-data
        topic: microbiome

---
