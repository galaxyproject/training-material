---
layout: learning-pathway
title: Bacterial comparative genomics
description: |
  This learning path aims to teach you how to analyze bacterial genomes, starting with genome retrieval from NCBI and GTDB. The workflow will include quality control (CheckM2), dereplication (dRep), and pangenomic analysis (PPanGolin). Functional and structural annotation will be performed, with a focus on identifying mobile genetic elements and antimicrobial resistance genes. Phylogenetic analysis will compare core genome trees and presence-absence pangenome trees to reveal complementary evolutionary insights.

draft: true
type: use
tags: [bacteria, comparative genomics, pangenome, phylogeny, microgalaxy]

editorial_board:
- bebatut
funding:
- ifb

pathway:
  - section: "Module: Dataset construction"
    description: |
      For comparative genomics, dataset construction is the first essential step. Here we explain how to build a dataset mixing public and private genomes or assemblies.
    tutorials:
      - name: bacterial-comparative-genomics-dataset-construction
        topic: genome-annotation

  - section: "Module: Genome quality control"
    description: |
      Before genome comparison, quality of genomes and assemblies needs to be checked.
    tutorials:
      - name: bacterial-genome-quality-control
        topic: genome-annotation

  - section: "Module: Structural and functional annotation"
    description: |
      Genomes and assemblies must can be annotated to detect genes with a focus on identifying mobile genetic elements and antimicrobial resistance genes.
    tutorials:
      - name: bacterial-genome-annotation
        topic: genome-annotation
      - name: amr-gene-detection
        topic: genome-annotation

  - section: "Module: Pangenomic analysis"
    description: |
      Genomes and assemblies must can be annotated to detect genes with a focus on identifying mobile genetic elements and antimicrobial resistance genes.
    tutorials:
      - name: bacterial-pangenomics
        topic: genome-annotation

  - section: "Module: Phylogenetic analysis"
    description: |
      Phylogenetic analysis allows visualization and compares core genome trees and presence-absence pangenome trees to reveal complementary evolutionary insights
    tutorials:
      - name: bacterial-comparative-genomics
        topic: evolution
---
