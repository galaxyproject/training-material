---
layout: event-track

title: Microbiome analysis
description: |
  This track introduces microbiome data analysis using Galaxy, covering both foundational and advanced workflows. 
  
  It begins with metabarcoding / amplicon-based analysis using DADA2, a widely adopted standard approach for processing and visualizing microbiome data, which can be directly applied to real-world datasets.

  The read-based metagenomics section introduces taxonomic classification using Kraken2 and MetaPhlAn, along with downstream diversity analysis. It covers how sequencing reads are assigned to microbial taxa using both k-mer–based and marker gene–based approaches, how classification results can be interpreted, and how community composition can be explored using diversity metrics and visualizations. Together, these steps provide a foundation for understanding microbial community structure directly from sequencing reads without the need for assembly.

  The Metagenome-Assembled Genomes (MAGs) section provides a comprehensive overview of MAG generation, developed in collaboration with the [FAIRyMAGs project](https://github.com/usegalaxy-eu/FAIRyMAGs). It includes all key steps: quality control, removal of host and contaminant reads, assembly, binning (demonstrated with four different tools), bin refinement, and functional annotation.

  A complete end-to-end workflow is provided as a foundation for analyzing your own samples. The track then progresses to more advanced topics, including metatranscriptomics and pathogen detection.

  You can start with the tutorials at your own pace. However, if you plan to work through the MAG workflow, we recommend allocating several days, as it is a complex and computationally intensive topic.

  If you need support, contact us via the Slack channel [gta_microbiome](https://gtnsmrgsbord.slack.com/channels/{{page.slack_channel}}).

  🚨 The workflows used for this track were tested thoroughly on the EU server. If you encounter any issues with execution on other servers, we recommend retrying the tutorial on usegalaxy.eu!

slack_channel: gta_microbiome

contributions:
    organisers:
        - paulzierep
        - bebatut
    instructors:
        - annasyme
        - bebatut
        - EngyNasr
        - GarethPrice-Aus
        - igormakunin
        - jennaj
        - bernt-matthias
        - paulzierep
        - plushz
        - RZ9082
        - SaimMomin12
        - pratikdjagtap
        - santamccloud
        - gdefazio
        - minamehr

program:
  - section: "Introduction"
    description: |
      Introduction to Microbiome Analysis.
    tutorials:
      - name: introduction
        topic: microbiome

  - section: "Metabarcoding"
    description: |
      Identify and analyze the diversity of species in environmental samples by amplifying and sequencing specific genetic markers.
    tutorials:
      - name: dada-16S
        topic: microbiome

  - section: "Read-based Metagenomics"
    description: |
     After basic sequence 
     Profile the taxonomic diversity using collective DNA from environmental samples. 🚨 The interactive tools Pavian and Phinch are currently not working on usegalaxy.org/org.au and .fr. Please use usegalaxy.eu for this step. 
    tutorials:
      - name: quality-control
        topic: sequence-analysis
      - name: host-removal
        topic: microbiome
      - name: taxonomic-profiling
        topic: microbiome
      - name: diversity
        topic: microbiome

  - section: "FAIRyMAGs presents: Metagenome-Assembled Genomes (MAGs) building"
    description: |
      Learn all the steps to build Metagenome-Assembled Genomes (MAGs). Start (optionally) with quality control and contamination and host read removal. Then follow assembly and binning, followed by functional annotation. Finally, the full end-to-end workflow can be used as a basis to build MAGs for your own samples.
    tutorials:
      - name: quality-control
        topic: sequence-analysis
      - name: host-removal
        topic: microbiome
      - name: metagenomics-assembly
        topic: microbiome
      - name: metagenomics-binning
        topic: microbiome
      - name: bacterial-genome-annotation
        topic: genome-annotation
      - name: mags-building
        topic: microbiome

  - section: "Metatranscriptomics"
    description: |
      Profile the taxonomic diversity and the functional processes from collective RNA from environmental samples
    tutorials:
      - name: metatranscriptomics
        topic: microbiome

  - section: "Pathogen detection"
    description: |
      Detect and track pathogens from metagenomic Nanopore sequencing
    tutorials:
      - name: pathogen-detection-from-nanopore-foodborne-data
        topic: microbiome

---