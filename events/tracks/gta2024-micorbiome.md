---
layout: event-track

title: Microbiome analysis
description: This track explains the capability to analyze microbiome data using Galaxy, covering essential tasks such as amplicon analysis and its visualization. It progresses to more advanced topics, including assembly and binning - the requirements for reconstructing Metagenome-Assembled Genomes (MAGs) from microbiome samples, and finally extends to complex analyses like metatranscriptome studies and pathogen detection.



contributions:
    organisers:
        - paulzierep
        - bebatut
    instructors:
        - bebatut
        - paulzierep
        - bernt-matthias
        - plushz
        - EngyNasr

program:
  - section: "Basics" 
    description: |
      Introduction to Microbiome Analysis. If you encounter any issue please ask us in this Slack channel. 
    tutorials:
      - name: introduction
        topic: microbiome
      - name: dada-16S
        topic: microbiome
      - type: custom
        name: "MGnify amplicon workflow"
        description: |
          Tutorial (will come soon)
  
  - section: "Genome recovery" 
    description: |
      Genome recovery. If you encounter any issue please ask us in this Slack channel. 
    tutorials:  
      - name: metagenomics-assembly
        topic: microbiome
      - name: metagenomics-binning
        topic: microbiome

  - section: "Read-based based analysis" 
    description: |
      Read-based based analysis. If you encounter any issue please ask us in this Slack channel. 
    tutorials: 
      - name: taxonomic-profiling
        topic: microbiome
      - type: custom
        name: "[Calculating α and β diversity from microbiome taxonomic data](https://github.com/galaxyproject/training-material/pull/4282)"
        description: |
          Tutorial (will come soon)
      - name: metatranscriptomics
        topic: microbiome

  - section: "Pathogen detection" 
    description: |
      Pathogen detection. If you encounter any issue please ask us in this Slack channel. 
    tutorials: 
      - name: pathogen-detection-from-nanopore-foodborne-data
        topic: microbiome

---
