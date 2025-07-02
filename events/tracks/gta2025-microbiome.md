---
layout: event-track

title: Microbiome analysis
description: This track explains the capability to analyze microbiome data using Galaxy, covering essential tasks such as amplicon analysis and its visualization. It progresses to more advanced topics, including assembly and binning - the requirements for reconstructing Metagenome-Assembled Genomes (MAGs) from microbiome samples, and finally extends to complex analyses like metatranscriptome studies and pathogen detection. Start with the tutorial at your own pace. If you need support contact us via the Slack Channel [gta_microbiome](https://gtnsmrgsbord.slack.com/channels/{{page.slack_channel}}). 🚨 The workflows used for this track where tested thoroughly on the EU server, if you encounter any issues with execution on other servers, we recommend to retry the tutorial on usegalaxy.eu!

# [#gta_microbiome](https://gtnsmrgsbord.slack.com/channels/{{page.slack_channel}}). Please note that the tutorials of this track where all tested successfully on useGalaxy.eu, therefore it is recommended to run all tools on this server !

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


program:
  - section: "Introduction"
    description: |
      Introduction to Microbiome Analysis.
    tutorials:
      - name: introduction
        topic: microbiome

  - section: "Metabarcoding"
    description: |
      Identify and analyze the diversity of species in environmental samples by amplifying and sequencing specific genetic markers
    tutorials:
      - name: dada-16S
        topic: microbiome

  - section: "Metagenomics"
    description: |
      Assemble, bin, and profile the taxonomic diversity using collective DNA from environmental samples. 🚨 The interactive tools Pavian and Phinch are currently not working on usegalaxy.org/org.au and .fr. Please use usegalaxy.eu for this step. 
    tutorials:
      - name: metagenomics-assembly
        topic: microbiome
      - name: metagenomics-binning
        topic: microbiome
      - name: taxonomic-profiling
        topic: microbiome
      - name: diversity
        topic: microbiome

  - section: "Metatranscriptomics"
    description: |
      Profile the taxonomic diversity and the functional processes from collective RNA from environmental samples
    tutorials:
      - name: metatranscriptomics
        topic: microbiome

#  - section: "Metaproteomics"
#    description: |
#      Generation and search of any metaproteomics database, verification and quantification of the proteins/peptides, statistical analysis of the quantified peptides
#    tutorials:
#      - name: clinical-mp-1-database-generation
#        topic: proteomics
#      - name: clinical-mp-2-discovery
#        topic: proteomics
#      - name: clinical-mp-3-verification
#        topic: proteomics
#      - name: clinical-mp-4-quantitation
#        topic: proteomics
#      - name: clinical-mp-5-data-interpretation
#        topic: proteomics

  - section: "Pathogen detection"
    description: |
      Detect and track pathogens from metagenomic Nanopore sequencing
    tutorials:
      - name: pathogen-detection-from-nanopore-foodborne-data
        topic: microbiome

---
