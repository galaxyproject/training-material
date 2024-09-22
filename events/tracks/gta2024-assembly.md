---
layout: event-track

title: Assembly
description: Learn all about Genome assembly. Start with the tutorial at your own pace. If you need support contact us via the Slack Channel [gta_assembly](https://gtnsmrgsbord.slack.com/archives/C07NGNT0DB5).

contributions:
    organisers:
        - delphine-l
    instructors:
        - delphine-l
        - SaimMomin12
        - bgruening
        - clsiguret
        - mschatz
        - GarethPrice-Aus
        - annasyme
        - trungN
        - igormakunin


program:
  - section: "Introduction" 
    description: |
      In this section, we will introduce what is a genome assembly, how it works, and the metrics to evaluate the quality of an assembly.  
    tutorials:
      - name: get-started-genome-assembly
        topic: assembly
      - name: ERGA-post-assembly-QC
        topic: assembly

  - section: "Genome Assembly" 
    description: |
      In this section, we will present genome assemblies for different types of organisms.  
    tutorials:
      - name: mrsa-illumina
        topic: assembly
      - name: vgp_workflow_training
        topic: assembly
      - name: metagenomics-assembly
        topic: assembly

  - section: "Non-nuclear genome assembly" 
    description: |
      When sequencing organisms with non-nuclear DNA, the data contain sequencing for both nuclear and non-nuclear DNA. In this section we will learn how to assemble organelles genome.  
    tutorials:
      - type: custom
        name: "Mitochondrial Genome Assembly"
        description: |
          Tutorial (will come soon)
      - name: chloroplast-assembly
        topic: assembly
        
  - section: "Assembly decontamination" 
    description: |
      In this section we will learn how remove sequences that do not belong to the organism you want to sequence.  
    tutorials:
      - type: custom
        name: "Assembly decontamination"
        description: |
          Tutorial (will come soon)
---



