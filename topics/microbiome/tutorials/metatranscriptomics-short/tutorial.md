---
layout: tutorial_hands_on

title: Metatranscriptomics analysis using microbiome RNA-seq data (short)
zenodo_link: https://zenodo.org/record/4776250
questions:
- How to analyze metatranscriptomics data?"
- What information can be extracted of metatranscriptomics data?
- How to assign taxa and function to the identified sequences?
objectives:
- Choose the best approach to analyze metatranscriptomics data
- Understand the functional microbiome characterization using metatranscriptomic results
- Understand where metatranscriptomics fits in 'multi-omic' analysis of microbiomes
- Visualise a community structure
level: Introductory
subtopic: metatranscriptomics
tags:
- metatranscriptomics
- microgalaxy
time_estimation: 3H
key_points:
- Metatranscriptomics data have the same QC profile that RNA-seq data
- A lot of metatranscriptomics sequences are identified as rRNA sequences
- With shotgun data, we can extract information about the studied community structure and also the functions realised by the community
- Metatranscriptomics data analyses are complex and must be careful done, specially when they are done without combination to metagenomics data analyses
contributors:
- pratikdjagtap
- subinamehta
- jraysajulga
- bebatut
- emmaleith
- pravs3683
- shiltemann
- paulzierep
- EngyNasr
redirect_from:
- /topics/metagenomics/tutorials/metatranscriptomics-short/tutorial
edam_ontology:
- topic_3941 # Metatranscriptomics
- topic_3697 # Microbial ecology
- topic_0637 # Taxonomy
- topic_1775 # Function analysis
- topic_0080 # Sequence analysis

recordings:
- youtube_id: HNYop3vLpoM
  date: '2023-05-17'
  speakers:
  - paulzierep
  captioners:
  - paulzierep
  length: 1H5M
  galaxy_version: '23.01'
- captioners:
  - EngyNasr
  - shiltemann
  date: '2021-02-15'
  galaxy_version: '21.01'
  length: 1H30M
  youtube_id: EMaos5u1_a8
  speakers:
  - pratikdjagtap
  - timothygriffin
  - subinamehta
  - shiltemann

---


{% include topics/microbiome/tutorials/metatranscriptomics/content.md short=true %}
