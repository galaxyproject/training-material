---
layout: learning-pathway
tags: [beginner]
type: use

editorial_board:
- shiltemann
- hexylena
- bebatut
funding:
- gallantries

title: Gallantries Grant - Intellectual Output 2 - Large-scale data analysis, and introduction to visualisation and data modelling


description: |
  This Learning Pathway collects the results of Intellectual Output 2 in the Gallantries Project

cover-image: ./shared/images/Gallantries_logo.png
cover-image-alt: "Gallantries logo with the carpentries wrench in galaxy 2 stripes 1 strip colour scheme."

priority: 5
draft: true

pathway:



  - section: "Year 1: Introduction to large-scale analyses in Galaxy"
    description: |
        Galaxy offers support for the analysis of large collections of data. This submodule will cover the upload, organisation, and analysis of such large sets of data and files. [SC2.1; SC1.3,5]
    tutorials:
      - name: upload-rules
        topic: galaxy-interface
      - name: upload-rules-advanced
        topic: galaxy-interface
      - name: ncbi-sarf
        topic: galaxy-interface
      - name: history-to-workflow
        topic: galaxy-interface
      - name: collections
        topic: galaxy-interface
      - name: workflow-automation
        topic: galaxy-interface
      - name: workflow-editor
        topic: galaxy-interface
      - name: workflow-parameters
        topic: galaxy-interface


  - section: "Year 1: Introduction to the human microbiome analyses"
    description: |
       The human microbiome consists of a community of thousands of species of microorganisms. Sequencing of this community is often performed to identify which species of microorganism are present. This aids in diagnostics and treatment of patients. [SC2.1-3,6; SC1.4,5]
    tutorials:
      - name: beer-data-analysis
        topic: metagenomics
      - name: nanopore-16S-metagenomics
        topic: metagenomics

  - section: "Year 1: Advanced microbiome analysis"
    description: |
       By using more complex sequencing techniques, it is possible to not only obtain information about which organisms are present in the microbiome, but also their activity. This can e.g. aid in identification of antibiotic resistance. This more complex sequencing requires more complex data analysis [SC2.1-4,6; SC1.4,5]
    tutorials:
      - name: pathogen-detection-from-nanopore-foodborne-data
        topic: metagenomics

  - section: "Year 2: Cancer Analysis"
    description: |
       The previous submodules focused on scaling up in terms of number of samples. This submodule will focus on scaling up in terms of complexity. Cancer is a disease of the genome, it is a multifaceted and heterogeneous disease. This leads to complex datasets and analysis pipelines [SC2.3,4; SC1.5]
    tutorials:
      - name: mapping-by-sequencing
        topic: variant-analysis
      # https://github.com/galaxyproject/training-material/pull/3802

  - section: "Year 2: Intro to machine learning"
    description: |
       Going beyond conventional statistics, many scientific data analyses benefit from machine learning techniques for modelling of datasets. This is widely used in biomedical domain. [SC2.4,5; SC1.4]
    tutorials:
      - name: intro-to-ml-with-r
        topic: statistics

  - section: "Year 2: Introduction to the Galaxy visualisation framework"
    description: |
       (This module was cancelled due to insufficiencies in the Galaxy Visualisation Framework.) Galaxy has many options for visualisation of scientific data. This module will cover how to use this framework to create and share visualisation. [SC2.2-3; SC1.1,3,6]
    tutorials: []

  - section: "Year 3: Visualisation of complex multidimensional data"
    description: |
       For advanced visualisation, tools such as Circos may be utilized where Galaxy’s basic visualisation framework does not suffice. [SC2.2-3; SC1.5]
    tutorials:
      - name: circos
        topic: visualisation
      - name: circos-microbial
        topic: visualisation

  - section: "Year 3: Introduction to Visualisation with R and Python"
    description: |
       When the available visualisation options do not suffice, custom plots and visualisations can be created using one of several extensive visualisation libraries available in R and Python. This module will cover the basics of using R and Python to create custom plots and visualisations. [SC2.3; SC1.1]
    tutorials:
      - name: data-manipulation-olympics-viz-r
        topic: data-science
      - name: python-plotting
        topic: data-science

---

Success Criteria:

- SC2.1) Large-scale data analyses and -handling. In this module, learners will gain competency in managing, organizing, and analysing large collections of datasets.
- SC2.2) Analysis of high-dimensional datasets. Real-world scientific studies often involve more complex datasets. For example, combining data from different experiments or timepoints. This more complex experimental setup translates to increased complexity in data analysis.
- SC2.3) Data visualisation. This module will cover the basics of data visualisation to aid with exploration, interpretation of complex datasets.
- SC2.4) Data modelling. This module will introduce learners to the basics data modelling techniques. This is often required for the identification of patterns in data required for e.g. classification.
- SC2.5) Machine learning. This module will also cover more advanced data modelling techniques such as machine learning.
- SC2.6) Reasoning about impact of computation on results. Many choices must be made during data analysis. This includes experimental design, choice of data analysis tools and their parameter settings, and external reference databases. Each of these choices will impact the results. Accurate interpretation of results is only possible with an understanding and awareness of the impact of these factors.

