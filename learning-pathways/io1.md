---
layout: learning-pathway
tags: [beginner]
type: use

editorial_board:
- fpsom
- shiltemann
- hexylena
funding:
- gallantries

title: Gallantries Grant - Intellectual Output 1 - Introduction to data analysis and -management, statistics, and coding

description: |
  This Learning Pathway collects the results of Intellectual Output 1 in the Gallantries Project

cover-image: ./shared/images/Gallantries_logo.png
cover-image-alt: "Gallantries logo with the carpentries wrench in galaxy 2 stripes 1 strip colour scheme."

priority: 5
draft: true

pathway:
  - section: "Year 1: Coding in Python"
    description: |
        Intro to Coding in Python. Covers variables, functions, and data structures [SC1.1,2]
    tutorials:
      - name: python-basics
        topic: data-science
      - name: python-advanced-np-pd
        topic: data-science

  - section: "Year 1: Coding in Python Modular (Avans)"
    description: |
        Intro to Coding in Python. Covers variables, functions, and data structures [SC1.1,2]

        In collaboration with Avans Hogeschool, an associated Partner we produced the following lessons
    tutorials:
      - topic: data-science
        name: python-math
      - topic: data-science
        name: python-functions
      - topic: data-science
        name: python-types
      - topic: data-science
        name: python-iterables
      - topic: data-science
        name: python-flow
      - topic: data-science
        name: python-loops
      - topic: data-science
        name: python-exceptions
      - topic: data-science
        name: python-files
      - topic: data-science
        name: python-basics-recap
      - topic: data-science
        name: python-glob
      - topic: data-science
        name: python-argparse
      - topic: data-science
        name: python-subprocess
      - topic: data-science
        name: python-venv
      - topic: data-science
        name: python-conda

  - section: "Year 1: Coding in R"
    description: |
        Intro to Coding in R. Covers variables, functions, and data structures [SC1.1,2]
    tutorials:
      - name: r-basics
        topic: data-science
      - name: r-advanced
        topic: data-science
      - name: r-dplyr
        topic: data-science

  - section: "Year 1: Intro to Command Line"
    description: |
        This submodule will cover the basics of the shell (variables, for loops), needed for data handling [SC1.1,2,6]
    tutorials:
      - name: cli-basics
        topic: data-science
      - name: cli-advanced
        topic: data-science
      - name: cli-bashcrawl
        topic: data-science
      - name: snakemake
        topic: data-science

  - section: "Year 1: Intro to Git and GitHub"
    description: |
        This submodule will cover the basics of research software development and sharing (committing, branching, forking, GitHub, etc.) [SC1.1,2,6]
    tutorials:
      - name: bash-git
        topic: data-science
      - name: git-cli
        topic: data-science
      - name: github-command-line-contribution
        topic: contributing
      - name: github-interface-contribution
        topic: contributing

  - section: "Year 2: Introduction to Genomics"
    description: |
       This submodule covers the biological background, as well as the technological concepts involved in genome sequencing, and their effects on downstream data analysis. [SC1.3,4,6]

  - section: "Year 2: Quality Control"
    description: |
       This submodule will cover the evaluation of the quality of datasets, and how to improve quality by a cyclic process of cleaning, trimming and filtering datasets and re-evaluating the quality. [SC1.3-5]
    tutorials:
      - name: quality-control
        topic: sequence-analysis

  - section: "Year 2: Mapping"
    description: |
       This submodule will cover the comparison of genome sequencing samples to a reference genome. The concept of reference data is relevant in many data analyses across life sciences; connecting to online databases and incorporating this data into an analysis. [SC1.3,4]
    tutorials:
      - name: mapping
        topic: sequence-analysis

  - section: "Year 3: Variant Analysis"
    description: |
       This submodule will cover the topic of variant calling; after mapping of sequences to the reference genome, the regions that are different from the reference genome (variants) must be determined, and evaluated for impact. As any two individuals will by definition show many differences, the challenge of distinguishing between healthy variation and potential disease-causing variants is one of the main challenges in variant calling. [SC1.3-5]
    tutorials:
      - name: bash-variant-calling
        topic: data-science

  - section: "Year 3: Transcriptomics"
    description: |
       DNA only describes the potential of the genome; which genes are actually active within the cell and impacting the health and function of the organism, is determined via transcriptomics (RNA sequencing). By integrating data from these two levels of analysis (DNA and RNA), a clearer picture of the state of the cell can be obtained. [SC1.3-5]
    tutorials:
      - name: rna-seq-bash-star-align
        topic: transcriptomics

---

In total, this module will form a course of around 10 days (± 2 days depending on exact analysis stories we identify). Some of these introductory submodules will build on existing training material available in the GTN or Carpentries (~15%).

Success Criteria:

- SC1.1) Basic coding skills. This module will cover the basics of the R and Python coding languages for novices. No coding experience will be assumed nor expected. Basic coding concepts will be introduced (variables, functions, data structures).
- SC1.2) Research software development. We will cover best practices for research software development. It will follow Open Science principles, and include topics such as collaborative code development (e.g. git), reproducible research, code review, and quality control.
- SC1.3) Familiarity with federated data analysis, management, and compute infrastructures. We will introduce the Galaxy platform, a user-friendly web-based analysis platform capable of distributing work across public/private clouds and High-Performance Computing (HPC) resources.
- SC1.4) Basic statistical analysis skills. This submodule will cover the basic concepts involved in statistical analysis of scientific data.
- SC1.5) Data acquisition and integration. Scientific data analyses often require interaction with external datasets. We will cover ways to retrieve data from online data sources, transform it to the required format, and integrate it into the analysis.
- SC1.6) Reproducibility and data sharing. A cornerstone of scientific research is reproducibility. We will cover how to effectively share data and analysis pipelines in order to make scientific results optimally reproducible.
