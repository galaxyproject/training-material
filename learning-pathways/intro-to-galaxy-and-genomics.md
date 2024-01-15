---
layout: learning-pathway
tags: [beginner]
type: use

editorial_board:
- shiltemann
funding:
- gallantries

title: Introduction to Galaxy and Sequence analysis
description: |
  This learning path aims to teach you the basics of Galaxy and analysis of sequencing data.
  You will learn how to use Galaxy for analysis, and will be guided through the most common
  first steps of any genome analysis; quality control and a mapping or assembly of your genomic
  sequences.

priority: 1

pathway:
  - section: "Module 1: Introduction to Galaxy"
    description: |
      Get a first look at the Galaxy platform for data analysis. We start with a
      short introduction (video slides & practical) to familiarize you with the Galaxy
      interface, and then proceed with a slightly longer introduction tutorials where
      you perform a first, very simple, analysis.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
      - name: galaxy-intro-101
        topic: introduction

  - section: "Module 2: Basics of Genome Sequence Analysis"
    description: When analysing sequencing data, you should always start with a quality control step to clean your data and make sure your data is good enough to answer your research question. After this step, you will often proceed with a mapping (alignment) or genome assembly step, depending on whether you have a reference genome to work with.
    tutorials:
      - name: quality-control
        topic: sequence-analysis
      - name: mapping
        topic: sequence-analysis
      - name: general-introduction
        topic: assembly
      - name: chloroplast-assembly
        topic: assembly

---

New to Galaxy and/or the field of genomics? Follow this learning path to get familiar with the basics!

