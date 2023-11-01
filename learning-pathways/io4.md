---
layout: learning-pathway
tags: [beginner]
type: use

editorial_board:
- yvanlebras
- bebatut
funding:
- gallantries

title: Gallantries Grant - Intellectual Output 4 - Data analysis and modelling for evidence and hypothesis generation and knowledge discovery

description: |
  This Learning Pathway collects the results of Intellectual Output 4 in the Gallantries Project

cover-image: ./shared/images/Gallantries_logo.png
cover-image-alt: "Gallantries logo with the carpentries wrench in galaxy 2 stripes 1 strip colour scheme."

priority: 5
draft: true

pathway:
  - section: "Year 1: Biodiversity data handling and visualisation"
    description: |
      learners will understand how to handle biodiversity data and analyse it, as well as elements of visualisation, identifying the optimal visualisation for a dataset. [SC1.1,SC1.4, SC2.1, SC2.3, SC4.1-3]
    tutorials:
      - name: PAMPA-toolsuite-tutorial
        topic: ecology
      - name: regionalGAM
        topic: ecology
      - name: biodiversity-data-exploration
        topic: ecology
      - name: gbif_cleaning
        topic: ecology

  - section: "Year 2: Metabarcoding and environmental DNA data analysis"
    description: |
      analysis of environmental DNA samples requires integrative analysis of highly diversified samples, and new techniques to scale with the data [SC1.4, SC1.5, SC2.1, SC3.1, SC4.1-4]
    tutorials:
      - name: Obitools-metabarcoding
        topic: ecology

  - section: "Year 3: Species distribution modeling"
    description: |
      As an application of data modeling, we will use species migration and biodiversity to teach learners how to build models for complex data and visualise the results. [SC1.1, SC2.4, SC4.1-4]
    tutorials:
      - name: species-distribution-modeling
        topic: ecology

---

Success Criteria:
- SC4.1) Statistical analysis. This will build on the basic statistics covered in IO1 to give a much better statistical comprehension often needed in more advanced analyses like modeling.
- SC4.2) Interactive data visualisation. For most cases, existing visualisations are sufficient, but knowing which visualisation is appropriate and why can be a key point often missed. Additionally sometimes analyses will require custom visualisation such as for geographic information system data.
- SC4.3) Hypothesis generation. When a researcher is handed a large pile of data, figuring out which questions to ask, and what the expected answer is, is the first step of good science.
- SC4.4) Advanced data modelling. Given a hypothesis for some data, a researcher should know how to model changes across some unknown variables, predicting into the future or filling in potential missing gaps in data.
