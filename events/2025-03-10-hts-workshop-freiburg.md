---
layout: event
title: "Workshop on high-throughput sequencing data analysis with Galaxy"

# external: "https://galaxyproject.org/events/2025-03-10-galaxy-workshop-freiburg/"

description: |
    This course introduces scientists to the data analysis platform Galaxy. The course is a beginner course; no programming skills are required.

date_start: 2025-03-10
date_end: 2025-03-14

cost: free
audience: Scientists with no or little Galaxy experience who want to analyse sequencing data.
contact_email: schneidd@informatik.uni-freiburg.de
# async: false
mode: onsite

registration:
  link: https://forms.gle/cC9ynBrypEWwNA9u8
  deadline: 2025-02-14

contributions:
  organisers: [Sch-Da, teresa-m]
  instructors: [pavanvidem, Nilchia, teresa-m, wm75, EngyNasr, paulzierep]
  funding: [eurosciencegateway, deNBI, nfdi4plants]

location:
  geo:
    lat: 47.993460
    lon: 7.845110
  name: University of Freiburg
  address: Werthmannstrasse 4
  postcode: 79104
  city: Freiburg
  country: Germany

infrastructure:
  tiaas: true    # tiaas = Training Infrastructure as a Service, and can be requested (for free) from all major Galaxies

  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU

tags:
- introduction
- sequence-analysis
- epigenetics
- transcriptomics
- variant-analysis
- microbiome

# Program of your course
# Add GTN tutorials by supplying the topic and tutorial name
program:
  - section: "Monday: Introduction and Quality Control"  # section title is optional
    description: We will do at least one coffee break in the morning, one in the afternoon, and 1h lunch break around noon.
    tutorials:
      - type: custom
        name: "Welcome"
        time: "09:15"
      - name: galaxy-intro-peaks2genes
        topic: introduction
      - name: quality-control
        topic: sequence-analysis
      - type: custom
        name: "End"
        time: "16:00"

  - section: "Tuesday: ChIP-Sequencing"
    description: We will do at least one coffee break in the morning, one in the afternoon, and 1h lunch break around noon.
    tutorials:
      - name: formation_of_super-structures_on_xi
        topic: epigenetics
        time: "09:15-17:00"

  - section: "Wednesday: RNA Sequencing"
    description: We will do at least one coffee break in the morning, one in the afternoon, and 1h lunch break around noon.
    tutorials:
      - name: ref-based
        topic: transcriptomics
        time: "09:15-17:00"


  - section: "Thursday: Variant Calling/Exome Sequencing"
    description:  We will do at least one coffee break in the morning, one in the afternoon, and 1h lunch break around noon.
    tutorials:
      - name: exome-seq
        topic: variant-analysis
        time: "09:15-17:00"


  - section: "Friday: Metagenomics"
    description:  We will do at least one coffee break in the morning, one in the afternoon, and 1h lunch break around noon.
    tutorials:
      - name: pathogen-detection-from-nanopore-foodborne-data
        topic: microbiome
        time: "09:15-16:00"

---
# Welcome to the Comprehensive Galaxy Workshop: From Introduction to Advanced Applications

Embark on a deep dive into the world of Galaxy, the leading platform for data-intensive biomedical research. This workshop is designed for researchers, students, and data analysts who wish to
harness the full potential of Galaxy in various genomic studies. Whether you're a beginner seeking to understand the basics or an intermediate user looking to refine your skills in specific applications,
this workshop offers valuable insights and hands-on experiences.

This comprehensive workshop will increase your proficiency in using Galaxy and enhance your ability to conduct sophisticated analyses in various fields of biological research.
Connect with experts and peers, gain practical skills, and take your research capabilities to the next level.

## Topics:

- Galaxy Introduction and Quality Control
- ChIP-seq Analysis
- RNA-seq Analysis
- Variant Calling
- Microbiome Analysis

## Learning Objectives:

- To get introduced to Galaxy as a data analysis platform and the GTN training material
- To learn how to use tools and databases in Galaxy
- To be able to run different tools and workflows

## Schedule

Please see the Program tab.

## Registaion

The number of places in this workshop is limited. The course takes place only in person (not online, not hybrid). There are no fees for the workshop but you need to care for your accommodation, travel costs and catering by yourself.


