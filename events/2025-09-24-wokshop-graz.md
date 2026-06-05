---
layout: event
title: "From data to discovery - Galaxy workshop at University of Graz"

description: |
    This course introduces scientists to the data analysis platform Galaxy. The course is an intermediate course; there is no requirement of any programming skills.

date_start: 2025-09-24
date_end: 2025-09-26

cost: free
audience: Scientist with no or little Galaxy experience who want to analyse high-throughput data.
contact_email: nilchia@informatik.uni-freiburg.de
# async: false
mode: onsite

registration:
  link: https://galaxyproject.org/events/2025-09-22-graz-workshop/
  deadline: 2025-09-24

contributions:
  organisers: [sebastian-preissl]
  instructors: [pavanvidem, Nilchia]
  funding: [mwk, deNBI, uni-graz]

location:
  geo:
    lat: 47.0805168
    lon: 15.4465704
  name: University of Graz
  address:  Mozartgasse 14
  postcode: 8010
  city: Graz
  country: Austria

infrastructure:
  tiaas: true    # tiaas = Training Infrastructure as a Service, and can be requested (for free) from all major Galaxies

  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU

tags:
- introduction
- sequence-analysis
- transcriptomics
- proteomics

# Program of your course
# Add GTN tutorials by supplying the topic and tutorial name
program:
  - section: "Wednesday: Welcome - Galaxy introduction - Quality control"
    description: We will do at least one coffee break in the morning, one in the afternoon, and a one-hour lunch break from 12:30 to 13:30.
    tutorials:
      - type: custom
        name: "Welcome"
        time: "10:00"
      - name: galaxy-intro-short
        topic: introduction
      - name: galaxy-intro-peaks2genes
        topic: introduction
      - type: custom
        name: "[RNA-seq introduction](https://docs.google.com/presentation/d/12A8oLXh0dPqdiXVYzqAZCmmBNGoWfmWLS3g4SDSzTXs/edit?slide=id.p#slide=id.p)"
        description: |
          [Lecture Slides](https://docs.google.com/presentation/d/12A8oLXh0dPqdiXVYzqAZCmmBNGoWfmWLS3g4SDSzTXs/edit?slide=id.p#slide=id.p)
      - type: custom
        name: "End"
        time: "17:00"


  - section: "Thursday: Transcriptomics"
    description:  We will do at least one coffee break in the morning, one in the afternoon, and a one-hour lunch break from 12:30 to 13:30.
    tutorials:
      - name: ref-based
        topic: transcriptomics
        time: "09:00-17:00"

  - section: "Friday: Proteomics and Metabolomics"
    description:  We will do at least one coffee break in the morning, one in the afternoon, and a one-hour lunch break from 12:30 to 13:30.
    tutorials:
      - type: custom
        name: "Beginning"
        time: "09:00"
      - name: introduction
        topic: proteomics
      - name: maxquant-msstats-dda-lfq
        topic: proteomics
      - name: introduction
        topic: metabolomics
      - name: lcms
        topic: metabolomics
      - type: custom
        name: "End"
        time: "17:00"

---
# Welcome to the Comprehensive Galaxy Workshop: From data to discovery - Galaxy workshop at University of Graz

Embark on a deep dive into the world of Galaxy, the leading platform for data-intensive biomedical research. This workshop is designed for researchers, students, and data analysts who wish to
harness the full potential of Galaxy in Transcriptomic, Proteomics and Metabolomics. 

This comprehensive workshop will increase your proficiency in using Galaxy and enhance your ability to conduct sophisticated analyses in various fields of biological research.
Connect with experts and peers, gain practical skills, and take your research capabilities to the next level.

## Topics

- Galaxy Introduction and Quality Control
- Transcriptomics
- Proteomics and Metabolomics

## Learning Objectives

- To get introduced to Galaxy as data analysis platform and the GTN training material
- To learn how to use use tools and databases in Galaxy
- To be able to run different tools and workflows

## Schedule

Please see the Program tab.

## Registration


The number of places in this workshop is limited. The course takes place only in person (not online, not hybrid). There are no fees for the workshop but you need to care for your accommodation, travel costs and catering by yourself.

If you have further questions, please contact Amirhossein Nilchi ([nilchia@informatik.uni-freiburg.de](mailto:nilchia@informatik.uni-freiburg.de))


