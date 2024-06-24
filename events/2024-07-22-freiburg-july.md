---
layout: event

# Status of this page, remove these lines when ready
draft: true  # will hide your event from the GTN events list
status: wip  # add a banner that warns that the contents of the page are still subject to change


# Description of your event
title: Workshop on high-throughput sequencing data analysis with Galaxy
description: |
  This course introduces scientists to the data analysis platform Galaxy. The course is a beginner course; there is no requirement of any programming skills.
cover-image:         # image for your corse, put in 'events/images' folder
cover-image-alt:     # supply alt text describing your image


# Practical Information
date_start: 2024-07-22
date_end: 2024-07-26 

cost: free 
audience: Scientist with no or little Galaxy experience who want to analyse sequencing data.
contact_email: erxleben@informatik.uni-freiburg.de
async: false 
mode: In-person

registration:
  link: https://docs.google.com/forms/d/e/1FAIpQLSeGHShGhMvFvK0Jf3TNn0xgSMVboabWiTPfP2s3L1iDM0qTzA/viewform
  deadline: 2024-05-29


# Location of the event
# For online events, just the 'name' is enough
location:
  name: University of Freiburg, Germany
  address: Werthmannstrasse 4
  city: Freiburg
  country: Germany
  #region: # optional
  postcode: 79104
  geo:
    lat: 51.9109324
    lon: 4.4680514


# People involved, organisers, speakers, funders, etc
# Must be defined in CONTRIBUTORS.yaml file
contributions:
  organisers:
  - erxleben 
  instructors:
  - erxleben
  - teresa-m
  funding:
  - dataPlant
  - EOSC Eurosice Gateway 
  - denbi


# Galaxy and other infrastructure that will be used for your event.
# This will be used to create the setup instructions for participants
infrastructure:
  tiaas: true    # tiaas = Training Infrastructure as a Service, and can be requested (for free) from all major Galaxies

  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU
      tiaas_link: https://usegalaxy.eu/join-training/smorgasbord3



  custom:  # optional, any other setup instructions you want to add to the "Setup" tab
    description: |
      Before joining the course, please make sure to:
        - Bring a laptop with at least 8GB of RAM.
        - Do the [Intro to Galaxy](http://training.galaxyproject.org/topics/introduction/tutorials/galaxy-intro-short/tutorial.html) tutorial if you are not yet familiar with Galaxy

# Program of your course
# Add GTN tutorials by supplying the topic and tutorial name
# For non-GTN sessions, add a "type:custom" session and description
program:
  - section: "Monday: Introduction"  # section title is optional
    description: |
      Welcome - Galaxy introduction - Quality control - From 9:15am to 4pm

  - section: "Tuesday: ChIP-Sequencing"
    description: |
      ChIP-Sequencing - from 9:15am to 5pm

  - section: "Wednesday: RNA-Sequencing"
    description: |
      RNA-Sequencing - from 9:15am to 5pm


  - section: "Thusday: Variant Calling/Exome Seqencing"
    description: |
      Variant Calling/Exome Seqencing - from 9:15am to 5pm


  - section: "Friday: Metagenomics"
    description: |
      Metagenomics/ Foodborne - from 9:15am to 5pm

---

# Welcome to the Comprehensive Galaxy Workshop: From Introduction to Advanced Applications

Embark on a deep dive into the world of Galaxy, the leading platform for data-intensive biomedical research. This workshop is designed for researchers, students, and data analysts who wish to harness the full potential of Galaxy in various genomic studies. Whether you're a beginner seeking to understand the basics or an intermediate user looking to refine your skills in specific applications, this workshop offers valuable insights and hands-on experiences.

Program Highlights:

    Galaxy Introduction and Quality Control: Start with the fundamentals of navigating Galaxy, understanding its core features, and utilizing its tools for data manipulation and analysis. Learn the importance of quality control and how to implement these practices to ensure the integrity and accuracy of your research data.

    ChIP-seq Analysis: Explore the techniques and tools available in Galaxy for Chromatin Immunoprecipitation Sequencing (ChIP-seq). This session will guide you through the process of analyzing protein-DNA interactions, essential for understanding regulatory networks and mechanisms in genomics.

    RNA-seq Analysis: Delve into RNA sequencing analysis using Galaxy. Learn how to interpret expression data, compare differential gene expression, and uncover the complexities of transcriptomics in a user-friendly environment.

    Variant Calling: Master the methods of variant calling with Galaxy. This part of the workshop focuses on detecting genetic variants from sequencing data, crucial for studies in genetics and personalized medicine.

    Microbiome Analysis: Gain expertise in microbiome data analysis. This session introduces tools and workflows to analyze microbial communities, helping you understand microbiome diversity and its implications on health and disease.

This comprehensive workshop will not only increase your proficiency in using Galaxy but also enhance your ability to conduct sophisticated analyses in various fields of biological research. Connect with experts and peers, gain practical skills, and take your research capabilities to the next level.
