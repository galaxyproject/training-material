---
layout: event

title: Galaxy Training Academy 2024
description: |
  The Galaxy Academy is a online training event for Beginners as well as learners who would like to improve their Galaxy data analysis skills. Over the course of one week, we will have a different topic and focus every day.

cover-image: events/images/galaxy-academy-logo.png
cover-image-alt: logo for the Galaxy Academy event consisting of a laptop surrounded by illustrations of DNA molecules

registration:
  link: https://forms.gle/cxzVatt7MAgiMX12A
  deadline: 2024-09-30
  open: true

date_start: 2024-10-07
date_end: 2024-10-11

cost: free
audience: Everyone who would like to get to know Galaxy and learn bioinformatics data analysis, whether you want to master a specific analysis or learn a new skill.
contact_email: academy@galaxyproject.org

async: true
mode: online

contributions:
    organisers:
        - erxleben
        - annasyme
        - nekrut
        - dannon
        - delphine-l
        - GarethPrice-Aus 
        - jdavcs
        - mschatz
        - natalie-wa
        - nakucher
        - shiltemann
        - teresa-m
    instructors:
        - ahmedhamidawan
        - erxleben
        - annasyme
        - anuprulez
        - abretaud
        - bebatut
        - bgruening
        - clsiguret
        - dannon
        - dianichj
        - deeptivarshney
        - delphine-l
        - elichad
        - EngyNasr
        - emmaustin20
        - GarethPrice-Aus 
        - igormakunin
        - jdavcs
        - lldelisle
        - bernt-matthias
        - foellmelanie
        - mschatz
        - natalie-wa
        - natefoo
        - paulzierep
        - pavanvidem
        - plushz
        - pratikdjagtap
        - RZ9082
        - rlibouba
        - SaimMomin12
        - stephanierobin
        - subinamehta
        - teresa-m
        - timothygriffin
        - tcollins2011
        - nomadscientist
        - wm75
    funding:
        - eurosciencegateway
        - biont
        - nfdi4plants
        - deNBI
        - by-covid
        - elixir-europe
        - mwk
        - abromics
        - ifb


location:
  name: online

infrastructure:
  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU
    - server: https://usegalaxy.org
      name: Galaxy Main
    - server: https://usegalaxy.org.au/
      name: Galaxy AU

  support:
    platform: Slack


program:
  - section: "Monday: Galaxy introduction"
    description: |
      Kick off the week with a hands-on introduction to Galaxy, covering everything from basic navigation and data manipulation to reproducing published analyses, quality control, and mapping sequences to a reference genome. Whether you're new to Galaxy or looking to strengthen your skills, today's sessions will equip you with the foundational knowledge needed for more advanced topics. 
      In the morning you can take part in the Icebreaker by joining us in the [Slack introduction channel](https://gtnsmrgsbord.slack.com/archives/C07NKAJ8THA). Or you can directly start with the tutorials. If you need support contact us via the [Slack introduction channel](https://gtnsmrgsbord.slack.com/archives/C07NKAJ8THA).
  - section: Start to get to know Galaxy
    subsection: true
    tutorials:
      - name: galaxy-intro-101
        topic: introduction
      - name: history
        topic: galaxy-interface
      - name: collections
        topic: galaxy-interface
      - name: data-manipulation-olympics
        topic: introduction
      - name: galaxy-intro-ngs-data-managment
        topic: introduction
  - section: Quick start or fresh up your Galaxy knowledge
    subsection: true
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
  - section: "Fundamentals of Sequences analysis"
    subsection: true
    tutorials:
      - name: quality-control
        topic: sequence-analysis
      - name: mapping
        topic: sequence-analysis

  - section: "Tuesday to Thursday: Pick a track"
    description: |
      Over the course of these three days, you can choose your preferred track and learn how to use Galaxy for data analysis in this research field. If you find multiple topics interesting, feel free to explore more than one track. Each track will guide you through the process, from basic to more advanced analyses, to accommodate learners of all levels.

      For assistance, you can access support through the Slack channel associated with each track module. Please note that while we strive to accommodate all time zones, responses to specific questions may take a bit longer depending on the availability of experts in your time zone.
    tracks:  # instead of tutorials, you can also define tracks, this will create a button per tracks that will lead to a different page
      - title: Proteomics
        link: events/tracks/gta2024-proteomics.md
      - title: Assembly
        link: events/tracks/gta2024-assembly.md
      - title: Transcriptomics
        link: events/tracks/gta2024-transcriptomics.md
      - title: Single Cell
        link: events/tracks/gta2024-single-cell.md
      - title: Microbiome
        link: events/tracks/gta2024-microbiome.md
      - title: Bacterial Genomics
        link: events/tracks/gta2024-bacterial-genomics.md
      - title: BY-COVID
        link: events/tracks/gta2024-bycovid.md
      - title: Machine Learning
        link: events/tracks/gta2024-ml.md

  - section: "Friday: Grab bag"
    description: |
      Can't get enough? Then please pick one of our FAIR tutorials or any of the tutorials of the GTN. Please be aware that only trainings that are part of the introduction day or a learning track have been tested on all instances for the event. The trainers present on Slack will do their best to help you if you have a problem and answer questions, but they may not be expert in the topic you selected.
      You can directly start with you prefered tutorial. If you need support contact us via the Slack Channel [#gta_friday-grab-bag](https://gtnsmrgsbord.slack.com/archives/C07N2A4HQ15).
  - section: Fair training
    subsection: true
    tutorials:
      - type: custom
        name: "[An overview of the RO-Crate concept and its implementations](https://gallantries.github.io/video-library/videos/ro-crates/intro/slides/)"
        description: |
          Lecture Video
      - name: ro-crate-intro
        topic: fair
      - type: custom
        name: "[Registering Galaxy workflows in WorkflowHub](https://gallantries.github.io/video-library/videos/ro-crates/workflowhub/tutorial/)"
        description: |
          Lecture Video
      - name: ro-crate-galaxy-best-practices
        topic: fair
---
# Welcome to the Galaxy Training Academy
Do you want to learn how to use Galaxy, a open source data analysis platform? Then you are at the right place. We offer here a 5-day Global Online and Asynchronous learning event.

On the first day you can make yourself familiar with the Galaxy platform. In the next days you can follow different tracks, please go to the program tab for more informaiton.

You can set your own pace on your learning journey using our provided self-learning materials. Next to the program, you will find Slack channels you can join to exchange with the trainers and other participants during the event. Here you will also find help if you have qustions or run into an issue during the training. We try to cover all time zones with helpers for each topic, but please be patient if you do not get an immediat response.

You only need a browser and an account on one of the galaxy instances registered for this event. Please have a look at the setup tab.
