---
layout: event

title: Galaxy Academy 2024
description: |
  The Galaxy Academy is a online training event for Beginners as well as learners who would like to improve their Galaxy data analysis skills. Over the course of one week, we will have a different topic and focus every day.


registration:
  link: TODO
  deadline: 2024-10-07
  open: false

date_start: 2024-10-07
date_end: 2024-10-11

cost: free
audience: Everyone who would like to get to know Galaxy and learn bioinformatics data analysis, whether you want to master a specific analysis or learn a new skill.
contact_email: academy@galaxyproject.org
async: true


contributions:
    organisers:
        - erxleben
        - annasyme
        - nekrut
        - dannon
        - delphine-l
        - jdavcs
        - natalie-wa
        - nakucher
        - shiltemann
        - teresa-m
    instructors:
        - erxleben
        - anuprulez
        - abretaud
        - bebatut
        - clsiguret
        - dianichj
        - deeptivarshney
        - delphine-l
        - EngyNasr
        - lldelisle
        - bernt-matthias
        - foellmelanie
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

program:
  - section: "Monday: Galaxy introduction"
    description: |
      Kick off the week with a hands-on introduction to Galaxy, covering everything from basic navigation and data manipulation to reproducing published analyses, quality control, and mapping sequences to a reference genome. Whether you're new to Galaxy or looking to strengthen your skills, today's sessions will equip you with the foundational knowledge needed for more advanced topics.
    tutorials:
      - type: custom
        name: "Start to get to know Galaxy "
      - name: galaxy-intro-101-everyone
        topic: introduction
      - name: data-manipulation-olympics
        topic: introduction
      - name: galaxy-reproduce
        topic: introduction
      - name: options-for-using-galaxy
        topic: introduction
      - type: custom
        name: "Fundamentals of Sequences analysis "
      - name: quality-control
        topic: sequence-analysis
      - name: mapping
        topic: sequence-analysis

  - section: "Tuesday to Thursday: Pick a track"
    description: Today you will learn all about transcriptomics
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
        link: events/tracks/gta2024-micorbiome.md
      - title: Bacterial Genomics
        link: events/tracks/gta2024-bacterial-genomeics.md
      - title: BY-COVID
        link: events/tracks/gta2024-bycovid.md
      - title: Machine Learning
        link: events/tracks/gta2024-ml.md

  - section: "Friday: Grab bag"
    description: |
      Can't get enough? Then please pick any of the tutorials of the GTN. Please be aware that only trainings that are part of the introduction day or a learning track have been tested on all instances for the event. The trainers present on Slack will do their best to help you if you have a problem and answer questions, but they may not be expert in the topic you selected.
---
# Welcome to the Galaxy Training Academy
Do you want to learn how to use Galaxy, a open source data analysis platform? Then you are at the right place. We offer here a 5-day Global Online and Asynchronous learning event.

On the first day you can make yourself familiar with the Galaxy platform. In the next days you can follow different tracks, please go to the program tab for more informaiton.

You can set your own pace on your learning journey using our provided self-learning materials. Next to the program, you will find Slack channels you can join to exchange with the trainers and other participants during the event. Here you will also find help if you have qustions or run into an issue during the training. We try to cover all time zones with helpers for each topic, but please be patient if you do not get an immediat response.

You only need a browser and an account on one of the galaxy instances registered for this event. Please have a look at the setup tab.
