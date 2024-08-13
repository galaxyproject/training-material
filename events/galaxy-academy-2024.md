---
layout: event

title: Galaxy Academy 2024
description: |
  The Galaxy Academy is a online training event for Beginners as well as learners who would like to improve there Galaxy data analysis skills. Over the course of one week we will have a different topic and focus every day. 

draft: true

date_start: 2024-10-07
date_end: 2024-10-11

cost: free 
audience: Everyone who would like to get to know Galaxy, learn bioinformatics data analysis, or master a specific new kind of analysis is welcome.
contact_email: academy@galaxyproject.org
async: true


contributions:
    organisers:
        - shiltemann
        - dannon
        - nakucher
        - jdavcs
        - nekrut
        - delphine-l
        - annasyme
        - teresa-m
        - erxleben
    instructors:
        - pratikdjagtap
        - foellmelanie
        - subinamehta
        - timothygriffin
        - wm75
        - bebatut
        - natefoo
        - teresa-m
        - erxleben
        - lldelisle
        - stephanierobin
        - rlibouba
        - abretaud
        - paulzierep
        - bernt-matthias
        - plushz
        - EngyNasr
        - clsiguret
        - pavanvidem
        - delphine-l
        - anuprulez
        - deeptivarshney
    funding:
        - gallantries


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
      On the first day, you will get to know Galaxy and some basics of sequence data analysis. Feel free to skip this day if you don't need this introduction. Please enjoy this introduction and prepare for the upcoming days.
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
      - title: Transciptomics
        link: events/tracks/gta2024-transcriptomics.md
      - title: Single Cell
        link: events/tracks/gta2024-single-cell.md
      - title: Microbime
        link: events/tracks/gta2024-micorbiome.md
      - title: Bacterial Genomics
        link: events/tracks/gta2024-bacterial-genomeics.md
      - title: BY-COVID
        link: events/tracks/gta2024-bycovid.md
      - title: Machine Learning
        link: events/tracks/gta2024-ml.md

  - section: "Friday: Grab bag"
    description: | 
      Can't get enough? Then please pick any of the tutorials of the GTN. Please be aware that we can only give consider tutorials of the first day or any that is part of a Track of the Galaxy Training Academ 2024.
---
# Wellcome to the Galaxy Training Academy
Do you want to learn how to use Galaxy, a open source data analysis platform. Than you are at the right place. We offer here a 5-day Global Online and Asynchronous learining event.

On the first day you can make your self familiar with the Galaxy platform. In the next days you can follow different tracks, please go to the program tab for more informaiton. 

You can follow our porvieded leraning pathes in your own past using our provieded self-learining materials. Next to the programm you will find Slack channels you can join to exchange with othe and other participants during the event. Here you will also find help if you have qustions or run into an issue during the training. We try to cover all time zoons with helpers for each topic, but pleses be pacient if you do not get an emediat respons.

You only need a browser and an account at a galaxy instance registerd for this event. Please have a look at the setup tab.