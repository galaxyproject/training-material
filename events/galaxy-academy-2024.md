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
    funding:
        - gallantries


location:
  name: online

program:
  - section: "Monday: Galaxy introduction"
    description: |
      On the first day, you will get to know Galaxy and some basics of sequence data analysis. Feel free to skip this day if you don't need this introduction. Please enjoy this introduction and prepare for the upcoming days.
    tutorials:
      - name: galaxy-intro-101-everyone
        topic: introduction
      - name: data-manipulation-olympics
        topic: introduction
      - name: galaxy-reproduce
        topic: introduction
      - name: options-for-using-galaxy
        topic: introduction
      - name: quality-control
        topic: sequence-analysis
      - name: mapping
        topic: sequence-analysis
 
  - section: "Tuesday to Thursday: Pick a track"
    description: Today you will learn all about transcriptomics
    tracks:  # instead of tutorials, you can also define tracks, this will create a button per tracks that will lead to a different page
      - title: Proteomics
        link: events/tracks/gta2024-proteomics.md

  - section: "Friday: Grab bag"
    description: | 
      Can't get enough? Then please pick whichever tutorial you like from the track of the previous days or from the list below.
    tutorials:
      - name: troubleshooting
        topic: admin

---

some text
