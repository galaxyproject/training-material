---
layout: event

# Status of this page, remove these lines when ready
draft: true  # will hide your event from the GTN events list
status: wip  # add a banner that warns that the contents of the page are still subject to change


# Description of your event
title: Galaxy Ecology training and collabroation fest
description: |
  As colleagues from CCAMLR / Antarctica conservation from Belgium are coming to Concarneau marine station to work on Galaxy Ecology use in their projects, we propose to create an open event, where colleagues from CCAMLR / Antarctica conservation from UK, New Zealand and others are interested to participate remotely, so people wanted to learn how to use and contribute to Galaxy Ecology European Platform during these 3 days.
cover-image:         # image for your corse, put in 'events/images' folder
cover-image-alt:     # supply alt text describing your image

# Tags to help users find your event. Tag with a topic id to have it show op on the topic community page
tags:
- ecology

# Practical Information
date_start: 2025-04-11
date_end: 2025-11-13 # 

cost: free # 
audience: This event is intended for any scientist interesting in ecology data analysis and notably biodiversity indicator production.
contact_email: yvan.le-bras@mnhn.fr
async: false
mode: hybrid

registration:
  link: contact Yvan
  # deadline: 2024-01-01
  # open: false # set this if registration is not open yet

# Location of the event
# For online events, just the 'name' is enough
location:
  name: Concarneau marine station and online   # can be e.g. "Online" for online
  address: quai de la croix
  city: Concarneau
  country: France
  #region: # optional
  postcode: 29900
  geo:
    lat: 47.86846
    lon: -3.91699


# People involved, organisers, speakers, funders, etc
# Must be defined in CONTRIBUTORS.yaml file
contributions:
  organisers:
  - yvanlebras
  instructors:
  - PaulineSGN
  - TuturBabar
  funding:
  - pndb


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
        - Do the [Intro to Galaxy](http://training.galaxyproject.org/topics/introduction/tutorials/galaxy-intro-short/tutorial.html) tutorial if you are not yet familiar with Galaxy

# Program of your course
# Add GTN tutorials by supplying the topic and tutorial name
# For non-GTN sessions, add a "type:custom" session and description
program:
  - section: "Tuesday: Introduction"  # section title is optional
    description: |
      Galaxy platform and Galaxy Ecology introduction
      Markdown formatted
    tutorials:
      - type: custom
        name: Galaxy & Galaxy Ecology introduction
        time: "10:30 - 12:30"
      - type: custom
        name: Lunch
        time: "12:30 - 14h:00"
      - type: custom
        name: How we can contribute code to Galaxy
        time: "14:00 - 15:00"
        description: |
          the "classical" way" (Pauline) & the "Jupyter notebook way" (Arthur)
      - type: custom
        name: Coffee break
        time: "15:00 - 15h:30"
      - type: custom
        name: collaboration fest / hackathon
        time: "15:30 - 18:00"

  - section: "Wednesday: Hack!"
      - type: custom
        name: collaboration fest / hackathon
        time: "9:30 - 18:00"

  - section: "thursday: Hack again!"
      - type: custom
        name: collaboration fest / hackathon
        time: "9:30 - 12:30"
---

Longer description of the course. This will be added to the overview page of your course.

You can add anything you like here, text, images, tables.

All in markdown format.
