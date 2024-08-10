---
layout: event

# Status of this page, remove these lines when ready
draft: true  # will hide your event from the GTN events list
status: wip  # add a banner that warns that the contents of the page are still subject to change


# Description of your event
title: My Training Event Title
description: |
  Short description of the event (one or two sentences).
  A longer description can be placed at the bottom of this document.
cover-image:         # image for your corse, put in 'events/images' folder
cover-image-alt:     # supply alt text describing your image


# Practical Information
date_start: 1970-04-01
date_end: 1970-04-02 # optional, if event is more than one day

cost: free # Or, e.g. 150 EUR
audience: This event is intended for PhD students interested in Genomics. A basic knowledge in R is useful but not required.
contact_email: organisers@example.com
async: false # if asynchronous, we will not display the time columns on the program
mode: online # In-person

registration:
  link: https://example.org
  deadline: 2024-01-01


# Location of the event
# For online events, just the 'name' is enough
location:
  name: Erasmus Medical Center    # can be e.g. "Online" for online
  address: Dr. Molewaterplein 40
  city: Rotterdam
  country: The Netherlands
  #region: # optional
  postcode: 3015 GD
  geo:
    lat: 51.9109324
    lon: 4.4680514


# People involved, organisers, speakers, funders, etc
# Must be defined in CONTRIBUTORS.yaml file
contributions:
  organisers:
  - shiltemann
  - hexylena
  instructors:
  - bebatut
  - fpsom
  funding:
  - gallantries


# Galaxy and other infrastructure that will be used for your event.
# This will be used to create the setup instructions for participants
infrastructure:
  tiaas: true    # tiaas = Training Infrastructure as a Service, and can be requested (for free) from all major Galaxies

  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU
      tiaas_link: https://usegalaxy.eu/join-training/smorgasbord3
    - server: https://usegalaxy.org
      name: Galaxy Main
      tiaas_link:

  support: # optional, remove if not using online support
     platform: Slack
     join_link:    # invite link; GTN Slack by default
     channel: my-event  # instructors can create channels on the GTN slack themselves.
     link: "https://gtnsmrgsbord.slack.com/archives/C032C2MRHAS" # will use the #general channel on GTN slack by default.

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
      Short description of the program in this section.
      Markdown formatted
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
        time: "09:00 - 10:00"
      - name: data-manipulation-olympics
        topic: introduction
        time: "10:00"
      - type: custom
        name: Coffee Break
        time: "10:30"
      - name: data-manipulation-olympics
        topic: introduction
        time: "10:45 - 12:00"
      - name: ansible-galaxy
        topic: admin
        time: "13:00"
      - type: custom
        name: Custom Session
        time: "16:00 - 17:00"
        description: |
          Description of the custom session here, in markdown, can add
          [links](https://example.com) if needed
      - type: custom
        name: Wrap-up & Drinks
        time: "17:00 - 18:00"
        description: |
          Time for some well-deserved drinks and socializing!

  - section: "Tuesday: An advanced look at .."
    description: |
      Short description of the program for this section.

  - section: "Track 1"
    subsection: true  # will treat this section as a subsection of the previous (i.e. smaller heading), useful to split day into tracks
    description: "you can further subdivide into multiple subsections/tracks as well"
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
        time: "09:00 - 10:00"
      - type: custom
        name: "Wrap-up"

  - section: "Track 2"
    subsection: true  # will treat this section as a subsection of the previous (i.e. smaller heading), useful to split day into tracks
    description: "blabla"
    tutorials:
      - name: mothur-miseq-sop-short
        topic: microbiome

  - section: "Wednesday: Track buttons example"
    description: Today we have two different tracks, click on one of the buttons below to view the track program.
    tracks:  # instead of tutorials, you can also define tracks, this will create a button per tracks that will lead to a different page
      - title: Track 1
        link: events/tracks/example-track1.md
      - title: Track 2
        link: events/tracks/example-track2.md
      - title: Ecology
        link: learning-pathways/intro-to-galaxy-and-ecology.md # can also link to learning pathways (or any other GTN page)
      - title: Climate
        link: https://usegalaxy.eu # external links also possible
      - title: Plants
        link: http://example.com
      - title: CYOA
      - title: Track 7

---

Longer description of the course. This will be added to the overview page of your course.

You can add anything you like here, text, images, tables.

All in markdown format.
