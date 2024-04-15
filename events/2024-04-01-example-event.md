---
layout: event
title: My Training Event Title
draft: true  # will hide your event from the GTN events list, remove once you are ready to announce your event

description: |
    Short one or two line description of the event

cover-image: assets/images/gat.png
cover-image-alt: GTN Logo on a spiral galaxy background with text galaxy admin training


contributions:
  organisers: # GTN contributors or funders, must be defined in CONTRIBUTORS.yaml
  - shiltemann
  - hexylena
  instructors:
  - bebatut
  - fpsom
  funding:
  - gallantries

tags: [Topic 1, Topic 2, 5-day course]

date_start: 2024-04-01
date_end: 2024-04-02 # optional, if event is more than one day

# Required, but minimally the Name field for online events
location:
  name: Bioinf Dept
  address: 42 E Main St.
  city: Reyjkjavik
  country: Iceland
  #region: # optional
  postcode: 912NM

cost: free # Or, e.g. 150 EUR, must be space separated, must include a currency in ISO 4217 format
audience: This event is intended for PhD students interested in Genomics. A basic knowledge in R is useful but not required.
contact_email: organisers@example.com
async: false # if asynchronous, we will not display the time columns on the program
mode: online # In-person?

registration:
  link: https://example.org
  deadline: 2024-01-01

feedback:
  link: https://example.org
  deadline: 2024-05-01


program:
  - section: "Monday: Introduction"
    description: Short description of the program.
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


---

Longer description of the course. This will be placed after the practical info, but before the course program.
