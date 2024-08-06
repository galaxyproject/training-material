---
layout: event-track

title: Example Track 1
description: In this track we focus on ...

program:
  - section: "Part 1: Topic"  # section title is optional
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

  - section: "Part 2: An advanced look at .."
    description: |
      Short description of the program for this section.


---
