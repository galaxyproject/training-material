---
# layout: learning-pathway  # (uncomment this line to activate it)

title: Title of your pathway
description: |
  Description of the pathway. What will be covered, what are the learning objectives, etc?
  Make this as thorough as possible, 1-2 paragraphs. This appears on the index page that
  lists all the learning paths, and at the top of the pathway page

type: use  # 'use' for science topics, or admin-dev or instructors
tags: [some, keywords, here ]
editorial_board:
- shiltemann

cover-image: assets/images/gat.png # optional, define your own cover image for the pathway index page, default is the GTN logo
cover-image-alt: alt text for this image

pathway:
  - section: "Module 1: Title"
    description: |
      description of the module. What will be covered, what should learners expect, etc.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
      - name: galaxy-intro-101
        topic: introduction

  - section: "Module 2: Title"
    description: |
      description of the tutorial
      will be shown under the section title
    tutorials:
      - name: quality-control
        topic: sequence-analysis
      - name: mapping
        topic: sequence-analysis
      - name: general-introduction
        topic: assembly
      - name: chloroplast-assembly
        topic: assembly
      - name: "My non-GTN session"
        external: true
        link: "https://example.com"
        type: hands_on  # or 'slides'


# you can make as many sections as you want, with as many tutorials as you want

---

You can put some extra information here. Markdown syntax can be used. This is shown after the description on the pathway page, but not on the cards on the index page.

