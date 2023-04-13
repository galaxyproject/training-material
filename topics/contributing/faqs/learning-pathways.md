---
title: Defining a Learning Pathway
box_type: hands_on
layout: faq
contributors: [shiltemann]
---

[Learning Pathways]({% link learning-pathways/ %}) are sets of tutorials curated by community experts to form a coherent set of lessons around a topic, building up knowledge step by step.

To define a learning pathway, create a file in the `learning-pathways/` folder. An example file is also given in this folder ([pathway-example.md]({% link learning-pathways/pathway-example.md %})). It should look something like this:

```
---
layout: learning-pathway

title: Title of your pathway
description: |
  Description of the pathway. What will be covered, what are the learning objectives, etc?
  Make this as thorough as possible, 1-2 paragraphs. This appears on the index page that
  lists all the learning paths, and at the top of the pathway page
tags: [some, keywords, here ]

cover-image: path/to/image.png # optional cover image, defaults to GTN logo
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

```

And that's it!

We are happy to receive contributions of learning pathways! Did you teach a workshop around a topic using GTN materials? Capture the program as a learning pathways for others to reuse!
