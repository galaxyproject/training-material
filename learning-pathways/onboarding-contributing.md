---
layout: learning-pathway
tags: [introduction, contributing, cofest]
type: instructors

priority: 1

editorial_board:
- nomadscientist
- PhilReedData

redirect_from:
- "/learning-pathways/building_tutorials"

title: Onboarding for new contributors
description: |
  This learning path aims to teach you how to contribute to GTN training materials. From editing typos to rendering test GTN sites, you will progress from being a user to a contributor.
  For support throughout these tutorials, join our Galaxy [training network on Matrix](https://matrix.to/#/#Galaxy-Training-Network_Lobby:gitter.im) to ask questions!

pathway:
  - section: "Module 0: Introduction to Galaxy"
    description: |
      If you've never used the Galaxy interface before, you will need to take our introductory tutorial first. Here, you will get a first look at the Galaxy platform for data analysis. You only need the short introduction to enable minor contributions, but for larger contributions, you can gain more familiarity through doing both tutorials.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
      - name: galaxy-intro-101
        topic: introduction

  - section: "Module 1: Edit a small error"
    description: |
      Make a minor correction to an existing training material. You will edit one page using GitHub, so find a typo and get ready!
    tutorials:
      - name: Making a minor correction to the GTN
        type: faq
        link: faqs/gtn/gtn_minor_correction

  - section: "Module 2: Add yourself as a contributor"
    description: |
      While you can make minor corrections using the web-based Github interface, for larger changes (and, indeed, for building new materials altogether), you will need to work more extensively with Github. When you make larger changes, you should be acknowledged for your work - we will therefore use this Github training to also add yourself to our contributors list, so that we can acknowledge you going forward!
    tutorials:
      - name: github-contribution
        topic: contributing
  - section: "Reference: Cool content"
    description: |
        If you want to add icons, images, cool learning boxes and more, you will need to use this next resource. Treat it as the Wikipedia of building training material - it's not meant to be a tutorial, but rather a resource to look up or scan through for ideas.
    tutorials:
      - name: create-new-tutorial-content
        topic: contributing

  - section: "Module 3: Make a larger change"
    description: |
      You will also need to decide on a larger change to make in Galaxy! You may already have ideas on what to fix; you can reference our FAQ for ideas on what to change; or contact a [Community of Practice](https://galaxyproject.org/community/sig/) to see if they have anything that needs doing. We will show you how to visualise those changes, and see how what you do will impact the materials. Finally, we show you had to add yourself as an editor, to acknowledge your contribution!
    tutorials:
      - name: How can I get started with contributing?
        type: faq
        link: faqs/gtn/contributors_getting_started
      - name: rendering_gtn
        topic: contributing
      - name: Using the new Contributions Annotation framework
        type: faq
        link: topics/contributing/faqs/new-contributions

# Ideally, we have a github-Desktop option as a CYOA in the github command line contribution
# Ideally, we then have an FAQ or tutorial on 'Reviewing GTN materials' here as the next Module.
# Where should these go? https://training.galaxyproject.org/training-material/topics/dev/faqs/contributing.html & https://training.galaxyproject.org/training-material/faqs/gtn/contributors_git_advanced.html?

  - section: "Module 4: Make an even bigger change!"
    description: |
     Tutorials often need their workflows updated to the latest set of tools. Or, you might find that you want to want to build something entirely new! Pick the tutorial that is right for you!
    tutorials:
      - name: updating_tutorial
        topic: contributing
      - name: create-new-tutorial
        topic: contributing
      - name: create-new-tutorial-slides
        topic: contributing

  - section: "The End!"
    description: |
      And now you're done! There are still loads of resources to help you improve your training [conceptually](../learning-pathways/train-the-trainers.html) or [structurally](../topics/contributing/).
---

New to making edits or creating training materials? Follow this learning pathway to get familiar with the basics!
