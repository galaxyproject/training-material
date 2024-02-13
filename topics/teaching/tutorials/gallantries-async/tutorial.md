---
layout: tutorial_hands_on

title: Asynchronous training
subtopic: practises
draft: true
time_estimation: 1h
questions:
- What is asynchronous training?
- How to deliver asynchronous training using Galaxy?
- What should be prepared for an asynchronous training event?
- How is it different than hybrid training?
objectives:
- Describe asynchronous training
- Organize an asynchronous training event using Galaxy
key_points:
- With hybrid training events, remote instructors pair up with on-site helpers to deliver 1 training across multiple sites simultaneously
- With asynchronous events, instead the contents are prepared ahead of time
- Asynchronous, online training events reduce costs and improve accessibility, even further than hybrid trainings
contributions:
  authorship:
    - bebatut
    - fpsom
    - hexylena
  funding:
    - gallantries
---

**asynchronous training events**, i.e. MOOC style courses pair pre-recorded and written lesson materials with online participants following along at their own pace and convenience.

With this model, we aim to bring training events to the trainees while **eliminate** environmental impact of instructor *and* participant travel.

# What is an asynchronous training event?

*To illustrate a typical Gallantries/Smörgåsbord event and the different roles (1 person can have several roles), we have decided to share a story with several fictional characters*

**Taylor** is a training coordinator of an international research consortium. They recently surveyed researchers in the consortium on their need in HTS data analysis. The results were overwhelming: over 60 scientists at the several locations would like to learn how to analyze their own RNA-seq data. To fulfill this demand, Taylor decided to organize a Smörgåsbord style event: a 3-days online workshop. Taylor is the **global organizer** of this event, in charge of finding the date, contact and coordinate with the teachers, advertise the event, find the instructors, etc. After the event, Taylor will aggregate all feedback from participants, helpers and hosts and share them.

Taylor recruited 2 **helpers** for each day, from amongst the participating incitations. One of the helpers, **Casey**,  has some previous experience with RNA-seq data analysis and Galaxy. Before the workshop, they went through the training material and tested it. They will help participants during the workshop when they are stuck or get different results than the instructors. They will also give direct feedback to the instructor about the pace, any possible local issues etc.

Taylor also recruited 4 **instructors**, **Farah** among them. Farah is a trained bioinformatician and experienced instructor, regularly giving training to scientists on HTS data analysis. She will give the introduction and the first steps of RNA-seq data analysis lessons. Farah is located in Germany and teaches from a room at her institute in front of her computer. **2 weeks** before an event she needs to record her lesson as she would normally teach it during a live training, and provide that video to Taylor. Taylor will generate captions (e.g. YouTube's auto-captions) to Farah who must review and manually correct any captioning mistakes, before returning those to Taylor.

**Alex** is a PhD student in molecular biology based in Greece. They would like to learn about RNA-seq data analysis to be able to analyze the data they generated. They heard about the workshop and wish to join online. As a **participant**, Alex will actively participate in the workshop by running their first RNA-seq data analysis given the instructor's instructions, ask for help from helpers when stuck, will raise their questions on the participant chat, and will give feedback using sticky notes and the dedicated feedback form.

# Cost of an event

With the asynchronous training model, the cost of organizing and participating in the workshop are essentially eliminated: the instructors do not need to travel to the venue, participants can learn wherever is convenient for them, no rooms nor coffee must be organised.

# Event Timeline

Here is an example event timeline we provided for a past Smörgåsbord event:

Date        | Time Left | What's happening
---         | ---       | ---
31 January  | 6 weeks   | You tell us which sessions you want to record
4 February  | 5 weeks   | The final schedule is published
17 February | 4 weeks   | GTN CoFest! Here we will help you create workflows + workflow tests for your tutorials.
28 February | 2 weeks   | **Training Videos Due**. Final Course Materials online
7 March     | 1 week    | Captions deadline & 30-second intro videos due
14 March    | 0         | It begins!

We also provided instructors with an instructors guideline page [tailored to the event](https://gallantries.github.io/video-library/instructors) and more important, a set of recording guidelines [like the Gallantries'](https://gallantries.github.io/video-library/events/smorgasbord2/recording.html) which helped instructors know what was best practice during video recording.


# Workshop Checklists

To help you organize an asynchronous event, we have created some checklists, by timing but also [by role](#checklists-by-role). Most of the items in these checklists are not specific to Gallantries/Smörgåsbord events.

## Before the workshop

### Global organizers

{% include topics/teaching/tutorials/gallantries-async/organizer-before.md %}

### Local helpers

{% include topics/teaching/tutorials/gallantries-async/helper-before.md %}

### Instructors

{% include topics/teaching/tutorials/gallantries-async/instructor-before.md %}

## During the workshop

- Chat online on Slack or another platform
- Enjoy!

### Global organizers

{% include topics/teaching/tutorials/gallantries-async/organizer-during.md %}

### Local helpers

{% include topics/teaching/tutorials/gallantries-async/helper-during.md %}

### Instructors

{% include topics/teaching/tutorials/gallantries-async/instructor-during.md %}

## After the workshop

Debrief
    Sticky notes collection from sites
    share experience back: feedback from sites (extract form from issue in github)
    form feedback from instructors + helpers

### Global organizers

{% include topics/teaching/tutorials/gallantries-async/organizer-after.md %}

### Instructors

{% include topics/teaching/tutorials/gallantries-async/instructor-after.md %}



# Checklists by role

## Global organizers

### Before the workshop

{% include topics/teaching/tutorials/gallantries-async/organizer-before.md %}

### During the workshop

{% include topics/teaching/tutorials/gallantries-async/organizer-during.md %}

### After the workshop

{% include topics/teaching/tutorials/gallantries-async/organizer-after.md %}

## Local helpers

### Before the workshop

{% include topics/teaching/tutorials/gallantries-async/helper-before.md %}

### During the workshop

{% include topics/teaching/tutorials/gallantries-async/helper-during.md %}

## Instructors

### Before the workshop

{% include topics/teaching/tutorials/gallantries-async/instructor-before.md %}

### During the workshop

{% include topics/teaching/tutorials/gallantries-async/instructor-during.md %}

### After the workshop

{% include topics/teaching/tutorials/gallantries-async/instructor-after.md %}
