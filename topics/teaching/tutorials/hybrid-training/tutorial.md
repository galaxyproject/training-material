---
layout: tutorial_hands_on

title: Hybrid training
subtopic: practises
enable: false
time_estimation: 1h
questions:
- What is hybrid training?
- How to deliver hybrid training using Galaxy?
- What should be prepared for an hybrid training event?
objectives:
- Describe hybrid training
- Organize an hybrid training event using Galaxy
key_points:
- With hybrid training events, remote instructors pair up with on-site helpers to deliver 1 training across multiple sites simultaneously
- Hybrid training events reduce costs and improve accessibility
contributions:
  authorship:
    - bebatut
    - fpsom
  funding:
    - erasmusplus
---

**Combine with Australia Biocommon paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008715**

**hybrid training events**, i.e. pairing-up on-site helpers with remote instructors across multiple sites simultaneously:

![Hybrid training diagram, a single remote instructor broadcasts training to multiple sites](./images/hybrid_training.png)

With this model, we aim to bring training events to the trainees while reducing the environmental impact of instructor travel.

# What is an hybrid training event?

*To illustrate a typical Gallantries event and the different roles (1 person can have several roles), we have decided to share a story with several fictional characters*

**Max** is a training coordinator of an international research consortium. They recently surveyed researchers in the consortium on their need in HTS data analysis. The results were overwhelming: over 60 scientists at the several locations would like to learn how to analyze their own RNA-seq data. To fulfill this demand, Max decided to organize a Gallantries event: a 3-days workshop with 4 different locations simultaneously (Greece, Estonia, France and Spain). Max is the **global organizer** of this event, in charge of finding the date, contact and coordinate with the local hosts, advertize the event, find the instructors, etc. After the event, Max will aggregate all feedback from participants, helpers and hosts and share them.

Max first contacted 4 **local hosts** (1 at each site):

One of them is **Imani** who is based in Greece. Imani is in charge of finding a suitable room (preferably with computers), checking the local setup, recruiting local helpers, advertizing the event locally to participants, and organizing the local catering and social event. Imani will also give the introduction and wrap-up on each day, and collecting feedback from participants.

Imani recruited 2 **local helpers** for each day. One of the helpers, **Casey**,  has some previous experience with RNA-seq data analysis and Galaxy. Before the workshop, they went through the training material and tested it. They will help participants during the workshop when they are stuck or get different results than the instructors. They will also give direct feedback to the instructor about the pace, any possible local issues etc.

Max also recruited 4 **instructors**, **Farah** among them. Farah is a trained bioinformatician and experienced instructor, regularly giving training to scientists on HTS data analysis. They will teach in the morning of the 2nd day (introduction and the first steps of RNA-seq data analysis). Farah is located in Germany and teaches from a room at their institute in front of their computer. During their session, they will adapt their pace given the feedback they receive from the local helpers, check the status of participants job on a dedicated page, and will also answer questions from participants written on chat.

**Alex** is a PhD student in molecular biology based in Greece. They would like to learn about RNA-seq data analysis to be able to analyze the data they generated. They heard about the workshop and join one site for the workshop close to their institute. As a **participant**, Alex will actively participate in the workshop by running their first RNA-seq data analysis given the instructor's instructions, ask for help from local helpers when stuck, will raise their questions on the participant chat, and will give feedback using sticky notes and the dedicated feedback form.

# Cost of an event

With the hybrid training model, the cost of organizing and participating in the workshop are minimized: the instructors do not need to travel to the venue and the event can be the closest as possible to participants (e.g. in their institute).

We recommend that local hosts organize drinks and coffee for the breaks and, if possible, lunch. These will be the major costs of hosting such an event.

To cover those costs, local hosts can ask for a small participation fee. This will also increase the number of registered participants showing up.

# Workshop Checklists

To help you organize an hybrid event, we have created some checklists, by timing but also [by role](#checklists-by-role). Most of the items in these checklists are not specific to Gallantries events.

## Before the workshop

### Global organizers

{% include {{ page.dir }}organizer-before.md %}

### Local hosts

{% include {{ page.dir }}host-before.md %}

### Local helpers

{% include {{ page.dir }}helper-before.md %}

### Instructors

{% include {{ page.dir }}instructor-before.md %}

## During the workshop

Templates for chat
How chat works
emphasize communication between helpers
put sticky notes there quickly

### Global organizers

{% include {{ page.dir }}organizer-during.md %}

### Local hosts

{% include {{ page.dir }}host-during.md %}

### Local helpers

{% include {{ page.dir }}helper-during.md %}

### Instructors

{% include {{ page.dir }}instructor-during.md %}

## After the workshop

Debrief
    Sticky notes collection from sites
    share experience back: feedback from sites (extract form from issue in github)
    form feedback from instructors + helpers

### Global organizers

{% include {{ page.dir }}organizer-after.md %}

### Local hosts

{% include {{ page.dir }}host-after.md %}

### Local helpers

{% include {{ page.dir }}helper-after.md %}

### Instructors

{% include {{ page.dir }}instructor-after.md %}

# Checklists by role

## Global organizers

### Before the workshop

{% include {{ page.dir }}organizer-before.md %}

### During the workshop

{% include {{ page.dir }}organizer-during.md %}

### After the workshop

{% include {{ page.dir }}organizer-after.md %}

## Local hosts

### Before the workshop

{% include {{ page.dir }}host-before.md %}

### During the workshop

{% include {{ page.dir }}host-during.md %}

### After the workshop

{% include {{ page.dir }}host-after.md %}

## Local helpers

### Before the workshop

{% include {{ page.dir }}helper-before.md %}

### During the workshop

{% include {{ page.dir }}helper-during.md %}

### After the workshop

{% include {{ page.dir }}helper-after.md %}

## Instructors

### Before the workshop

{% include {{ page.dir }}instructor-before.md %}

### During the workshop

{% include {{ page.dir }}instructor-during.md %}

### After the workshop

{% include {{ page.dir }}instructor-after.md %}
