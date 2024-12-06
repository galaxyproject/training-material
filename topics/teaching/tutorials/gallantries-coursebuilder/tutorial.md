---
layout: tutorial_hands_on

title: Course Builder
subtopic: prepare
time_estimation: 10m
questions:
  - Is this service appropriate for my event?
objectives:
  - Identify if it is appropriate
key_points:
  - Can make hosting Gallantries/Smörgåsbord style courses easier
contributions:
  authorship:
  - hexylena
  funding:
  - gallantries
requirements:
  - type: "internal"
    topic_name: contributing
tags:
  - cyoa
---

As part of the [Gallantries Project](https://gallantries.github.io) we developed a small "Course Builder" which can help you build a schedule of pre-existing videos from the GTN Video library to use with your courses.

> <agenda-title></agenda-title>
>
> In this tutorial, we will see:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Try the Course Builder

You can try the [course builder](https://gallantries.github.io/course-builder/) online. It works as follows:

1. The "Welcome!" Tab

   Pick some modules from the left hand menu by clicking on them. They will turn blue when added to your event.

   > <hands-on-title>Add some modules</hands-on-title>
   > Add the following modules:
   >
   > Basics:
   > 1. Setup
   > 1. Code of Conduct
   > 1. Feedback
   > Sessions:
   > 1. Quality Control
   > 1. Mapping
   > Galaxy Intro:
   > 1. A Very Short Introduction to Galaxy
   > Transcriptomics:
   > 1. Introduction to Transcriptomics
   > 1. Reference based RNA-Seq data analysis
   {: .hands_on}

2. The "Session" Tab

   This tab allows you to organise your modules into sessions. Do you want to split up the materials by day? Or morning/afternoon? Or by week if you're running a longer course?
   
   > <hands-on-title>Create Sessions and Organise Content</hands-on-title>
   > 1. Create a new section
   >
   >    name: Galaxy Intro
   >    description: "Morning Session"
   >
   > 1. Create a new section
   >
   >    name: RNA-Seq
   >    description: "Afternoon Session"
   >
   > 1. Add Very Short Introduction to Galaxy to the first session by dragging and dropping it.
   > 1. Add RNA Seq modules to the second session
   {: .hands_on}

3. The "Configure Event" Tab

   This tab lets you add metadata to the event, it builds on the data in the video-library so you may not find your name in the list.
   
   > <hands-on-title>Update event metadata</hands-on-title>
   > 1. Change the event title and description
   > 1. Change the event dates
   > 1. Select some trainers
   > 1. Select some contacts (can be randomly, or yourself if you are present)
   > 1. Select some affiliated institutions
   {: .hands_on}

3. The "Export" Tab

   If you've made it this far you're ready to contribute this file to the Video Library which will provide a workshop webpage for you!

   > <hands-on-title>Create the event webpage</hands-on-title>
   > 1. Fork the [video-library repo](https://github.com/gallantries/video-library)
   > 1. Add your file to `docs/events/<your-event-identifier>/index.md`
   > 1. Create a pull request
   >    - If any instructors/institutes/organisers were missing from your
   >      selections, then write the missing information in the Pull Request
   >      and the maintainers will help you add it in the correct place.
   {: .hands_on}
