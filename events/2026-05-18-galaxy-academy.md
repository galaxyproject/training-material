---
layout: event
redirect_from:
  - events/2025-05-18-galaxy-academy

title: Galaxy Training Academy 2026
description: |
  The Galaxy Training Academy is a self-paced online training event for beginners and advanced learners who want to improve their data analysis skills in Galaxy and/or in popular fields in bioinformatics.
  Over the course of one week, we offer a diverse selection of learning tracks for you.
  <button id="program-button" class="btn btn-info" onclick="$('#program-tab').tab('show');">Start the Course!</button>


registration:
  link: https://forms.gle/L8ATuTQuzgsrP1iw9
  deadline: 2026-05-16
  open: true


date_start: 2026-05-18
date_end: 2026-05-22


cost: free
audience: Everyone who would like to get to know Galaxy and learn bioinformatics data analysis, whether you want to master a specific analysis or learn a new skill.
contact_email: academy@galaxyproject.org

async: true
mode: online

tags:
- microbiome
- single-cell
- proteomics
- introduction
- galaxy-interface
- assembly
- statistics
- variant-analysis
- climate
- humanities
- ecology

contributions:
    organisers:
        - delphine-l
        - scottcain
        - natalie-wa
        - shiltemann
    instructors:
        - ahmedhamidawan
        - annasyme
        - annefou
        - anuprulez
        - abretaud
        - nekrut
        - dadrasarmin
        - Nilchia
        - bebatut
        - bgruening
        - clsiguret
        - Sch-Da
        - dannon
        - dianichj
        - deeptivarshney
        - delphine-l
        - elifsu-simula
        - elichad
        - EngyNasr
        - emmaustin20
        - evenmm
        - GarethPrice-Aus
        - hrhotz
        - hvelab
        - igormakunin
        - khaled196
        - jdavcs
        - j34ni
        - jennaj
        - jhahnfeld
        - jochenblom
        - lisanna
        - lfenske-93
        - meltemktn
        - bernt-matthias
        - PfisterMaxJLU
        - foellmelanie
        - meltemktn
        - mschatz
        - hujambo-dunia
        - mirandaembl
        - natalie-wa
        - natefoo
        - Oliver_Rupp
        - oschwengers
        - pauldg
        - PaulineSGN
        - paulzierep
        - pavanvidem
        - plushz
        - poterlowicz-lab
        - pratikdjagtap
        - RZ9082
        - rlibouba
        - reyhaneh-tavakoli
        - SaimMomin12
        - sanjaysrikakulam
        - scottcain
        - silviadg87
        - stephanierobin
        - subinamehta
        - teresa-m
        - timothygriffin
        - tcollins2011
        - nomadscientist
        - wm75
        - yvanlebras
    funding:
        - eurosciencegateway
        - biont
        - nfdi4plants
        - deNBI
        - elixir-europe
        - mwk
        - abromics
        - ifb
        - FAIR2Adapt
        - pndb
        - dataterra
        - gaiadata

location:
  name: online

infrastructure:
  tiaas: false
  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU
    - server: https://usegalaxy.org
      name: Galaxy Main
    - server: https://usegalaxy.org.au/
      name: Galaxy AU
  support:
    platform: Slack


program:
  - section: "Welcome!"
    description: When you are ready to start, just work your way through the program below by watching videos and/or following the text-based tutorials. There are no live sessions, so you can determine your own schedule. A large team of instructors is available on Slack to answer your questions 24/7! Enjoy!


  - section: "Start Here: Course Introduction"
    tutorials:
      - name: Welcome & Course logistics
        type: custom
        description: |
          [<i class="fas fa-video" aria-hidden="true"></i> Video: Welcome to the course!](https://youtu.be/rOCR6Rm5i20)
      - type: custom
        name: Daily Icebreaker
        description: |
          Whenever you start each day, start by answering the daily ice breaker! Post your answers in the Slack  [#social](https://gtnsmrgsbord.slack.com/channels/social) channel.

          **Monday:** Introduce yourself and let us know what you are hoping to learn this week!

          **Tuesday:** Let us know what track are following today and one thing you are excited to learn.

          **Wednesday:** Midweek check-in! What’s one thing you’ve learned so far?

          **Thursday:** What’s one question you still have about Galaxy or your track?

          **Friday:** What’s your GTA2026 win of the week? Big or small, we want to celebrate it!

          Post your answers to Slack [#social](https://gtnsmrgsbord.slack.com/channels/social) channel. *(See the setup tab for instructions for joining Slack)*.

  - section: Daily Meet the Experts
    description: Each day, we will release a new video in this special series. Read more about the speakers in our [news post](/training-material/news/2026/05/11/gta-meet-experts.html)
    tutorials:
      - type: custom
        name: "**Monday:** Delphine Larivière"
        description: |
          Watch the [video](https://youtu.be/90tlBmTTFFY)! (Will release on Monday)
      - type: custom
        name: "**Tuesday:** Aaron Quinlan"
        description: |
           Watch the [video](https://youtu.be/t1_VoudLyxY)! (Will release on Tuesday)
      - type: custom
        name: "**Wednesday:** Stephanie Hicks"
        description: |
           Watch the [video](https://youtu.be/HHsSaZJhhzU)! (Will release on Wednesday)
      - type: custom
        name: "**Thursday:** Timothy Griffin"
        description: |
           Watch the [video](https://youtu.be/ZTqG14VZcaY)! (Will release on Thursday)
      - type: custom
        name: "**Friday:** Gareth Price"
        description: |
           Watch the [video](https://youtu.be/MeazAZ3IHDw)! (Will release on Friday)


  - section: "Monday: Galaxy introduction"
    description: |
      Kick off the week with a hands-on introduction to Galaxy, covering everything from basic navigation and data manipulation to reproducing published analyses, quality control, and mapping sequences to a reference genome. Whether you're new to Galaxy or looking to strengthen your skills, today's sessions will equip you with the foundational knowledge needed for more advanced topics.

      In the morning you can take part in the Icebreaker by joining us in the [Slack #social channel](https://gtnsmrgsbord.slack.com/channels/social)

      **Need help with these tutorials?** Ask your questions via the Slack [#gta_introduction](https://gtnsmrgsbord.slack.com/channels/gta_introduction) channel.

  - section: "New to Galaxy? First introduction to Galaxy"
    subsection: true
    description: Pick **one** of the following two tutorials (depending on your interest) if you are new to Galaxy.
    tutorials:
      - name: galaxy-intro-101
        topic: introduction
      - name: galaxy-intro-101-everyone
        topic: introduction

  - section: "Used Galaxy Before? Quick refresh of Galaxy"
    subsection: true
    description: Have you used Galaxy before but like a little refresher course? We suggest this shorter introduction tutorial.
    tutorials:
      - name: galaxy-intro-short
        topic: introduction

  - section: Go a bit further with Galaxy
    subsection: true
    description: Dive into advanced features of Galaxy such as histories, collections, data manipulation and data logistics.
    tutorials:
      - name: history
        topic: galaxy-interface
      - name: collections
        topic: galaxy-interface
      - name: data-manipulation-olympics
        topic: introduction
      - name: galaxy-intro-ngs-data-managment
        topic: introduction


  - section: "Fundamentals of Sequence analysis"
    subsection: true
    description: If you are planning to follow any genomics themed modules this week (such as assembly, transcriptomics, microbiome, single cell, variant analysis), we suggest you follow these introductory tutorials today.
    tutorials:
      - name: quality-control
        topic: sequence-analysis
      - name: mapping
        topic: sequence-analysis


  - section: "Tuesday to Thursday: Pick your track(s)"
    description: |
      Over the course of these three days, you can choose your preferred track and learn how to use Galaxy for data analysis in this research field. If you find multiple topics interesting, feel free to explore more than one track. Each track will guide you through the process, from basic to more advanced analyses, to accommodate learners of all levels.

      For assistance, you can access support through the Slack channel associated with each track module. Please note that while we strive to accommodate all time zones, responses to specific questions may take a bit longer depending on the availability of experts in your time zone.
    tracks:
      - title: Assembly
        link: events/tracks/gta2026-assembly.md
      - title: Climate
        link: events/tracks/gta2026-climate.md
      - title: Digital Humanities
        link: events/tracks/gta2026-dh.md
      - title: Ecology
        link: events/tracks/gta2026-ecology.md
      - title: From Zero to Hero with Python
        link: events/tracks/gta2026-bioNT.md
      - title: Machine Learning
        link: events/tracks/gta2026-ml.md
      - title: Microbiome
        link: events/tracks/gta2026-microbiome.md
      - title: Proteomics
        link: events/tracks/gta2026-proteomics.md
      - title: Single Cell
        link: events/tracks/gta2026-single-cell.md
      - title: Transcriptomics
        link: events/tracks/gta2026-transcriptomics.md
      - title: Variant Analysis
        link: events/tracks/gta2026-variant-analysis.md


  - section: "Friday: Grab bag"
    description: |
      Can't get enough? You can explore the rest of the [Galaxy Training Network material](/training-material/) and follow any tutorial that interests you!

      The trainers will still be present on Slack for you if you have a problem and answer questions. However, they may not be experts in the topic you selected.

      **Need more time to finish your module?** Of course you can continue with that as well today!

      **Need help with these tutorials?** Ask your questions via Slack!



---
# Welcome to the Galaxy Training Academy 2026
Do you want to learn how to use Galaxy, an open source data analysis platform? Then you are at the right place. We offer here a **5-day Global Online and Asynchronous learning event**.


**Content of the event**

We provide you with training materials which you can study at your own pace and on your own time throughout the week. Have a look at our [program](https://training.galaxyproject.org/training-material/events/2025-05-18-galaxy-academy.html#program).

<button id="program-button" class="btn btn-info" onclick="$('#program-tab').tab('show');">Program</button>

On the first day, you can make yourself familiar with the Galaxy platform. In the next days, you can learn about Proteomics, Assembly, Transcriptomics, Single Cell, Microbiome or Machine Learning data analysis in Galaxy.

**New for 2026**

We are pleased to announce that in addition to the learning tracks, we will also have a "Meet the Expert" video for each day of the course. See the announcement page for it [for more information](/training-material/news/2026/05/11/gta-meet-experts.html).

**Course format**

The event is asynchronous which means you can set your own pace on your learning journey using our provided self-learning materials. You will start the event by yourself at your preferred time.
Don't worry, asynchronous does not mean that you are alone! The trainers and your fellow participants will be present in **Slack chat** for questions and discussions. We will be available at anytime during the event across all time zones.

**How to participate**

Register throught the link at the top of the page.
You only need a browser and an account on one of the Galaxy instances registered for this event.

**How to get help**

You will not be alone! If you ever need help, you can contact one of our many trainers worldwide via **Slack chat**. Next to the program, you will find Slack channels you can join to exchange with the trainers and other participants during the event. Here you will also find help if you have questions or run into an issue during the training. We try to cover all time zones with helpers for each topic, but please be patient if you do not get an immediate response.

**Certificates**

You can get a certificate for each track you completed. Here are the steps to get your certificate:
- In any order:
  -  Fill the anonymous feedback survey for the Galaxy Training Academy: [Feedback Survey](https://gta-feedback.galaxyproject.org/). You will get a unique token that will be use to claim your certificates. Save it preciously, you won't be able to find it again once you close the page. If you loose the token. fill the survey again. 
  - Fill the form(s) for the track(s) you completed.
- Claim the certificates with your email adress and your token: [Claim Certificates](https://gta-feedback.galaxyproject.org/claim-certificates). You can claim the certificates as many times as you want.

You will be able to claim your certificates until June 21st. 

<!-- ## Do you want to join the GTA as a trainer?
Please fill out our [Google Form](https://docs.google.com/forms/d/e/1FAIpQLSdyNhUhWtk2W6nxDwpxao4aRgV9921_PGGyQpGIaRbZynHnAQ/viewform?usp=header) to indicate in what capacity you would like to help.


 -->
