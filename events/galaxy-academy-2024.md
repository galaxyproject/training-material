---
layout: event

title: Galaxy Training Academy 2024
description: |
  The Galaxy Academy is a self-paced online training event for beginners as well as learners who would like to improve their Galaxy data analysis skills. Over the course of one week, we will have a different topic and focus every day.

  <button id="program-button" class="btn btn-info" onclick="$('#program-tab').tab('show');">Start the Course!</button>

cover-image: events/images/galaxy-academy-logo.png
cover-image-alt: logo for the Galaxy Academy event consisting of a laptop surrounded by illustrations of DNA molecules

registration:
  link: https://forms.gle/cxzVatt7MAgiMX12A
  deadline: 2024-10-04
  open: true

date_start: 2024-10-07
date_end: 2024-10-11

cost: free
audience: Everyone who would like to get to know Galaxy and learn bioinformatics data analysis, whether you want to master a specific analysis or learn a new skill.
contact_email: academy@galaxyproject.org

async: true
mode: online

contributions:
    organisers:
        - teresa-m
        - natalie-wa
        - nakucher
        - erxleben
        - annasyme
        - nekrut
        - dannon
        - delphine-l
        - GarethPrice-Aus
        - jdavcs
        - mschatz
        - shiltemann
    instructors:
        - ahmedhamidawan
        - erxleben
        - annasyme
        - anuprulez
        - abretaud
        - bebatut
        - bgruening
        - clsiguret
        - dannon
        - dianichj
        - deeptivarshney
        - delphine-l
        - elichad
        - EngyNasr
        - emmaustin20
        - GarethPrice-Aus
        - igormakunin
        - jdavcs
        - lldelisle
        - bernt-matthias
        - foellmelanie
        - mschatz
        - natalie-wa
        - natefoo
        - paulzierep
        - pavanvidem
        - plushz
        - pratikdjagtap
        - RZ9082
        - rlibouba
        - SaimMomin12
        - stephanierobin
        - subinamehta
        - teresa-m
        - timothygriffin
        - tcollins2011
        - nomadscientist
        - wm75
    funding:
        - eurosciencegateway
        - biont
        - nfdi4plants
        - deNBI
        - by-covid
        - elixir-europe
        - mwk
        - abromics
        - ifb


location:
  name: online

infrastructure:
  tiaas: true
  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU
      tiaas_link: https://usegalaxy.eu/join-training/gta2024
    - server: https://usegalaxy.org
      name: Galaxy Main
    - server: https://usegalaxy.org.au/
      name: Galaxy AU

  support:
    platform: Slack


program:
  - section: "Monday: Introduction"
    description: When you are ready to start, just work your way through the program below by watching videos and/or following the text-based tutorials. There are no live sessions, so you can determine your own schedule. A large team of instructors is available on Slack to answer your questions 24/7! Enjoy!
  - section: Course introduction
    subsection: true
    tutorials:
      - name: Welcome & Course logistics
        type: custom
        description: |
          [<i class="fas fa-video" aria-hidden="true"></i> Video: Welcome to the course!](https://youtu.be/OyMpSNEDyEA)
      - type: custom
        name: Daily Icebreaker
        description: |
          **Please take a moment to introduce yourself and tell us one fun fact about yourself!**

          Post your answer to Slack [#social](https://gtnsmrgsbord.slack.com/channels/social) channel. *(See the setup tab for instructions for joining Slack)*.


  - section: "Galaxy introduction"
    subsection: true
    description: |
      Kick off the week with a hands-on introduction to Galaxy, covering everything from basic navigation and data manipulation to reproducing published analyses, quality control, and mapping sequences to a reference genome. Whether you're new to Galaxy or looking to strengthen your skills, today's sessions will equip you with the foundational knowledge needed for more advanced topics.
      In the morning you can take part in the Icebreaker by joining us in the [Slack introduction channel](https://gtnsmrgsbord.slack.com/archives/C07NKAJ8THA). Or you can directly start with the tutorials.

      **Need help with these tutorials?** Ask your questions via the [Slack introduction channel](https://gtnsmrgsbord.slack.com/archives/C07NKAJ8THA).
  - section: Start to get to know Galaxy
    subsection: true
    tutorials:
      - name: galaxy-intro-101
        topic: introduction
      - name: history
        topic: galaxy-interface
      - name: collections
        topic: galaxy-interface
      - name: data-manipulation-olympics
        topic: introduction
      - name: galaxy-intro-ngs-data-managment
        topic: introduction
  - section: Quick start or fresh up your Galaxy knowledge
    subsection: true
    tutorials:
      - name: galaxy-intro-short
        topic: introduction
  - section: "Fundamentals of Sequences analysis"
    subsection: true
    tutorials:
      - name: quality-control
        topic: sequence-analysis
      - name: mapping
        topic: sequence-analysis

  - section: "Tuesday to Thursday: Choose your own Adventure!"
    tutorials:
      - type: custom
        name: Daily Icebreakers
        description: |
          **Tuesday:** For today's ice breaker, we would love to know one weird fact that you know for no reason.
          **Wednesday:** We would love to hear where you find inspiration. Maybe you find inspirtaiton through nature or maybe you have a prominant role model in your life—we'd love to learn more about you!
          **Thursday:** if you could take any one movie prop from a movie set, what would it be?

          Post your answers each day to Slack [#social](https://gtnsmrgsbord.slack.com/channels/social) channel.


  - section: "Pick a track"
    subsection: true
    description: |
      Over the course of these three days, you can choose your preferred track and learn how to use Galaxy for data analysis in this research field. If you find multiple topics interesting, feel free to explore more than one track. Each track will guide you through the process, from basic to more advanced analyses, to accommodate learners of all levels.

      For assistance, you can access support through the Slack channel associated with each track module. Please note that while we strive to accommodate all time zones, responses to specific questions may take a bit longer depending on the availability of experts in your time zone.
    tracks:  # instead of tutorials, you can also define tracks, this will create a button per tracks that will lead to a different page
      - title: Proteomics
        link: events/tracks/gta2024-proteomics.md
      - title: Assembly
        link: events/tracks/gta2024-assembly.md
      - title: Transcriptomics
        link: events/tracks/gta2024-transcriptomics.md
      - title: Single Cell
        link: events/tracks/gta2024-single-cell.md
      - title: Microbiome
        link: events/tracks/gta2024-microbiome.md
      - title: Bacterial Genomics
        link: events/tracks/gta2024-bacterial-genomics.md
      - title: BY-COVID
        link: events/tracks/gta2024-bycovid.md
      - title: Machine Learning
        link: events/tracks/gta2024-ml.md

  - section: "Friday: Grab bag"
    description: |
      Can't get enough? Then please pick one of our FAIR tutorials below, or any of the tutorials of the [GTN](https://training.galaxyproject.org). Please be aware that only trainings that are part of the introduction day or a learning track have been tested on all instances for the event. The trainers present on Slack will do their best to help you if you have a problem and answer questions, but they may not be expert in the topic you selected.
      You can directly start with you prefered tutorial.

      **Need help with these tutorials?** Ask your questions via the Slack Channel [#gta_friday-grab-bag](https://gtnsmrgsbord.slack.com/archives/C07N2A4HQ15).
  - section: Fair training
    subsection: true
    tutorials:
      - type: custom
        name: Daily Icebreaker
        description: |
          **For the last ice breaker, we would love to know what the most interesting or exciting thing you learned this week is!!**

          Post your answer on Slack [#social](https://gtnsmrgsbord.slack.com/channels/social) channel
      - type: custom
        name: "[An overview of the RO-Crate concept and its implementations](https://gallantries.github.io/video-library/videos/ro-crates/intro/slides/)"
        description: |
          [<i class="fas fa-video" aria-hidden="true"></i> Lecture Video](https://gallantries.github.io/video-library/videos/ro-crates/intro/slides/)
      - name: ro-crate-intro
        topic: fair
      - type: custom
        name: "[Registering Galaxy workflows in WorkflowHub](https://gallantries.github.io/video-library/videos/ro-crates/workflowhub/tutorial/)"
        description: |
          [<i class="fas fa-video" aria-hidden="true"></i> Lecture Video](https://gallantries.github.io/video-library/videos/ro-crates/workflowhub/tutorial/)
      - name: ro-crate-galaxy-best-practices
        topic: fair

  - section: Explore the GTN!
    subsection: true
    description: |
      Feel free to try any other tutorials you see on the GTN, and we will do our best to guide you!
    tracks:
      - title: Go to the GTN!
        link: https://training.galaxyproject.org

---
# Welcome to the Galaxy Training Academy
Do you want to learn how to use Galaxy, a open source data analysis platform? Then you are at the right place. We offer here a **5-day Global Online and Asynchronous learning event**.

We provide you with training materials which you can study at your own pace and on your own time throughout the week. Don't worry, asynchronous does not mean that you are alone! If you ever need help, you can contact one of our many trainers worldwide via **Slack chat** .


On the first day you can make yourself familiar with the Galaxy platform. In the next days you can follow different tracks, please go to the program tab for more informaiton.

You can set your own pace on your learning journey using our provided self-learning materials. Next to the program, you will find Slack channels you can join to exchange with the trainers and other participants during the event. Here you will also find help if you have qustions or run into an issue during the training. We try to cover all time zones with helpers for each topic, but please be patient if you do not get an immediat response.

You only need a browser and an account on one of the galaxy instances registered for this event. Please have a look at the setup tab.

**Certificates:** You can earn a certificate for the tracks you attended by answering a short questionnaire to verify your participation.
