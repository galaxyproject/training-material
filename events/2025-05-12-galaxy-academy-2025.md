---
layout: event

title: Galaxy Training Academy 2025
description: |
  The Galaxy Training Academy is a self-paced online training event for beginners and advanced learners who want to improve their Galaxy data analysis skills.
  Over the course of one week, we offer a diverse selection of learning track for you.

# <button id="program-button" class="btn btn-info" onclick="$('#program-tab').tab('show');">Start the Course!</button>

#cover-image: <div style="background:white">
#              <img src="./events/images/Galaxy_GTA2025_Transparent.png">
#              </div>
cover-image: events/images/Galaxy_GTA2025.png
cover-image-alt: logo for the Galaxy Training Academy Event 2025

registration:
  link: https://forms.gle/vSfms7xDVk4r7dhx8
  deadline: 2025-05-16
  open: true

date_start: 2025-05-12
date_end: 2025-05-16

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

contributions:
    organisers:
        - delphine-l
        - teresa-m
        - scottcain
        - natalie-wa
        - shiltemann
        - dianichj
        - dadrasarmin
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
        - paulzierep
        - pavanvidem
        - plushz
        - poterlowicz-lab
        - pratikdjagtap
        - RZ9082
        - rlibouba
        - reytakop
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

location:
  name: online

infrastructure:
  tiaas: true
  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU
      # tiaas_link:
    - server: https://usegalaxy.org
      name: Galaxy Main
    - server: https://usegalaxy.org.au/
      name: Galaxy AU
    - server: https://usegalaxy.fr
      name: Galaxy FR
      tiaas_link: https://usegalaxy.fr/join-training/gta2025
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
          [<i class="fas fa-video" aria-hidden="true"></i> Video: Welcome to the course!](https://youtu.be/P80S51m6GLk)
      - type: custom
        name: Daily Icebreaker
        description: |
          **Please take a moment to introduce yourself and tell us is your favorite season and why!**

          Post your answer to Slack [#social](https://gtnsmrgsbord.slack.com/channels/social) channel. *(See the setup tab for instructions for joining Slack)*.

  - section: "Galaxy introduction"
    subsection: true
    description: |
      Kick off the week with a hands-on introduction to Galaxy, covering everything from basic navigation and data manipulation to reproducing published analyses, quality control, and mapping sequences to a reference genome. Whether you're new to Galaxy or looking to strengthen your skills, today's sessions will equip you with the foundational knowledge needed for more advanced topics. In the morning you can take part in the Icebreaker by joining us in the [Slack introduction channel](https://gtnsmrgsbord.slack.com/channels/gta_introduction).  Or you can directly start with the tutorials. **Need help with these tutorials?** Ask your questions via the [Slack introduction channel](https://gtnsmrgsbord.slack.com/channels/gta_introduction).
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
          **Tuesday:** For today's icebreaker, we would love to know what is your favourite holiday in your culture and why?

          **Wednesday:** We would love to hear what would you do if you knew you couldn't fail?

          **Thursday:** Tell us what is the best piece of feedback you have ever received?

          Post your answers each day to Slack [#social](https://gtnsmrgsbord.slack.com/channels/social) channel.

  - section: "Tuesday to Thursday: Pick a track"
    subsection: true
    description: |
      Over the course of these three days, you can choose your preferred track and learn how to use Galaxy for data analysis in this research field. If you find multiple topics interesting, feel free to explore more than one track. Each track will guide you through the process, from basic to more advanced analyses, to accommodate learners of all levels.

      For assistance, you can access support through the Slack channel associated with each track module. Please note that while we strive to accommodate all time zones, responses to specific questions may take a bit longer depending on the availability of experts in your time zone.
    tracks:  # Instead of tutorials, you can also define tracks, this will create a button per track that will lead to a different page
      - title: Proteomics
        link: events/tracks/gta2025-proteomics.md
      - title: Assembly
        link: events/tracks/gta2025-assembly.md
      - title: Transcriptomics
        link: events/tracks/gta2025-transcriptomics.md
      - title: Single Cell
        link: events/tracks/gta2025-single-cell.md
      - title: Microbiome
        link: events/tracks/gta2025-microbiome.md
      - title: Climate
        link: events/tracks/gta2025-climate.md
      - title: Machine Learning
        link: events/tracks/gta2025-ml.md
      - title: From Zero to Hero with Python
        link: events/tracks/gta2025-bioNT.md
      - title: Variant Analysis
        link: events/tracks/gta2025-variant-analysis.md

  - section: "Friday: Grab bag"
    description: |
      Can't get enough? Then please pick one of our FAIR or Plant Galaxy tutorials below or any of the tutorials within one of the tracks. The trainers present on Slack will do their best to help you
      if you have a problem and answer questions, but they may not be experts in the topic you selected.
      You can directly start with your preferred tutorial.

      **Need help with these tutorials?** Ask your questions via the Slack Channel [#gta_friday-grab-bag](https://gtnsmrgsbord.slack.com/channels/gta_friday-grab-bag).
    tutorials:
      - type: custom
        name: Daily Icebreaker
        description: |
          **For the last icebreaker, we would love to know what the most interesting or exciting thing you learned this week is!!**

          Post your answer on Slack [#social](https://gtnsmrgsbord.slack.com/channels/social) channel

  - section: Fair training
    subsection: true
    tutorials:
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

  - section: Plant Galaxy
    subsection: true
    tutorials:
      - name: tapscan-streptophyte-algae
        topic: sequence-analysis



---
# Welcome to the Galaxy Training Academy 2025
Do you want to learn how to use Galaxy, an open source data analysis platform? Then you are at the right place. We offer here a **5-day Global Online and Asynchronous learning event**.

**Content of the event**

We provide you with training materials which you can study at your own pace and on your own time throughout the week. Have a look at our [program](https://training.galaxyproject.org/training-material/events/2025-05-12-galaxy-academy-2025.html#program).

<button id="program-button" class="btn btn-info" onclick="$('#program-tab').tab('show');">Program</button>

On the first day, you can make yourself familiar with the Galaxy platform. In the next days, you can learn about Proteomics, Assembly, Transcriptomics, Single Cell, Microbiome or Machine Learning data analysis in Galaxy. On the last day, the Friday we additionally offer a ***FAIR training*** program to you.

**Course format**

The event is asynchronous which means you can set your own pace on your learning journey using our provided self-learning materials. You will start the event by yourself at your preferred time.
Don't worry, asynchronous does not mean that you are alone! If you ever need help, you can contact one of our many trainers worldwide via **Slack chat**.

**How to participate**

First, you will need to register. **We will open registrations on the 3rd of March.**
You only need a browser and an account on one of the Galaxy instances registered for this event. Please have a look at the [setup page](https://training.galaxyproject.org/training-material/events/2025-05-12-galaxy-academy-2025.html#setup).

<button id="program-button" class="btn btn-info" onclick="$('#setup-tab').tab('show');">Setup</button>



**How to get help**

You will not be alone! If you ever need help, you can contact one of our many trainers worldwide via **Slack chat**. Next to the program, you will find Slack channels you can join to exchange with the trainers and other participants during the event. Here you will also find help if you have questions or run into an issue during the training. We try to cover all time zones with helpers for each topic, but please be patient if you do not get an immediate response.


**Certificates**

You will be able to obtain a certificate by the end of the event. More information is coming soon.



## Do you want to join the GTA as a trainer?
Please fill out our [form](https://forms.gle/V9QqSDNg2UmQaDHy8) to indicate in what capacity you would like to help.

