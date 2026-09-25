---
layout: event
google_form_id: 1789990830
status: wip
title: Workflow4Experimenters 2027
description: 'Analyze your data with Galaxy and the Workflow4Metabolomics infrastructure!
  During this one-week-long on-site course (preceded by half-a-week of webinars),
  you will learn how to use the W4M-promoted Galaxy tools for metabolomics, and analyze
  your own LC-MS, GC-MS or NMR data through tutoring sessions. '
cover-image: W4E_short_logo.jpg
cover-image-alt: W4E official logo, which is a green square with W4E written in black and white

tags:
- metabolomics

external: https://workflow4metabolomics.github.io/website/W4E/w4e2027.html

contributions:
  organisers:
  - workflow4metabolomics
  - bidieme
  - cdalle036
  - cdelporte
  - cecilecanlet
  - dcenteno
  - fsouard
  - hechth
  - isabelle802
  - lecorguille
  - melpetera
  - mtremblayfr
  - RJMW
  - sylvainchereau
  - yguitton
  instructors:
  - cecilecanlet
  - isabelle802
  - melpetera
  - mtremblayfr
  - RJMW
  - yguitton
  funding:
  - rfmf
  - metabohub
  - ifb

date_start: 2027-04-05
date_end: 2027-04-30

cost: 1300€ for academic and 2500€ for private institution (to cover expenses for trainers, organization, materials and meals); special discount of 100€ is offered for [RFMF](https://www.rfmf.fr/) members 2027.
audience: This event is intended for anyone interested in learning how to process metabolomic data in Galaxy, including PhD students, technicians, engineers, post-docs, scientists. As the school is based on a bring-your-own-data format, best benefit is achieved when participants do have data to bring to the training at the time of the school.
contact_email: workflow4metabolomics@proton.me
async: false
mode: onsite

registration:
  link: to be opened soon (october)
  deadline: 2027-02-06
  open: false

location:
  name: Pole Numerique Rennes Beaulieu (PNRB)
  address: 263 Av. du General Leclerc
  city: Rennes
  country: France
  postcode: 35042
  geo:
    lat: 48.1143556
    lon: -1.6405555


infrastructure:
  tiaas: false

  servers:
    - server: https://workflow4metabolomics.usegalaxy.fr/
      name: "Galaxy FR - W4M subdomain"
      tiaas_link: 

  custom:
    description: |
      Before joining the in-person part of the course, please make sure to:
        - Bring a laptop with wifi feature and at least an web navigator installed.
        - Answer your tutors' e-mails to ensure your data is ready for tutoring sessions.
        - Have a look at the homework given during the webinar session (in particular the suggested GTN tutorials listed bellow).
        - Galaxy Homework: [A short introduction to Galaxy](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-short/tutorial.html)
        - MS pre-processing Homework: [Mass spectrometry: LC-MS preprocessing with XCMS ](https://training.galaxyproject.org/training-material/topics/metabolomics/tutorials/lcms-preprocessing/tutorial.html)
        - Data processing Homework: [Mass spectrometry: LC-MS data processing ](https://training.galaxyproject.org/training-material/topics/metabolomics/tutorials/lcms-dataprocessing/tutorial.html)

program:
  - section: "Webinar Monday"
    description: Introduction to Galaxy
    tutorials:
      - type: custom
        name: Introduction to W4E2027
        time: "9:30 - 10:00"
        description: Welcome
      - type: custom
        name: Galaxy basics
        time: "10:00 - 11:30"
        description: "Galaxy session: initiation and data management - by M. Petera"
      - type: custom
        name: Know your data
        time: "11:30 - 12:00"
        description: "A short overview of what MS data is - by I. Schmitz"

  - section: "Webinar Tuesday"
    description: "Data processing - Paralel session for MS and NMR"
    tutorials:
      - type: custom
        name: MS session
        time: "9:30 - 11:30"
        description: "MS data pre-processing"
      - type: custom
        name: RMN session
        time: "9:30 - 11:30"
        description: "NMR data pre-processing"
      - type: custom
        name: Data processing
        time: "11:30 - 12:00"
        description: Filtering, signal drift and batch effect correction, quality control...

  - section: "Webinar Wednesday"
    description: Statistics
    tutorials:
      - type: custom
        name: Statistics
        time: "9:30 - 12:00"
        description: Yes, only statistics the whole morning

  - section: "Webinar Thursday"
    description: "Annotation - Paralel session for MS and NMR"
    tutorials:
      - type: custom
        name: MS session
        time: "9:30 - 12:00"
        description: MS data annotation
      - type: custom
        name: RMN session
        time: "9:30 - 12:00"
        description: NMR data annotation

  - section: "Webinar Friday"
    description: Advanced session from an invited speaker
    tutorials:
      - type: custom
        name: To be determined based on participant votes
        time: "9:30 - 12:00"
        description: To be determined based on participant votes

  - section: "In-person week: Tutoring session!"
    description:  |
      Tutoring on your own data, with also time for some well-deserved drinks and socializing!
      Also includes a keynote presentation and an advanced session (to be defined).

---

Analyze your data with Galaxy and the Workflow4Metabolomics infrastructure!

The Workflow4Experimenters 2027 session will take place in April 2027. 
During this one-week course (entirely in English), participants will learn how to use the W4M infrastructure and analyze their own LC-MS, GC-MS, or NMR data.

For this new session, we continue with the acclaimed two-step format:

- From 5th to 9th of April: online theoretical sessions (methods and tools). This session is *only accessible to person who will be present in Rennes for the in-person session.*
- From April, 26th to 30th: tutoring on your own data, at [PNRB](https://formation-continue.univ-rennes.fr/pole-sciences-et-ingenierie) ([Rennes University](https://goo.gl/maps/9AN5Kq24mCNBiTwi7), France).


**The pre-registration is soon to be opened. As soon as it opens, you will find here the link to the pre-registration form.**

**Scientific comitee**: C. Delporte (ULB, Bruxelles), C. Dalle (U.DAB IRBA, Brétigny-sur-Orge), Y. Guitton (Laberca Oniris/INRAE, Nantes), D. Centeno & M. Pétéra (PFEM INRAE, Clermont-Ferrand), G. Le Corguillé (Abims, Roscoff), B. Diémé (PFEM Université Clermont Auvergne), F. Souard (ULB, Bruxelles & Université de Grenoble), C. Canlet, M. Tremblay-Franco (Toxalim INRAE, Toulouse), I. Schmitz (PBS, CNRS, Rouen), S. Chéreau (INRAE, Rennes), H. Hecht (RECETOX, Brno), R.Weber (University of Birmingham)

**Contact**: workflow4metabolomics@proton.me

**Sponsor**: W4M is jointly developed and maintained by the French Bioinformatics Infrastructure (IFB, ELIXIR-FR), the French Infrastructure for Metabolomics and Fluxomics (MetaboHUB) and the Réseau Francophone de Métabolomique et Fluxomique (RFMF).
