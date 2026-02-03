---
layout: event
title: "Galaxy Beyond Basics: Mastering Workflows, Automation, and Scalability"

# external: "https://galaxyproject.org/events/2025-03-10-galaxy-workshop-freiburg/"

description: |
    Join us for an **intensive, week-long, in-person training** designed to elevate your Galaxy expertise to new heights. This workshop is tailored for **data scientists, advanced Galaxy users, and team leaders** who need to **scale, automate, and publish** their data analysis workflows for **batch processing and production-level applications**.

cover-image: events/images/Galaxy_AGT2026.png
cover-image-alt: logo for the Advanced Galaxy Training

date_start: 2026-10-12
date_end: 2026-10-16


cost: 850-900 euros (To be confirmed)
audience: |
    This course is designed for **data scientists, advanced Galaxy users, and team leaders** who need to scale up their data analysis workflows for **batch processing, automation, or production-level applications**. Whether you're looking to optimize workflows, integrate external tools, or enhance reproducibility and scalability, this training will provide the skills to leverage Galaxy effectively. The sessions will be conducted in **French**, while the training materials (slides) will be in **English**.
contact_email: berenice.batut@france-bioinformatique.fr

# async: false
mode: onsite

registration:
  link: 
  deadline: 2026-05-31
  open: false

tags:
- galaxy-interface
- workflow
- command-line
- advanced

contributions:
    organisers:
        - abretaud
        - bebatut
        - lecorguille
        - scorreard
        - Marie59
        - lleroi
        - loraine-gueguen
        - clsiguret
        - rlibouba
        - fmareuil
    instructors:
        - abretaud
        - bebatut
        - lecorguille
        - scorreard
        - Marie59
        - lleroi
        - rlibouba
        - fmareuil
    #testing:
    funding:
        - ifb

location:
  name: Institut Pasteur
  address: 25-28 Rue du Dr Roux
  postcode: 75015
  city: Paris
  country: France
  geo:
    lat: 48.8408253
    lon: 2.3081039

infrastructure:
  tiaas: false
  servers:
    - server: https://usegalaxy.eu
      name: Galaxy EU
    - server: https://usegalaxy.fr
      name: Galaxy FR


program:
  - section: "Monday: Introduction & Workflow development" 
    description: |
      This half-day session starts with a **welcome, icebreaker, and participant introductions** to set a collaborative tone. It continues with a **brief overview of Galaxy** and its workflow capabilities, followed by an **introduction to workflow development**. Participants will **design clean, efficient workflows**, **customize them using parameters**, and **generate user-friendly workflow reports**—blending theory with hands-on practice.
    tutorials:
      - type: custom
        name: Welcome
        time: "14:00"
      - type: custom
        name: Icebreaking
        time: "14:15"
      - type: custom
        name: Participant presentation
        time: "14:30"
      - type: custom
        name: Galaxy overview
        time: "14:50"
      - type: custom
        name: Break
        time: "15:00"
      - type: custom
        name: Galaxy Workflow introduction
        time: "15:30"
      - name: workflow-editor
        topic: galaxy-interface
        time: "15:40"
      - name: workflow-parameters
        topic: galaxy-interface
        time: "16:10"
      - name: workflow-reports
        topic: galaxy-interface
        time: "16:25"
      - type: custom
        name: End of the day
        time: "17:00"

  - section: "Tuesday: Workflow FAIRification, documentation and export"
    description: |
      The second day begins with a **recap of the previous day's key concepts**, followed by an introduction to **UseGalaxy.fr** and its unique features. Participants will then dive into **Workflow FAIRification**, learning to annotate workflows with metadata, apply best practices for consistency, and implement tests to ensure reliability—culminating in publishing workflows to **WorkflowHub** and **Dockstore** via the IWC. The session continues with **Workflow Documentation**, where attendees will create high-resolution workflow images and develop interactive tutorials using a "Choose Your Own Tutorial" approach. Finally, the day covers **Workflow Export**, focusing on creating **RO-Crates** for reproducibility, evaluating their completeness, and submitting workflows to **LifeMonitor** for performance tracking—equipping participants with the skills to share, document, and monitor their workflows effectively.
    tutorials:
      - type: custom
        name: Recap of the previous day
        time: "9:00"
      - type: custom
        name: UseGalaxy.fr
        time: "9:15"
      - name: workflow-fairification
        topic: galaxy-interface
        time: "9:30"
      - type: custom
        name: Break
        time: "10:30"
      - name: workflow-fairification
        topic: galaxy-interface
        time: "11:00"
      - type: custom
        name: Lunch
        time: "12:00"
      - name: workflow-posters
        topic: galaxy-interface
        time: "13:30"
      - name: create-new-tutorial
        topic: contributing
        time: "14:00"
      - type: custom
        name: Break
        time: "15:00"
      - type: custom
        name: Overview of the RO-Crate concept and its implementations
        time: "15:30"
      - name: ro-crate-in-galaxy
        topic: fair
        time: "16:00"
      - name: ro-crate-submitting-life-monitor
        topic: fair
        time: "16:30"
      - type: custom
        name: End of the day
        time: "17:00"

  - section: "Wednesday: Scaling workflow and Galaxy using command-line"
    description: |
      The third day starts with a **recap of the previous sessions**, followed by inspiring examples of **large-scale projects using Galaxy** to showcase real-world applications. Participants will then explore **Workflow Scaling**, learning to **execute workflows from the command line** using Planemo and **automate batch processing** with shell scripts, while analyzing performance for efficiency. The session will also cover **Scaling Galaxy Use with BioBlend**, where attendees will **programmatically interact with Galaxy** using Python, design scripts for batch workflow execution, and evaluate scalability compared to manual processes—empowering them to optimize workflows for large-scale data analysis. The day will end with explanation of the "Bring Your Own Work".
    tutorials:
      - type: custom
        name: Recap of the previous day
        time: "9:00"
      - type: custom
        name: Big projects using Galaxy
        time: "9:15"
      - name: workflow-automation
        topic: galaxy-interface
        time: "9:30"
      - type: custom
        name: Break
        time: "10:30"
      - name: workflow-automation
        topic: galaxy-interface
        time: "11:00"
      - type: custom
        name: Lunch
        time: "12:00"
      - name: bioblend-api
        topic: dev
        time: "13:30"
      - type: custom
        name: Break
        time: "15:00"
      - name: bioblend-api
        topic: dev
        time: "15:30"
      - type: custom
        time: "16:00"
        name: Django to Galaxy?
        time: "16:30"
      - type: custom
        name: End of the day
        time: "17:00"
        
  - section: "Thursday: Bring Your Own Work (BYOW)"
    description: |
      The fourth day is dedicated to **"Bring Your Own Work"**, offering participants hands-on time to **apply the skills and tools** learned throughout the training to their own projects. With guidance from trainers, attendees will **refine their workflows**, **troubleshoot challenges**, and **implement solutions** using their personal data. The session encourages collaboration and peer support, allowing participants to **document their progress**, optimize their workflows, and showcase their results—ensuring they leave with practical, actionable outcomes for their research.
    tutorials:
      - type: custom
        name: Group work
        time: "9:00"
      - type: custom
        name: Break
        time: "10:15"
      - type: custom
        name: Group work
        time: "10:45"
      - type: custom
        name: Lunch
        time: "12:00"
      - type: custom
        name: Group work
        time: "13:30"
      - type: custom
        name: Break
        time: "15:00"
      - type: custom
        name: Group work
        time: "15:30"
      - type: custom
        name: Bring Your Own Work group recap
        time: "16:30"
      - type: custom
        name: End of the day
        time: "17:00"

  - section: "Friday: Storage, data, general recap, and supplementary exercises"
    description: |
      The final half-day begins with a **recap of the previous day's progress**, followed by an introduction to **"Bring Your Own Storage"**, exploring how to integrate personal or institutional storage solutions with Galaxy. Participants will then learn about **managing databases in Galaxy** and the **IDC (Intergalactic Data Commission) effort**, gaining insights into efficient data organization and accessibility. The session concludes with a **general recap**, supplementary exercises to reinforce key concepts, and **feedback and closing remarks**, ensuring participants leave with a comprehensive understanding and resources for continued success.
    tutorials:
      - type: custom
        name: Recap of the previous day
        time: "9:00"
      - type: custom
        name: Bring Your Own storage
        time: "9:30"
      - type: custom
        name: Managing databases in Galaxy and the IDC effort
        time: "9:45"
      - type: custom
        name: Break
        time: "10:00"
      - type: custom
        name: General recap & Supplementary exercises
        time: "10:30"
      - type: custom
        name: Feedback
        time: "11:30"
      - type: custom
        name: Closing remarks
        time: "11:45"
      - type: custom
        name: End of the workshop
        time: "12:00"


---

Over five days, you'll embark on a **comprehensive journey** through Galaxy's advanced capabilities:

- **Monday: Introduction & Workflow Development**

  Start with a **welcome and icebreaker** to foster collaboration, followed by a **brief overview of Galaxy** and its workflow features. Dive into **hands-on workflow development**, where you'll learn to design **clean, efficient workflows**, customize them with parameters, and generate **user-friendly workflow reports**—combining theory with practical application.

- **Tuesday: Workflow FAIRification, Documentation, and Export**

  Begin with a **recap of Day 1**, then explore **UseGalaxy.fr** and its unique features. Learn to **annotate workflows with metadata**, apply best practices for **FAIR compliance**, and implement tests to ensure reliability. Publish your workflows to **WorkflowHub and Dockstore** via the IWC. Develop **high-resolution workflow visualizations** and create **interactive tutorials** using a "Choose Your Own Tutorial" approach. Finally, master **workflow export** by creating **RO-Crates** for reproducibility and submitting workflows to **LifeMonitor** for performance tracking.

- **Wednesday: Scaling Workflows & Galaxy Using Command-Line and API**

  Start with a **recap and real-world examples** of large-scale Galaxy projects. Learn to **execute workflows from the command line** using Planemo, **automate batch processing** with shell scripts, and analyze performance for efficiency. Discover how to **scale Galaxy use with BioBlend**, designing Python scripts for batch workflow execution and evaluating scalability. The day concludes with an introduction to the **"Bring Your Own Work"** session.

- **Thursday: Bring Your Own Work (BYOW)**

  Dedicate the day to **applying your new skills** to your own projects. With guidance from trainers, **refine your workflows**, **troubleshoot challenges**, and **implement solutions** using your personal data. Collaborate with peers, **document your progress**, and optimize your workflows to leave with **actionable results** for your research.

- **Friday: Storage, Data Management, Recap, and Closing**

  The final half-day begins with a **recap of the week's progress**, followed by a session on **"Bring Your Own Storage"**, exploring how to integrate personal or institutional storage with Galaxy. Learn about **managing databases in Galaxy** and the **IDC (Intergalactic Data Commission) effort** for efficient data organization. The workshop concludes with a **general recap**, supplementary exercises, and **feedback and closing remarks**, ensuring you leave with a **comprehensive understanding** and resources for continued success.

This training will be conducted in **French**, while the materials (slides) will be in **English**.

## Learning Objectives

At the end of the workshop, you will be able to:

- **Workflow development**
  - **Understand** the key aspects of workflows by **identifying** their core components and purpose.
  - **Create** clean, non-repetitive workflows by **applying** best practices for process design.
  - **Use** workflow parameters to **customize** and **optimize** workflows for specific tasks.
  - **Generate** user-friendly workflow reports to **display** workflow results in a structured way.
- **Workflow FAIRyfication**
  - **Annotate** a Galaxy workflow with essential metadata to ensure it is **findable** and **reusable**.
  - **Apply** best practices to data analysis workflows to improve **consistency** and **interoperability**.
  - **Implement** robust tests to **validate** workflow reliability and accuracy.
  - **Publish** a Galaxy workflow on WorkflowHub and Dockstore via its integration into the IWC, demonstrating enhanced **findability, accessibility, interroperability and usability** for the scientific community.
- **Workflow Documentation**
  - **Design** a high-resolution workflow image optimized for documentation and presentations.
  - **Develop** a hands-on tutorial with a "Choose Your Own Tutorial" approach, including:
    - A **step-by-step tutorial** with skeleton generation from the workflow.
    - A **real-time tutorial** that runs and explains the workflow interactively.
  - **Produce** a final documentation package that includes both tutorial formats and high-resolution visuals.
- **Workflow Export**
  - **Apply** the process of **creating a Galaxy Workflow Run RO-Crate** by packaging a workflow with its metadata, inputs, and outputs, ensuring it is reproducible and FAIR-compliant.
  - **Evaluate** the completeness and accuracy of a Galaxy Workflow Run RO-Crate by reviewing its structure, metadata, and included files for adherence to best practices.
  - **Submit** a workflow to **LifeMonitor**, **analyzing** the platform's feedback to assess workflow performance and improve its reliability for future use.
- **Workflow Scaling using command-line**
  - **Execute** workflows from the command line using the **Planemo `run` subcommand**, demonstrating the ability to run and monitor workflows outside the Galaxy interface.
  - **Develop** simple shell scripts to **automate the execution of multiple workflows** concurrently or sequentially, optimizing efficiency and scalability.
  - **Analyze** the performance and resource usage of workflows run via shell scripts, **evaluating** the effectiveness of scaling strategies for large-scale data processing.
- **Scaling Galaxy Use with the API and BioBlend**
  - **Utilize** the **BioBlend library** to programmatically interact with Galaxy, **executing** workflows, managing datasets, and automating repetitive tasks.
  - **Design** a Python script using BioBlend to **scale Galaxy workflows** for batch processing, ensuring efficient resource use and reproducibility.
  - **Evaluate** the performance and scalability of workflows executed via BioBlend, **comparing** results with manual Galaxy interactions to identify improvements.
- **"Bring Your Own Work"**
  - **Apply** the concepts and tools learned during the training to **develop or refine your own workflows** using your personal data, with guidance from trainers.
  - **Troubleshoot** challenges in your workflow or data analysis, **implementing** solutions with the support of trainers and peers.
  - **Demonstrate** progress in your project by **documenting** your workflow, results, and any optimizations made during the sessions.

## Requirements

- Prior knowledge and experience using Galaxy
- Prior knowledge and experience using command line
- Fluent in French (materials will be in English and discussions will happen in French)

## Registration

This **in-person-only workshop** is limited to **20 participants** to ensure personalized attention and hands-on support. The course will proceed with a **minimum of 10 registrations**. Registration will open in **April 2026**—stay tuned for updates!




