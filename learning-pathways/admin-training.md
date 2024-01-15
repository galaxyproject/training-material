---
layout: learning-pathway
title: Admin Training Course
type: admin-dev
description: |
  Learn how to setup, configure, and maintain your own Galaxy server. This learning pathway
  will guide you through all the steps required to setup your own Galaxy server with Ansible,
  configuring it, and making it production ready.

cover-image: assets/images/gat.png
cover-image-alt: GTN Logo on a spiral galaxy background with text galaxy admin training
editorial_board:
- hexylena
- natefoo
- slugger70

tags: [Galaxy administrators, 5-day course]

pathway:
  - section: "Monday: Setting up Galaxy with Ansible"
    description: This module covers getting a Galaxy server setup with Ansible, a server you will develop furhter in the rest of the modules
    tutorials:
      - name: introduction
        topic: admin
      - name: ansible
        topic: admin
      - name: ansible-galaxy
        topic: admin

  - section: "Tuesday: Making the server useful"
    description: |
      Here we pivot to focus on making our server useful; adding tools and data,
      configuring quotas and authentication
    tutorials:
      - name: backup-cleanup
        topic: admin
      - name: cvmfs
        topic: admin
      - name: reference-genomes
        topic: admin
      - name: apptainer
        topic: admin
      - name: tool-management
        topic: admin
      - name: users-groups-quotas
        topic: admin
      - name: data-library
        topic: admin
      - name: bioblend-api
        topic: dev

  - section: "Wednesday: Clusters"
    description:
    tutorials:
      - name: connect-to-compute-cluster
        topic: admin
      - name: job_conf.xml
        external: true
        link: "https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/config/sample/job_conf.xml.sample_advanced"
        #type: slides ## hands_on
      - name: job-destinations
        topic: admin
      - name: pulsar
        topic: admin

  - section: "Thursday: Expanding"
    description: ""
    tutorials:
      - name: celery
        topic: admin
      - name: tiaas
        topic: admin
      - name: reports
        topic: admin
      - name: gxadmin
        topic: admin
      - name: monitoring
        topic: admin
      - name: sentry
        topic: admin

  - section: "Friday: Grab bag"
    description: |
      Here we have some additional topics, some of which are not admin related.
      Please feel free to pick and choose the tutorials that are interesting for you.
    tutorials:
      - name: troubleshooting
        topic: admin
      - name: tool-integration
        topic: dev
      - name: processing-many-samples-at-once
        topic: galaxy-interface
      - name: upload-rules
        topic: galaxy-interface
      - name: create-new-tutorial
        topic: contributing
      - name: interactive-tools
        topic: admin

---

This learning path covers all the topics usually taught during the regular 5-day admin
training workshops. Each module represents what is usually roughly taught in one day during
the events.
