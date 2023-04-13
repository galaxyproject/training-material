---
layout: learning-pathway
title: Admin Training Course
description: |
  Learn how to setup, configure, and maintain your own Galaxy server. This learning pathway
  will guide you through all the steps required to setup your own Galaxy server with Ansible,
  configuring it, and making it production ready.

tags: [Galaxy administrators, 5-day course]

pathway:
  - section: "Module 1: Setting up Galaxy with Ansible"
    description: This module covers getting a Galaxy server setup with Ansible, a server you will develop furhter in the rest of the modules
    tutorials:
      - name: introduction
        topic: admin
      - name: ansible
        topic: admin
      - name: ansible-galaxy
        topic: admin
      - name: database
        topic: admin
      - name: uwsgi
        topic: admin
      - name: systemd-supervisor
        topic: admin
      - name: production
        topic: admin
      - name: toolshed
        topic: admin
      - name: management
        topic: admin

  - section: "Module 2: Making the server useful"
    description: |
      Here we pivot to focus on making our server useful; adding tools and data,
      configuring quotas and authentication
    tutorials:
      - name: users-groups-quotas
        topic: admin
      - name: external-auth
        topic: admin
      - name: cvmfs
        topic: admin
      - name: connect-to-compute-cluster
        topic: admin

  - section: "Module 3 "
    description: blabla
    tutorials:
      - name: bioblend-api
        topic: dev
      - name: heterogeneous-compute
        topic: admin
      - name: object-store
        topic: admin
      - name: reports
        topic: admin
      - name: gxadmin
        topic: admin
      - name: monitoring
        topic: admin

  - section: "Module 4"
    description: ""
    tutorials:
      - name: maintenance
        topic: admin
      - name: tiaas
        topic: admin
      - name: job-metrics
        topic: admin
      - name: interactive-tools
        topic: admin
      - name: jenkins
        topic: admin
      - name: advanced-galaxy-customisation
        topic: admin

  - section: "Module 5: Grab bag"
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

---

This learning path covers all the topics usually taught during the regular 5-day admin
training workshops. Each module represents what is usually roughly taught in one day during
the events.
