---
layout: learning-pathway
title: Tool development for a nice & shiny subdomain
type: admin-dev
description: |
  Discover Galaxy's communities and learn how to  create your subdomain and enrich it by writing, testing and submiting your tools on Galaxy. This learning pathway
  will guide you through all the steps required to build a tool for Galaxy with Planemo for batch tools and how write an interactive tool.
cover-image: assets/images/galaxy_subdomain.png
cover-image-alt: Image of a researcher or developer on a computer thinking of building a community.
editorial_board:
- Marie59

tags: [subdomain, community, tool development, 3-day course]

pathway:

  - section: "Day 1: Set up your subdomain for your community"
    description: This first part explains how to discover the already existing communities (to avoid replication), how to build your subdomain, and finally how to set up your community
    tutorials:
      - name: sig_define
        topic: community
      - name: subdomain
        topic: admin
      - name: sig_create
        topic: community

  - section: "Day 2: Build a batch tool"
    description: This module covers getting your package on Conda, a local Galaxy instance with Planemo, write a Galaxy tool, publish it, and make it visible on a Galaxy server.
    tutorials:
      - name: tool-from-scratch
        topic: dev
      - name: tools_subdomains
        topic: community
      - name: community-tool-table
        topic: community

  - section: "Day 3: Build an interactive tool"
    description: |
      Here we go through how to build a docker image and write the correct wrapper for your interactive tool, and then again make it visible on a Galaxy server.
    tutorials:
      - name: interactive-tools
        topic: dev
      - name: tools_subdomains
        topic: community



---

This learning path covers topics on how to build your subdomain, to enrich it, and create a dynamic community.
