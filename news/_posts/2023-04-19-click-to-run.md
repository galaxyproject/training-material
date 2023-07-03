---
title: "New Feature: Click-to-run Workflows"
contributions:
  authorship: [hexylena]
tags: [new feature, gtn]
layout: news
tutorial: topics/metagenomics/tutorials/mothur-miseq-sop-short/workflows/
---

The GTN has implemented "click to run" workflows! One click, will get you into Galaxy with a workflow imported and ready to run.
If you are on either usegalaxy.org or usegalaxy.eu, we have implemented easy buttons which you can click that will direct you straight to your Galaxy server.

Clicking a single button will bring you directly into the associated Galaxy interface, with the workflow loaded into your account, and ready to run.

We have additionally added a 'generic' method by which you can launch workflows on your preferred server via a redirect.

## How does the magic work?

In the backend the GTN has implemented very limited support for the GA4GH TRS API. Galaxy has support for loading and running workflows directly from compatible endpoints which we leverage to enable this cool feature. You can read more in the [GTN API]({% link api/index.html %}).
