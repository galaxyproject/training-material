---
title: "New Feature: Prometheus Metrics endpoint"
contributions:
  authorship: [hexylena]
tags: [new feature, gtn]
layout: news
---

[Prometheus](https://prometheus.io/) is an increasingly popular monitoring solution whereby sites expose a `/metrics` endpoint, and then a central Prometheus server can scrape it for data which is then plotted or visualised or alerted upon. We at the GTN are now testing out a prometheus endpoint for our completely static site which we'll use to expose some maybe useful statistics about the GTN that others can then scrape and display in ways that are useful for them.

You can see the current GTN metrics endpoint at [/metrics](https://training.galaxyproject.org/metrics)
