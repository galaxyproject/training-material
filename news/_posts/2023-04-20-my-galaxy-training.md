---
title: "New Feature: my.galaxy.training"
contributions:
  authorship: [hexylena]
tags: [new feature, gtn]
layout: news
---

The GTN has set up a very simple "redirection" service based on [my.home-assistant.io](https://my.home-assistant.io/) which Helena discovered after reading some Home Assistant documentation and saw a really neat link which led to her own internal home assistant.

So we have built a similar service for the GTN:

1. You write a URL like [https://my.galaxy.training/?path=/workflows/list](https://my.galaxy.training/?path=/workflows/list)
2. A user visits this URL
3. They're asked for which Galaxy server is their 'home' server
4. They will be redirected to that path, on that server, e.g. https://usegalaxy.eu/workflows/list

We will be using this to implement a 'button' that users can click to directly import a workflow on their home server, and as a fallback for the "Click to Run" tools enabled by Tutorial Mode.

Additionally there is a [my.gat.galaxy.training](https://my.gat.galaxy.training) which provides the same service for the Galaxy Admin Training workshop series. There we will use this feature to link admins directly to the VM they are working on when they're working through the training and we need them to e.g. see something in their `/admin` page or in their tool data tables.

## Potential Uses

Regardless of which server a user will use, in the cloud or locally, you can now link to:

- Tools
- Galaxy UI pages:
  - Workflow editor
  - Workflow import
  - TRS import
  - Visualisations
  - etc!

## Try them out

Both of these services can be used from today onwards!

- [my.galaxy.training](https://my.galaxy.training/)
- [my.gat.galaxy.training](https://my.gat.galaxy.training/)

We have already updated some Galaxy Admin training materials to use `my.gat.galaxy.training`.
