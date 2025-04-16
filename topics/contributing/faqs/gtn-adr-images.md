---
title: "GTN ADR: Image Storage"
area: adr
box_type: tip
layout: faq
contributors: [hexylena, shiltemann, bebatut]
---

[FAQ: What is an ADR?]({% link topics/contributing/faqs/gtn-adr.md %})

## Context and Problem Statement

Contributors to the GTN have image and occasionally datasets they wish to include in the GTN. These datasets are generally quite small (kilobytes) but, are necessary for the understanding of a tutorial.

## Decision Drivers

* We prioritise contributor UX very highly, we cannot ask them to learn multiple systems. Git + Markdown is already enough.
* We wish to be able to sufficiently serve the website offline, with just a clone.

## Considered Options

* Storage in git directly
* In another system (e.g. S3)
* Allowing linked images anywhere on the internet.

## Decision Outcome

Chosen option: "Storage in git directly", because it is the simplest solution that meets our requirements, and doesn't require development we cannot fund, and doesn't risk dead links over time.

### Consequences

* Good, because it is simple and doesn't require additional development.
* Bad, because it will permanently inflate the size of the repository, and it will never decrease. (We can offset this with 

## Pros and Cons of the Options

### Storage in S3

* Good, because it's cheap and well known.
* Bad, because we would need to build a way for users to upload images as part of a GTN tutorial development, and then link to them in markdown.
* Bad, because then the website would not be hostable offline.

### Hotlinking

* Good, because it's easy for contributors
* Bad, because unnecessary impact on someone else's bandwidth
* Bad, because the links will rot over time, images and tutorials will not be able to be followed.
