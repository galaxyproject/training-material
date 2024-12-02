---
title: Credit where it's due: GTN Reviewers in the spotlight
layout: news
tags:
- gtn infrastructure
- new feature
- automation
contributions:
  authorship:
  - hexylena
  - nomadscientist
  infrastructure:
  - hexylena
cover: news/images/reviewing.png
coveralt: A screenshot of the GTN's short introduction to Galaxy tutorial. There are two authors and two editors, but now shown is a new reviewers with 13 individuals, some overlapping with editors and authors.
---

We would like to recognise and thank all of the reviewers who have contributed to the GTN tutorials; your efforts are greatly appreciated, and we are grateful for your contributions to the GTN community. Today, we are highlighting your efforts on every single learning material across the GTN.

@gtn:nomadscientist requested the ability to annotate reviewers of a tutorial. What a great idea! We needed a way to give credit to the reviewers who have contributed to the GTN tutorials, as it is usually a somewhat thankless job, there is not a lot of visibility in reviewing code even though it is an incredibly valuable step in the process of developing (e-)learning materials. We quickly [implemented support](https://github.com/galaxyproject/training-material/commit/bce05249d60f571e72c4508da63433b68d243b59) for that contribution role in the GTN infrastructure. Everyone can now manually annotate when a colleague or coworker reviews a tutorial outside of GitHub (as not everyone is familiar or comfortable reviewing learning materials there!)

However given our extensive automation, the we took that one step further! The GTN has recently implemented a new automation that collects metadata about every pull request that is merged into the GTN. This metadata includes the reviewers of learning materials, so of course we can automatically annotate this on every single material within our codebase, leading to our updated headers including up to dozens of previously uncredited reviewers per tutorial.

Thank you all for your hard work and dedication to the GTN community! 
