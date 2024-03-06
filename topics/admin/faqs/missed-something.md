---
title: Got lost along the way?
area: utilities
box_type: comment
layout: faq
contributors: [hexylena]
---

{% assign prevStep = include.step | plus: -1 %}
If you missed any steps, you can compare against the [reference files](https://github.com/hexylena/git-gat/tree/step-{{ include.step }}), or see what changed since [the previous tutorial](https://github.com/hexylena/git-gat/compare/step-{{ prevStep }}...step-{{ include.step }}#files_bucket).

If you're using `git` to track your progress, remember to add your changes and commit with a good commit message!
