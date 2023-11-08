---
title: How do I know what I can do with a role? What variables are available?
area: ansible
box_type: tip
layout: faq
contributors: [hexylena, nsoranzo, slugger70]
---

You don't. There is no standard way for reporting this, but well written roles by trusted authors (e.g. geerlingguy, galaxyproject) do it properly and write all of the variables in the README file of the repository. We try to pick sensible roles for you in this course, but, in real life it may not be that simple.

So, definitely check there first, but if they aren't there, then you'll need to read through `defaults/` and `tasks/` and `templates/` to figure out what the role does and how you can control and modify it to accomplish your goals.
