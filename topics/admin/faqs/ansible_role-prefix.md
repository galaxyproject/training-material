---
title: What is the difference between the roles with `role:` prefix and without?
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

The bare role name is just simplified syntax for the roles, you could equally specifiy `role: <name>` every time but it's only necessary if you want to set additional variables like `become_user`
