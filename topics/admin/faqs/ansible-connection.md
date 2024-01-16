---
title: Variable connection
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

When the playbook runs, as part of the setup, it collects any variables that are set. For a playbook affecting a group of hosts named `my_hosts`, it checks many different places for variables, including "group_vars/my_hosts.yml". If there are variables there, they're added to the collection of current variables. It also checks "group_vars/all.yml" (for the built-in host group `all`). There is a precedence order, but then these variables are available for roles and tasks to consume.
