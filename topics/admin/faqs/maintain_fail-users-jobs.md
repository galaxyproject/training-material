---
title: Failing all jobs from a specific user
area: ansible
box_type: tip
layout: faq
contributors: [hexylena, bgruening]
---

This command will let you quickly fail every job from the user 'service-account' (replace with your preferred user)

```
gxadmin tsvquery jobs --user=service-account --nonterminal | awk '{print $1}' |  xargs -I {} -n 1 gxadmin mutate fail-job {} --commit
```
