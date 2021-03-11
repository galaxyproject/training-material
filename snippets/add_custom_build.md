---
title: Adding a custom database/build (dbkey)
description: Galaxy may have several reference genomes built-in, but you can also create your own.
area: datasets
box_type: tip
layout: faq
---

- In the top menu bar, go to the *User*, and select *Custom Builds*
- Choose a name for your reference build `{{ include.name }}`
- Choose a dbkey for your reference build `{{ include.dbkey }}`
- Under **Definition**, select the option `FASTA-file from history`
- Under **FASTA-file**, select your fasta file `{{ include.fasta }}`
- Click the **Save** button
