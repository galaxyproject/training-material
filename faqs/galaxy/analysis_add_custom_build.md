---
title: Adding a custom database/build (dbkey)
description: Galaxy may have several reference genomes built-in, but you can also create your own.
area: analysis
box_type: tip
layout: faq
contributors: [shiltemann]
---

- In the top menu bar, go to the *User*, and select *Custom Builds*
- Choose a name for your reference build {% if include.name %}`{{ include.name }}`{% endif %}
- Choose a dbkey for your reference build {% if include.dbkey %}`{{ include.dbkey }}`{% endif %}
- Under **Definition**, select the option `FASTA-file from history`
- Under **FASTA-file**, select your fasta file {% if include.fasta %}`{{ include.fasta }}`{% endif %}
- Click the **Save** button
