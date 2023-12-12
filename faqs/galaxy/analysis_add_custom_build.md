---
title: Adding a custom database/build (dbkey)
description: Galaxy may have several reference genomes built-in, but you can also create your own.
area: analysis
box_type: tip
layout: faq
contributors: [shiltemann]
---

- Navigate to the History that contains your fasta for the reference genome
- Standarize the [fasta format]({% link faqs/galaxy/reference_genomes_custom_genomes.md %})
- In the top menu bar, go to **User -> Preferences -> Manage Custom Builds**
- Create a unique **Name** for your reference build {% if include.name %}`{{ include.name }}`{% endif %}
- Create a unique **Database** (dbkey) for your reference build {% if include.dbkey %}`{{ include.dbkey }}`{% endif %}
- Under **Definition**, select the option `FASTA-file from history`
- Under **FASTA-file**, select your fasta file {% if include.fasta %}`{{ include.fasta }}`{% endif %}
- Click the **Save** button
