---
title: Adding a custom database/build (dbkey)
description: Galaxy may have several reference genomes built-in, but you can also create your own.
area: analysis
box_type: tip
layout: faq
contributors: [shiltemann]
---

- Navigate to the History that contains your fasta for the reference genome
- Standardize the [fasta format]({% link faqs/galaxy/reference_genomes_custom_genomes.md %}), and wait until the fasta dataset is green
- Click your username in the masthead, then **Preferences -> Manage Custom Builds** (the page only offers fasta files that were in your current history when it opened)
- Create a unique **Name** for your reference build {% if include.name %}`{{ include.name }}`{% endif %}
- Create a unique **Key** (the dbkey) for your reference build {% if include.dbkey %}`{{ include.dbkey }}`{% endif %}
- Under **Definition**, select the option `FASTA-file from history`
- Under **FASTA-file**, select your fasta file {% if include.fasta %}`{{ include.fasta }}`{% endif %}
- Click the **Save** button
- Wait a minute and reload the page; the build is listed once Galaxy has counted its chromosomes
