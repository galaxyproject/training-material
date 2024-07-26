---
title: "Tutorial Feature: Easier launching of WorkflowHub & Dockstore Workflows"
contributions:
  authorship: [hexylena]
  funding: [by-covid]
tags: [feature update, gtn, tutorials]
layout: news
---

While most GTN tutorials include their associated workflow directly in the GTN, some may wish to write a tutorial about running workflows from WorkflowHub, e.g. leveraging the [IWC workflow collection](https://workflowhub.eu/projects/33).

We've added two snippets which make that easier than ever:

{% snippet faqs/galaxy/workflows_run_wfh.md title="gromacs-mmgbsa/main" wfhub_id="248" version="4" %}

And one for Dockstore:

{% snippet faqs/galaxy/workflows_run_ds.md title="K-mer Profiling HiFi" dockstore_id="github.com/iwc-workflows/kmer-profiling-hifi-VGP1/main" version="v0.1.5" %}

And one for workflows hosted natively in the GTN:

{% snippet faqs/galaxy/workflows_run_trs.md path="topics/assembly/tutorials/largegenome/workflows/Galaxy-Workflow-Data_QC.ga" title="Galaxy Workflow Data QC" %}

You can read more about these snippets and how to use them yourself in your tutorials in the [GTN Contribution Guide section on Workflows]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %}#workflows).
