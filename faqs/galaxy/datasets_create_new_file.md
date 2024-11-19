---
title: Creating a new file
description: Galaxy allows you to create new files from the upload menu. You can supply the contents of the file.
area: datasets
box_type: tip
layout: faq
contributors: [bebatut,shiltemann,wm75,hexylena,simonbray]
optional_parameters:
  name: The name of the dataset
  format: The format of the dataset
  genome: The genome of the dataset
  convertspaces: Ask the user to convert spaces to tabs in the settings
examples:
  Creating a new file: {}
  Creating a specific file:
    name: "SARS-CoV-2 feature mapping"
    format: "tabular"
    convertspaces: "true"
---

* Click {% icon galaxy-upload %} **Upload Data** at the top of the tool panel
* Select {% icon galaxy-wf-edit %} **Paste/Fetch Data** at the bottom
* Paste the file contents into the text field
{% if include.name %}
* Change the dataset name from "New File" to `{{ include.name }}`
{% endif %}
{%- if include.format -%}
* Change **Type** from "Auto-detect" to `{{ include.format }}`
{%- endif -%}
{%- if include.genome -%}
* Change **Genome** to `{{ include.genome }}`
{%- endif -%}
{%- if include.convertspaces -%}
* From the Settings menu ({% icon galaxy-gear %}) select **Convert spaces to tabs**
{%- endif -%}
* Press **Start** and **Close** the window
