---
title: Creating a new file
description: Galaxy allows you to create new files from the upload menu. You can supply the contents of the file.
area: datasets
box_type: tip
layout: faq
contributors: [bebatut,shiltemann,wm75,hexylena,simonbray]
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
