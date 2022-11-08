---
title: Creating a new file
description: Galaxy allows you to create new files from the upload menu. You can supply the contents of the file.
area: datasets
box_type: tip
layout: faq
contributors: [bebatut,shiltemann,wm75,hexylena,simonbray]
---

* Open the Galaxy Upload Manager
* Select **Paste/Fetch Data**
* Paste the file contents into the text field
{% if include.name %}
* Change the dataset name from "New File" to `{{ include.name }}`
{% endif %}
{% if include.format %}
* Change **Type** from "Auto-detect" to `{{ include.format }}`
{% endif %}
{% if include.genome %}
* Change **Genome** to `{{ include.genome }}`
{% endif %}
{% if include.convertspaces %} * From the Settings menu ({% icon galaxy-gear %}) select **Convert spaces to tabs** {% endif %}
* Press **Start** and **Close** the window
