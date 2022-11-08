---
title: Changing the datatype
description: Galaxy will try to autodetect the datatype of your files, but you may need to manually set this occasionally.
area: datasets
box_type: tip
layout: faq
contributors: [bebatut,nsoranzo,hexylena,shiltemann,ajadi-abiola]
---

* Click on the {% icon galaxy-pencil %} **pencil icon** for the dataset to edit its attributes
* In the central panel, click on the {% icon galaxy-chart-select-data %} **Datatypes** tab on the top
* Select {% if include.datatype %}`{{ include.datatype }}`{% else %} your desired datatype {% endif %}
  - tip: you can start typing the datatype into the field to filter the dropdown menu
* Click the **Save** button
