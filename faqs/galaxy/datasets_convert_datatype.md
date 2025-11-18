---
title: Converting the file format
description: Some datasets can be transformed into a different format. Galaxy has some built-in file conversion options depending on the type of data you have.
area: datasets
box_type: tip
layout: faq
contributors: [shiltemann,hexylena,lldelisle,scottcain]
optional_parameters:
  conversion: The datatype to convert the dataset to
examples:
  Convert to BAM:
    conversion: BAM
---

* Click on the {% icon galaxy-pencil %} pencil icon for the dataset to edit its attributes.
* In the central panel, click {% icon galaxy-chart-select-data %} Datatypes tab on the top.
* In the {% icon galaxy-gear %} Convert to Datatype section, select {% if include.conversion %}`{{ include.conversion }}`{% else %} your desired datatype {% endif %} from "Target datatype" dropdown.
* Click the **Create Dataset** button to start the conversion.
