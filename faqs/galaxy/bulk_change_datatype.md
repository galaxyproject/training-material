---
title: Changing the data type of multiple datasets or collections
description: This will set the data type for multiple datasets, collections, and the datasets therein. Does not change the datasets themselves.
area: collections
box_type: tip
layout: faq
contributors: [kostrykin]
---

1. Click on {% icon galaxy-selector %} **Select Items** at the top of the history panel ![Select Items button]({% link topics/galaxy-interface/images/historyItemControls.png %})
2. Check {% if include.datasets_description %}{{ include.datasets_description }}{% else %}all the datasets and collections in your history you would like to change the data type for{% endif %}
3. Click **{% if include.n %}{{ include.n }}{% else %}n{% endif %} of N selected** and choose **Change data type** 
4. In the drop-down field, select {% if include.datatype %}`{{ include.datatype }}`{% else %} your desired data type {% endif %}
  - tip: you can start typing the data type into the field to filter the dropdown menu
4. Click the **OK** button
