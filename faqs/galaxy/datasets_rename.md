---
title: Renaming a dataset
area: datasets
box_type: tip
layout: faq
contributors: [bebatut,nsoranzo,shiltemann,wm75,hexylena]
optional_parameters:
  name: The new name for the dataset
examples:
  Without prescribing a name: {}
  Suggest the user rename the dataset genomes.fa:
    name: genomes.fa
---

- Click on the {% icon galaxy-pencil %} **pencil icon** for the dataset to edit its attributes
- In the central panel, change the **Name** field {% if include.name %} to `{{ include.name }}` {% endif %}
- Click the **Save** button
