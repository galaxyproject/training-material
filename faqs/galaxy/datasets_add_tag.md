---
title: Adding a tag
description: Tags can help you to better organize your history and track datasets.
area: datasets
layout: faq
box_type: tip
contributors: [bebatut,wm75,hexylena,shiltemann]
---

1. Click on the dataset to expand it
2. Click on **Add Tags** {% icon galaxy-tags %}
3. Add a tag {% if include.tag %}named `{{include.tag}}` {% else %} starting with `#` {% endif %}
   - Tags starting with `#` will be automatically propagated to the outputs of tools using this dataset.
4. Press <kbd>Enter</kbd>
5. Check that the tag appears below the dataset name

