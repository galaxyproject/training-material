---
title: Copy a dataset to a new history
description: Sometimes you may want to use a dataset in multiple histories. You do not need to re-upload the data, but you can copy datasets from one history to another.
area: histories
box_type: tip
layout: faq
contributors: [lecorguille,shiltemann,hexylena]
---

1. Click on the {% icon galaxy-gear %} icon (**History options**) on the top of the history panel
2. Click on **Copy Dataset**
3. Select the desired files
{% if include.history_name %}
4. "New history name:" `{{ include.history_name }}`
{% else %}
4. Give a relevant name to the "New history"
{% endif %}
5. Click on the new history name in the green box that have just appear to switch to this history
