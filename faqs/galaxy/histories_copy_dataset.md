---
title: Copy a dataset between histories
description: Sometimes you may want to use a dataset in multiple histories. You do not need to re-upload the data, but you can copy datasets from one history to another.
area: histories
box_type: tip
layout: faq
contributors: [lecorguille,shiltemann,hexylena,bebatut,lldelisle]
---

There 3 ways to copy datasets between histories

1. From the original history

    1. Click on the {% icon galaxy-gear %} icon which is on the top of the list of datasets in the history panel
    2. Click on **Copy Datasets**
    3. Select the desired files
    {% if include.history_name %}
    4. "New history named:" `{{ include.history_name }}`
    {% else %}
    4. Give a relevant name to the "New history"
    {% endif %}
    5. Validate by 'Copy History Items'
    5. Click on the new history name in the green box that have just appear to switch to this history

2. Using the {% icon galaxy-columns %} **Show Histories Side-by-Side**

    1. Click on the {% icon galaxy-dropdown %} dropdown arrow top right of the history panel (**History options**)
    2. Click on {% icon galaxy-columns %} **Show Histories Side-by-Side**
    3. If your target history is not present
        1. Click on 'Select histories'
        2. Click on your target history
        3. Validate by 'Change Selected'
    3. Drag the dataset to copy from its original history
    4. Drop it in the target history

3. From the target history

    1. Click on **User** in the top bar
    2. Click on **Datasets**
    3. Search for the dataset to copy
    4. Click on its name
    5. Click on **Copy to current History**
