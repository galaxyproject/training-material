---
title: Changing database/build (dbkey)
description: You can tell Galaxy which dbkey (e.g. reference genome) your dataset is associated with. This may be used by tools to automatically use the correct settings.
area: datasets
box_type: tip
layout: faq
contributors: [shiltemann,hexylena]
---

- Click on the {% icon galaxy-pencil %} **pencil icon** for the dataset to edit its attributes
- In the central panel, change the **Database/Build** field
- Select your desired database key from the dropdown list{% if include.dbkey %}: `{{ include.dbkey }}`{% endif %}
- Click the **Save** button
