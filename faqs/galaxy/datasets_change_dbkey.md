---
title: Changing database/build (dbkey)
description: You can tell Galaxy which dbkey (e.g. reference genome) your dataset is associated with. This may be used by tools to automatically use the correct settings.
area: datasets
box_type: tip
layout: faq
contributors: [shiltemann,hexylena,nekrut]
---

- Click the desired dataset's name to expand it.
- Click on the "?" next to database indicator:

  ![UI for changing dbkey]({% link shared/images/datasets_dbkey.svg %})

- In the central panel, change the **Database/Build** field
- Select your desired database key from the dropdown list{% if include.dbkey %}: `{{ include.dbkey }}`{% endif %}
- Click the **Save** button
