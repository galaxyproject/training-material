---
title: Adding a tag to a collection
area: collections
box_type: tip
layout: faq
contributors: [shiltemann, hexylena, delphine-l]
---

1. Click on the collection in your history to view it
2. Click on **Edit** {% icon galaxy-pencil %} next to the collection name at the top of the history panel
3. Click on **Add Tags** {% icon galaxy-tags %}
4. Add a tag {% if include.tag %} named `{{include.tag}}`{% else %} starting with `#` {% endif %}
   - Tags starting with `#` will be automatically propagated to the outputs any tools using this dataset.
5. Click **Save** {% icon galaxy-save %}
6. Check that the tag appears below the collection name

