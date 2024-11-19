---
title: Making an element collapsible in a report
description: If you have extraneous information you might want to let a user collapse it.
area: reports
box_type: tip
layout: faq
contributors: [hexylena, guerler]
---

This applies to any GalaxyMarkdown elements, i.e. the things you've clicked in the left panel to embed in your Workflow Report or Page

By adding a `collapse=""` attribute to a markdown element, you can make it collapsible. Whatever you put in the quotes will be the title of the collapsible box.

````markdown
```
history_dataset_type(history_dataset_id=3108c91feeb505da, collapse="[TITLE]")
```
````
