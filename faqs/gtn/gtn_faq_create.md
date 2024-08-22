---
title: Creating a GTN FAQ
description: An explanation on how to create a GTN faq
area: contributors
layout: faq
box_type: tip
contributors: [hexylena]
description: "If you have a snippet of knowledge that is reusable, we recommend you to share with the GTN community, and we encourage you to create an FAQ for it!"
---

If you have a snippet of knowledge that is reusable, we recommend you to share with the GTN community, and we encourage you to create an FAQ for it!

**Creating the FAQ: The Easy Way**

[Fill out this Google Form](https://forms.gle/2JVMfd1AgtenZPvv9). Every day our bot will import the FAQs submitted via this Google Form, and we will process them, perhaps requesting small changes, so we recommend that you have a GitHub account already.

**For Advanced Users**

Have a look at the existing FAQs in the [`faqs/galaxy/` folder](https://github.com/galaxyproject/training-material/tree/main/faqs/galaxy) of the GTN repository for some examples.

A news post is a markdown file that looks as follows:

{% raw %}
```markdown
---
title: Finding Datasets
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


- To review all active Datasets in your account, go to **User > Datasets**.

Notes:
- Logging out of Galaxy while the Upload tool is still loading data can cause uploads to abort. This is most likely to occur when a dataset is loaded by browsing local files.
- If you have more than one browser window open, each with a different Galaxy History loaded, the Upload tool will load data into the most recently used history.
- Click on refresh icon {% icon galaxy-refresh %} at the top of the History panel to display the current active History with the datasets.
```
{% endraw %}
