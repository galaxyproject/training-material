---
title: Changing the datatype of a collection
description: This will set the datatype for all files in your collection. Does not change the files themselves.
area: collections
box_type: tip
layout: faq
contributors: [shiltemann, hexylena]
---

1. Click on **Edit** {% icon galaxy-pencil %} next to the collection name in your history
2. In the central panel, click on the {% icon galaxy-chart-select-data %} **Datatypes** tab on the top
3. Under **new type**, select {% if include.datatype %}`{{ include.datatype }}`{% else %} your desired datatype {% endif %}
  - tip: you can start typing the datatype into the field to filter the dropdown menu
4. Click the **Save** button


**Cannot find the feature?**

If you are on a smaller Galaxy server, i.e. not one of the large (multi)national public servers, you may not be able to find this operation, and there is no indication it is missing or why it is disabled.

Galaxy has recently started putting [more features behind a setting and deployment configuration](https://docs.galaxyproject.org/en/master/admin/production.html#use-celery-for-asynchronous-tasks) that needs to be enabled by the server administrator.
Your administrator will need to deploy Celery and potentially additionally flower and redis to their stack to enable changing the datatype of a collection. Consider sending your Galaxy administrator the link to the [simpler deployment option]({% link topics/admin/tutorials/celeryless/tutorial.md %}) or more complex [GTN tutorial for setting up redis and flower]({% link topics/admin/tutorials/celery/tutorial.md %}).
