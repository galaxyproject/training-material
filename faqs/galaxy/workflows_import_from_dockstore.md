---
title: Import workflows from DockStore
area: workflows
box_type: tip
layout: faq
contributors: [nekrut]
---

[Dockstore](https://dockstore.org/) is a free and open source platform for sharing reusable and scalable analytical tools and workflows.

1. Ensure that you are logged in to your Galaxy account.
1. {% unless include.dockstore_id %}Go to [DockStore](https://dockstore.org).
2. Select any Galaxy workflow you want to import.{% else %}Go to [DockStore](https://dockstore.org/workflows/{{ include.dockstore_id }}){% endunless %}
3. Click on "Galaxy" dropdown within the "Launch with" panel located in the upper right corner.
4. Select a galaxy instance you want to launch this workflow with.
5. You will be redirected to Galaxy and presented with a list of workflow versions.
6. Click the version you want (usually the latest labelled as "main")
7. You are done!

The following short video walks you through this uncomplicated procedure:

{% include _includes/youtube.html id="K2wFrSLFpa0" title="Importing from Dockstore" %}
