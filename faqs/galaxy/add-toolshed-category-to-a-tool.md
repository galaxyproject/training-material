---
title: Add Toolshed category to a tool
layout: faq
area: tools
box_type: tip
google_form_id: 1730825288
contributors:
- nomadscientist
---
1. Find the target tool in the [Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/repository). Note -

the easiest way to do this from the Galaxy interface is to (A) search for the tool,

then (B) select the drop-down menu **See tool in toolshed**. 2. Follow the **Development respository** URL. 3. Go to the `.shed.yml` file. 4. In the `categories:` metadata section, add your Toolshed category (which must correspond to those already in the [Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/repository).

Example format: ``` categories:

 - Single Cell

 - Spatial Omics

 - Transcriptomics ```
