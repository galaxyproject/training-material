---
title: Add Toolshed category to a tool
layout: faq
area: tools
box_type: tip
google_form_id: 1730825288
contributors:
- nomadscientist
---

1. Find the target tool in the [Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/repository).
   
   Note: the easiest way to do this from the Galaxy interface is to (A) search for the tool and open it, then (B) select **See in Tool Shed** from the drop-down menu at the top of the tool form (only shown for tools installed from a Toolshed).
2. Follow the **Development repository** URL.
3. Go to the `.shed.yml` file.
4. In the `categories:` metadata section, add your Toolshed category (which must correspond to those already in the [Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/repository)).

Example format: 

```yaml
categories:
 - Single Cell
 - Spatial Omics
 - Transcriptomics
```
