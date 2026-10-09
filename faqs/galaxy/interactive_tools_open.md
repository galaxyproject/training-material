---
title: Open interactive tool
area: interactive tools
box_type: tip
layout: faq
contributors: [mgramm1, shiltemann, Sch-Da]
optional_parameters:
  tool: The name of the interactive tool
examples:
  Basic example:
    tool: "Rstudio"
---

1. Go to **Interactive Tools** on the **Activity Bar** on the left
2. If your **Activity Bar** does not show the **Interactive Tools** icon, make sure you are logged in first. 
3. Wait for {{ include.tool | default: "your interactive tool" }} to be in the **running** state (Job Info).
4. You might need to wait a bit until your interactive tool is loaded and ready for you.
5. Click on {{ include.tool | default: "the tool name" }}
   - Wait until the {% icon external-link %} icon next to the tool name appears; otherwise, you will run into an error message. 
   - Clicking on the {% icon external-link %} icon will open it in a new tab
   
