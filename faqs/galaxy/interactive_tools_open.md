---
title: Open interactive tool
area: interactive tools
box_type: tip
layout: faq
contributors: [mgramm1, shiltemann]
optional_parameters:
  tool: The name of the interactive tool
examples:
  Basic example:
    tool: "Rstudio"
---

1. Go to **Interactive Tools** on the **Activity Bar** on the left
2. Wait for {{ include.tool | default: "your interactive tool" }} to be in the **running** state (Job Info)
3. Click on {{ include.tool | default: "the tool name" }}
   - Clicking on the {% icon external-link %} icon will open it in a new tab
