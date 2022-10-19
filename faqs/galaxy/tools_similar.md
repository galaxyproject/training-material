---
title: Multipile similar tools available
area: tools
box_type: tip
layout: faq
contributors: [shiltemann]
---

Sometimes there are multiple tools with very similar names. If the parameters in the tutorial don't match with what you see in Galaxy, please try the following:

1. Use **Tutorial Mode** {% icon curriculum %} in Galaxy, and click on the blue tool button in the tutorial to automatically open the correct tool and version (not available for all tutorials yet)

   {% snippet faqs/galaxy/tutorial_mode.md %}

2. Check that the **entire tool name** matches what you see in the tutorial
   {% if include.toolname %}
   - In this case, make sure the full tool name is: **"{{ include.toolname }}"**
   {% endif %}
