---
title: Importing a workflow
area: workflows
box_type: tip
layout: faq
contributors: [shiltemann,mblue9,hexylena,wm75]
---

1. Click on {% icon galaxy-workflows-activity %} *Workflows* in the Galaxy activity bar (on the left side of the screen, or in the top menu bar of older Galaxy instances). You will see a list of all your workflows
2. Click on {% icon galaxy-upload %} **Import** at the top-right of the screen
3. {% if include.import_url %}Paste the following URL into the box labelled *"Archived Workflow URL"*: `{{ include.import_url }}`{% else %}Provide your workflow
  - Option 1: Paste the URL of the workflow into the box labelled *"Archived Workflow URL"*
  - Option 2: Upload the workflow file in the box labelled *"Archived Workflow File"*{% endif %}
4. Click the **Import workflow** button

Below is a short video demonstrating how to import a workflow from GitHub using this procedure:

{% include _includes/youtube.html id="FXlkreylhww" title="Importing a workflow from URL" %}
