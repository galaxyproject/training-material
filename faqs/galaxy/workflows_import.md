---
title: Importing a workflow
area: workflows
box_type: tip
layout: faq
contributors: [shiltemann,mblue9,hexylena]
---

- Click on *Workflow* on the top menu bar of Galaxy. You will see a list of all your workflows.
- Click on {% icon galaxy-upload %} **Import** at the top-right of the screen
- {% if include.import_url %}Paste the following URL into the box labelled *"Archived Workflow URL"*: `{{ include.import_url }}`{% else %}Provide your workflow
  - Option 1: Paste the URL of the workflow into the box labelled *"Archived Workflow URL"*
  - Option 2: Upload the workflow file in the box labelled *"Archived Workflow File"*{% endif %}
- Click the **Import workflow** button

Below is a short video demonstrating how to import a workflow from GitHub using this procedure:

{% include _includes/youtube.html id="hoP36Te5wko" title="Importing a workflow from URL" %}
