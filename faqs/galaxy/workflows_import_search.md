---
title: Importing a workflow using the search
area: workflows
box_type: tip
layout: faq
contributors: [bebatut]
---

- Click on *Workflow* on the top menu bar of Galaxy. You will see a list of all your workflows.
- Click on the {% icon galaxy-upload %} **Import** icon at the top-right of the screen
- Click on **search form** in **Import a Workflow from Configured GA4GH Tool Registry Servers (e.g. Dockstore)**
{% if include.trs_server %}
- Select *"TRS Server"*: `{{ include.trs_server }}`
{% else %}
- Select the relevant TRS Server
{% endif %}
{% if include.search_query %}
- Type `{{ include.search_query }}` in the search query
{% else %}
- Type the query
{% endif %}
{% if include.workflow_name %}
- Expand the workflow named `{{ include.workflow_name }}`
{% else %}
- Expand the correct workflow
{% endif %}
- Click on the wanted version

    The workflow will be imported in your workflows
