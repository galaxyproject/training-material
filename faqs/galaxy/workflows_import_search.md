---
title: Importing a workflow using the search
area: workflows
box_type: tip
layout: faq
contributors: [bebatut,wm75]
---

1. Click on *Workflow* in the top menu bar of Galaxy. You will see a list of all your workflows.
2. Click on the {% icon galaxy-upload %} **Import** icon at the top-right of the screen
3. On the new page, select the **GA4GH servers** tab, and configure the **GA4GH Tool Registry Server (TRS) Workflow Search** interface as follows:
   1. *"TRS Server"*: {% if include.trs_server %}**{{ include.trs_server }}**{% else %}the TRS Server you want to search on (Dockstore or Workflowhub){% endif %}
   2. {% if include.search_query %}*"search query"*: `{{ include.search_query }}`
      {% else %}Type in the *search query*{% endif %}
   3. {% if include.workflow_name %}Expand the workflow named `{{ include.workflow_name }}`
      {% else %}Expand the correct workflow{% endif %} by clicking on it
   4. Select the version you would like to {% icon galaxy-upload %} import

The workflow will be imported to your list of workflows. Note that it will also carry a little green check mark next to its name, which indicates that this is an original workflow version imported from a TRS server. If you ever modify the workflow with Galaxy's workflow editor, it will lose this indicator.
