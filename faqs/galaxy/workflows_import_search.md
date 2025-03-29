---
title: Importing a workflow using the Tool Registry Server (TRS) search
area: workflows
box_type: tip
layout: faq
contributors: [bebatut,wm75]
redirect_from: [/faqs/galaxy/workflows_import_from_workflowhub]
---

1. Click on {% icon galaxy-workflows-activity %} *Workflows* in the Galaxy activity bar (on the left side of the screen, or in the top menu bar of older Galaxy instances). You will see a list of all your workflows
2. Click on {% icon galaxy-upload %} **Import** at the top-right of the screen
3. On the new page, select the **GA4GH servers** tab, and configure the **GA4GH Tool Registry Server (TRS) Workflow Search** interface as follows:
   1. *"TRS Server"*: {% if include.trs_server %}**{{ include.trs_server }}**{% else %}the TRS Server you want to search on (Dockstore or workflowhub.eu){% endif %}
   2. {% if include.search_query %}*"search query"*: `{{ include.search_query }}`
      {% else %}Type in the *search query*{% endif %}

      {% unless include.disable_image %}![galaxy TRS workflow search field, name:vgp is entered in the search bar, and five different workflows all labelled VGP are listed]({% link topics/assembly/images/vgp_assembly/workflow_list.png %}){% endunless %}
   3. {% if include.workflow_name %}Expand the workflow named `{{ include.workflow_name }}`
      {% else %}Expand the correct workflow{% endif %} by clicking on it
   4. Select the version you would like to {% icon galaxy-upload %} import

The workflow will be imported to your list of workflows. Note that it will also carry a little blue-white shield icon next to its name, which indicates that this is an original workflow version imported from a TRS server. If you ever modify the workflow with Galaxy's workflow editor, it will lose this indicator.

{% unless include.disable_video %}
Below is a short video showing the entire uncomplicated procedure:

{% include _includes/youtube.html id="hoP36Te5wko" title="Importing via search from WorkflowHub" %}
{% endunless %}
