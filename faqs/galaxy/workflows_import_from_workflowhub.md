---
title: Import workflows from WorkflowHub
area: workflows
box_type: tip
layout: faq
contributors: [gallardoalba, abueg]
---

[WorkflowHub](https://workflowhub.eu/) is a workflow management system which allows workflows to be FAIR (Findable, Accessible, Interoperable, and Reusable), citable, have managed metadata profiles, and be openly available for review and analytics.

> <hands-on-title> Import a workflow </hands-on-title>
>
> 1. Click on the **Workflow** menu, located in the top bar.
> 2. Click on the **Import** button, located in the right corner.
> 3. In the section "Import a Workflow from Configured GA4GH Tool Registry Servers (e.g. Dockstore)", click on **Search form**.
> 
> 4. In the **TRS Server: *workflowhub.eu*** menu you should type {% if include.filter %}`{{ include.filter }}`{% else %}your query.{% endif %}
>    ![galaxy TRS workflow search field, name:vgp is entered in the search bar, and five different workflows all labelled VGP are listed]({% link topics/assembly/images/vgp_assembly/workflow_list.png)
> 5. Click on the desired workflow, and finally select the latest available version.
{: .hands_on}

After that, the imported workflows will appear in the main workflow menu. In order to initialize the workflow, we just need to click in the {% icon workflow-run %} **Run workflow** icon.
