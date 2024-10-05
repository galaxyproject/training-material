---
title: Ensuring Workflows meet Best Practices
area: workflows
box_type: tip
layout: faq
contributors: [hexylena, elichad]
---

When you are editing a workflow, there are a number of additional steps you can take to ensure that it is a Best Practice workflow and will be more reusable.

1. Open a workflow for editing
1. In the workflow menu bar, you'll find the {% icon galaxy-wf-options %} **Workflow Options** dropdown menu.
1. Click on it and select {% icon galaxy-wf-best-practices %} **Best Practices** from the dropdown menu.

   ![screenshot showing the best practices menu item in the gear dropdown.]({% link faqs/galaxy/images/best-practices1.png %})

1. This will take you to a new side panel, which allows you to investigate and correct any issues with your workflow.

   ![screenshot showing the best practices side panel. several issues are raised like a missing annotation with a link to add that, and non-optional inputs that are unconnected. Additionally several items already have green checks like the workflow defining creator information and a license.]({% link faqs/galaxy/images/best-practices2.png %})

The Galaxy community also has a [guide on best practices for maintaining workflows](https://planemo.readthedocs.io/en/latest/best_practices_workflows.html). This guide includes the best practices from the Galaxy workflow panel, plus:

* adding tests to the workflow
* publishing the workflow on GitHub, a public GitLab server, or another public version-controlled repository
* registering the workflow with a workflow registry such as WorkflowHub or Dockstore
