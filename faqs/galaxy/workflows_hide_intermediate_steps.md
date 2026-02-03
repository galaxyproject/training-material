---
title: Hiding intermediate steps
area: workflows
box_type: tip
layout: faq
contributors: [shiltemann,hexylena]
---

When a workflow is executed, the user is usually primarily interested in the final product and not in all intermediate steps.
By default all the outputs of a workflow will be shown, but we can explicitly tell Galaxy which outputs to show and which to hide for a given workflow.
This behaviour is controlled by clicking on the eye icon on for all dataset that should be hidden:

![Asterisk for `out_file1` in the `Cut` tool]({% link shared/images/workflow_editor_hide_output.png %})

