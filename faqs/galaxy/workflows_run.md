---
title: Running a workflow
area: workflows
box_type: tip
layout: faq
contributors: [shiltemann,hexylena,wm75, Marie59,katherine-d21,bebatut, Sch-Da]
---

1. Click on {% icon galaxy-workflows-activity %} **Workflows** on the **Activity Bar** on the left.
2. At the top of the resulting page you will have the option to switch between the *My workflows*, *Workflows shared with me* and *Public workflows* tabs.
3. {% if include.tab %}Select the tab *{{ include.tab }}*{% else %}Select the tab you want to see all workflows in that category{% endif %}
4. {% if include.name %}Search for the workflow named `{{ include.name }}`{% else %}Search for your desired workflow{% endif %}.

    ![Select workflow]({% link topics/climate/images/bgc_calib/bgc_workflow.png %}){:width="15%"}

5. Click on the workflow name: a pop-up window opens with a preview of the workflow.
6. To run it directly: click {% icon workflow-run %} **Run** (top-right). This will take you to the *workflow run form*.
7. Configure the workflow
   - **Send results to a new history**: if enabled, will send the results to a new history instead of your current active history. You can provide a name for the new history here as well.
   - **Re-use jobs with identical parameters**. This will check if any identical jobs have already been run before, and save compute time and energy by re-using the previous results. Great to use if you previously ran (part of) this workflow on the same data already.
8. Set the workflow parameters (e.g. input data)
   {% if include.parameters %} - {{include.parameters}}{% endif %}
{% unless include.tab == "My workflows" %}
{{include}}
9. **Recommended**: click **Import** (left of Run) to make your own local copy under *Workflows / My Workflows*.
{% endunless %}
10. Click on {% icon workflow-run %} **Run Workflow** in the upper right corner to start the workflow.
11. You will now see the *workflow invocation* page showing the progress of your workflow. This page can always be accessed via {% icon galaxy-panelview %} **Workflow invocations** on the Activity bar (left-hand menu).
    - If you sent the results to a new history, you can view this history by clicking on the {% icon galaxy-histories-activity %} history link in the top left corner of the workflow invocation page.

