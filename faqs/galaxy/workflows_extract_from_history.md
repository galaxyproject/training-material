---
title: Extracting a workflow from your history
description: Galaxy can automatically create a workflow based on the analysis you have performed in a history. This means that once you have done an analysis manually once, you can easily extract a workflow to repeat it on different data.
area: workflows
box_type: tip
layout: faq
contributors: [shiltemann,hexylena,nsoranzo]
---


1. **Clean up** your history: remove any failed (red) jobs from your history by clicking on the {% icon galaxy-cross %} button.

   This will make the creation of the workflow easier.

2. Click on {% icon galaxy-gear %} (**History options**) at the top of your history panel and select **Extract workflow**.

   ![`Extract Workflow` entry in the history options menu]({% link topics/introduction/images/history_menu_extract_workflow.png %})

   The central panel will show the content of the history in reverse order (oldest on top), and you will be able to choose which steps to include in the workflow.

3. Replace the **Workflow name** to something more descriptive.

4. **Rename** each workflow input in the boxes at the top of the second column.

5. If there are any steps that shouldn't be included in the workflow, you can **uncheck** them in the first column of boxes.

6. Click on the **Create Workflow** button near the top.

   You will get a message that the workflow was created.

