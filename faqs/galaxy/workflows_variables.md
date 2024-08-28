---
title: Rename outputs in a workflow using workflow variables
area: workflows
box_type: tip
layout: faq
contributors: [delphine-l]
---

Workflow text inputs can be used as parameters for tools but also in the workflow editor for renaming files. 

# Create a workflow variable

1. Open the workflow editor
2. Click on **Inputs** in the toolbar on the left
3. Click on **Simple inputs used for Workflow logic**
4. Click on the created input box in the workflow editor
5. Label your input. Example: `Species Name`

# Use a workflow variable to rename a dataset

To use the workflow parameter for renaming a dataset, use the syntax: `${Parameter Name}`. For example, to rename a dataset using my species name: 
1. Open the workflow editor
2. Click on the tool in the workflow to get the details of the tool on the right-hand side of the screen
3. Scroll down and click on **Configure Output**
4. In **Rename Dataset**, enter the new dataset name: `Tool run on ${Species Name}`