---
layout: tutorial_hands_on

title: Updating a tutorial
subtopic: getting-started
priority: 6

questions:
- How can I update a tutorial in the GTN?

objectives:
- Implement tutorial update on existing GTN material

time_estimation: 1H
key_points:
- The GTN has a number of fantastic features that make it a cutting-edge resource.
- Jumping into existing materials can be daunting, but by following these steps, you can be a contribute to this vibrant community!

contributions:
  authorship:
    - nomadscientist

requirements:
-
  type: "internal"
  topic_name: community
  tutorials:
      - github-interface-contribution

---

Here, we provide a clear set of instructions for updating a GTN tutorial that uses the Galaxy interface for performing analysis. The GTN, and Galaxy, have a number of features to make it as reproducible and shareable as possible, so navigating what needs to be done for an overall update may feel daunting - but not anymore!

We encourage you to pick a tutorial to try this on, so you can use this tutorial side-by-side.

# Phase 1: Find a tutorial with an outdated workflow

{% include _includes/cyoa-choices.html option1="No" option2="Yes" default="No"
       text="Do you already have a tutorial and workflow you want to update?" %}

<div class="No" markdown="1">

> <hands-on-title>Check for an outdated workflow</hands-on-title>
>
> 1. **Find a {% icon hands_on %} Hands-on tutorial** in your training topic of interest (for example, Single-cell)
> 2. **Select the {% icon workflow %} workflow** from the header
> 3. **Import the workflow** to your server
> 4. Go to the **Workflow** menu and select {% icon galaxy-wf-edit %} **Edit workflow**
> 5. **Click through the tools** in the workflow and check the {% icon tool-versions %} to see if the tools are outdated.
{: .hands_on}

If the workflow has tools that are up to date (or very close!), great! That tutorial does not need updating! Try another one!

</div>

<div class="Yes" markdown="2">
> <hands-on-title>Import the workflow</hands-on-title>
>
> 1. **Navigate to your target {% icon hands_on %} Hands-on tutorial**
> 2. **Select the {% icon workflow %} workflow** from the header
> 3. **Import the workflow** to your server
> 4. Go to the **Workflow** menu and select {% icon galaxy-wf-edit %} **Edit workflow**
{: .hands_on}

</div>

# Phase 2: Check that nobody else is working on this

Always good to check! If you're not sure how, skip to the next step!

# Phase 3: Update the workflow

> <hands-on-title>Update the workflow</hands-on-title>
>
> 1. Select {% icon galaxy-wf-options %} **Workflow Options** from the **Workflow Editor**
> 2. Select {% icon upgrade_workflow %} **Upgrade workflow** to automatically update the tools.
> 3. Address any issues that may arise.
> 4. **Save** the workflow.
{: .hands_on}

# Phase 4: Test & fix the workflow

> <hands-on-title>Test the workflow</hands-on-title>
>
> 1. Import the input datasets from the tutorial you are upgrading (follow the instructions on the tutorial itself for this).
> 2. Run your updated workflow on the datasets.
> 3. Address any issues that may arise, and **NOTE DOWN** all changes.
> 4. Add yourself as an author of the workflow.
> 5. **Save** the final workflow.
> 6. **Save** the final (*answer key*) history (with the correctly run workflow).
{: .hands_on}

# Phase 5: Update the tutorial

> <hands-on-title>Update the tools in the tutorial</hands-on-title>
>
> 1. **Update tool versions**: Each {% icon hands_on %} *Hands-on* step in the tutorial will likely have tool (s) with versions in the text. Update these versions / links to be equivalent to your updated workflow.
> 2. **Update tool instructions**: (Simultaneously) update the text to address any differences arising during the update, i.e. new parameters to set or other changes.
> 3. **Update images**: Wherever images are used to show tool outputs, these will need updating. Use your final history for this.
> 4. **Update text**: Wherever results are referenced in the tutorial text, update these numbers (and possibly interpretation) to reflect the new answer key history.
> 5. **Update header**: In the metadata at the beginning of the tutorial, there are likely answer key history links. Update this by making your answer key history public and putting the link there.
> 6. **Update workflow file**: In the tutorial folder in the `training-material` repository, there is a subfolder titled `workflows`. Download your updated workflow as a file and deposit that file there.
> 7. Check if there are any links to workflows/histories in the tutorial text, and if so, update them.
> 7. Add your name as a *contributor* to the tutorial as an **editor** in the metadata.
{: .hands_on}

You may find this tutorial to be a helpful reference in the Markdown content of a GTN tutorial: [Creating content in Markdown]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %})

# Phase 6: Make a pull request and (optional) update workflow testing

At this point, you're welcome to make your Pull Request for the updated tutorial. However, it is likely to fail linting if the workflow testing is not also updated. This can be tricky, so we'd rather you make the Pull Request with all the work you did than let this stop you.

But if you can make the workflow tests, that would be amazing!
Here's the tutorial: [Creating workflow tests using Planemo]({% link faqs/gtn/gtn_workflow_testing.md %})

{% icon congratulations %} Congratulations! You've made it to the end! Thank you so much for contributing to the sustainability of Galaxy tutorials!
