---
layout: tutorial_hands_on

title: Updating tool versions in a tutorial
subtopic: getting-started
priority: 6

questions:
- How can I update the tool versions in a tutorial in the GTN?
- What else do these updates impact, and how do I update that for consistency?

objectives:
- Implement tutorial tool version update on existing GTN material

time_estimation: 1H
key_points:
- The GTN has a number of fantastic features that make it a cutting-edge resource.
- Jumping into existing materials can be daunting, but by following these steps, you can be a contribute to this vibrant community!

contributions:
  authorship:
    - nomadscientist
  editing:
    - wee-snufkin

requirements:
-
  type: "internal"
  topic_name: contributing
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
</div>
<div class="Yes" markdown="1">
</div>
> 2. **Select the {% icon workflow %} workflow** from the header
> 3. **Import the workflow** to your favourite Galaxy server
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 4. Go to the **Workflow** menu and select {% icon galaxy-wf-edit %} **Edit workflow**
> 5. **Click through the tools** in the workflow and check the {% icon tool-versions %} to see if the tools are outdated.
{: .hands_on}

If the workflow has tools that are up to date (or very close!), great! That tutorial does not need updating! Try another one!




> <hands-on-title>Import the workflow</hands-on-title>
>
> 1. **Navigate to your target {% icon hands_on %} Hands-on tutorial**
> 2. **Select the {% icon workflow %} workflow** from the header
> 3. **Import the workflow** to your favourite Galaxy server
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 4. Go to the **Workflow** menu and select {% icon galaxy-wf-edit %} **Edit workflow**
{: .hands_on}


# Phase 2: Check that nobody else is working on this

It's always a good idea to check, just in case!

You can:
 - Send a message to the [GTN Matrix Channel](https://matrix.to/#/#Galaxy-Training-Network_Lobby:gitter.im) is your quickest way forward.
 - Search through the [GTN Github Repository](https://github.com/galaxyproject/training-material) for existing draft Pull Requests.
 - Check with [individual communities](https://galaxyproject.org/community/sig/), who may have their own method of tracking. For example, the [🖖🏾Single-cell & sPatial Omics Community](https://galaxyproject.org/community/sig/singlecell/) have a shared [Click-Up board](https://sharing.clickup.com/9015477668/b/h/5-90152810734-2/557452707486fef) at the time of writing.

# Phase 3: Update the workflow

Now, you will update the workflow to using the latest {% icon tool-versions %} tool versions.

> <hands-on-title>Update the workflow</hands-on-title>
>
> 1. Select {% icon galaxy-wf-options %} **Workflow Options** from the **Workflow Editor**
> 2. Select {% icon upgrade_workflow %} **Upgrade workflow** to automatically update the tools.
> 3. Address any issues that may arise.
> 4. **Save** the workflow.
{: .hands_on}


# Phase 4: Test & fix the workflow

It's crucial to test the workflow, as often times the outputs will be different due to the new tool versions. It can also transpire that the newer tool versions lead to errors, either because they don't work or because you need to change parameter settings that were previously unavailable or not required. Testing is key!

> <hands-on-title>Test the workflow</hands-on-title>
>
> 1. Import the input datasets from the tutorial you are upgrading (follow the instructions on the tutorial itself for this).
> 2. Run your updated workflow on the input datasets.
> 3. Address any issues that may arise, and **NOTE DOWN** all changes.
> 4. Ensure the workflow meets workflow best practices using the {% icon galaxy-wf-best-practices %} **Best Practices** button and add your name as a contributor, if not there already.
>
>    {% snippet faqs/galaxy/workflows_best_practices.md %}
>
> 5. Add yourself as an author of the workflow.
> 6. **Save** the `updated workflow`.
> 7. Run your `updated workflow` on the input datasets in a fresh history.
> 7. **Save** this history as an `answer key history` & make your history publicly available in Published Histories.
>
>    {% snippet faqs/galaxy/histories_sharing.md %}
>
{: .hands_on}

# Phase 5: Update the tutorial

> <warning-title>If you run out of time</warning-title>
> You have already completed a large chunk of work (`updated workflow`+ `answer key history`) to get here, and we don't want to lose it!
> Updating the tutorial to match the `updated workflow` can be a separate contribution. So if you get to this point and run out of time, please:
> 1. Create an **issue** on the GTN Github Repository and include shareable links to your `updated workflow` & `answer key history`
> 2. Message on the GTN Matrix channel with a link to your issue and explaining which tutorial you updated.
> If, however, you are able to finish the task yourself, please read on!
{: .warning}

Note that you will need to have done either the [Contributing to the GTN Github tutorial using command line]({% link topics/contributing/tutorials/github-command-line-contribution/tutorial.md %}) or the [Contributing to the GTN Github tutorial using Github Desktop]({% link topics/contributing/tutorials/github-interface-contribution/tutorial.md %}). It's helpful as well to understand the folder structure in the GTN, particularly how images go in an image folder either in a tutorial or in the parent folder. Each `topic` has roughly the following structure:

Each topic has the following structure:

```
├── README.md
├── metadata.yaml
├── images
├── docker
│   ├── Dockerfile
├── slides
│   ├── index.html
├── tutorials
│   ├── tutorial1
│   │   ├── tutorial.md
│   │   ├── slides.html
│   │   ├── data-library.yaml
│   │   ├── workflows
│   │   │   ├── workflow.ga
```

The `tutorial.md` is what you'll be editing, however you will also at the end upload a `workflow.ga` file and likely some image files.

> <hands-on-title>Update the tools in the tutorial</hands-on-title>
>
> 1. **Update tool versions**: Each {% icon hands_on %} *Hands-on* step in the tutorial will likely have tool (s) with versions in the text. Update these versions / links to be equivalent to your `updated workflow`.
> 2. **Update tool instructions**: (Simultaneously) update the text to address any differences arising during the update, i.e. new parameters to set or other changes.
> 3. **Update images**: Wherever images are used to show tool outputs, these will need updating. Use your final `answer key history` for this.
> 4. **Update text**: Wherever results are referenced in the tutorial text, update these numbers (and possibly interpretation) to reflect the new `answer key history`.
> 5. **Update header**: In the metadata at the beginning of the tutorial, there are likely `answer key history` links. Update this by adding the link to your `answer key history`.
> 6. **Update workflow file**: In the tutorial folder in the `training-material` repository, there is a subfolder titled `workflows`. Download your `updated workflow` as a file and deposit that file there.
> 7. Check if there are any links to workflows/histories in the tutorial text, and if so, update them.
> 7. Add your name as a *contributor* to the tutorial as an **editor** in the metadata.
{: .hands_on}

You may find this tutorial to be a helpful reference in the Markdown content of a GTN tutorial: [Creating content in Markdown]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %})

# Phase 6: Make a pull request and (optional) update workflow testing

At this point, you're welcome to make your Pull Request for the updated tutorial. However, it is likely to fail linting if the workflow testing is not also updated. This can be tricky, so we'd rather you make the Pull Request with all the work you did than let this stop you.

But if you can make the workflow tests, that would be amazing!
Here's the tutorial: [Creating workflow tests using Planemo]({% link faqs/gtn/gtn_workflow_testing.md %})


{% icon congratulations %} Congratulations! You've made it to the end! Thank you so much for contributing to the sustainability of Galaxy tutorials!
