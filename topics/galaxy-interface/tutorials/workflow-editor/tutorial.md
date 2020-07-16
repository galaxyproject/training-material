---
layout: tutorial_hands_on

title: Crearing, Editing and Importing Galaxy Workflows
zenodo_link: ''
questions:
- How can you construct Galaxy workflows from scratch?
- How can you label outputs
- How can you include workflows in workflows?
- How can you tag workflows?
- How can you manage tool versions?
- How can you manage workflow versions?
objectives:
- Understand key aspects of workflows
- Create clean, non-repetitive workflows
time_estimation: 30m
key_points:
- Creating powerful workflows is easy in the workflow editor.
contributors:
- mvdbeek

---


# Introduction
{:.no_toc}


Workflows are a powerful feature in Galaxy that allow you to link multiple steps of complex analysis.
In this tutorial we will demonstrate how to use the Workflow Editor to construct multiple variants of a simple workflow. Note that these workflows are meant to illustrate different concepts. Not all Workflows require using all of the features described below, but we hope this tutorial will inspire you to make your analysis tasks more efficient.

> ### {% icon tip %} Tip: Create Workflows from existing histories
>  Workflows can be extrated from existing histories. You can find a tutorial for this [here]({{ site.baseurl }}/topics/galaxy-interface/tutorials/history-to-workflow/tutorial.html)
>
{: .tip}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Creating a new workflow

> ### {% icon hands_on %} Hands-on: Create a new Workflow
>
> - Click on **Workflow** in the top panel of the Galaxy page.
> - On the top right you will 2 buttons: **Create** and **Import**.
> - To create a new Workflow click on **Create**
> - Enter a Name and Annotation for your workflow and click **Save**
> - The Workflow Editor will open with a new, empty Workflow loaded
> ![A new empty workflow](../../images/workflow_editor_new_workflow.png "A new empty workflow")
{: .hands_on}

On the left hand side of the editor you see the available tools in the tool
panel. The center panel (or "canvas") holds the workflow layout. Steps will
appear in the center panel. On the right you see the attributes of the
Workflow, such as name, version, annotation and tags. Depending on the context
the contents of the right panel will change, but you can always return to these
attributes by clicking on the **Edit Attributes** button (the Pencil icon on
the upper right).

# Workflow Steps

Workflows logically connect a collection of steps. Possible step types are
currently workflow inputs, tools, and workflows.

## Input steps

You may have noticed that an additional section called "Inputs" has appeared in
the tool panel.  There are 3 input types, "Input dataset", "Input dataset
collection" and "Simple inputs used for workflow logic". Insert an input
dataset or dataset collection for each possible input to your workflow.
"Simple inputs used for workflow logic" allow the definition of parameters that
users can or should change when running your workflow. Please check out the
[Using Workflow Parameters]({{ site.baseurl
}}/topics/galaxy-interface/tutorials/workflow-parameters/tutorial.html) for a
detailed description of how to use these.


> ### {% icon hands_on %} Hands-on: Insert a dataset input and a dataset collection input
>
> 1. Expand the "Inputs" section in the tool panel and click on "Input dataset" to create a new dataset input.
> 2. Click on the new input dataset in the center panel. Set the folling parameter on the right side:
>    - *Label*: `A simple text input`
> 3. Click on *input dataset collection* in the "Inputs" section in the tool panel.
> 4. Click on the new input dataset collection in the center panel
>    - *Label : `A Dataset Collection`
>
> > ### {% icon comment %} Optional Input Datasets & Formats
> > Tools may have optional dataset inputs. If your workflow should use
> > optional datasets, you can set **Optional** to `Yes`. Doing this allows you
> > to connect such an input only to Tool inputs that are optional. You can
> > also restrict the format of an input dataset or input dataset collection.
> > This serves as documentation and prevents selection of incompatible datasets.
> {: .comment}
>
{: .hands_on}

## Inserting Tools

We're going to paste one text file side by side with each element of the input collection
wit Paste1.
We're then going to collapse the output back to a single dataset, and then split it out again.


## Connections

Connections can be made by clicking on an output terminal and dragging the
cursor from to an input terminal.  Input terminals that are compatible with the
current output are highlighted in green and, while input terminals that can't
be connected to are highlighted in Orange. When dragging an incompatible output
over an input a small textbox appears mentioning the reason why a connection
cannot be made.  A valid connection can be made if the format of an output is
allowed as input. A simple text file output for instance cannot be used when
the input requires a binary format. If a dataset collection is required but the
output of a node is a single dataset you will see the message "Cannot attach a
data output to a collection input". If an output of a step is connected to
another input one cannot change the input dataset to a dataset collection. In
order to connect inputs in such a case all outputs of the step must be
disconnected. Connections can be removed by hovering over an input terminal and
clicking.

## Labelling steps

Steps can be labelled. A click on a step will open the steps' setting on the right side. Any label will immeditately appear in the center panel as well.

## Labelling outputs and Post Job Actions

Open a step and scroll to the "Configure Output:" section on the right side of the editor. Here you can set a Label. Outputs with a label can be used as outputs in a subworkflow. You will also be able to set an output name for the dataset, to change the datatype and to add or remove tags.


# Workflow versions

Every time a workflow is saved a new version is created, so that you can go back and forth between new and old versions of a workflow. Click on the pencil symbol to bring up the workflow attributes. You can freely select different versions. You can change an old version of a workflow, and when you save it it will become the newest version.

# Importing Workflows

# Conclusion
{:.no_toc}

Workflows are cool!
