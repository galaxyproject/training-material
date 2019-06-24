---
layout: tutorial_hands_on

title: 'Workflows: Using Workflow Parameters'
zenodo_link: ''
questions:
- What are Workflow Parameters
- How can I define and use Workflow Parameters
- How can I read Parameters from Datasets
objectives:
- Learn how to use Workflow Parameters to improve your Workflows
time_estimation: "10m"
key_points:
- Use Workflow Parameters to make your Workflows more versatile
contributors:
- mvdbeek

---


# Introduction
{:.no_toc}

Workflows are a powerful feature in Galaxy that allow you to chain multiple steps of an analysis together.
To make a workflow reusable with slightly different settings you can define and use workflow parameters.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Add workflow parameters to a workflow

We will import a simple workflow and then demonstrate how workflow parameters
can be added.

> ### {% icon hands_on %} Hands-on: Import integer workflow
>
> 1. Click on Workflow in the top menu
> 2. Click on the {% icon galaxy-upload %} import workflow button on the right
> 3. Enter the following URL into the "Archived Workflow" URL box
>
>    ```
>    https://github.com/galaxyproject/training-material/raw/master/topics/galaxy-ui/tutorials/workflow-parameters/workflows/goenrichment-workflow.ga
>    ```
> 4. Click on import workflow
{: .hands_on}

# Add an integer workflow parameter

This workflow takes an input dataset and outputs the first 10 lines of that input dataset.
The number of lines that should be output can be set using the `Number of lines` parameter
of the select first tool.

Instead of selecting a specific number of lines we can choose to insert a workflow parameter
that allows the user to set this parameter when running the workflow.
To do that click on the `Select first` tool in the workflow and find the `Number of lines`
parameter in the `Details` panel on the right hand side of the workflow.
If you hover over the second symbol you will see a small text box appearing that says
"Add connection to module". Click on this symbol and you will see that a new input with the text `Number of lines` has
appeared in the `Select first` tool in the editor.

To insert the workflow parameter click on Inputs in the tool panel on the right side of the workflow editor
and select `Simple inputs used for workflow logic`. This will insert a new box into the workflow editor.

If you click on the Input parameter box in the editor you will see in the `Details panel` that the `Parameter type`
is set to `Text`. If you try to connect this parameter to the `Number of lines` parameter you will see that
the noodle turns orange and that you cannot create a connection. If you change the parameter type to
`Integer` on the right hand side you will be able to connect the workflow parameter to the `Number of lines` parameter.
Make sure to also set a descriptive label and annotation, and save the workflow by clicking on the `save` icon.
If you run the workflow you will now see the new parameter that needs to be set before the workflow can be run.


![Animation of simple integer workflow parameter](../../images/workflow_integer_param.gif "Integer workflow parameter")

> ### {% icon tip %} Tip: You can use a single parameter multiple times
>
> * You can connect a parameter to multiple steps or to multiple parameters within a step
{: .tip}

# Create a text parameter from multiple parameters

It is often necessary to create a text parameter that is flanked by additional text.
Take for instance a regular expression that finds `foo` in a string and replaces it
with `foobar` as in `s/(foo)/\1bar/`.
The Galaxy tool ``Regex Find And Replace`` can be used to tun such a regular expression. The find part of the regex (`s/(foo)/`) can be defined in the `Find regex` parameter and the replacement part can be entered in the `Replacement` parameter.
If we want to make the `foo` part of the regular expression configurable we can
compose this text parameter using the `Compose text parameter value` tool.

> ### {% icon hands_on %} Hands-on: Compose a text parameter
> 1. **Add a data input to the workflow**
> 2. **Add a `Simple inputs used for workflow logic` to the workflow**
>    - {% icon param-select %} *"Parameter type"*: `Text`
> 3. **Add the `Compose text parameter value` to the workflow** {% icon tool %}
>     - Add three repeats
>     - In the first repeat:
>       - {% icon param-text %} *"Enter text that should be part of the computed value"*: `/(`
>     - In the second repeat:
>       - {% icon param-text %} *"Enter text that should be part of the computed value"*: Leave empty and click "Add connection to module".
>     - In the third repeat:
>       - {% icon param-text %} *"Enter text that should be part of the computed value"*: `)/`.
>     - Connect the `Simple inputs used for workflow logic` to the `Compose text parameter value` parameter input
> 4. **Add the `Regex Find And Replace` to the workflow**
>     - Click on `Insert Check`
>     - Click on `Add connection to module` for the `Find Regex` parameter
>     - *"Replacement"*: `\1bar`
>     - Connect the output of the  `Compose text parameter value` tool to the `Find Regex` parameter.
>     - Connect the data input to the `Select lines from` input of the `Regex Find And Replace` tool
{: .hands_on}

Now upload a text dataset with the contents `wunder`
If you run this workflow on the dataset and you select `wunder` as the newly defined parameter  you will see 2 new datasets in your history. The first dataset has the data type ``expression.json`` and contains the composed parameter value `(wunder)`, the second dataset will contain the output of the `Regex Find And Replace` step. A click on the `i` button will show the used parameters. You will see that the `Find Regex` parameter will contain the values that you entered in the workflow run form. If you look at the dataset content you will see it is `wunderbar`.

# Read a parameter from a dataset

Often times it is necessary to calculate a parameter in one step of a workflow and then to use it in another step of the same workflow. This can be accomplished by reading the parameter from a dataset (As long as it is a text, integer, float, boolean or color parameter).
In this example we will construct a workflow where we calculate the sum of all values in a dataset and then divide the values in this dataset by the sum calculated in the previous step.

> ### {% icon hands_on %} Hands-on: Construct Workflow with Parameters read from a dataset
>
> 1. **Add a data input to a new workflow** {% icon tool %}
> 2. **Add the `Datamash` tool to the workflow and connect the data input** {% icon tool %}
>     - {% icon param-repeat %} **Operation to perform in each group**
>       - {% icon param-select %} *"Type"*: `sum`
>       - {% icon param-text %} *"On column"*: `1`
> 3. **Add the `Parse parameter value` tool to the workflow** {% icon tool %}
>     - *"Select type of parameter to parse"*: integer
>     - Connect `Datamash` output to `Parse parameter value` input
> 4.  **Add the `Compose text parameter value` tool to the workflow** {% icon tool %}
>     - Click `Insert Repat`
>     - In the first repeat:
>       - {% icon param-select %} *"Choose the type of parameter for this field"*: Text Parameter
>       - {% icon param-text %} *"Enter text that should be part of the computed value"*: c1/
>     - In the second repeat:
>       - {% icon param-select %} *"Choose the type of parameter for this field"*: Integer Parameter
>       - *"Enter integer that should be part of the computed value"*: Click on `Add connection to module`
>       - Connect `Parse parameter value` output to input
> 5. **Add the `Compute an expression on every row ` tool to the workflow**
>     - *"Add expression as a new column to"*: Click on `Add connection to module`
>     - Connect data input to tool
>     - Connect `Compose text parameter value` output parameter to `Add expression` parameter input
{: .hands_on}

If run on a tabular dataset this workflow will produce a new dataset, where the last column
will be the result of dividing the value in the first column by the sum of all values in the
first column.

> ### {% icon tip %} Tip: You can try many different parameter values at once
>
> * Often times you need to try a couple of different parameter values and pick the best one.
>   You can create or compute a dataset for each parameter you would like to try and run the
>   `Parse parameter value` tool on it, running all downstream tools once for each parameter.
{: .tip}

# Conclusion
{:.no_toc}

Galaxy Workflows chain together different steps of an analysis. To make your workflows
more useful to your colleagues you can add workflow parameters. Sometimes a parameter is not
known in advance, but can be calculated as part of the workflow. Now you know how to read these
parameters from datasets!
