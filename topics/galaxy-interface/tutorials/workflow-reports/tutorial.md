---
layout: tutorial_hands_on

title: Workflow Reports
tags:
- workflows
zenodo_link: ''
questions:
- What are workflow reports?
- How can I view a workflow report?
- How can I define a workflow report?
objectives:
- Understanding, viewing and creating workflow reports
time_estimation: 30m
key_points:
- Workflow reports help you display the most important results of a workflow in an organized fashion.
contributors:
- shiltemann
level: Intermediate
subtopic: interface
---


# Introduction
{:.no_toc}

Workflows are a powerful Galaxy feature that allows you to scale up your analysis by performing an end-to-end analysis with a single click of a button. In order to aid interpretation of workflow results, *workflow reports* may be configured to combine and display the most important analysis results in a single, customizable view.

This is especially useful if you are configuring a Galaxy workflow to share with others. Not everybody is familiar with Galaxy, and having all the important results shown on a single page can be very useful.

This tutorial will guide you through the steps of defining such a report for your workflow, and how to view workflow reports after running the workflow.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Import an example workflow

For this tutorial, we will use the workflow from the [Galaxy 101 for everyone tutorial]({% link topics/introduction/tutorials/galaxy-intro-101-everyone/tutorial.md %}). If you have not done this tutorial yet, the only thing you need to know is that this is a workflow that takes as input a table of data about different species of iris plants, this table is subsequently sorted and filtered, and some plots are made. The specifics of the workflow are not important for this tutorial, only that it outputs a number of different kinds of outputs (images, tables, etc).

We will start by importing this workflow into your Galaxy account:

> ### {% icon hands_on %} Hands-on: Import the workflow
>
> 1. **Import the workflow** into Galaxy
>    - Copy the URL (e.g. via right-click) of [this workflow]({{ site.baseurl }}{{ page.dir }}workflows/galaxy-101-everyone.ga) or download it to your computer.
>    - Import the workflow into Galaxy
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
{: .hands_on}


# Run the workflow and view the default report

Galaxy will produce a default report for any workflow. This default report shows the workflow inputs, outputs, and a description of the workflow on a single web page. You will usually want to customize this report yourself, but it provides a good starting point.

Let's run the workflow and view the default report.


> ### {% icon hands_on %} Hands-on: Run the workflow
>
> 1. {% tool [Import](upload1) %} the file `iris.csv` via link
>
>    ```
>    https://zenodo.org/record/1319069/files/iris.csv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 2. Run **GTN Training: Galaxy 101 For Everyone** {% icon workflow %} using the following parameters:
>    - *"Send results to a new history"*: `No`
>    - {% icon param-file %} *"1: Iris Dataset""*: the `iris.csv` file we just uploaded
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
> 3. **View the workflow outputs** {% icon galaxy-eye %} once the workflow has completed
>    - The workflow produces several text and tabular outputs, and two plot (image) outputs
>
{: .hands_on}

After the workflow has completed, we can access the workflow report. The report does not appear in your history, but can be accessed from the **User -> Workflow Invocations** menu on the top bar. An invocation of a workflow means one run (execution) of the workflow.


> ### {% icon hands_on %} Hands-on: View the default workflow report
>
> 1. In the top menu bar, go to **User -> Workflow Invocations**
>
> 2. Our latest workflow run should be listed at the top.
>    - Click on it to expand it:
>
>    ![screenshot of the workflow invocations menu, with our latest invocation at the top](./images/invocations-list.png)
>
> 3. Click **View Report** in the expanded view of the workflow invocation.
>
> 4. You should see a page like this. It contains:
>    - The input file
>    - The (text-based) output files
>    - A summary of the workflow itself
>
>    ![screenshot of the beginning of the default workflow report](./images/report-default.png)
>
{: .hands_on}



# Customize the workflow report

This is a great start, but we might want to customize this report to fit our needs.


> ### {% icon hands_on %} Hands-on: Open the workflow report editor
>
> 1. Open the workflow in the **workflow editor**
>
>    {% snippet faqs/galaxy/workflows_edit.md %}
>
> 2. Click on **Edit Report** {% icon galaxy-wf-edit %} in the top-right of the screen
>
> 3. You should see something like the image below, you will find
>    - Text editor in the center, with the default report specified in [Markdown format](https://www.markdownguide.org/getting-started/)
>    - A list of components that can be added to the report in the left-hand panel
>
>    ![screenshot of the workflow report edit interface. Markdown-formatted editor in the center, with a list of components that can be added to the report in the left-hand panel](./images/report-edit.png)
> 4. Scroll down the report and look at all the components
>    - notice that there is no plot image output shown, even though we know that was created, we will add this to the report later.
>
> 5. To **edit this report**, we can edit the markdown directly. For example, let's
>     - change the title of the report to `# Iris Analysis`
>     - add a line of introduction text for whoever will read the report:
>
>        <pre class="highlight"><code><del># Workflow Execution Report </del>
>        # Iris analysis
>        Below you will find the results for the plant analysis.</code></pre>
>
>    > ### {% icon tip %} Tip: Markdown format
>    >
>    > The report is specified in [Markdown format](https://www.markdownguide.org/getting-started/), this is a simple markup language that
>    > Some basics of the Markdown syntax can be found in [this cheatsheet](https://www.markdownguide.org/cheat-sheet/)
>    {: .tip}
>
> 6. Let's play around with some components we can add via the left-hand panel
>    - Under **Miscellaneous**, select "**Galaxy version** as text" and "**Current Time** as text"
>    - You will see bits of Markdown are added to your report
>    - You can add some text around these parts as well
>    - Make sure the beginning of your report looks something like this:
>
>      ````markdown
>      # Iris Analysis
>      Below are the results for the Iris analysis workflow.
>
>      This workflow was run on:
>
>      ```galaxy
>      generate_time()
>      ```
>
>      With Galaxy version:
>
>      ```galaxy
>      generate_galaxy_version()
>      ```
>      ````
>
> 7. Let's try to add the missing plot outputs as well:
>    - On the left-hand panel, under **Insert Objects**, in the **History** section, choose **Image**
>    - You should see a list of outputs to insert into the report:
>      ![menu for inserting an image output into the report](./images/report-add-image-options.png)
>    - Hmmm, no obvious options to insert the plot outputs. We will need to **label the outputs** in our workflow first, before we can use them here.
>
> 8. But before we do that, let's save our changes and run the workflow again to view their effects.
>    - Click on {% icon galaxy-cross %} (**Return to Workflow**) in the top-right of the screen.
>    - Click on {% icon galaxy-save %} (**Save Workflow**) to save our changes to the report.
>
> 9. **Run the workflow** again
>    - Select `iris.csv` as the input
>
>    {% snippet faqs/galaxy/workflows_run.md %}
>
> 10. **View the new workflow report**, you should see your changes, something like:
>
>     {% snippet faqs/galaxy/workflows_report_view.md %}
>
>     ![screenshot of the workflow report with our edits included](./images/report-edit1.png)
>
> In the next section, we will add labels to our workflow outputs to more easily add them to our workflow report
>
{: .hands_on}



## Add labels to workflow outputs

As you saw in the previous step, we might need to edit the workflow to add labels outputs so we can easily distinguish between the outputs when adding them to the report. This is especially useful for large workflows with many outputs.

> ### {% icon hands_on %} Hands-on: Add output labels to the workflow
>
> 1. Open the workflow in the **workflow editor**
>
>    {% snippet faqs/galaxy/workflows_edit.md %}
>
> 2. Click on one of the **Scatterplot** {% icon tool %} boxes
>    - On the right-hand panel, you should see the settings for the tool
>    - Scroll to the bottom and find the **Configure Output: ..** sections
>    - From the box we can see that output1 (the first) is a `png` output, and output2 is a `pdf` output of the plot.
>    - Let's use the png output for our report.
>    - Also take note of the plot title, since we have 2 runs of the scatterplot tool in this workflow. One plot is about the petals of the iris, and one about the sepals (in this screenshot it is the sepal plot)
>
>    ![](./images/workflow-add-labels1.png)
>
> 3. **Add an output label** for the `png` output
>    - Click on "Configure output: output1"
>    - Add a descriptive label (e.g. "Sepal plot (PNG)")
>    - Notice that the label on the box changs as well
>    ![](./images/workflow-add-labels2.png)
>
> 4. Repeat this process and **add an output label** for the other plot as well
>    - Label it something like "Petal plot (PNG)"
>
> 5. Click on the **Unique** {% icon tool %} toolbox
>    - This tool gives a list of all the unique Iris species found in the dataset
>    - Let's add this to our report as well
>    - **Add a label** to this output (e.g. "Iris Species")
>
> 6. {% icon galaxy-save %} **Save** the workflow
>
{: .hands_on}


Now that we have added our output labels, let's go back to our report editor and add these outputs






# Share or publish your workflow report
