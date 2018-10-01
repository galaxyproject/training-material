---
layout: tutorial_hands_on
topic_name: genome_annotation
tutorial_name: functional-annotation-workflow
---

# Functional Annotation Workflow

Upon completion of the structural annotation, the genome is ready to be functionally annotated.

> ### Agenda
>
> 1. Workflow
> 2. Completion
>
{: .agenda}

# Workflow

To begin, access Galaxy ([CPT Public Galaxy](https://cpt.tamu.edu/galaxy-pub), [CPT TAMU Galaxy](https://cpt.tamu.edu/galaxy/)]; the data  must be fetched from Apollo into Galaxy. Using the search bar at the top of the Tool panel on the left, enter “Retrieve data.” Click on the hyperlink of the same name underneath CPT: Get Data, and the tool with parameters to adjust will appear. Using the *Organism* drop-down menu, select the structurally annotated phage, then click “Execute.”

![](../../images/functional-annotation-workflow-screenshots/1_retrieve_data_tool.png)

Upon execution, if everything was in the proper place, a message in a green box will appear to inform of the successful execution.

![](../../images/functional-annotation-workflow-screenshots/2_retrieve_data_success_message.png)

Now that the data has been retrieved, it is ready to be run in the functional annotation workflow. At the top of the screen, click on the “Shared Data” drop-down list and select “Workflows.”

![](../../images/functional-annotation-workflow-screenshots/3_shared_data_workflow.png)

In this collection of workflows, look for “PAP 201# Functional (v#.#), where # is the largest number indicating the most recent version of this workflow. Click on the drop-down menu for the most recent functional workflow, and select “Import.”

![](../../images/functional-annotation-workflow-screenshots/4_import_functional_workflow.png)

> ### {% icon tip %} Note that…
> The above image may not precisely reflect the current functional workflow. It is possible the functional annotation workflow has been updated since the creation of this page.
{: .tip}

A message in a green box should appear, indicating a successfully imported workflow. Clicking on the “… starting using this workflow …” hyperlink will bring the user to their imported workflows on Galaxy.

![](../../images/functional-annotation-workflow-screenshots/5_import_functional_workflow_success.png)

Alternatively, clicking on “Workflows” at the top of the Galaxy page (next to the “Shared Data” drop-down menu) will direct the user to the same page. Find the functional annotation workflow that was just imported and click on the drop-down menu; select “Run.”

![](../../images/functional-annotation-workflow-screenshots/6_imported_workflows_run.png)

This will yield a list of parameters to run the functional annotation workflow. Specifically, it is important to ensure that datasets are associated with their counterpart in these parameters.
> 1. Genome Sequence should contain the “#. Sequence(s) from Apollo” set (where # varies dependent on their place in the current History).
> 2. Apollo Organism JSON File should contain the “#. Metadata from Apollo” set (where # varies dependent on their place in the current History).
> 3. Annotation Set should contain the “#. Annotation from Apollo” set (where # varies dependent on their place in the current History).

![](../../images/functional-annotation-workflow-screenshots/7_workflow_parameters.png)

When the proper parameters have been set, select “Run workflow” at the top or bottom of the page. A message in a green box will appear, indicating a successful invocation of the functional annotation workflow.

![](../../images/functional-annotation-workflow-screenshots/8_successful_workflow_execution.png)

This workflow includes multiple computationally-intensive steps. With high server load, it may take up to **several days** for the workflow to complete, as the workflow is executing a large number of analyses on the behalf of the user:

> * BLAST against numerous databases
> * InterProScan
> * phage spanning search tools
> * various other analyses

> ### {% icon tip %} Note that…
> Check back on this workflow periodically; if any of the queued jobs have failed (the datasets in the History column have turned red), click on the dataset. Report a bug by clicking on the bug icon.
>
>![](../../images/functional-annotation-workflow-screenshots/9_report_bug.png)
>
> Some individual jobs (E.G.: BLAST and InterProScan) may remain yellow (“running”) for many hours.
{: .tip}

# Completion

Once all the datasets and tools have completed, then functional annotation within Apollo my begin.

<!-- LINK ANNOTATION INFORMATIONAL UPON COMPLETION. -->

Conclusion about the technical key points. And then relation between the techniques and the biological question to end with a global view.
