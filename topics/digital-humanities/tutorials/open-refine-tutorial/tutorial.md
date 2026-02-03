---
layout: tutorial_hands_on
title: OpenRefine Tutorial for researching cultural data
level: Introductory
zenodo_link: 'https://doi.org/10.5281/zenodo.17047254'
questions:
- How to use OpenRefine in Galaxy to clean your data?
- How to use a workflow in Galaxy to extract and visualise information from your data?
objectives:
- Start OpenRefine as an Interactive Tool in Galaxy
- Use OpenRefine to clean your data (remove duplicates, separate multiple values from the same field, etc.)
- Export your cleaned data from OpenRefine to Galaxy
- Use a pre-existing workflow in Galaxy to extract specific information and visualise your findings
time_estimation: 2H
key_points:
- You can use OpenRefine online through Galaxy.
- OpenRefine allows you to work interactively with messy data.
- Galaxy allows you to run workflows with your data.
- With Galaxy, you can visualise your data in various ways.
contributions:
  authorship:
    - dianichj
    - dadrasarmin
    - Sch-Da
  reviewing:
    - Sch-Da
  funding:
    - nfdi4culture
requirements:
  - type: internal
    topic_name: digital-humanities
    tutorials:
      - introduction_to_dh
answer_histories:
  - label: "UseGalaxy.eu"
    history: https://usegalaxy.eu/u/armin.dadras/h/visualise-amount-of-objects-in-museum-collection
    date: 2025-09-19
---
This tutorial shows how to use **OpenRefine** in Galaxy to clean and visualize data from the **humanities and social sciences**. It has two parts:
- **Introduction to OpenRefine**, based on {% cite Hooland_2013 %} and adapted for Galaxy.
- **Introduction to running Galaxy workflows** to visualize cleaned data and extract specific information.


## What is OpenRefine?

**OpenRefine** is a free, open-source “data wrangler” built for messy, heterogeneous, evolving datasets. It imports common formats (CSV/TSV, Excel, JSON, XML) and domain-specific ones used across GLAM (Galleries, Libraries, Archives and Museums) and official statistics (MARC, RDF serializations, PC-Axis).

It is **non-destructive**—OpenRefine does not alter your source files, but works on copies and saves projects locally. Facets and filters enable you to audit categories, surface outliers, and triage inconsistencies without writing code. Its **clustering** tools consolidate near-duplicates using both key-collision methods (fingerprint, n-gram, phonetic) and edit-distance/nearest-neighbour methods (Levenshtein, PPM) so you can standardize names and places at scale while keeping human oversight.

For enrichment, OpenRefine speaks the **Reconciliation API** to match local values to external authorities (e.g. **Wikidata**, **ROR**) and, optionally, retrieve richer metadata. Transformations—both point-and-click and **GREL** formulas—are recorded as a stepwise, undoable history that you can export as JSON and re-apply to other datasets, enabling reproducible cleaning and easy peer review. Finished tables can be exported cleanly to **CSV/TSV**, ODS/XLS(X), SQL statements, templated JSON, Google Sheets, or can be exported back to Galaxy.

## From Cleaning to Analysis in Galaxy

Once your dataset has been cleaned with OpenRefine, you often want to analyze it further or visualize specific aspects. This is where **Galaxy Workflows** become essential: they enable you to build reproducible pipelines that operate on your curated data, transitioning from one-off cleaning to structured analysis.

## What are Galaxy Workflows?

**Galaxy Workflows** are structured, step-by-step pipelines that you build and run entirely in the browser—either extracted from a recorded analysis *history* or assembled in the visual editor. They can be annotated, shared, published, imported, and rerun, making them ideal for teaching, collaboration, and reproducible research.

A captured analysis is easy to share: export the workflow as JSON (**`.ga`**: tools, parameters, and Input/Output) or export a provenance-rich run as a **[Workflow Run RO-Crate](https://www.researchobject.org/workflow-run-crate/)** bundling the definition with inputs, outputs, and invocation metadata. This lowers the barrier to entry (no local installs; web User-Interface (UI) with pre-installed tools and substantial compute) while preserving best practices (histories track tool versions and parameters; workflows are easily re-applied to new data).

For findability and credit, the community uses **[WorkflowHub](https://workflowhub.eu/)**—a curated registry that supports multiple workflow technologies (including Galaxy) and promotes **FAIR** principles; it offers Spaces/Teams, permissions, versioning, and **DOIs via DataCite**, with metadata linking to identifiers like **[ORCID](https://orcid.org/)** so contributions enter scholarly knowledge graphs and are properly acknowledged.

In practice, you can iterate on a workflow in a familiar Graphic-User-Interface (GUI), export the exact definition or a run package, and deposit it where peers can discover, reuse, review, and cite it—closing the loop between simple authoring and robust scholarly dissemination.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Hands on: Get the data

We will work with a slightly adapted dataset from the **[Powerhouse Museum](https://powerhouse.com.au/)** (Australia’s largest museum group) containing a metadata collection. The museum shared the dataset online before giving API access to its collection. We slightly adapted the dataset and uploaded it to Zenodo for long-term reuse. The tabular file (**36.4 MB**) includes **14 columns** for **75,811** objects, released under a **[Creative Commons Attribution Share Alike (CCASA) license](http://creativecommons.org/licenses/by-nc/2.5/au/)**. We will answer three questions: From **which category** does the museum have the most objects? From **which year** does the museum have the most objects? And **what objects does the museum have from that year?**

**Why this dataset?** It is credible, openly published, and realistically messy—ideal for practising problems scholars encounter at scale. Records include a **Categories** field populated from the **Powerhouse Museum Object Names Thesaurus (PONT)**, a controlled vocabulary reflecting Australian usage. The tutorial deliberately surfaces common quality issues—blank values that are actually stray whitespace, duplicate rows, and multi-valued cells separated by the pipe character `|` (including edge cases where **double pipes** `||` inflate row counts)—so we can practice systematic inspection before any analysis. Without careful atomization and clustering, these irregularities would bias statistics, visualizations, and downstream reconciliation.

We suggest that you download the data from the Zenodo record as explained below. This helps us with the reproducibility of the results.

> <hands-on-title>Upload your data</hands-on-title>
>
> 1. Create a new history for this tutorial and name it "Powerhouse Museum — OpenRefine"
> 2. Import the files from [Zenodo]({{page.zenodo_link}}):
>
>    ```
>    https://zenodo.org/records/17047254/files/phm_collection_adapted.tsv
>    https://zenodo.org/records/17047254/files/stopwords-en.txt
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Ensure that the datatype of "phm_collection_adapted" is "tsv". Otherwise, use convert datatype.
> 4. Verify that the datatype of "stopwords-en" is "txt". If not, convert the datatype.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

In this first part, we will focus on working with the metadata from the Powerhouse Museum. 

# Use OpenRefine to explore and clean your dataset

Access OpenRefine as an interactive tool in Galaxy and explore your data.

## Start OpenRefine

> <hands-on-title>Opening the dataset with OpenRefine</hands-on-title>
>
> 1. Open the tool {% tool [OpenRefine](interactive_tool_openrefine) %}: Working with messy data
>    - *"Input file in tabular format"*:  `openrefine-phm-collection.tsv`
>
> 2. Click on "Run Tool".
>
>    ![OpenRefine tool interface in Galaxy](images/openrefine.png)
>
> 3. After around 30 seconds, a red dot appears over the interactive tools section on the left panel. Click on "interactive tools". A new window opens. Make sure to wait until you see the symbol with an arrow pointing outside the box, which indicates that you can start OpenRefine in a new tab. Now you can open OpenRefine by clicking on its name.
>
>    ![Open OpenRefine tool as an Interactive tool](images/interactive_tools.png)
>
> 4. Here, you can see the OpenRefine GUI. Click on `Open Project`.
>
>    ![Open OpenRefine interface](images/openrefine_interface.png)
>
> 5. Click on `Galaxy file`. If the file does not appear, you may have started OpenRefine before it was fully loaded. Retry steps 3 and 4, and the file should be visible.
>
>    ![Open OpenRefine Open Project as an input](images/openrefine_open_project.png)
>
> 6. You can see the data that has been loaded for you.
>
>    ![Open OpenRefine GUI](images/openrefine_gui.png)
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many rows does this table have?
>
> > <solution-title></solution-title>
> >
> > 1. 75809
> >
> {: .solution}
{: .question}

Great, now that the dataset is in OpenRefine, we can start cleaning it.

## Remove duplicates

In large datasets, errors are common. Some basic cleaning exercises can help enhance the data quality. One of those steps is to remove duplicate entries. 

> <hands-on-title>Removing duplicates</hands-on-title>
>
> 1. Click on the triangle on the left of `Record ID`.
>
>    ![Sort Record ID](images/sort.png)
>
> 2. Click on `Sort...`.
>
> 3. Select `numbers` and click on `OK`.
>
>    ![Sort Record ID options](images/sort2.png)
>
> 4. Above the table, click on `Sort` and select `Reorder rows permanently`.
>
>    ![Sort Record ID reorder permanently](images/sort3.png)
>
> 5. Click on the triangle left of the `Record ID` column. Hover over `Edit cells` and select `Blank down`.
>
>    ![Blank down Record ID](images/sort4.png)
>
> 6. Click on the triangle left of the `Record ID` column. Hover over `Facet`, then move your mouse to `Customized facets` and select `Facet by blank (null or empty string)`.
>
>    ![Facet by blank Record ID](images/sort5.png)
>
> 7. On the left, a new option appears under `Facet/Filter` with the title `Record ID`. It shows two choices, `true` and `false`. Click on `true`.
>
>    ![Facet by blank true Record ID](images/sort6.png)
>
> 8. Click on the triangle to the left of the column called `All`. Hover over `Edit rows`, and select `remove matching rows`.
>
>    ![Remove matching rows Record ID](images/deduplicate.png)
>
> 9. Close the `Record ID` under `Facet/Filter` by clicking on the cross (x) to see all rows.
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many rows have been removed?
>
> > <solution-title></solution-title>
> >
> > 1. 84
> >
> {: .solution}
{: .question}

The dataset no longer contains duplicates based on the Record ID. However, we need to perform further cleaning to enhance the dataset.

## Use GREL

There are many ways to manipulate your dataset in OpenRefine. One of them is the Google Refine Expression Language (GREL). With the help of GREL, you can, for example, create custom facets or add columns by fetching URLs. We will use it to find and replace errors. For more information, refer to the [GREL documentation](https://openrefine.org/docs/manual/expressions).

Take a look at the `Categories` column of your dataset. Most objects were attributed to various categories, separated by "\|". However, several fields contain "\|\|" instead of "\|" as a separator. We want to unify those.

> <hands-on-title>Find and replace typos using GREL</hands-on-title>
>
> To remove the occurance of double pipe "\|\|" from the file we can do the following:
> 1. Click on the triangle on the left of `Categories` and select `Text filter`.
> 2. On the left, using the `Facet/Filter` section, search for the occurrence of \| and \|\|. There are 71061 rows with "\|" and 9 rows with "\|\|". We would like to remove these nine lines, as they were added by mistake.
> 3. Click on the triangle on the left of `Categories`, hover over `Edit cells`, and click on `Transform...`.
> 4. In the new window, use the following text `value.replace('||', '|')` as "Expression" and click on `OK`.
>
>    ![Custom text transform on column Categories](images/filter_grel3.png)
>
>    The expression replaces \|\| with \|. If you search for the occurrence of \|\| again, you will no longer get any results. 
>
>    Many different categories describe the object. You may notice duplicates categorising the same object twice.
>    We also want to remove those to ensure we only have unique categories that describe a single object.
>
> 5. Click on the triangle on the left of `Categories`, hover over `Edit cells`, and click on `Transform...`.
>
>    ![Edit cells Categories](images/filter_grel.png)
>
>    ![Transform Categories](images/filter_grel2.png)
>
> 6. In the new window, use the following text `value.split('|').uniques().join('|')` as "Expression" and click on `OK`.
>
{: .hands_on}

These expressions split categories at the pipe separator and join the unique ones within this column. As a result, duplicate categories for one object are deleted.

> <question-title></question-title>
>
> 1. How many cells had duplicated categories?
>
> > <solution-title></solution-title>
> >
> > 1. 1,668
> >
> {: .solution}
{: .question}

## Atomization

Once the duplicate records have been removed, we can examine the content of the "Categories" column more closely. Different categories are separated from each other by pipe (\|). 
Each entry can be assigned to more than one category. To leverage those keywords, the values in the Categories column must be split into individual cells using the pipe character.

> <hands-on-title>Atomization</hands-on-title>
>
> 1. Click on the triangle on the left of `Categories`, hover over `Edit cells`, and click on `Split multi-valued cells...`.
>
>    ![Atomization of Categories](images/split_multi_valued_cells.png)
>
> 2. Define the `Separator` as `\|` (pipe). Click on `OK`.
>
{: .hands_on}

Are you ready for a little challenge? Let's investigate the categories column of the museum items.

> <question-title></question-title>
>
> 1. How many rows do you have after atomizing the categories column?
> 2. How many entries do not have any category?
>
> > <solution-title></solution-title>
> >
> > 1. 168,476
> > 2. Click on the triangle on the left of `Categories` and hover over `facet` and move your mouse over `Customized facets`, and click on `Facet by blank (null or empty string)`. The `true` value for blank entries is 447.
> >    
> >   ![Facet Blank of atomized Categories](images/facet_categories_blank.png)
> >
> {: .solution}
{: .question}

## Faceting

Now that the `Categories` field is cleaned, we can check the occurrence of categories with various facets.

> <hands-on-title>Faceting</hands-on-title>
>
> 1. Click on the triangle on the left of `Categories`, hover over `Facet`, and click on `Text facet`.
> 2. On the left panel, it mentions the total number of choices. The default value of `count limit` is low for this dataset, and we should increase it. Click on `Set choice count limit`.
>
>    ![Text faceting of atomized Categories](images/text_facet.png)
>
> 3. Enter `5000` as the new limit and click on `Ok`.
>
>    ![Increasing the limit of text facetring](images/text_facet2.png)
>
> 4. Now, you see all categories. Click on `count` to see the categories sorted in descending order.
>
>    ![Text faceting of atomized Categories sorted by count](images/text_facet3.png)
>
{: .hands_on}

You can now see, from which category the museum has the most objects, one of our initial questions about the dataset.

> <question-title></question-title>
>
> 1. What are the top 3 categories? How many items are associated with each of them?
>
> > <solution-title></solution-title>
> >
> > 1. Numismatics (8011), Ceramics (7389), and Clothing and Dress (7279).
> >    Congratulations, you have just answered our first question: from which category does the museum have the most objects?
> >    It is numismatic objects, meaning coins. This makes a lot of sense; coins have a long history and convey a lot of information. They are therefore very interesting for researchers.
> >    Moreover, they are robust and compact, making them durable and relatively easy for museums to store.
> >
> {: .solution}
{: .question}


## Clustering

The clustering allows you to solve issues regarding case inconsistencies, incoherent use of either the singular or plural form, and simple spelling mistakes. We apply those to the object categories for the next step of cleaning.

> <hands-on-title>Clustering of similar categories</hands-on-title>
>
> 1. Click on the `Cluster` button on the left in the `Facet/Filter` tab.
> 2. Use `Key collision` as clustering method. Change the Keying function to `n-Gram fingerprint` and change the n-Gram size to `3`.
>
>    ![Cluster and edit column Categories](images/cluster.png)
>
> 3. Click on the `Cluster` button in the middle window.
>
>    ![Clustered and merged similar Categories](images/cluster2.png)
>
> 4. Here, you can see different suggestions from OpenRefine to cluster different categories and merge them into one. In our tutorial, we merge all the suggestions by clicking on `Select all` and then clicking on `Merge selected and re-cluster`.
>
> 5. Now, you can close the clustering window by clicking on `close`.
>
>    Be careful with clustering! Some settings are very aggressive, so you might end up clustering values that do not belong together!
>
>    Now that the different categories have been clustered individually, we can reassemble them in the respective object single cell.
>    
> 6. Click the Categories triangle and hover over the `Edit cells` and click on `Join multi-valued cells`.
> 7. Choose the pipe character (`|`) as a separator and click on `OK`.
>
>    ![Join multi-valued cells on Categories](images/join.png)
>    
> The rows now look like before, with a multi-valued Categories field.
>
{: .hands_on}

You have now successfully split, cleaned and re-joined the various categories of objects in the museum's metadata! Congratulations.
As before the splitting of columns, we are now back to 75725 rows.
When you are satisfied with your data, choose whether to export the dataset to your Galaxy history or download it directly to your computer.

## Exporting your data back to Galaxy

Exporting your data back to Galaxy allows you to analyse or visualise it with further tools in the platform.
But OpenRefine also allows you to export your operation history, detailing all the steps you took in JSON format.
This way, you can import it later and reproduce the exact same analysis. To do so:

> <hands-on-title>Exporting the OpenRefine history</hands-on-title>
> 
> 1. Click on `Undo/Redo` on the left panel.
> 2. Click on `Extract...`.
>
>    ![Extract OpenRefine](images/extract_tasks.png)
>
> 3. Click on the steps that you want to extract. Here, we selected everything.
> 4. Click on `Export`. Give your file a name to save it on your computer.
>
> ![Extract OpenRefine GUI](images/extract_tasks2.png)
>
{: .hands_on}

However, you will also ensure that you save your data. You can download your cleaned dataset in various formats, such as CSV, TSV, and Excel, within OpenRefine. 
For further analysis or visualisation, we suggest you export it to your Galaxy history.

> <hands-on-title>Exporting the results and history</hands-on-title>
>
> 1. Click on `Export` at the top of the table.
> 2. Select `Galaxy exporter`. Wait a few seconds. In a new page, you will see a text as follows: "Dataset has been exported to Galaxy, please close this tab". When you see this, you can close that tab.
>
>    ![Export results of OpenRefine](images/export_results3.png)
>
> 3. You can now close the OpenRefine interactive tool. For that, go to your history with the orange item `OpenRefine on data [and a number]`. This is your interactive tool. Click "OK" on the small square (it says "Stop this interactive tool" when you mouse over it). You do not need it for your next steps.
> 4. You can find a new dataset (named something like "openrefine-Galaxy file.tsv") in your Galaxy History (with a green background). It contains your cleaned dataset for further analysis.
> 5. You can click on the eye icon ({% icon galaxy-eye %}) and investigate the table.
>
>    ![Cleaned dataset](images/dataset_cleaned.png)
>
{: .hands_on}

Awesome work! However, you may recall that we still have two unanswered questions about our data: From which year does the museum have the most objects? And what objects does the museum have from that year?

# Run a Galaxy Workflow on your cleaned data

Congratulations, you have successfully cleaned your data and improved its quality!
But what can you do with it now?
This depends on your research objectives. For us, it is interesting to extract further information from the data.
To make it easier for you, we have created a workflow that links all the tools needed for this analysis.
We wanted to know which year the museum had the most objects and what they were.
You can follow along and answer those questions with us, or explore the Galaxy tools on your own, to adapt the analysis to your needs.
In this case, be sure to check out our other tutorials, particularly the introductory ones.

## How to find and run existing workflows

> <hands-on-title>Run a Galaxy workflow on your dataset</hands-on-title>
>
> There are different ways to import or create a workflow to Galaxy. For example, you can import a workflow from the registered workflows on [WorkflowHub](https://workflowhub.eu/) which is a registry for describing, sharing, and publishing scientific computational workflows. To do that, you have to navigate to the [WorkflowHub](https://workflowhub.eu/) and find the workflow of interest. In this tutorial, we are working [with this workflow](https://workflowhub.eu/workflows/1884?version=1). When you open the link to this workflow on WorkflowHub, you see the following page:
>
> ![Workflowhub page](images/workflowhub.png)
>
> Please click on the `Run on Galaxy` button on top right. After doing this, you will be redirected to your Galaxy account and see the workflow automatically in your middle panel as follows:
>
> ![Workflow imported to Galaxy](images/workflow.png)
>
> This is one way of importing a workflow to your account.
> 
> Let's assume you have done this and imported the workflow to your Galaxy account.
> 1. You can find all workflows available to you by clicking on the Workflows Icon ({% icon galaxy-workflows-activity %}) on the left panel.
>
>    ![Workflows button](images/workflows.png)
>
> 2. Then, you can select and run different workflows (if you have any workflows in your account). Here, let's click on the Run button ({% icon workflow-run %}) of the workflow we provided to you in this tutorial.
>
>    ![Select this workflow](images/select_workflow.png)
>
> 3. Determine the inputs as follows:
>    
>    Input: `openrefine-Galaxy file.tsv`—This is the file you cleaned in OpenRefine.
>    
>    stop_words_english: `stop_words_english.txt`, which is the file we provided to you in this tutorial.
>
>    ![Determine the inputs of the workflow](images/workflow_inputs.png)
>
> 4. Click on the `Run Workflow` button at the top.
> 5. You can follow the stages of different jobs (computational tasks). They will be created, scheduled, executed, and completed. When everything is green, your workflow has run fully and the results are ready.
>
>    ![Overview of the workflow](images/workflow_overview.png)
>
{: .hands_on}

What can you see here? To follow along, we made all substeps of the task available as outputs. 
To answer our question of which year most elements in the museum derive from, we first cut the column of production time from the table.
You can see this in the file `Cut on Data (Number)`.
From this, we filter only the dates that derive from specific years, not year ranges. (See `Filter Tabular on Data (Number)`.) 
You can click on the arrow in the circle button of this dataset (`Run job again`) to see what exact input was used to exclude year ranges.
Regular expressions help clean remaining inconsistencies in the dataset. (Dataset: `Column Regex Find And Replace on Data (Number)`)
Sorting the production date in descending order, as done in the dataset `Sort on data (lowest Number)`, reveals that one faulty dataset, which is supposed to have been created in 2041, is part of the table. 
We remove it in the next step with the tool `Remove beginning`.

The tool **Datamash** allows for summarising how many elements arrived at the museum in each year.  (Dataset: `Datamash on data (Number)`.)
After we apply this tool, the dataset is no longer 7738 lines long, but only 259. 
This is because the number of times each year appeared in the table was summed up in a second column.
Sorting in ascending order (`Sort on data (Number)`) shows a chronological dataset with the earliest entries in the beginning and the most recent entries at the end of the table.
 We can easily visualise this in a (particularly crowded) bar chart directly within Galaxy. (`Bar chart on data (Number)`)
But this is not the most optimal view to show us which year most objects derive from.

To determine from which year most objects originate, we use another sorting order (`Sort on Data (highest number)`). 

> <question-title></question-title>
>
> 1. From what year does the museum have the most objects?
>
> > <solution-title></solution-title>
> >
> > 1. The dataset `Sort on data (highest number)` shows the number of objects by year, sorted from most to least. 288 items are noted for the year 1969 in the first row. This is the year from which the museum has the most (clearly datable) objects.
> >
> {: .solution}
{: .question}

The next four steps parse this year as a conditional statement step by step. (`Select first on data (Number)`, `Cut on data (Number)`, `Parse parameter value on data (Number)` and `Compose text parameter value`.)  
This means, even if you upload another dataset, the highest number is always selected and taken as an input for the next steps.

Based on this input, which is determined by the year with the highest input, `Search in textfiles on data (Number)` searches for object descriptions from the 288 objects of the most prominent year. 
The table is very rich in information, but not that easy to digest.
To make the table more accessible, we create a word cloud of the object descriptions with the offered stop word list.
If you click on the stop word list we provided, you see what "fill words" are excluded from the word cloud. In essence, only words conveying meaning remain.
This helps us quickly determine what kinds of objects the museum has from this popular year.
The dataset `Word cloud image` shows that most objects from the museum are negatives from Davis Mist, a famous Australian photographer, who created them that year and donated them to the museum.

![Word cloud of objects' descriptions](images/display_1969.png)

# Conclusion

Congratulations! You used OpenRefine to clean your data and ran a workflow from Galaxy with your results! You now know how to perform basic steps in Galaxy, run OpenRefine as an interactive tool, and transfer data from Galaxy to OpenRefine and back. On the way, you have learned basic data cleaning techniques, such as facetting, to enhance the quality of your data. To extract further information from the cleaned data, running a pre-designed workflow showed you a glimpse into Galaxy. Of course, you can always conduct your own analysis using the tools most useful to you, instead.
