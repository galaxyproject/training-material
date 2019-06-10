---
layout: tutorial_hands_on

title: Group tags for complex experimental designs
zenodo_link: ''
questions:
- What are group tags?
- How can I use group tags to perform multi-factor analyses with collections
objectives:
- Learn how to set group tags
- Learn how to select group tags in tools
time_estimation: 10m
key_points:
- Group tags allow complex analyses without reshaping or unhiding datasets in a collection
contributors:
- mvdbeek

---


# Introduction
{:.no_toc}

Advanced uses of Galaxy often require the use of dataset collections,
which can contain between one and tens of thousands of datasets.
Grouping datasets in this way has numerous advantages:
  - It is easy to represent a single collection in the History
  - Dataset names ("Element Identifiers") are immutable and preserved
  - Collections can be split and nested in arbitrary ways

While collections can be split in any way, doing so for multi-factor analysis
quickly becomes cumbersome and messy. An alternative is to label collection
elements with special group tags. These tags can be displayed in the Tool form,
allowing users to select subsets of collections.

This tutorial outlines how to set and use group tags with the DESeq2 tool.
For a more detailed description and background for differential expression
testing see the [Reference-based RNA-Seq data analysis]({{ site.baseurl }}/topics/transcriptomics/tutorials/ref-based/tutorial.html).


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Setting group tags using the apply rules tool

There are three ways to set group tags:
  - Using the rule builder / apply_rules tool
  - Using the "Tag elements from file" tool
  - Manually adding dataset tags with the prefix `group:`

We will use the first two methods in this tutorial.
The second method works at any step during the analysis.

## Set group tags during upload

> ### {% icon hands_on %} Hands-on: Set group tags during upload
>
> 1. Create a new history for this tutorial
>
>    {% include snippets/create_new_history.md %}
>
> 2. Open the Galaxy Upload Manager ({% icon galaxy-upload %} on the top-right of the tool panel)
> 3. Click on **Rule-based** on the top
>
>    ![Rule-based upload](../../images/group-tags/rule_upload.png)
>
>    As you can see in this dialog, data can be selected from a history dataset or pasted in directly
>
> 4. Set **Upload data as:** to `Collection(s)`
> 5. Paste the following links into the text box
>
>    ```
>    https://zenodo.org/record/1185122/files/GSM461176_untreat_single.counts
>    https://zenodo.org/record/1185122/files/GSM461177_untreat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461178_untreat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461179_treat_single.counts
>    https://zenodo.org/record/1185122/files/GSM461180_treat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461181_treat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461182_untreat_single.counts
>    ```
>
> 6. Click **Build**
> 7. We will add a regex that creates 3 new columns with accession, treatment and library type:
>    - Click on the **Columnn** button and then **Using a Regular Expression**
>    - Select **Create columns matching expression groups**
>    - Paste `.*(GSM.*)_(.*)_(.*).counts` in *"Regular Expression"*
>    - Set *"Number of Groups"* to 3
>    - Click on **Apply**
>
>       We should have now a table with 4 columns: link, sample name, treatment, sequencing type
>
>    - Click on **Rules** and then **Add / Modify Column Definitions**
>    - Click on **Add Definition** and select:
>      - *"URL"*: Column A
>      - *"List Identifiers"*: Column B
>      - *"Group Tags"*: Columns C and D
>    - Click **Apply**
>    - Enter a name for the new collection
>    - Click **Build**
>
> 8. Expand the generated collection and the files in it and check their names and tags
>
>    ![Group tags in Galaxy UI are prefixed with group:](../../images/group-tags/group-tags.png)
{: .hands_on}

## Set group tags using the "Tag elements from file" tool

We now want to add group tags using the "Tag elements from file" tool.

> ### {% icon hands_on %} Hands-on: Upload and create a collection
>
> 1. Create a new history for this tutorial
> 2. Import the following files
>
>    ```
>    https://zenodo.org/record/1185122/files/GSM461176_untreat_single.counts
>    https://zenodo.org/record/1185122/files/GSM461177_untreat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461178_untreat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461179_treat_single.counts
>    https://zenodo.org/record/1185122/files/GSM461180_treat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461181_treat_paired.counts
>    https://zenodo.org/record/1185122/files/GSM461182_untreat_single.counts
>    ```
>
>    {% include snippets/import_via_link.md %}
>
> 3. Click on the {% icon galaxy-selector %} icon (**Operations on multiple datasets**)
> 4. Check all new datasets
> 5. Click on **For all selected...** and then **Build Dataset List**
> 6. Enter a name for the new collection and click **Create list**
{: .hands_on}

We have now a collection with our files. We can now either upload a tabular file containing the element identifiers
and the tags we want to apply, or we can extract the element identifiers and extract the tags using a Regular Expression.
We will do the latter.

> ### {% icon hands_on %} Hands-on: Set group tags using the "Tag elements from file" tool
> 1. **Extract element identifiers** {% icon tool %}
>      - {% icon param-collection %} *"Dataset collection"*: created collection
> 2. **Replace Text in entire line** {% icon tool %}
>      - {% icon param-file %} *"File to process"*: output of **Extract element identifiers** {% icon tool %}
>      - In *"Replacement"*:
>         - In *"1: Replacement"*
>            - *"Find pattern"*: `(.*)_(.*)_(.*).counts`
>            - *"Replace with"*: `\1_\2_\3.counts\tgroup:\2\tgroup:\3`
>
>     This step add an additional columns that can be used with the ``Tag elements from file`` tool
>
> 3. Change the datatype to `tabular`
>
>    {% include snippets/change_datatype.md datatype="tabular" %}
>
> 4. **Tag elements from file** {% icon tool %}
>      - {% icon param-collection %} *"Input Collection"*: created collection
>      - {% icon param-collection %} *"Tag collection elements according to this file"*: output of **Replace Text** {% icon tool %}
{: .hands_on}

You should now have a properly tagged collection of tabular files that can be used in DESeq2.

# Using group tags in tool, e.g. DESeq2

DESeq2 has two modes for specifying factors. One can either
select datasets corresponding to factors, or use group tags
to specify factors. We will use the grop tags present in
our collection to specify factors.

The tool interface will prompt you with the group tags that are available for your inputs:

![Group tags in the tool UI](../../images/group-tags/tool-ui.png)

> ### {% icon hands_on %} Hands-on: Running **DESeq2** with group tags
>
> 1. **DESeq2** {% icon tool %} with the following parameters:
>    - *"how"*: `Select group tags corresponding to levels`
>       - {% icon param-collection %} *"Count file(s) collection"*: Generated collection
>       - In *"Factor"*:
>          - In "1: Factor"
>              - *"Specify a factor name"*: `Treatment`
>              - In *"Factor level"*:
>                  - In *"1: Factor level"*:
>                      - *"Specify a factor level"*: `treat`
>                      - *"Select groups that correspond to this factor level"*: `Tags: treat`
>                  - In *"2: Factor level"*:
>                      - *"Specify a factor level"*: `untreat`
>                      - *"Select groups that correspond to this factor level"*: `Tags: untreat`
>          - {% icon param-repeat %} Click on *"Insert Factor"* (not on "Insert Factor level")
>          - In "2: Factor"
>              - "Specify a factor name" to `Sequencing`
>              - In *"Factor level"*:
>                  - In *"1: Factor level"*:
>                      - *"Specify a factor level"*: `paired`
>                      - *"Select groups that correspond to this factor level"*: `Tags: paired`
>                  - In *"2: Factor level"*:
>                      - *"Specify a factor level"*: `single`
>                      - *"Select groups that correspond to this factor level"*: `Tags: single`
>    - *"Files have header?"*: `No`
>    - *"Output normalized counts table"*: `Yes`
{: .hands_on}

# Conclusion
{:.no_toc}

We can select a subset of Collections using the special group tag.
