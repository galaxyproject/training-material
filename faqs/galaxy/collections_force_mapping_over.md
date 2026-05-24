---
title: Force a tool to process collection jobs one by one instead of using the collection as input
area: collections
box_type: tip
layout: faq
contributors: [paulzierep, Smkingsize, scottcain]
---

## To map over or not to map over

Galaxy tools can either **use a collection as a single input** or **map over** a collection to process each element individually.
Whether a collection is consumed as a whole or mapped over depends entirely on **how the tool is designed**.
For example, {% tool [fastp](toolshed.g2.bx.psu.edu/repos/iuc/fastp/fastp/1.0.1+galaxy3) %} processes each FASTQ file individually. When you supply a collection to fastp, Galaxy indicates that:
`The supplied input will be mapped over this tool.`

![Screenshot showing fastp mapping over.]({% link faqs/galaxy/images/fastp_mapping_over_collection.png %})

Other tools can process an entire collection at once. For example {% tool [collection_column_join](toolshed.g2.bx.psu.edu/repos/iuc/collection_column_join/collection_column_join/0.0.3) %} operates on all files in the collection together.

![Screenshot showing column join without mapping over.]({% link faqs/galaxy/images/column_join_no_mapping_over.png %})


## The catch

Some tools *allow* a collection as a single input, but you may want them to process **each element one-by-one** instead.
A common example is {% tool [metaspades](toolshed.g2.bx.psu.edu/repos/nml/metaspades/metaspades/4.2.0+galaxy0) %}:
If you supply a collection of FASTA files, the tool will treat them as a single dataset and perform **co-assembly**, which is its default behavior.
However, in many workflows you want to assemble each sample **individually**, not all together.

Because the tool form does not offer an option to switch this behavior, you can **force mapping-over** by creating a **nested list**.


## Solution - Create a nested collection

Convert your original collection into a **collection of collections** (`list:list:`).
This forces any tool - including ones that normally process the whole collection - to run on **each subcollection individually**.

### Solution for Galaxy Server with Version > 25.1: Create a nested collection using Nest collection (where available)

Use the {% tool [Nest collection](__NEST__) %}  tool to convert your original collection into a **collection of collections**.

### Solution for Galaxy Server with Version < 25.1 or if nest collection isn't available: Create a nested collection using APPLY_RULES

Use the {% tool [Apply rules](__APPLY_RULES__) %}  tool to convert your original collection into a **collection of collections**.

1. Open {% tool [Apply rules](__APPLY_RULES__) %}
2. Select your collection
3. Click **Edit**

Adding another nesting level to create `list:list` requires a few changes:

1. From **Column** menu select `Basename of Path of URL`
   - From Column A
   - Click **Apply**
   This creates a new Column B with the same list identifier as Column A.
2. From **Rules** menu select `Add / Modify Column Definitions`
   - Click `Add Definition` button and select `List Identifier(s)`
     - *"Select a column"*: `A`
     - *"... Assign Another Column"*: `B`
   - Click Apply
3. Click **Save**
4. Click **Run Tool**

The rule logic in a workflow editor to follow for paired reads in the Apply rules tool looks like this:

![Screenshot showing logic to create the new collection.]({% link faqs/galaxy/images/nested_collection_for_mapping_over.png %})

