---
title: Flatten a list of list of paired datasets into a list of paired datasets
area: rule-builder
box_type: tip
layout: faq
contributors: [hexylena]
---

Sometimes you find yourself with a `list:list:paired`, i.e. a collection of collection of paired end data, and you really want a `list:paired`, a flatter collection of paired end data. This is easy to resolve with {% tool [Apply rules](__APPLY_RULES__) %}:

1. Open {% tool [Apply rules](__APPLY_RULES__) %} 
2. Select your collection
3. Click **Edit**

You'll now be in the Apply rules editing interface. There are three columns (if it's a `list:list:paired`)

1. The outermost `list` identifier(s)
2. The next `list` identifier(s)
3. The paired-end indicator

Flattening this top level list, so it's just a `list:paired` requires a few changes:

1. From **Column** menu select `Concatenate Columns`
   - *"From Column"*: `A`
   - *"From Column"*: `B`
   This creates a column with the top list identifier, and the inner list identifier, which will be our new list identifier for the flattened list. 
2. From **Rules** menu select `Add / Modify Column Definitions`
   - Click `Add Definition` button and select `Paired-end Indicator`
     - *"Paired-end Indicator"*: `C`
   - Click `Add Definition` button and select `List Identifier(s)`
     - *"List Identifier(s)"*: `D`
   - Click Apply
3. Click **Save**
4. Click **Execute**

The tool will execute and reshape your list, congratulations!
