---
title: Identifying and formatting Tabular Datasets
area: datatypes
description: Format help for Tabular/BED/Interval Datasets
layout: faq
box_type: tip
contributors: [jennaj, beachyesh]
---


A Tabular datatype is human readable and has tabs separating data columns. Please note that tabular data is different from comma separated data (.csv) and the common datatypes are: `.bed`, `.gtf`, `.interval`, or `.txt`.
1. Click the pencil icon {% icon galaxy-pencil %} to reach the **_Edit Attributes_** form.
   1. Change the datatype (3rd tab) and save.
   2. Label columns (1st tab) and save.
   3. Metadata will be assigned, then the dataset can be used.
2. If the required input is a BED or Interval datatype, adjusting (``.tab`` → ``.bed``, ``.tab`` → ``.interval``) maybe possible using a combination of **_Text Manipulation_** tools, to create a dataset that matches required specifications.
3. Some tools require that BED format be followed, even if the datatype Interval (with less strict column ordering) is accepted on the tool form.
   - These tools will fail, if they are run with malformed BED datasets or non-specific column assignments.
   - **Solution**: reorganize the data to be in BED format and rerun.
