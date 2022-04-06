---
title: Indentifying and Formatting Tabular Dataset 
area: Datatypes
description: Format help for Tabular/BED/Interval Datasets
layout: faq          
---
 

- A Tabular datatype is human readable and has tabs separating data columns.
  - Note: tabular data is different from comma seperated data (.csv)
1. Common tabular datatypes are `.bed`, `.gtf`, `.interval`, or `.txt`.
2. The datatype metadata attribute can be directly reassigned to tabular format data.
3. Click the {% icon galaxy-pencil %} ✏️ icon to reach the **_Edit Attributes_** form. 
   - Change the datatype (3rd tab) and save
   - label columns (1st tab) and save
   -  Metadata will be assigned, then the dataset can be used.
4. If the required input is a BED or Interval datatype, adjusting (``.tab`` → ``.bed``, ``.tab`` → ``.interval``) maybe possible using a combination of **_Text Manipulation_** tools, to create a dataset that matches required specifications.
5. Some tools require that BED format be followed, even if the datatype Interval (with less strict column ordering) is accepted on the tool form.
   - These tools will fail, if they are run with malformed BED datasets or non-specific column assignments.
   - **Solution** is to reorganize the data to be in BED format and rerun. 
