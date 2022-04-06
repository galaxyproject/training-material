---
title: Datatype Tabular
area: Datatypes
description: Format help for Tabular/BED/Interval Datasets
layout: faq          
---
 

- A Tabular datatype is human readable and has tabs separating data columns.
  - Note: tabular data is different from comma seperated data (.csv)
- Common tabular datatypes are `.bed`, `.gtf`, `.interval`, or `.txt`.
- The datatype metadata attribute can be directly reassigned to tabular format data.
- Click the {% icon galaxy-pencil %} icon to reach the **_Edit Attributes_** form. 
  - Change the datatype (3rd tab) and save, then label columns (1st tab) and save.  Metadata will be assigned, then the dataset can be used.
- If the required input is a BED or Interval datatype, adjusting (.tab → .bed, .tab → .interval) maybe possible using a combination of **_Text Manipulation_** tools, to create a dataset that matches the BED or Interval datatype specifications.

#### Tips
- Some tools require that BED format is followed, even if the datatype Interval (with less strict column ordering) is accepted on the tool form.
- These tools will fail if run with malformed BED datasets or non-specific column assignments.

**Solution**: Reorganize the data to be in BED format and rerun. 
