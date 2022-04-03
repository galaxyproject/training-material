---
title: Datatype Tabular
area: Getting inputs right    # FAQs will be grouped by these areas on the FAQ page
description: Format help for Tabular/BED/Interval Datasets
layout: faq          
contributors: [beachyesh] 
---

[Back to Support Hub](https://galaxyproject.org/support/)

Related Topics 
- [Fastq Manipulation and Quality Control](https://galaxyproject.org/tutorials/ngs/#fastq-manipulation-and-quality-control)  
- [How to format fastq data for tools that require .fastqsanger format?](https://galaxyproject.org/support/fastqsanger/)
- [Understanding compressed fastq data (fastq.gz)](https://galaxyproject.org/support/compressed-fastq/)
- [Common datatypes explained](https://galaxyproject.org/learn/datatypes/)
- [Input datatype misassignment and errors](https://galaxyproject.org/support/tool-error/)
 
### Help for Tabular Datasets

#### Tabular or Interval or BED or GFF or TXT or ???
- A [Tabluar](https://galaxyproject.org/learn/datatypes/#tabular) datatype is any that is human readable and has tabs seperating data columns.
 - Note: tabular data is different from comma seperated data (.csv)
- Common tabular datatypes are .bed, .gtf, .interval, or .txt.
- The datatype metadata attribute can often be directly reassigned to tabular format data.
- Click the {% icon galaxy-pencil %} icon to reach the **_Edit Attributes_** form. In the center panel, using tabs to navigate, change the datatype (3rd tab) and save, then label columns (1st tab) and save. Metadata will assign, then the dataset can be used.
- If the required input is a [BED](https://galaxyproject.org/learn/datatypes/#bed) or [Interval](https://galaxyproject.org/learn/datatypes/#interval) datatype, adjusting (.tab → .bed, .tab → .interval) may be possible using a combination of **_Text Manipulation_** tools, to create a [dataset](https://galaxyproject.org/learn/managing-datasets/) that matches the BED or Interval [datatype](https://galaxyproject.org/learn/datatypes/) specifications.

#### Tips
- Some tools require that BED format is followed, even if the datatype Interval (with less strict column ordering) is accepted on the tool form.
- These tools will fail if run with malformed BED datasets or non-specific column assignments.
- Solution: Reorganize the data to be in BED format and rerun. [More error troubleshooting help here](https://galaxyproject.org/support/tool-error/).
