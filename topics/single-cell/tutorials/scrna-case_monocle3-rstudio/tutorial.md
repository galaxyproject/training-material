---
layout: tutorial_hands_on

title: 'Trajectory Analysis: Monocle3 in RStudio'
subtopic: single-cell-CS
priority: 6
zenodo_link: 'https://zenodo.org/record/7455590'

questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?

objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives

time_estimation: 1H

key_points:
- The take-home messages
- They will appear at the end of the tutorial

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin
        - scrna-case_alevin-combine-datasets
        - scrna-case_basic-pipeline
        - scrna-case_JUPYTER-trajectories
        - scrna-case_monocle3-trajectories
        
tags:
- single-cell
- trajectory-analysis
- paper-replication

contributions:
  authorship:
    - wee-snufkin

notebook:
  language: r
  snippet: topics/single-cell/tutorials/scrna-case_monocle3-rstudio/preamble.md
---

```r
cells_annotated <- file.choose() 
genes_annotated <- file.choose() 
expression_annotated <- file.choose() 
```
# Working in RStudio

## Uploading files

If you are working in RStudio Cloud, you have to download the generated files from your history first. To do so, just click on the {% icon galaxy-save %} save icon for `Cell metadata (obs)`, `Gene metadata (var)` and `Expression matrix`. Then, return to the RStudio and click on ‘Upload’ button in the right bottom window toolbar and choose already downloaded files to upload. You should now see all three filed in this window. You might want to rename the files to make their names shorter.

![Screenshot of Files tab in RStudio, highlighting 'Upload' and 'Rename' buttons and listing three uploaded and renamed files: 'cell_metadata', 'gene_metadata', 'expression_matrix'.](../../images/scrna-casestudy-monocle/r_files_tab.png "The view of the Files tab with uploaded files and highlighted relevant buttons")

If you are 


## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

><hands-on-title> Data upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from
>
{: .hands_on}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
