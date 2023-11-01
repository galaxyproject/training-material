---
layout: tutorial_hands_on

title: 'Cite-Seq Tool Data Processing into RStudio Visualization (Cite-Seq, Seurat, R)'
subtopic: single-cell-CS-code
priority: 2

questions:
- How can I use Seurat's Cite-Seq capabilities?
- How can I visualize and interpret multimodal data in Seurat?

objectives:
- Learn to use Galaxy's Seurat tool with Cite-seq capabilities to create a Seurat Object 
- Understand the parameters of the Seurat tool 
- Move between Galaxy and RStudio to holistically explore Cite-Seq data

time_estimation: 3H

key_points:
- Being able to switch between Galaxy and RStudio when analyzing datasets can be useful when looking to adjust default parameters within Seurat's functions and workflow.
- Seurat in RStudio gives more flexibility and specificity of analyses, but Galaxy offers great reproducibility and ease of analysis.
- Beginning to learn the syntax and use of R will expand your access to bioinformatic analyses. 

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        
tags:
- Cite-Seq
- RStudio

contributions:
  authorship:
    - Camila-goclowski

notebook:
  language: r
  snippet: topics/single-cell/tutorials/scCiteSeq-Tool-to-RStudio/preamble.md
---

{% snippet topics/single-cell/faqs/notebook_warning.md %}

