---
title: On the Scanpy PlotEmbed step, my object doesn’t have Il2ra or Cd8b1 or Cd8a etc.
box_type: question
layout: faq
contributors: [nomadscientist, mtekman]
redirect_from:
  - /topics/transcriptomics/tutorials/scrna-case_basic-pipeline/faqs/plotembed_results_missing
---

Check your `Anndata` object - it should be `7874 x 14832`, i.e. 7874 cells x 14832 genes. Is it actually 2000 genes only (i.e. and therefore missing the above markers)? You may have selected to remove genes at the Scanpy FindVariableGenes step (last toggle, ‘Remove genes not marked as highly variable’ < Select NO.) (Most likely you did this correctly the first time, but later in investigating how many got marked as highly variable, may have run this tool again and removed the nonvariable ones. We’ve updated the text to more clearly prevent this, but you may have gotten caught out!)




