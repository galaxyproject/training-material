---
title: "New Feature: Automatic RMarkdown"
contributions:
  authorship: [hexylena]
  funding: [erasmusplus, avans-atgm]
tags: [gtn infrastructure, new feature]
cover: topics/data-science/images/rstudio/r-preview-output.png
coveralt: Image showing content from a tutorial, rendered as an rmarkdown html via knitting. A table of contents appears on the left, and code and outputs on the right.
layout: news
---

Further building on the work in the [automatic Jupyter Notebook]({% link news/_posts/2021-09-24-jupyter.md %}), we've now re-written the Jupyter export to be faster, and more importantly added support for R and RMarkdown! Check out an [example material]({% link topics/data-science/tutorials/r-dplyr/tutorial.md %}). Here we take the content of the tutorial, written like normal GTN markdown, and we automatically convert it to Jupyter Notebooks and now RMarkdown documents! [Check out the documentation]({% link topics/contributing/tutorials/create-new-tutorial-content/tutorial.md %}#automatic-jupyter-notebooks--rmarkdown) on how to setup your tutorials to support this.
