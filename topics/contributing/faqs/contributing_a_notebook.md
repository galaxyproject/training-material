---
title: Contributing a Jupyter Notebook to the GTN
area: notebooks
box_type: tip
layout: faq
contributors: [hexylena, nsoranzo]
---

Problem: I have a notebook that I'd like to add to the GTN.

Solution: While we do not support directly adding notebooks to the GTN, as all of our notebooks are generated from the tutorial Markdown files, there is an alternative path! Instead you can:

1. Install [`jupytext`](https://pypi.org/project/jupytext/)
2. Use it to convert the ipynb file into a Markdown file (`jupytext notebook.ipynb --to markdown`)
3. Add this Markdown file to the GTN
4. Fix any missing header metadata

Then the GTN's infrastructure will automatically convert that Markdown file directly to a notebook on deployment. This approach has the advantage that Markdown files are more diff-friendly than ipynb, making it much easier to review updates to a tutorial.
