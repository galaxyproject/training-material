---
title: "Update: Embeddable UseGalaxy Workflow List now includes searches WorkflowHub.eu"
contributions:
  authorship: [hexylena]
  testing: [paulzierep]
tags: [feature update, gtn]
layout: news
---

Based on [a request from Paul](https://github.com/galaxyproject/training-material/issues/4494), more advanced querying was needed.

As such we've added a couple of alternative query parameters that you may use:

Query Parameter | Example Values | Interpretation
--- | --- | ---
`?query=` | `?query=single-cell` | This will search for the text `singlecell` anywhere in the conjoined fields of title, tags, and authors. This is similar to `all` but will look for exact phrases including spaces, rather than splitting up queries word-by-word.
`?all=` | `?all=longreads+microbiome` | Each term, separated by a `+` (which is interpreted as a ` ` space character in URL parsing), must appear somewhere in those fields. If any term is missing, that workflow will not be included
`?any=` | `?any=longreads+shortread` | As long as one of these terms is present, the result will be shown. So long or short read workflows will match this query. (Assuming they're tagged properly!)
`?none=` | `?none=testing+training` | If you want to exclude one or more terms, you can list them here.

Try some examples here:
- [microbiome tutorials, without nanopore]({% link workflows/embed.html %}?all=microbiome&none=nanopore)
- [single-cell tutorials using the 'old' text query]({% link workflows/embed.html %}?query=single-cell)
- [single-cell tutorials using the new ?all query (it's the same.)]({% link workflows/embed.html %}?all=single-cell)
- [Just Paul's workflows]({% link workflows/embed.html %}?all=paulzierep)

If you have any ideas for new features or improvements, please [open an issue on GitHub](https://github.com/galaxyproject/training-material/issues/)
