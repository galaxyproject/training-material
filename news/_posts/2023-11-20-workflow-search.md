---
title: "Update: Workflow List now searches WorkflowHub.eu, advanced query syntax"
contributions:
  authorship: [hexylena]
  testing: [paulzierep, wm75]
tags: [feature update, gtn]
layout: news
---

## WorkflowHub

We have now added support for [WorkflowHub.eu](https://workflowhub.eu) in our [cross-galaxy workflow search]({% link workflows/list.html %}) interface that lets you find workflows from around the universe. The support for WorkflowHub helps us showcase all of the best-practice workflows currently available there! These can all be imported and run directly in Galaxy.

<iframe src="{% link workflows/embed.html %}?all=vgp+Delphine+Lariviere" height="300px" width="100%" class="gtn-embed"></iframe>

The "Load in Galaxy" button uses our [my.galaxy.training]({{ site.baseurl }}/news/2023/04/20/my-galaxy-training.html) service to let you choose which Galaxy server you're redirected to. These links will work for any recent Galaxy server.

## Querying

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
