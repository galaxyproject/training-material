---
title: Enhancing tabular dataset previews in reports/pages
description: There are lots of fun advanced features!
area: reports
box_type: tip
layout: faq
contributors: [jmchilton, hexylena]
---

There are a number of options, specifically for tabular data, that can allow it to render more nicely in your workflow reports and pages and anywhere that GalaxyMarkdown is used.

- `title` to give your table a title
- `footer` allows you to caption your table
- `show_column_headers=false` to hide the column headers
- `compact=true` to make the table show up more inline, hiding that it was embedded from a Galaxy dataset.

The existing `history_dataset_display` directive displays the dataset name and some useful context at the expense of potentially breaking the flow of the document

> <code-in-title>Galaxy Markdown</code-in-title>
> ````markdown
> ```galaxy
> history_dataset_display(history_dataset_id=1e8ab44153008be8) 
> ```
> ````
{: .code-in}

> <code-out-title>Example Screenshot</code-out-title>
> ![a tabular dataset rendered, it has a title and a download button and sortable columns]({% link faqs/galaxy/images/report_table_history_dataset_display.png %})
{: .code-out}


The existing `history_dataset_embedded` directive was implemented to try to inline results more and make the results more readable within a more... curated document. It is dispatches on tabular types and puts the results in a table but the table doesn't have a lot of options. 

> <code-in-title>Galaxy Markdown</code-in-title>
> ````markdown
> ```galaxy
> history_dataset_embedded(history_dataset_id=1e8ab44153008be8) 
> ```
> ````
{: .code-in}

> <code-out-title>Example Screenshot</code-out-title>
> ![the same as before but no title nor download button. just a rendered table with sortable columns]({% link faqs/galaxy/images/report_table_history_dataset_embedded.png %})
{: .code-out}


The `history_dataset_as_table` directive mirrors the `history_dataset_as_image` directive: it tries harder to coerce the data into a table and provides new table—specific options. The first of these is "show_column_headers` which defaults to `true`.



> <code-in-title>Galaxy Markdown</code-in-title>
> ````markdown
> ```galaxy
> history_dataset_as_table(history_dataset_id=1e8ab44153008be8,show_column_headers=false)
> ```
> ````
{: .code-in}

> <code-out-title>Example Screenshot</code-out-title>
> ![the same as before but no title nor download button nor column headers]({% link faqs/galaxy/images/report_table_history_dataset_as_table.png %})
{: .code-out}


There is also a `compact` option. This provides a much more inline experience for tabular datasets:

> <code-in-title>Galaxy Markdown</code-in-title>
> ````markdown
> ```galaxy
> history_dataset_as_table(history_dataset_id=1e8ab44153008be8,show_column_headers=false,compact=true)
> ```
> ````
{: .code-in}

> <code-out-title>Example Screenshot</code-out-title>
> ![again the same screenshot, no table metadata, and now it lacks the small margin around it.]({% link faqs/galaxy/images/report_table_history_dataset_compact.png %})
{: .code-out}


Figures in general should have titles and legends — so there is the "title" and "footer" options also.

> <code-in-title>Galaxy Markdown</code-in-title>
> ````markdown
> ```galaxy
> history_dataset_as_table(history_dataset_id=1e8ab44153008be8,show_column_headers=false,title='Binding Site Results',footer='Here is a very good figure caption for this table.')
> ```
> ````
{: .code-in}

> <code-out-title>Example Screenshot</code-out-title>
> ![the same table with now a tasteful title and small caption below it describing that the author would write a caption if he knew what a binding site was.]({% link faqs/galaxy/images/report_table_history_dataset_title.png %})
{: .code-out}

