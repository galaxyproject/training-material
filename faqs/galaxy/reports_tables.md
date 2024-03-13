---
title: Enhancing tabular dataset previews in reports/pages
description: There are lots of fun advanced features!
area: reports
box_type: tip
layout: faq
contributors: [jmchilton, hexylena]
---

The existing `history_dataset_display` directive displays the dataset name and some useful context at the expense of potentially breaking the flow of the document

> > <code-in-title>Galaxy Markdown</code-in-title>
> > ````markdown
> > ```galaxy
> > history_dataset_display(history_dataset_id=1e8ab44153008be8) 
> > ```
> > ````
> {: .code-in}
>
> > <code-out-title>Example Screenshot</code-out-title>
> > ![a tabular dataset rendered, it has a title and a download button and sortable columns]({% link images/report_table_history_dataset_display.png %})
> {: .code-out}
{: .code-2col}




The existing `history_dataset_embedded` directive was implemented to try to inline results more and make the results more readable within a more... curated document. It is dispatches on tabular types and puts the results in a table but the table doesn't have a lot of options. 

> > <code-in-title>Galaxy Markdown</code-in-title>
> > ````markdown
> > ```galaxy
> > history_dataset_embedded(history_dataset_id=1e8ab44153008be8) 
> > ```
> > ````
> {: .code-in}
>
> > <code-out-title>Example Screenshot</code-out-title>
> > ![the same as before but no title nor download button. just a rendered table with sortable columns]({% link images/report_table_history_dataset_embedded.png %})
> {: .code-out}
{: .code-2col}

The `history_dataset_as_table` directive mirrors the `history_dataset_as_image` directive: it tries harder to coerce the data into a table and provides new table—specific options. The first of these is "show_column_headers` which defaults to `true`.



> > <code-in-title>Galaxy Markdown</code-in-title>
> > ````markdown
> > ```galaxy
> > history_dataset_as_table(history_dataset_id=1e8ab44153008be8,show_column_headers=false)
> > ```
> > ````
> {: .code-in}
>
> > <code-out-title>Example Screenshot</code-out-title>
> > ![the same as before but no title nor download button nor column headers]({% link images/report_table_history_dataset_as_table.png %})
> {: .code-out}
{: .code-2col}

There is also a `compact` option. This provides a much more inline experience for tabular datasets:



> > <code-in-title>Galaxy Markdown</code-in-title>
> > ````markdown
> > ```galaxy
> > history_dataset_as_table(history_dataset_id=1e8ab44153008be8,show_column_headers=false,compact=true)
> > ```
> > ````
> {: .code-in}
>
> > <code-out-title>Example Screenshot</code-out-title>
> > ![again the same screenshot, no table metadata, and now it lacks the small margin around it.]({% link images/report_table_history_dataset_compact.png %})
> {: .code-out}
{: .code-2col}

Figures in general should have titles and legends — so there is the "title" and "footer" options also.

> > <code-in-title>Galaxy Markdown</code-in-title>
> > ````markdown
> > ```galaxy
> > history_dataset_as_table(history_dataset_id=1e8ab44153008be8,show_column_headers=false,title='Binding Site Results',footer='Here is a very good figure caption for this table.')
> > ```
> > ````
> {: .code-in}
>
> > <code-out-title>Example Screenshot</code-out-title>
> > ![the same table with now a tasteful title and small caption below it describing that the author would write a caption if he knew what a binding site was.]({% link images/report_table_history_dataset_title.png %})
> {: .code-out}
{: .code-2col}
