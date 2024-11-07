---
title: Supporting Tutorial Mode (GTN-in-Galaxy) in a tutorial
area: contributors
layout: faq
box_type: tip
contributors: [shiltemann]
---

GTN tutorials can be viewed directly withing Galaxy, we call this *Tutorial mode* ([read news post]({% link faqs/galaxy/tutorial_mode.md %}))

In this mode, tool names in hands-on boxes become clickable, directly opening the tool in Galaxy, at the right version.

To enable this feature, a bit of metadata needs to be added to the tool names in hands-on boxes as follows:

{% raw %}
```
{% tool [Name](Toolshed ID) %}
```
{% endraw %}

For example:

{% raw %}
```
{% tool [bedtools intersect intervals](toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_intersectbed/2.30.0+galaxy1) %}
```
{% endraw %}

To find the toolshed ID of a tool:

1. Open the tool in Galaxy
2. Click on the {% icon dropdown %} options menu (dropdown icon)
3. Select **Copy Tool ID**

  ![screenshot of menu with toolshed link for a tool]({% link faqs/gtn/images/tool_toolshed_link.png %})

Example of a hands-on box using this feature:

{% raw %}
```
> <hands-on-title> Counting SNPs </hands-on-title>
>
> 1. {% tool [Datamash](toolshed.g2.bx.psu.edu/repos/iuc/datamash_ops/datamash_ops/1.8+galaxy0) %} (operations on tabular data):
>
>    - *"Input tabular dataset"*: select the output dataset from **bedtools intersect intervals** {% icon tool %}
>    - *"Group by fields"*: `Column: 4` (the column with the exon IDs)
>
{: .hands_on}
```
{% endraw %}
