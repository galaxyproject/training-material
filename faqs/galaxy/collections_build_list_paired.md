---
title: Creating a paired collection
area: collections
box_type: tip
layout: faq
contributors: [shiltemann, hexylena]
---


* Click on **Operations on multiple datasets** (check box icon) at the top of the history panel ![Operations on multiple datasets button]({% link topics/galaxy-interface/images/historyItemControls.png %})
* Check all the datasets in your history you would like to include
* Click **For all selected..** and choose **Build List of Dataset Pairs**

* Change the filter text in the text box next to *unpaired forward* to
{% if include.fw_filter %}`{{ include.fw_filter }}`{% else %}a common selector for the forward reads{% endif %}
* Change the filter text in the text box next to *unpaired reverse* to
{% if include.rv_filter %}`{{ include.rv_filter }}`{% else %}a common selector for the reverse reads{% endif %}
* Click **Pair these datasets** for each valid *forward* and *reverse* pair, or **Auto-pair** if all pairs look arranged as intended
* Enter a name for your collection
* Click **Create List** to build your collection
* Click on the checkmark icon at the top of your history again
