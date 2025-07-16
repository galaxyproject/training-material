---
title: Creating a paired collection
area: collections
box_type: tip
layout: faq
contributors: [shiltemann, hexylena,lldelisle, pavanvidem]
---


* Click on {% icon galaxy-selector %} **Select Items** at the top of the history panel ![Select Items button]({% link topics/galaxy-interface/images/historyItemControls.png %})
* Check {% if include.datasets_description %}{{ include.datasets_description }}{% else %}all the datasets in your history you would like to include{% endif %}
* Click **{% if include.n %}{{ include.n }}{% else %}n{% endif %} of N selected** and choose **Advanced Build List**

  ![build paired collection menu item]({% link topics/galaxy-interface/images/buildList.png %}){:width="15%"}

* You are in collection building wizard. Choose **Flat List** and click 'Next' button at the right bottom corner.

  ![collection building wizard paired list]({% link topics/galaxy-interface/images/paired_list_selector.png %}){:width="15%"}

* Check and configure auto-pairing. Commonly matepairs have suffix `_1 ` and `_2` or `_R1` and `_R2`. Click on 'Next' at the bottom.

  ![edit and build a paired list collection]({% link topics/galaxy-interface/images/paired_list_edit.png %}){:width="15%"}

* Edit the  **List Identifier** as required.
* Enter a name for your collection
* Click **Build** to build your collection
* Click on the checkmark icon at the top of your history again