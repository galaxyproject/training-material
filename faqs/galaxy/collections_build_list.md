---
title: Creating a dataset collection
area: collections
box_type: tip
layout: faq
contributors: [shiltemann, hexylena, pavanvidem]
---

* Click on {% icon galaxy-selector %} **Select Items** at the top of the history panel ![Select Items button]({% link topics/galaxy-interface/images/historyItemControls.png %})
* Check {% if include.datasets_description %}{{ include.datasets_description }}{% else %}all the datasets in your history you would like to include{% endif %}
* Click **{% if include.n %}{{ include.n }}{% else %}n{% endif %} of N selected** and choose **Advanced Build List**

  ![build list collection menu item]({% link topics/galaxy-interface/images/buildList.png %}){:width="15%"}

* You are in collection building wizard. Choose **Flat List** and click 'Next' button at the right bottom corner.

  ![collection building wizard flat list]({% link topics/galaxy-interface/images/flat_list_selector.png %}){:width="15%"}

* Double clcik on the file names to edit. For example, remove file extensions or common prefix/suffixes to cleanup the names.

  ![edit and build a list collection]({% link topics/galaxy-interface/images/flat_list_edit.png %}){:width="15%"}

* Enter a name for your collection {% if include.name %} to {{ include.name }} {% endif %}
* Click **Build** to build your collection
* Click on the checkmark icon at the top of your history again

