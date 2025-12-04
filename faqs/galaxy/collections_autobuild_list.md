---
title: Creating a dataset collection with autobuild
area: collections
box_type: tip
layout: faq
contributors: [Marie59]
---

* Click on {% icon galaxy-selector %} **Select Items** at the top of the history panel ![Select Items button]({% link topics/climate/images/bgc_calib/collection_emptyhist.png %})
* Check {% if include.datasets_description %}{{ include.datasets_description }}{% else %}all the datasets in your history you would like to include{% endif %}
* Click **{% if include.n %}{{ include.n }}{% else %}n{% endif %} of N selected** and choose **Auto build List**

  ![Collection building with autobuild]({% link topics/climate/images/bgc_calib/build_collection.png %}){:width="15%"}

* Enter a name for your collection {% if include.name %} to {{ include.name }} {% endif %}
* Turn off **Remove file extension**

  ![Put a name and remove extension]({% link topics/climate/images/bgc_calib/collection_hist.png %}){:width="15%"}

* Click **Build** to build your collection
* Click on the checkmark icon at the top of your history again

Once the collection is created, all files turn green. You can limit visible files using the eye icons in the history panel.
