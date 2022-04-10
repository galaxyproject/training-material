---
title: Data Privacy Features 
area: History
box_type: tip
layout: faq
contributors: [desmond-jasper]
---

Privacy controls are only enabled if desired. Otherwise, datasets remain public but unlisted. 
This means that a dataset you've created is virtually invisible until you publish a link to it. 

below are steps to setting private Histories

1. Changing the privacy settings of individual dataset.
-  click on the **dataset** name for a dropdown.
-  clicking the '**pencil - {% icon galaxy-pencil %} icon**
- move on the **Permissions** tab.
- On the permission tab is two input tab
- On the second input with a label of **access** 
-  Search for the name of the user to grant permission
-  click on **save permision** 
-[watch animation](https://galaxyproject.org/learn/privacy-features/set-perm.gif)
- *note: **Adding additional roles to the 'access' permission along with your "private role" does not do what you may expect. 
Since roles are always logically ANDed together, only you will be able to access the dataset, since only you are a member of your "private role"**.*

2. Make all datasets in the current history private.
- Open the History Options galaxy-gear menu {% icon galaxy-gear %}at the top of your history panel
- click the **Make Private** option in the dropdown menu available 
- sets the default settings for all new datasets in this history to private.
- [watch animation](https://galaxyproject.org/learn/privacy-features/this-hist-priv-perm.gif)

3.Set the default privacy settings for new histories
-  click **user** button on top of the main channel for a dropdown{% icon galaxy-dropdown %} 
-  click on the *prefrences* under the dropdown {% icon galaxy-dropdown %}
-  select **Set Dataset Permissions for New Histories**icon {% icon cofest %}
- add a permission and click **save permission** 
- [watch animation](https://galaxyproject.org/learn/privacy-features/new-hist-perm.gif)
- *note: Changes made here will only affect histories created after these settings have been stored.*
