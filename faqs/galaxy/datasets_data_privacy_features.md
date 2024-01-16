---
title: How to set Data Privacy Features?
area: histories
box_type: tip
layout: faq
contributors: [jennaj,desmond-jasper]
---

Privacy controls are only enabled if desired. Otherwise, datasets by defaults remain private and unlisted in Galaxy. This means that a dataset you've created is virtually invisible until you publish a link to it.

Below are three optional steps to setting private Histories, a user can make use of any of the options below depending on what the user want to achieve:

1. Changing the privacy settings of individual dataset.

   - Click on the **dataset** name for a dropdown.
   - Clicking the '**pencil - {% icon galaxy-pencil %} icon**
   - Move on the **Permissions** tab.
   - On the permission tab is two input tab
   - On the second input with a label of **access**
   - Search for the name of the user to grant permission
   - Click on **save permission**

   ![gif of the process described above, in Galaxy](https://galaxyproject.org/learn/privacy-features/set-perm.gif)

   *Note*: Adding additional roles to the 'access' permission along with your "private role" does not do what you may expect. Since roles are always logically added together, only you will be able to access the dataset, since only you are a member of your "private role".

2. Make all datasets in the current history private.

   - Open the History Options galaxy-gear menu {% icon galaxy-gear %} at the top of your history panel
   - Click the **Make Private** option in the dropdown menu available
   - Sets the default settings for all new datasets in this history to private.

   ![gif of the process described above, in Galaxy](https://galaxyproject.org/learn/privacy-features/this-hist-priv-perm.gif)

3. Set the default privacy settings for new histories

   - Click **user** button on top of the main channel for a dropdown {% icon galaxy-dropdown %}
   - Click on the *preferences* under the dropdown {% icon galaxy-dropdown %}
   - Select **Set Dataset Permissions for New Histories** icon {% icon cofest %}
   - Add a permission and click **save permission**

   ![gif of the process described above, in Galaxy](https://galaxyproject.org/learn/privacy-features/new-hist-perm.gif)

   *Note*: Changes made here will only affect histories created after these settings have been stored.
