---
title: Importing data from Sierra LIMS
layout: faq
area: data upload
box_type: tip
google_form_id: 1730224450
contributors: [hexhowells]
---
This section will guide you through generating external links to your data stored in the Sierra LIMS system to be downloaded directly into Galaxy.

1. Go to the [Sierra portal](https://www.bioinformatics.babraham.ac.uk/sierra/sierra.pl) and login to your account.
2. Click on the **Sample ID** of the sample you want to download data from.
3. Click on the **Edit Sample Details** button.
4. At the bottom of the page there will be an input box for creating a link, enter a description for the link in the **Reason for link** section, and click **Create link**. This will reload the page and add a new link to the sample under **Authorised links to this sample**.
5. Go back to the sample page or click on the hyperlink called **link** to take you back.
6. In the **Results** section select the lane you want to access your data from.
7. The bottom of the page, under the **Links** section, will now contain a list of `wget` commands with links for accessing all the files within that sample/lane.
8. Since this list is for `wget` commands, you need to extract out the links from the command. You can copy the link in the first set of double quotes for each line and {% icon galaxy-wf-edit %} **Paste/Fetch Data** them directly into Galaxy to download the files.
