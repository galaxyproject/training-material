---
title: Transfer entire histories from one Galaxy server to another
area: histories 
box_type: tip
layout: faq
contributors: [jennaj, AnomalyCodes]
---
1. Click on {% icon galaxy-gear %} in the history panel of the *sender* Galaxy server
2. Click on **Export to File**
3. Select either exporting history **to a link** or **to a remote file**
4. Click on the link text to generate a new archive for the history *if* exporting to a link
5. Wait for the link to generate
6. Copy the link address or click on the generated link to download the history archive
7. Click on **User** on the top menu of the *receiver* Galaxy server
8. Click on **Histories** to view saved histories
9. Click on **Import history** in the grey button on the top right
10. Select the appropriate importing method based on the choices made in steps 3 and 6
    - Choose **Export URL from another galaxy instance** if link address was copied in step 6
    - Select **Upload local file from your computer** if history archive was downloaded in step 6
    - Choose **Select a remote file** if history was exported to a remote file in step 3
11. Click the link text to check out your histories if import is successful


If history being transferred is too large, you may:
1. Click on {% icon galaxy-gear %} in the history panel of the *sender* Galaxy server
2. Click **Copy Datasets** to move just the important datasets into a new history
3. Create the archive from that smaller history

