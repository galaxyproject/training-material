---
title: Missing or Incomplete Dataset metadata 
area: datasets   
description:  Finding and Correcting Metadata
layout: faq         
---

How to notice missing Dataset Metadata:
- Dataset will not be downloaded when using the disk icon
- Tools error when using a previously successfully used specific dataset
- Tools error with a message that ends with: ``OSError: [Errno 2] No such file or directory.``

Solution:
Click on the the dataset's {% icon galaxy-pencil %} ✏️ to reach the _Edit Attributes_ forms and click either:
- **Directly reset metadata** 
  - Find the tab for the metadata you want to change, make the change, and save.
- **Autodetect metadata**
  - Click on the _Auto-detect button_. The dataset will turn yellow in the history while the job is processing.




