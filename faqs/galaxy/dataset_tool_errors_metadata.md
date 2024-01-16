---
title: How to find and correct tool errors related to Metadata?
area: troubleshooting
description:  Finding and Correcting Metadata
layout: faq
box_type: tip
contributors: [jennaj, beachyesh]
---

Tools can error when the wrong dataset attributes (metadata) are assigned. Some of these wrong assignments may be:
 - Tool outputs, which are automatically assigned without user action.
 - Incorrect autodetection of datatypes, which need manual modification.
 - Undetected attributes, which require user action (example: assigning database to newly uploaded data).

How to notice missing Dataset Metadata:
- Dataset will not be downloaded when using the disk icon {% icon galaxy-save %}.
- Tools error when using a previously successfully used specific dataset.
- Tools error with a message that ends with: ``OSError: [Errno 2] No such file or directory``.

Solution:

Click on the dataset's pencil icon {% icon galaxy-pencil %} to reach the _Edit Attributes_ forms and do one of the following as applies:
- **Directly reset metadata**
  - Find the tab for the metadata you want to change, make the change, and save.
- **Autodetect metadata**
  - Click on the _Auto-detect button_. The dataset will turn yellow in the history while the job is processing.
