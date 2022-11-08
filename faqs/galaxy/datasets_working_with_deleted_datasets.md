---
title: Working with deleted datasets
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

Deleted datasets and histories can be recovered by users as they are retained in Galaxy for a time period set by the instance administrator. Deleted datasets can be undeleted or permanently deleted within a History. _Links to show/hide deleted (and hidden) datasets are at the top of the History panel._

- To review or adjust an individual dataset:
    1. Click on the name to expand it.
    2. If it is only deleted, but not permanently deleted, you'll see a message with links to recover or to purge.
        - Click on _Undelete_ it to recover the dataset, making it active and accessible to tools again.
        - Click on _Permanently remove it from disk_ to purge the dataset and remove it from the account quota calculation.

- To review or adjust multiple datasets in batch:
    1. Click on the checked box icon {% icon galaxy-selector %} near the top right of the history panel to switch into "Operations on Multiple Datasets" mode.
    2. Accordingly for each individual dataset, choose the selection box. Check the datasets you want to modify and choose your option (show, hide, delete, undelete, purge, and group datasets).
