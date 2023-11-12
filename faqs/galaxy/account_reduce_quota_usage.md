---
redirect_from: [/faqs/galaxy/reduce_quota_usage]
title: How can I reduce quota usage while still retaining prior work (data, tools, methods)?
area: account
layout: faq
box_type: tip
contributors: [jennaj, bernandez]
---
* [**Download**]({% link faqs/galaxy/datasets_download_datasets.md %}) Datasets as individual files or entire Histories as an archive. Then *purge* them from the public server.
* **Transfer/Move** [Datasets]({% link faqs/galaxy/datasets_moving_datasets_between_galaxy_servers.md %}) or [Histories]({% link faqs/galaxy/histories_transfer_entire_histories_from_one_galaxy_server_to_another.md %}) to another Galaxy server, including your own Galaxy. Then purge.
* [**Copy**]({% link faqs/galaxy/histories_copy_dataset.md %}) your most important Datasets into a new/other History (inputs, results), then purge the original full History.
* [**Extract**]({% link faqs/galaxy/workflows_extract_from_history.md %}) a **Workflow** from the History, then purge it.
* **Back-up your work**. It is a best practice to download an archive of your FULL original Histories periodically, even those still in use, as a backup.

**Resources** Much discussion about all of the above options can be found at the [Galaxy Help forum](https://help.galaxyproject.org/).
