---
title: How can I Reduce quota usage while still retaining prior work (data, tools, methods)?
area: Account
layout: faq
box_type: tip
contributors: [jennaj, bernandez]
---


* *There are several options to preserve data/work while also reducing your quota usage.* Some examples are below.
* **Download** Datasets as individual files or entire Histories as an archive. Then *purge* them from the public server.
* **Copy** Datasets or Histories to another Galaxy server, including your own Galaxy. Then purge.
* **Copy** your most important Datasets into a new/other History (inputs, results), then purge the original full History. This can be a quicker alternative if the orignal History is very large or complex to review/purge dataset-by-dataset.
* Extract a **Workflow** from the History, then purge it. Workflows do not consume quota space as they do not contain any datasets, only tools/parameters. IF you first saved back the original inputs and final result datasets, the analysis can be rerun to recreate the intermediate datasets and to compare the new result with your original. Be aware that minor differences are expected over time, as tools and Galaxy itself are updated.
* **Back-up your work**. It is a best practice to download an archive of your FULL original Histories periodically, even those still in use, as a backup.
* **In short, keep active work on the server and archive/download older work + backups of current work.** You can always reupload prior downloaded Datasets or Histories later. Or, rerun an analysis with a Workflow and the original inputs.
* **Resources** Much discussion about all of the above options can be found at the Galaxy Help forum https://help.galaxyproject.org/ (in context with difference use cases). The [Support FAQs](/support/) and [Admin Docs](https://docs.galaxyproject.org/) are more places to find help. If any option is unclear or you run into problems, ask for more help! [https://galaxyproject.org/support/#help-resources](/support/#help-resources)
