---
title: Understanding job statuses
description: Job statuses will help you understand the stages of your work.
area: datasets
layout: faq
box_type: tip
contributors: [jennaj, garimavs]
---

Compare the color of your datasets to these job processing stages.

- **Grey**: The job is queued. Allow this to complete!
- **Yellow**: The job is executing. Allow this to complete! 
- **Green**: The job has completed successfully.
- **Red**: The job has failed. Check your inputs and parameters with Help examples and GTN tutorials. Scroll to the bottom of the tool form to find these.
- **Light Blue**: The job is paused. This indicates either an input has a problem or that you have exceeded the disk quota set by the administrator of the Galaxy instance you are working on.
- **Grey, Yellow, Grey again**: The job is waiting to run due to admin re-run or an automatic fail-over to a longer-running cluster.

{% icon galaxy-info %}  Don't lose your queue placement! **It is _essential_ to allow queued jobs to remain queued, and to never interrupt an executing job.** If you delete/re-run jobs, they are added back to the _end of the queue again_. 

Related FAQs
- [Troubleshooting errors]({% link faqs/galaxy/analysis_troubleshooting.md %})
- [My jobs aren't running!]({% link faqs/galaxy/analysis_jobs_not_running.md %})
- [Extended Help for Differential Expression Analysis Tools]({% link faqs/galaxy/analysis_differential_expression_help.md %})
