---
title: Understanding job statuses
description: Job statuses will help you understand the stages of your work.
area: datasets
layout: faq
box_type: tip
contributors: [jennaj, garimavs]
---

The following job statuses will help you better understand the working stage of the process.

- **Green**: The job was completed successfully.
- **Yellow**: The job is executing. Allow this to complete! Should they run longer, they will fail with a "wall-time" error and turn _red_.
- **Grey**: The job is being evaluated to run (new dataset) or is queued. Allow this to complete.
- **Red**: The job has failed.
- **Light Blue**: The job is paused. This indicates either an input has a problem or that you have exceeded the disk quota set by the administrator of the Galaxy instance you are working on.
- **Grey, Yellow, Grey again**: The job is waiting to run due to admin re-run or an automatic fail-over to a longer-running cluster.
- **Bright blue with moving arrow**: May be found in earlier Galaxy versions. Applies to the "Get Data → Upload File" tool only - the upload job is queuing or running.

**It is _essential_ to allow queued jobs to remain queued and not delete/re-run them.**
