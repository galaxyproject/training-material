---
title: Get the workflow invocation
area: workflows
box_type: tip
layout: faq
contributors: [bebatut]
---

- Go to the workflow invocations page
    - Before Galaxy 24.0: Go to **User** > **Workflow Invocations**
    - In Galaxy 24.0: Go to **Data** > **Workflow Invocations**
    - Above Galaxy 24.1: Go to **Workflow Invocation** in the activity bar on the left

- Open the most recent item
- Find the invocation id:
    - Below 24.0, you can get it here:

        ![The image depicts a user interface from a computational workflow or job management system. At the top left, there is a "View Report" button, which likely allows users to access a detailed report of the job or workflow. The interface features two horizontal progress bars: the first indicates that all 4 steps have been successfully scheduled, and the second shows that all 3 jobs are complete. A "Download BioCompute Object" button is present, presumably for downloading a standardized description of the bioinformatics protocol or workflow. The interface includes collapsible sections labeled "Inputs," "Outputs," and "Steps," which likely contain detailed information about the inputs used, outputs generated, and steps involved in the workflow, respectively. At the top right, an invoice number is displayed: "Invoice: 6e32d21c3708b2b6." This interface is part of a system used for managing and monitoring computational jobs, particularly in bioinformatics or computational biology contexts.]({% link faqs/galaxy/images/workflow-invocation-below-24-0.png %})

    - Above Galaxy 24.1 (activity bar), you can find the workflow invocation id from the URL. For example, `https://usegalaxy.org/workflows/invocations/be5c48c113145dd5` means that the workflow invocation id is `be5c48c113145dd5`.
