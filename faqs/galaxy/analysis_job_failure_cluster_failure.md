---
title: Understanding 'canceled by admin' or a cluster failure in job failure
area: analysis
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

The initial error message could be: 
`This job failed because it was cancelled by an administrator.`
`Please click the bug icon to report this problem if you need help.`
Or
`job info:`
`Remote job server indicated a problem running or monitoring this job.`

- Causes:
    - Server or cluster error.
    - Less frequently, input problems are a factor.

- Solutions:
    - Try at least one rerun. Server/cluster errors like this are usually transient. 
    - Review the Solution section for the Input problems above.
    - If after any corrections, the job still fails, please report the technical issue [following the extended issue guidelines](https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_reporting_issues.html).