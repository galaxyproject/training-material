---
redirect_from: [/faqs/galaxy/analysis_job_failure_cluster_failure]
title: Understanding 'canceled by admin' or cluster failure error messages
area: troubleshooting
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

The initial error message could be:
<pre>
This job failed because it was cancelled by an administrator.
Please click the bug icon to report this problem if you need help.
</pre>

Or
<pre>
job info:
Remote job server indicated a problem running or monitoring this job.
</pre>

- Causes:
    - Server or cluster error.
    - Less frequently, input problems are a factor.

- Solutions:
    - Try at least one rerun. Server/cluster errors like this are usually transient.
    - Review the Solutions section of the [Understanding input error messages](https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_job_failure_input_problem.html) FAQ.
    - If after any corrections, the job still fails, please report the technical issue [following the extended issue guidelines](https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_reporting_issues.html).
