---
redirect_from: [/faqs/galaxy/analysis_job_failure_walltime]
title: Understanding walltime error messages
area: troubleshooting
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

The full error message will be reported as below, and can be found by clicking on the bug icon for a failed job run (red dataset):
<pre>
job info:
This job was terminated because it ran longer than the maximum allowed job run time.
Please click the bug icon to report this problem if you need help.
</pre>
Or sometimes,
<pre>
job stderr:
slurmstepd: error: *** JOB XXXX ON XXXX CANCELLED AT 2019-XX-XXTXX:XX:XX DUE TO TIME LIMIT ***

job info:
Remote job server indicated a problem running or monitoring this job.
</pre>

- Causes:
    - The job execution time exceeded the "wall-time" on the cluster node that ran the job.
    - The server may be undergoing maintenance.
    - Very often input problems also cause this same error.
- Solutions:
    - Try at least one rerun.
    - Check the server homepage for banners or notices. Selected servers also post to the [Galaxy status page](https://status.galaxyproject.org/).
    - Review the Solutions section of the [Understanding input error messages]({% link faqs/galaxy/troubleshooting_input_problem.md %}) FAQ.
    - Your data may actually be too large to process at a public Galaxy server. Alternatives include setting up a [private Galaxy server]({% link faqs/gtn/galaxy_usage.md %}).
