---
redirect_from: [/faqs/galaxy/analysis_job_failure_excess_memory]
title: Understanding 'exceeds memory allocation' error messages
area: troubleshooting
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

The error message to be displayed are as follows:
<pre>
job info:
This job was terminated because it used more memory than it was allocated.
Please click the bug icon to report this problem if you need help.
</pre>

Or
<pre>
stderr:
Fatal error: Exit code 1 ()
slurmstepd: error: Detected 1 oom-kill event(s) in step XXXXXXX.batch cgroup.
</pre>

Sometimes this message may appear at the bottom
<pre>
job stderr:
slurmstepd: error: Detected 1 oom-kill event(s) in step XXXXXXX.batch cgroup.
</pre>

In rare cases when the memory quota is exceeded very quickly, an error message such as the following can appear
<pre>
job stderr:
Fatal error: Exit code 1 ()
Traceback (most recent call last):
(other lines)
Memory Error
</pre>

**Note:** Job runtime memory is different from the amount of free storage space (quota) in an account.

- Causes:
    - The job ran out of memory while executing on the cluster node that ran the job.
    - The most common reasons for this error are input and tool parameters problems that must be adjusted/corrected.
- Solutions:
    - Try at least one rerun to execute the job on a different cluster node.
    - Review the Solutions section of the [Understanding input error messages](https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_job_failure_input_problem.html) FAQ.
    - Your data may actually be too large to process at a public Galaxy server. Alternatives include setting up a [private Galaxy server](https://training.galaxyproject.org/training-material/faqs/gtn/galaxy_usage.html).
