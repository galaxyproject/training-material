---
title: Understanding 'exceeds memory allocation' in job failure
area: analysis
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

The error message to be displayed are as follows:
`job info:`
`This job was terminated because it used more memory than it was allocated.`
`Please click the bug icon to report this problem if you need help.`
Or
`stderr:`
`Fatal error: Exit code 1 ()`
`slurmstepd: error: Detected 1 oom-kill event(s) in step XXXXXXX.batch cgroup.`
Sometimes this message may appear at the bottom
`job stderr:`
`slurmstepd: error: Detected 1 oom-kill event(s) in step XXXXXXX.batch cgroup.`
In rare cases when the memory quota is exceeded very quickly, an error message such as the following can appear
`job stderr:`
`Fatal error: Exit code 1 ()`
`Traceback (most recent call last):`
`(other lines)`
`Memory Error`
 
**Note:** Job runtime memory is different from the amount of free storage space (quota) in an account.

- Causes:
    - The job ran out of memory while executing on the cluster node that ran the job.
    - The most common reasons for this error are input and tool parameters problems that must be adjusted/corrected.
- Solutions:
    - Try at least one rerun to execute the job on a different cluster node.
    - Review the Solution section for the Input problems above.
    - Your data may actually be too large to process at a public Galaxy server. Alternatives include setting up a private Galaxy server. 