---
redirect_from:
- /faqs/galaxy/analysis_job_failure_value_error
title: Understanding ValueError error messages
area: troubleshooting
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

The full error is usually a longer message seen only after clicking on the bug icon or by reviewing the job details `stderr`.

How to do both is covered in the Troubleshooting errors [FAQ]({% link faqs/galaxy/analysis_troubleshooting.md %}).

<pre>
stderr
...
Many lines of text, may include parameters
...
...
ValueError: invalid literal for int() with base 10: some-sequence-read-name
</pre>

- Causes:
    - MACS2 produces this error the first time it is run. MACS is not the only tool that can produce this issue, but it is the most common.
- Solutions:
    - Try at least one rerun.
    - MACS/2 is not capable of interpreting sequence read names with spaces included. Try following these two:
        - Remove unmapped reads from the SAM dataset. There are several filtering tools in the groups **SAMTools** and **Picard** that can do this.
        - Convert the SAM input to BAM format with the tool **SAMtools: SAM-to-BAM**. When compressed input is given to MACS, the spaces are no longer an issue.