---
title: What information should I include when reporting a problem?
area: troubleshooting
layout: faq
box_type: tip
contributors: [hexylena]
---

Writing bug reports is a good skill to have as bioinformaticians, and a key point is that you should include enough information from the first message to help the process of resolving your issue more efficient and a better experience for everyone.

**What to include**

1. Which commands did you run, precisely, we want details. Which flags did you set?
2. Which server(s) did you run those commands on?
3. What account/username did you use?
4. Where did it go wrong?
5. What were the stdout/stderr of the tool that failed? Include the text.
6. Did you try any workarounds? What results did those produce?
7. (If relevant) screenshot(s) that show exactly the problem, if it cannot be described in text. Is there a details panel you could include too?
8. If there are job IDs, please include them as text so administrators don't have to manually transcribe the job ID in your picture.

It makes the process of answering 'bug reports' much smoother for us, as we will have to ask you these questions anyway. If you provide this information from the start, we can get straight to answering your question!

**What does a GOOD bug report look like?**

The people who provide support for Galaxy are largely volunteers in this community, so try and provide as much information up front to avoid wasting their time:

> I encountered an issue: I was working on (this server> and trying to run (tool)+(version number) but all of the output files were empty. My username is jane-doe.
>
> Here is everything that I know:
>
> - The dataset is green, the job did not fail
> - This is the standard output/error of the tool that I found in the information page (insert it here)
> - I have read it but I do not understand what X/Y means.
> - The job ID from the output information page is 123123abdef.
> - I tried re-running the job and changing parameter Z but it did not change the result.
>
> Could you help me?
{: .quote}
