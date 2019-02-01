---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: galaxy-cpt
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the tutorial (the questions and the objectives). It is nice also to have a scheme to sum up the pipeline used during the tutorial. The idea is to give to trainees insight into the content of the tutorial and the (theoretical and technical) key concepts they will learn.

> ### {% icon question %} One of my datasets turned red, what do I do?
>    > ### {% icon solution %} Solution
>    > * Sometimes tools fail and the output dataset turns red. A red dataset means that the tool has terminated with an error of some kind. If it was part of a **workflow** and downstream steps were waiting on the failed dataset, the workflow will not continue, and you may see one or more of the queued datasets remain in the paused (light blue) state.  Tools may return errors for a number of reasons, some of which the user can correct.
>    > * A common reason that tools fail is that the user specified an incorrect input dataset or tool parameters. Expand the dataset by clicking on the name. Click on the circular arrow "re-run" button to bring up the tool's original run paramters. Double-check that you selected the correct dataset(s) as input, and that you set any other parameters approriately, and try executing the tool again.
>    > * Click on the bug icon ![](../../images/galaxy-faq-screenshots/1_report_bug.png), and read the message that appears - it may give you a clue. [CPT IT staff](http://cpt.tamu.edu) will respond as soon as possible. **Do not rerun the job** unless you are changing some parameter to attempt to fix the problem. Otherwise, the history will be clogged with unnecessary, stalled jobs.
>    > * Examples of common errors include:
>    >    > * *User did not enter common name* - this is where you need to choose the phage name from the drop-down box, or type in a new name for your new organism (unique to Apollo).
>    >    > * *Wrong input* - triple check the file that the tool requires as input. If you have the wrong file type, search the list of tools for a converter.
>    >    > * *Wrong tool used* - read the blurb below each tool to make sure it does what you think it should.
> {: .solution}
{: .question}

> ### {% icon question %} I am seeing an error message in galaxy; how do I report it?
>    > ### {% icon solution %} Solution
>    > * If it isn’t a job-related bug, take a screenshots and [follow the directions here.](https://cpt.tamu.edu/computer-resources/github-repo-list/)
> {: .solution}
{: .question}

> ### {% icon question %} Nothing is working!
>    > ### {% icon solution %} Solution
>    > * Check to make sure you are logged in.
>    > * Check your internet connection.
>    > * Try logging into an incognito window.
> {: .solution}
{: .question}
