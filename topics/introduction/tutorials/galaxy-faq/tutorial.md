---
layout: tutorial_hands_on
topic_name: introduction
tutorial_name: galaxy-faq
---

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

> ### {% icon question %} I am seeing an error message in Galaxy; how do I report it?
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

> ### {% icon question %} Structural workflow stalled at second to last step; Create or Update Organism job failed.
>    > ### {% icon solution %} Solution
>    > ![](../../images/galaxy-faq-screenshots/2_job_failed_structural_workflow.png)
>    >
>    > This gives an error that says this in the preview and in the bug report:
>    >
>    > ![](../../images/galaxy-faq-screenshots/3_error_report_failed_structural.png)
>    >
>    > This is because you failed to provide a name for the organism before running the workflow. To fix this, rerun the tool. Enter the appropriate organism name. Also, check *Yes* for "Resume dependencies from this job."
>    >
>    > ![](../../images/galaxy-faq-screenshots/4_rerun_tool_adjustments.png)
> {: .solution}
{: .question}

> ### {% icon question %} Converting GFF3 of a genome with frameshifted genes into Genbank format files.
>    > ### {% icon solution %} Solution
>    > This usually fails when the genes (particularly the tape measure chaperones) in Apollo were not annotated with the correct attributes [put in link to frameshift tutorial]. It can be corrected by adding the attribute tag 'frameshift' and value 'a' to both the genes, then retrieving the data again and re-running the GFF33 to Genbank tool. The error will appear like this in the Galaxy history:
>    >
>    > ![](../../images/galaxy-faq-screenshots/5_gff3_genbank_error.png)
>    >
>    > The Galaxy bug report looks like this:
>    >
>    > ![](../../images/galaxy-faq-screenshots/6_gff3_genbank_bug_report.png)
>    >
>    > Go back to the genome in Apollo and edit the information for BOTH pieces of the frameshifted protein:
>    >
>    > ![](../../images/galaxy-faq-screenshots/7_edit_frameshift_product_information.png)
>    >
>    > If the attributes were not added (or did not save):
>    >
>    > ![](../../images/galaxy-faq-screenshots/8_attributes_before.png)
>    >
>    > The correct entry:
>    >
>    > ![](../../images/galaxy-faq-screenshots/9_attributes_corrected.png)
>    >
>    > Now the tool successfully completes, and the frameshift features are properly merged:
>    >
>    > ![](../../images/galaxy-faq-screenshots/10_tool_successful_features_merged.png)
> {: .solution}
{: .question}

