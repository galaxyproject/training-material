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

> ### {% icon question %} Functional workflow failed at FASTA translate and TMHMM steps, or at steps involving GFF3 as input files.
>    > ### {% icon solution %} Solution
>    > The common cause for this problem is that the user picked the incorrect FASTA file or the GFF3 file as the input for the workflow.  The newly updated functional workflow (version v2020.07 and up) takes the "Annotation and Sequence" file (the GFF3 annotation combined with the DNA FASTA file) directly retrieved from Apollo as the input file. The workflow has a built-in step to split the combined file into two separate files, "GFF3 annotation from Apollo" and "FASTA sequence from Apollo", for downstream analysis steps.  If your functional workflow does not start with the "Annotation and Sequence" file directly retrieved from Apollo, or you are re-running certain step and you manually select the wrong GFF3 (or FASTA) file as input, your analysis will fail.  
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

> ### {% icon question %} PhageTerm broke on my data!
>    > ### {% icon solution %} Solution
>    > When running a [PhageTerm](https://cpt.tamu.edu/galaxy/root?tool_id=PhageTerm) job, the job fails. We do not have a good fix for this problem, as we didn't write the PhageTerm tool (it is published [here(https://www.ncbi.nlm.nih.gov/pubmed/28811656)). A couple common errors, and a few possible solutions are given below.
>    >
>    > Error 1:
>    >
>    > ![](../../images/galaxy-faq-screenshots/11_phageterm_error_1.png)
>    >
>    > ![](../../images/galaxy-faq-screenshots/12_phageterm_error_1_description.png)
>    >
>    > Error 2:
>    >
>    > ![](../../images/galaxy-faq-screenshots/13_phageterm_error_2.png)
>    >
>    > ![](../../images/galaxy-faq-screenshots/14_phageterm_error_2_description.png)
>    >
>    > Error 3 (related to using host genome in the optional input):
>    >
>    > ![](../../images/galaxy-faq-screenshots/15_phagerterm_error_3.png)
>    >
>    > **Possible solutions**:
>    > 1. Re-run [PhageTerm](https://cpt.tamu.edu/galaxy/root?tool_id=PhageTerm) in Galaxy, omitting an input for the optional parameters, including the corresponding paired-end reads and the host genome.
>    > 2. Re-open the genome in another place ([Re-open FASTA sequence](https://cpt.tamu.edu/galaxy/root?tool_id=edu.tamu.cpt.fasta.reopen)). Try running the job again.
>    > 3. It is possible that some phage genome data simply cannot be run through this tool in its current state. We are in communication with the authors to attempt to remedy the situation. If you have encountered a reproducible error when running PhageTerm on a particular dataset, please send us the information (cory.maughmer@tamu.edu) so that the developers can develop a fix for that specific problem.
> {: .solution}
{: .question}
> ### {% icon question %} Why won't my workflow run?
>    > ### {% icon solution %} Solution
>    > Sometimes when a worflow will not run, it is because tools in the Workflow may have new versions. **It is always a good practice to check "Workflows" section under "Shared Data" tab on the main page to get the newest version number and the last updated date info for public workflow, and re-import the newest workflow for your analysis**. If you run a workflow with tools that have had new parameters added, those must be manually set (even if they have a preset default value) the first time it is viewed in the workflow edit view. The workflow needs to be re-saved. Alternatively, for published workflows, they can be re-imported. 
>    >
>    > Here are some useful examples of warnings Galaxy will give in this case:
>    >
>    > Example 1: Galaxy warning messages indicating that tools have changed
>    >
>    > ![](../../images/galaxy-faq-screenshots/16_outdated_tool_warning.png)
>    >
>    >
>    > Example 2: Editor warnings upon the opening of a workflow with tools that have changed. The values mention here may need to be manually set prior to saving.
>    >
>    > ![](../../images/galaxy-faq-screenshots/17_editor_warnings.png)
>    >
>    >
>    > Example 3: If error messages are ignored and the user tries to run a workflow even with errors, Galaxy may give a "Workflow failed" message similar to this, which gives the first clue as to which tools are preventing the workflow from running.
>    >
>    > ![](../../images/galaxy-faq-screenshots/18_workflow_failed.png)
>    >
>    >
> {: .solution}
{: .question}
> ### {% icon question %} Why do my top protein or nucleotide hit lists contain duplicate entries? 
>    > **For Example:**
>    > ![](../../images/galaxy-faq-screenshots/19_top_protein_list_duplicates.png)
>    > ### {% icon solution %} Solution
>    > The duplicate entries are caused when processing of top hits from the BLASTp job separately counts organisms that have a unique accession, or multiple TaxIDs. For organisms with a representative genome in NCBI's RefSeq collection, this will result in duplicate organisms with identical TaxIDs, but unique accessions.(**Figure 1**) For some organisms with many representative genomes in the database, they will have been assigned multiple TaxIDs, each with a unique accession. Both these cases will result in what appear to be duplicates in the Top hits list. 
>    >
>    > Figure 1:
>    >
>    > ![](../../images/galaxy-faq-screenshots/20_top_hits.png)
>    >
>    > **Suggested Action:**
>    >
>    >Verify that entries are in fact, representing the same organism. The number of top hits displayed in the output list can be manually adjusted  when running the relatedness tool.
>    > 
>    >
> {: .solution}
{: .question}
> ### {% icon question %}Why is the dice score in the Top BLASTn results greater than 1?
>    > **For Example:**
>    > ![](../../images/galaxy-faq-screenshots/21_dice_score.png)
>    > ### {% icon solution %} Cause
>    > This happens when the aligned length is greater than the total length of query and/or subject length like in Line 1 of the above image.
>    >
>    > BLAST results contain full-length high-scoring pairs (HSPs), usually because the query is identical to the subject, plus additional duplicate smaller HSPs. The dice score calculation is (2 * # identical matches)/((genome 1 length) + (genome 2 length)). When adding up all the HSPs from a subject sequence, the additional smaller sequences push the value to slightly greater than one.
>    >
>    >
> {: .solution}
{: .question}
> ### {% icon question %} Why does the output from the GFF3 feature extraction tool have question marks up to the correct sequence length for each gene instead of the actual nucleotide sequence?
>    > 
>    > ### {% icon solution %} Cause
>    > The names of the GFF3 and the FASTA do not match. Thus, the downstream tools will run, but not reference the input fasta file to 'extract' features of the length called for by the gff3 file.
>    >
>    > **Suggested Action:**
>    >
>    > Change the names to match. Use the Rename tool for both GFF3 and FASTA files, or the FASTA sequence renamer tool.
>    >
>    > **Note:**
>    > The same symptom, problem, and solution occurs with ShineFind and other similar tools.
>    >
>    >
> {: .solution}
{: .question}
> ### {% icon question %} I need to know the parameters I used for an analysis, but I don't remember what I did. How can I reconstruct it?
>    > ### {% icon solution %} Solution
>    > **Suggested Action:**
>    > Check the parameters for the jobs in your history where the analysis was done. Dates and parameters can be checked by clicking on a dataset title to expand its options, then checking **'View details'** and/or **'Run this job again'**. This can be done by clicking on the circled i icon and/or cycling arrows icon respectively.
>    >
>    > ![](../../images/galaxy-faq-screenshots/22_galaxy_icons.PNG)
>    >
>    > **Note:**
>    > If you used an older published workflow from the CPT and you did not import it, email an admin for help determining what parameters or tools have been changed since then.
>    > 
> {: .solution}
{: .question}
