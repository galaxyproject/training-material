---
title: Troubleshooting Custom Genome fasta
area: reference genomes    
box_type: tip        
layout: faq        
contributors: [jennaj, Nurzhamalyrys] 
---


If a custom genome/transcriptome/exome dataset is producing errors, double check the format and that the chromosome identifiers between **ALL** inputs. Clicking on the bug icon {% icon galaxy-bug %} will often provide a description of the problem. This does not automatically submit a bug report, and it is not always necessary to do so, but it is a good way to get some information about [why a job is failing](https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_troubleshooting.html). 

- Custom genome not assigned as FASTA format

   - **Symptoms include**: Dataset not included in custom genome "From history" pull down menu on tool forms.
   - **Solution**: Check datatype assigned to dataset and assign **fasta** format.
   - **How**: Click on the dataset's pencil icon {% icon galaxy-pencil %} to reach the "Edit Attributes" form, and in the Datatypes tab > [redetect the datatype]({% link faqs/galaxy/datasets_detect_datatype.md %}).
   - If `fasta` is not assigned, there is a format problem to correct.
  
 - Incomplete Custom genome file load

   - **Symptoms include**: Tool errors result the first time you use the Custom genome.
   - **Solution**: Use **Text Manipulation → Select last lines from a dataset** to check last 10 lines to see if file is truncated.
   - **How**: Reload the dataset (switch to [FTP](https://galaxyproject.org/ftp-upload/) if not using already). Check your FTP client logs to make sure the load is complete.

- Extra spaces, extra lines, inconsistent line wrapping, or any deviation from strict FASTA format

   - **Symptoms include: RNA-seq tools (Cufflinks, Cuffcompare, Cuffmerge, Cuffdiff)** fails with error `Error: sequence lines in a FASTA record must have the same length!.`
   - **Solution**: File tested and corrected locally then re-upload or test/fix within Galaxy, then re-run.
   - **How**:
     - **Quick re-formatting** Run the tool through the tool **NormalizeFasta** using the options to wrap sequence lines at 80 bases and to trim the title line at the first whitespace.
     - **Optional Detailed re-formatting** Start with **FASTA manipulation → FASTA Width formatter** with a value between 40-80 (60 is common) to reformat wrapping. Next, use Filter and Sort → Select with ">" to examine identifiers. Use a combination of **Convert Formats → FASTA-to-Tabular, Text Manipulation** tools, then **Tabular-to-FASTA** to correct.
     - **With either of the above**, finish by using **Filter and Sort → Select** with `^\w*$` to search for empty lines (use "NOT matching" to remove these lines and output a properly format fasta dataset).
   
- Inconsistent line wrapping, common if merging chromosomes from various Genbank records (e.g. primary chroms with mito)

   - **Symptoms include**: Tools (**SAMTools, Extract Genomic DNA**, but rarely alignment tools) may complain about unexpected line lengths/missing identifiers. Or they may just fail for what appears to be a cluster error.
   - **Solution**: File tested and corrected locally then re-upload or test/fix within Galaxy.
   - **How**: Use **NormalizeFasta** using the options to wrap sequence lines at 80 bases and to trim the title line at the first whitespace. Finish by using **Filter and Sort → Select** with `^\w*$` to search for empty lines (use "NOT matching" to remove these lines and output a properly format fasta dataset).
   
- Unsorted fasta genome file

   - **Symptoms include**: Tools such as **Extract Genomic DNA** report problems with sequence lengths.
   - **Solution**: First try sorting and re-formatting in Galaxy then re-run.
   - **How**: To sort, follow instructions for [Sorting]({% link faqs/galaxy/reference_genomes_sorting_reference_genome.md %}) a Custom Genome.

- Identifier and Description in ">" title lines used inconsistently by tools in the same analysis

   - **Symptoms include**: Will generally manifest as a false genome-mismatch problem.
   - **Solution**: Remove the description content and re-run all tools/workflows that used this input. Mapping tools will usually not fail, but downstream tools will. When this comes up, it usually means that an analysis needs to be started over from the mapping step to correct the problems. No one enjoys redoing this work. Avoid the problems by formatting the genome, by double checking that the same reference genome was used for all steps, and by making certain the 'identifiers' are a match between all planned inputs (including reference annotation such as GTF data) **before using your custom genome**.
   - **How**: To drop the title line description content, use **NormalizeFasta** using the options to wrap sequence lines at 80 bases and to trim the title line at the first whitespace. Next, double check that the chromosome identifiers are an exact match between all inputs.
   
- Unassigned database

   - **Symptoms include**: Tools report that no build is available for the assigned reference genome.
   - **Solution**: This occurs with tools that require an assigned database metadata attribute. **SAMTools and Picard** often require this assignment.
   - **How**: Create a [Custom Build]({% link faqs/galaxy/analysis_add_custom_build.md %}) and assign it to the dataset.
