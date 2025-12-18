---
layout: tutorial_hands_on

title: 'Remove contamination and host reads'
zenodo_link: https://zenodo.org/records/17829290
questions:
- What preprocessing steps are required to obtain cleaned reads for downstream analysis?
- How can we identify and remove contaminant or host-derived reads from raw sequencing data?
objectives:
- Identify reads originating from contaminants or host genomes.
- Remove those reads to produce high-quality, clean metagenomic data suitable for downstream analyses.
time_estimation: 1H
key_points:
- Identifying and removing contaminant and host reads is a critical preprocessing step in metagenomic workflows.
- Clean reads improve the accuracy of downstream assembly, binning, and taxonomic profiling.
contributions:
  authorship:
    - minamehr    
    - bebatut
---


Metagenomic sequencing captures all DNA present in a sample, including the **microbial community**, **host DNA**, and potential **environmental or external contaminants** (such as human DNA introduced during sample handling or sequencing).  
Before performing taxonomic or functional analysis, it is essential to remove these host-drived or contaminant reads to avoid misleading downstream interpretations.

In this tutorial, we will learn how to identify and remove host or contaminant reads using Galaxy.
We will:
- Map raw reads to a **host reference genome** using Bowtie2 and extract unmapped reads.
- Repeat the process with unmapped reads against a **human reference genome** to remove potential human contamination.
- Generate a final set of **clean, non-host reads** suitable for downstream analyses such as assembly, binning, or profiling.

To test and illustrate the process, we use bee gut metagenome samples and focus on filtering out both host (bee) sequences and potential human contamination from the reads. 


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## Prepare Galaxy and data
Any analysis should get its own Galaxy history. So let's start by creating a new one:

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

Now, we need to import the data

> <hands-on-title>Import datasets</hands-on-title>
>
> 1. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library:
>
>    ```text
>    {{ page.zenodo_link }}/files/SRR24759598_1.fastq.gz
>    {{ page.zenodo_link }}/files/SRR24759598_2.fastq.gz
>    {{ page.zenodo_link }}/files/SRR24759616_1.fastq.gz
>    {{ page.zenodo_link }}/files/SRR24759616_2.fastq.gz
>    ```
>    
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 2. Create a paired collection named 'Raw reads' that includes both pairs.
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
{: .hands_on}

## Map reads to a host genome with Bowtie2
To remove host contamination, we start by mapping the reads to the host genome using Bowtie2 to detect and remove host-derived sequences.

> <hands-on-title>Remove host reads</hands-on-title>
>
> 1. {% tool [Bowtie2](toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.5.3+galaxy1) %} with the following parameters:
>    - *"Is this single or paired library"*: `Paired-end Dataset Collection`
>        - {% icon param-collection %} *"FASTQ Paired Dataset"*: `Input reads`
>        - *"Write unaligned reads (in fastq format) to separate file(s)"*: `Yes`
>        - *"Do you want to set paired-end options?"*: `Yes`
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>        - *"Select reference genome"*: `A. mellifera genome (apiMel3, Baylor HGSC Amel_3.0)`
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1: Default setting only`
>    - *"Do you want to tweak SAM/BAM Options?"*: `No`
>    - *"Save the bowtie2 mapping statistics to the history"*: `Yes`
>
> 2. Run the tool. The outputs will include:
>    - Mapping statistics report (`bowtie2.log`)
>    - Unaligned (unmapped) forward and reverse reads
>
> 3. These unmapped reads represent sequences **not belonging to the host** and will be used in the next step. 
>
> > <comment-title>Tip</comment-title>
> > Host reference genomes vary depending on the study organism. You can upload a FASTA file of your host genome if it is not available as a built-in index.
> {: .comment}
>
{: .hands_on}

## Re-pair unmapped reads
We now combine the unmapped forward and reverse reads into a new paired-end dataset for further processing.

> <hands-on-title> Combine unmapped forward and reverse reads into a paired collection </hands-on-title>
>
> 1. {% tool [Zip collections](__ZIP_COLLECTION__) %} with the following parameters:
>    - {% icon param-file %} *"Input 1"*: `output_unaligned_reads_l` 
>    - {% icon param-file %} *"Input 2"*: `output_unaligned_reads_r` 
>
>    2. This step creates a new paired-end collection that represents all reads **not aligned to the host genome**. Rename this output collection to 'Reads without bee reads'
>
>    > <comment-title> Note </comment-title>
>    >
>    > Zipping restores the normal paired-end structure, which is required for downstream tools or for rerunning the workflow on another reference.
>    {: .comment}
>
{: .hands_on}

## Summarize mapping statistics
Once the host mapping is complete, we use MultiQC to summarize and visualize the mapping statistics, helping us assess how many reads were removed and how many remain.

> <hands-on-title> Evaluate host read removal results </hands-on-title>
>
> 1. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.27+galaxy3) %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `Bowtie 2`
>                - {% icon param-file %} *"Output of Bowtie 2"*: `mapping_stats` (output of **Bowtie2** {% icon tool %})
>    - *"Report title"*: `Host Removal`
>
> 2. Run the tool and open the generated HTML report.
> 3. Review the mapping percentage, number of reads aligned, and number of unmapped reads.
>
>    > <comment-title>Tip</comment-title>
>    >
>    > Low mapping percentages in the report confirm that most host reads were successfully removed.
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What percentage of reads mapped to the bee reference genome?  
> 2. What percentage of reads were removed?
>
> > <solution-title></solution-title>
> >
> > 1. 2.8% for SRR24759598 and 1.1% for SRR24759616. 
> > 2. 2.8% for SRR24759598 and 1.1% for SRR24759616.
> >
> {: .solution}
>
{: .question}

## Remove potential human contamination
After removing host reads, we can run the **same workflow again** to eliminate possible **human contamination** that may remain in the dataset.

> <hands-on-title> Rerun the workflow using the human genome as reference</hands-on-title>
>
> 1. Use the unmapped reads **Reads without bee reads** as the input for this second run.
> 2. In the **Bowtie2** step:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>    - *"Select reference genome"*: `Human (GRCh38)`
>    - Keep all other parameters the same as in the first run.
>
> 3. Continue through the **Zip collections** and **MultiQC** steps as before.
> 4. The output of this second run represents your **final cleaned reads**, free from both host and human sequences. Rename this collection to 'Clean metagenomic reads'
>
>    > <comment-title>Note</comment-title>
>    > Rerunning the same workflow maintains reproducibility. 
>    > Only the reference genome and the input data change between the two runs.
>    {: .comment}
>
{: .hands_on}

> <question-title>Verify your final dataset</question-title>
>
> 1. How does the mapping percentage differ between the host and human filtering runs?
> 2. How many reads have been removed?
>
> > <solution-title></solution-title>
> >
> > 1. The second run (against the human genome) usually removes only a small number of additional reads.
> > 2. 0.1% in both samples.
> >
> {: .solution}
>
{: .question}


# Conclusion

In this tutorial, you learned how to:

- Identify and remove reads originating from host or contaminant genomes using **Bowtie2**.  
- Combine unmapped forward and reverse reads into a paired collection for reuse.  
- Summarize mapping statistics and verify host-read removal using **MultiQC**.  
- Rerun the same workflow with a human reference genome to remove residual human contamination.  

The resulting **clean reads** are now ready for downstream metagenomic analyses such as:
- **Assembly**   
- **Binning**   
- **Functional or taxonomic profiling** 

These preprocessing steps are essential to ensure accurate microbial community reconstruction without interference from host DNA.