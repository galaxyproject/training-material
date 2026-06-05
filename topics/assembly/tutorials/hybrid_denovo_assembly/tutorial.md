---
layout: tutorial_hands_on
title: Hybrid genome assembly - Nanopore and Illumina
zenodo_link: https://doi.org/10.5281/zenodo.15756327
tags:
  - assembly
  - nanopore
  - illumina
  - denovo
  - genomics
questions:
  - How do long- and short-read assembly methods differ?
objectives:
  - Understand how Nanopore and Illumina reads can be used together to produce a high-quality assembly
  - Become familiar with genome assembly and polishing programs
  - Learn how to assess the quality of a genome assembly, regardless of whether a reference genome is present or absent
  - Be able to assemble an unknown, previously undocumented genome to high-quality using Nanopore and Illumina reads!
time_estimation: 2h
key_points:
  - Performing de novo genome assembly with a combination of long and short reads results in a high quality genome
  - There are multiple methods for performing hybrid de novo assembly
  - One option is to assemble the long reads using `Flye` and then polish with the short reads using `Pilon`
  - Alternatively, `Unicycler` creates an assembly graph from the short reads and then untangles the graph with the long reads
contributions:
  authorship:
    - GraceAHall
    - tflowers15
  funding:
    - unimelb
    - melbournebioinformatics
    - AustralianBioCommons

---


*Assemble a genome! Learn how to create and assess genome assemblies using the powerful combination of Nanopore and Illumina reads.*

This tutorial explores how long and short read data can be combined to produce a high-quality 'finished' bacterial genome sequence. Termed 'hybrid assembly', we will use read data produced from two different sequencing platforms, Illumina (short read) and Oxford Nanopore Technologies (long read), to reconstruct a bacterial genome sequence.

In this tutorial, we will perform 'de novo assembly'. De novo assembly is the process of assembling a genome from scratch using only the sequenced reads as input - no reference genome is used.  This approach is common practise when working with microorganisms, and has seen increasing use for eukaryotes (including humans) in recent times.  

Using short read data (Illumina) alone for de novo assembly will produce a complete genome, but in pieces (commonly called a 'draft genome'). For the genome to be assembled into a single chromosome (plus a sequence for each plasmid), reads would need to be longer than the longest repeated element on the genome (usually ~7,000 base pairs, Note: Illumina reads are 350 base maximum). Draft bacterial genome sequences are cheap to produce (less than AUD\$60) and useful (>300,000 draft *Salmonella enterica* genome sequences published at NCBI [https://www.ncbi.nlm.nih.gov/pathogens/organisms/](https://www.ncbi.nlm.nih.gov/pathogens/organisms/)), but sometimes you need a high-quality 'finished' bacterial genome sequence.  There are <1,000 'finished' or 'closed' *Salmonella enterica* genome sequences.

In these cases, long reads can be used together with short reads to produce a high-quality assembly. Nanopore long reads (commonly >40,000 bases) can fully span repeats, and reveal how all of the genome fragments should be arranged. Long reads currently have higher error rates than short reads, so the combination of technologies is particularly powerful. Long reads provide information on the genome structure, and short reads provide high base-level accuracy.  

Combining read data from the long and short read sequencing platforms allows the production of a complete genome sequence with very few sequence errors, but the cost of the read data is about AUD$1,000 to produce the sequence. Understandably, we usually produce a draft genome sequence with very few sequence errors using the Illumina sequencing platform.

Nanopore sequencing technology is rapidly improving, expect the cost difference to reduce!!

- **Data:** Nanopore reads, Illumina reads, bacterial organism (*Bacillus subtilis*) reference genome
- **Tools:** Flye, Pilon, Unicycler, Quast, Busco
- **Pipeline:** Hybrid de novo genome assembly - Nanopore draft Illumina polishing
- **Pipeline:** Hybrid de novo genome assembly - Unicycler

This Galaxy tutorial based on material from the [Hybrid genome assembly - Nanopore and Illumina](https://www.melbournebioinformatics.org.au/tutorials/tutorials/hybrid_assembly/nanopore_assembly/) tutorial and workshop created by [Melbourne Bioinformatics](https://mdhs.unimelb.edu.au/melbournebioinformatics) at the University of Melbourne.


> <agenda-title></agenda-title>
> 
> In this tutorial we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background

## How do we produce the genomic DNA for a bacterial isolate?

Traditional *in vitro* culture techniques are important. Take a sample (e.g. a swab specimen from an infected sore) and streak a 'loopful' on to solid growth medium that suppoprts the growth of the bacteria. **Technology from the time of Louis Pasteur!**

Mixtures of bacterial types can be sequenced e.g. prepare genomic DNA from environmental samples containing bacteria - water, soil, faecal samples etc. (Whole Metagenome Sequencing)

![image of bacterial growth of a sample](../../images/denovo_assembly/sample_bacterial_growth.png)

One colony contains $$10^7 – 10^8$$ cells. The genomic DNA extracted from one colony is enough for Illumina sequencing. Larger amounts of genomic DNA are required for Nanopore sequencing.

## Shotgun sequencing - Illumina Sequencing Library

Genomic DNA is prepared for sequencing by fragmenting/shearing: multiple copies of Chromosome + plasmid  --> ~500 bp fragments

Note: Nanopore sequencing - there is usually no need to shear the genomic DNA as specialist methods are used to minimise shearing during DNA preparation. For Nanopore sequencing, the longer the DNA fragments, the better!

# Nanopore draft assembly, Illumina polishing

In this section, you will use `Flye` to create a draft genome assembly from Nanopore reads. We will perform assembly, then assess the quality of our assembly using two tools: `Quast` and `Busco`.

![workflow diagram for hybrid assembly starting with nanopore read](../../images/denovo_assembly/nanopore_illumina_hybrid_assembly.png)

## Upload data

Let's start with uploading the data.

> <hands-on-title>Import the data</hands-on-title>
> 
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Import from [Zenodo](https://zenodo.org/uploads/15756328):
>   
>    - FASTQ file with illumina forward reads: `illumina_reads_1.fastq`
>    - FASTQ file with illumina reverse reads: `illumina_reads_2.fastq`
>    - FASTQ file with nanopore reads: `nanopore_reads.fastq`
>    - FASTA file with reference genome: `reference_genome.fasta`
> 
>    ```
>    https://zenodo.org/records/15756328/files/illumina_reads_1.fastq
>    https://zenodo.org/records/15756328/files/illumina_reads_2.fastq
>    https://zenodo.org/records/15756328/files/nanopore_reads.fastq
>    https://zenodo.org/records/15756328/files/reference_genome.fasta
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
> 
{: .hands_on}


## A baseline for "high-quality" assemblies

To begin, we will identify what a high-quality assembly looks like.

When running assembly tools, we want to check the quality of assemblies we produce.
It is paramount that genome assemblies are high-quality for them to be useful. To get a baseline for what is considered a "high-quality" assembly, we will first run a common assembly QC tool - `Busco` - on a published genome similar to the organism we are working with today.

In your imported files, you should see a `reference_genome.fasta` item.
This is the published genome we will compare against.

### QC with Busco

BUSCO analysis uses the presence, absence, or fragmentation of key genes in an assembly to determine its quality.

BUSCO genes are specifically selected for each taxonomic clade, and represent a group of genes that each organism in the clade is expected to possess. At higher clades, 'housekeeping genes' are the only members, while at more refined taxa such as order or family, lineage-specific genes can also be used.

We expect the reference genome to have all of these genes. When running `Busco`, we expect it to find most (if not all) of these in the assembly.

Find and select the `Busco` tool in the tools panel using the search bar.
In this tutorial, we know our organism is within the 'Bacillales' order.

> <hands-on-title>Run Busco</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.8.0+galaxy1) %}:
>    - *"Sequences to analyse"*: `reference_genome.fasta`
>    - *"Mode"*:
>      - *"Select a gene predictor"*: `Augustus`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>      -  *"Lineage"*: `Bacillales`
>    - *"Which outputs should be generated"*: `Short summary text`
>   
> 2. ***Rename*** the `short summary` output to `Busco on Reference`.
> 
> 3. View the `Busco on Reference` output:
>    - It may look something like this:
>      ![output busco report for reference assembly](../../images/denovo_assembly/busco_reference_assembly.png)
> 
{: .hands_on}

> <question-title></question-title>
> 
> What can you conclude from the BUSCO report?
> 
> > <solution-title></solution-title>
> > 
> > - It seems that `Busco` could find almost all expected genes in the reference genome assembly.
> > - By looking at the results, we see that we have 449 / 450 Complete BUSCOs, and one Fragmented BUSCO.
> > - This will form the baseline for the `Busco` QC results expected of a high-quality genome assembly.
> > - From here, we will use our input DNA sequence data to assemble the genome of the sequenced organism, and will compare the QC results to that of the published `reference_genome.fasta` assembly.
> > 
> {: .solution}
> 
{: .question}


## Draft assembly with Flye + Nanopore reads

Our first assembly will use the long-read data to create a draft genome, then the short-read data to "polish" (improve) the draft into a better assembly.

We will start by using a long-read assembly tool called `Flye` to create an assembly using the Nanopore long-read data.

Search for `Flye` in the tools panel search bar and select the tool.

> <hands-on-title>Run Flye</hands-on-title>
>
> 1. {% tool [Flye](toolshed.g2.bx.psu.edu/repos/bgruening/flye/flye/2.9.6+galaxy0) %}:
>    - *"Input reads"*: `nanopore_reads.fastq`
>    - Leave all other parameters as default
>    - Scroll down and run Flye by clicking the blue `Run tool` button at the bottom of the page.
> 
> 2. View output:
>    - `Flye` produces a number of outputs. We only need the 'consensus' fasta file. You can delete the other outputs.
> 
> 3. ***Rename*** the `Flye on data XX: consensus` output to `Flye: Assembly`
> 
{: .hands_on}
 

## Assessing the Flye draft assembly quality

### Busco

We need to check if our assembly is good quality or not. It is paramount that genome assemblies are high-quality for them to be useful.

> <hands-on-title>Run Busco</hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.8.0+galaxy1) %}:
>    - *"Sequences to analyse"*: `Flye: Assembly`
>    - *"Mode"*:
>      - *"Select a gene predictor"*: `Augustus`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>      -  *"Lineage"*: `Bacillales`
>    - *"Which outputs should be generated"*: `Short summary text`
> 
> 2. ***Rename*** the `short summary` output to `Busco on Flye: Assembly`.
> 
> 3. View the `Busco on Flye: Assembly` output:
>    - It may look something like this:
>      ![output busco report for draft assembly](../../images/denovo_assembly/busco_draft_assembly.png)
>    - The `full table` is also useful. It gives a detailed list of the genes we are searching for, and information about whether they were missing, fragmented, or complete in our assembly.
>      ![output busco table for draft assembly](../../images/denovo_assembly/busco_table_draft_assembly.png)
> 
{: .hands_on}


> <question-title></question-title>
> 
> How does this compare with the reference genome BUSCO report?
> 
> > <solution-title></solution-title>
> > 
> > We see that many genes are fragmented or missing. Our draft genome assembly isn't as good as the reference genome ***yet***. This is because we have so far only used the Nanopore long-read sequences, which have higher base-level error rates than short reads. In the next steps we will use Illumina short reads to correct for errors in the assembled Nanopore reads.
> > 
> {: .solution}
> 
{: .question}


### Quast

Aside from `Busco`, we can use another method to perform assembly QC. `Quast` allows us to compare two assemblies to determine their similarity.

Although the organism we sequenced may be different, we can use `Quast` to compare our assembly with the provided reference genome to see how similar they are on the individual base level.

Search for the `Quast` tool in the tools panel.

> <hands-on-title>Run Quast</hands-on-title>
>
> 1. {% tool [Quast](toolshed.g2.bx.psu.edu/repos/iuc/quast/quast/5.3.0+galaxy0) %}:
>    - *"Assembly mode"*: `Individual assembly`
>    - *"Use customized names for the input files?"*: `No, use dataset names`
>    - *"Contigs/scaffolds file"*: `Flye: Assembly`
>    - *"Type of assembly?"*: `Genome`
>    - *"Use a reference genome?"*: `Yes`
>      - *"Reference genome"*: `reference_genome.fasta`
>    - *"Output files"*: Select `HTML reports`, `Tabular reports`
> 
> 2. ***Rename*** the `tabular report` output to `Quast on Flye: Assembly`.
> 
> 3. View output:
>    - `Quast` will produce a HTML report summarising it's results.
>    - Open the report. It may look something like this:
>      ![output quast report for reference assembly](../../images/denovo_assembly/quast_draft_assembly.png)
>    - Note the:
>      - Genome fraction (%)
>      - \# mismatches per 100 kbp
>      - \# indels per 100 kbp
>      - \# contigs information
> 
{: .hands_on}

> <question-title></question-title>
> 
> What can you conclude about our draft assembly?
> 
> > <solution-title></solution-title>
> > 
> > - From the output of `Quast`, our draft assembly seems to have good coverage and not too many contigs.
> > - Unfortunately, the mismatch / indel rate is quite high.
> > - Although we don't expect our organism to be ***identical*** to the supplied reference, we would expect fewer mismatches and indels, as the provided reference genome is ***very*** similar to organism which was sequenced.
> > 
> {: .solution}
> 
{: .question}


## Assembly Polishing with Pilon

We should be able improve our draft assembly with the Illumina reads available and correct some of these errors.

This process involves two steps. We will first align the Illumina reads to our draft assembly, then supply the mapping information to `Pilon`, which will use this alignment information to error-correct our assembly.

Illumina reads have much higher per-base accuracy than Nanopore reads. We will map the Illumina reads to our draft assembly using a short-read aligner called `BWA-MEM`. Then we can give `Pilon` this alignment file to polish our draft assembly.

### Map Illumina reads to draft assembly

Search for `Map with BWA-MEM` in the tools panel.

> <hands-on-title>Map with BWA-MEM</hands-on-title>
>
> 1. {% tool [Map with BWA-MEM](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.19) %}:
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a genome from history and build index`
>    - *"Use the following dataset as the reference sequence"*: `Flye: Assembly`
>    - *"Single or Paired-end reads"*: `Paired `
>    - *"Select first set of reads"*: `illumina_reads_1.fastq`
>    - *"Select second set of reads"*: `illumina_reads_2.fastq`
> 
> 2. View output:
>    - The output will be a BAM file (Binary Alignment Map). This is tabular data recording information about how reads were aligned to the draft assembly.
>    - We can now use this output BAM file as an input to `Pilon`.
> 
> 3. ***Rename*** the output to `Flye: Short read alignments`
> 
{: .hands_on}


### Polish assembly with Pilon

Search for `Pilon` in the tools panel.

> <hands-on-title>Run Pilon</hands-on-title>
>
> 1. {% tool [Pilon](toolshed.g2.bx.psu.edu/repos/iuc/pilon/pilon/1.20.1) %}:
>    - *"Source for reference genome used for BAM alignments"*: `Use a genome from History`
>      - *"Select a reference genome"*: `Flye: Assembly`
>    - *"Type automatically determined by pilon"*: 
>      - *"Input BAM file"*: `Flye: Short read alignments`
>    - *"Variant calling mode"*: `No`
> 
> 2. View output:
>    - `Pilon` gives a single output file - the polished assembly.
> 
> 3. ***Rename*** the output to `Flye: Polished assembly`
> 
{: .hands_on}



### Compare draft and polished assemblies

We are now interested to see how much `Pilon` improved our draft assembly.

> <hands-on-title>Run Quast</hands-on-title>
>
> Run `Quast` as before with the new `Flye: Polished assembly` data
> 
> 1. Select the history item for our initial `Quast` job, then click the rerun {% icon dataset-rerun %} button. This will load the settings used for the previous `Quast` job.
> 
> 2. Change the *"Contigs/scaffolds file"* input to `Flye: Polished assembly`.
> 
> 3. Click the `Run Tool` button to submit the job.
>
> 4. ***Rename*** the `tabular report` output to `Quast on Flye: Polished assembly`.
> 
> 5. After `Quast` has finished, open the HTML report. 
>    - Make note of `# mismatches per 100 kbp` and `# indels per 100 kbp`.
> 
{: .hands_on}

> <question-title></question-title>
> 
> Has our assembly improved?
>
> > <solution-title></solution-title>
> > 
> > Yes, we now have lower mismatches and indels per 100kbp.
> > 
> {: .solution}
> 
{: .question}

> <hands-on-title>Run Busco</hands-on-title>
>
> Run `Busco` as before with the new `Flye: Polished assembly` data
> 
> 1. Select the history item for our initial `Busco` job, then click the rerun {% icon dataset-rerun %} button. This will load the settings used for the previous `Busco` job.
> 
> 2. Change the *"Sequences to analyse"* input to `Flye: Polished assembly`.
> 
> 3. Click the `Run Tool` button to submit the job.
>
> 4. ***Rename*** the `short summary` output to `Busco on Flye: Polished assembly`.
> 
> 5. After `Busco` has finished, open the `Busco on Flye: Polished assembly` output.
> 
{: .hands_on}

> <question-title></question-title>
> 
> Have we identified more expected genes?
>
> > <solution-title></solution-title>
> > 
> > Yes, we now identified 423 complete BUSCO genes.
> > 
> {: .solution}
> 
{: .question}


## Section Summary

All going well, the polished assembly should be much higher quality than our draft.

The per-base accuracy of our assembly contigs should have markedly improved. This is reflected in the lower mismatches and indels per 100kbp reported by `Quast`, and the higher number of complete BUSCO genes. Our contiguity and coverage (as measured by the genome fraction (%) statistic reported by `Quast`) may not show the same level of improvement, as the polishing step is mainly aimed at improving per-base contig accuracy.

Our next step is to use a purpose-built hybrid de novo assembly tool, and compare its performance with our sequential draft + polishing approach.


## Section Questions

> <question-title></question-title>
> 
> Which read set - short or long - was used to create our draft?
>
> > <solution-title></solution-title>
> > 
> > Long reads (Nanopore) were used to create the draft. Long reads allow excellent re-creation of the proper structure of the genome, and adequately handle repeat regions. The drawback of long reads is a higher error rate of the technology compared to short reads. This results in more mismatches and indels.
> > 
> {: .solution}
> 
{: .question}

> <question-title></question-title>
> 
> How was the draft polished?
>
> > <solution-title></solution-title>
> > 
> > Illumina reads have higher per-base accuracy than Nanopore. Illumina reads were aligned to the draft assembly, then `Pilon` used this alignment information to improve locations with errors in the assembly.
> > 
> {: .solution}
> 
{: .question}

> <question-title></question-title>
> 
> How does `Quast` inform on assembly quality?
>
> > <solution-title></solution-title>
> > 
> > `Quast` shows summary information about the assembly contigs. If a reference genome is given, it informs the genome fraction (how much of the reference is covered by the assembly), if any genomic regions appear duplicated, and error information including the rate of mismatches and indels.
> > 
> {: .solution}
> 
{: .question}

> <question-title></question-title>
> 
> How does `Busco` inform on assembly quality?
>
> > <solution-title></solution-title>
> > 
> > `Busco` does not use a reference genome to compare. It attempts to locate key genes which should be present in the assembly, and reports whether it could/could not find those genes. If a key gene is found, it reports whether the gene was fragmented (errors) or complete.
> > 
> {: .solution}
> 
{: .question}



## Purpose-built hybrid assembly tool - Unicycler

In this section, we will use a purpose-built tool called `Unicycler` to perform hybrid assembly.

`Unicycler` uses our Nanopore and Illumina read sets together as input, and returns an assembly. Once we have created the assembly, we will assess its quality using `Quast` and `Busco` and compare with our previous polished assembly. We will also perform BUSCO analysis on the supplied reference genome itself, to record a baseline for our theoretical best `Busco` report.

![workflow diagram for hybrid assembly using unicycler](../../images/denovo_assembly/unicycler_hybrid_assembly.png)

### Hybrid de novo assembly with Unicycler

`Unicycler` performs assembly in the opposite manner to our previous approach. Illumina reads are used to create an assembly graph, then Nanopore reads are used to disentangle problems in the graph. The Nanopore reads serve to bridge Illumina contigs, and to reveal how the contigs are arranged sequentially in the genome.

### Run Unicycler

Find `Unicycler` in the tools panel. It is listed as `Create assemblies with Unicycler`

Run `Unicycler` using the Nanopore and Illumina read sets.

> <hands-on-title>Run Unicycler</hands-on-title>
>
> 1. {% tool [Unicycler](toolshed.g2.bx.psu.edu/repos/iuc/unicycler/unicycler/0.5.1+galaxy0) %}:
>    - *"Paired or Single end data?"*: `Paired `
>    - *"Select first set of reads"*: `illumina_reads_1.fastq`
>    - *"Select second set of reads"*: `illumina_reads_2.fastq`
>    - *"Select long reads. If there are no long reads, leave this empty"*: `nanopore_reads.fastq`
> 
> 2. View Output
>    - `Unicycler` will output three files 
>       - The assembly 
>       - An assembly graph 
>       - SPAdes graphs 
>    - We are interested in the `Final Assembly` output, which is the assembly as a `fasta` file. 
> 
> 3. ***Rename*** the `Final Assembly` output to `Unicycler: Assembly`.
> 
> > <comment-title></comment-title>
> >
> > If `nanopore_reads.fastq` does not appear in the dropdown list, its datatype needs to be changed. 
> >
> > - Click the pencil {% icon galaxy-pencil %} icon next to `nanopore_reads.fastq` in the history panel.
> > - Then select the `Datatypes` tab.
> > - For the parameter `New Type`: `fastqsanger`
> > - Leave all else default and execute the program.
> > 
> {: .comment}
> 
{: .hands_on}


### Comparing Unicycler assembly to Nanopore + Illumina polished assembly

`Busco` and `Quast` can be used again to assess our `Unicycler` assembly. As a purpose-built tool, it generally produces much better assemblies than our sequential approach. This is reflected as (`Quast`) a lower number of contigs, lower mismatches and indels per 100kb, and (`Busco`) greater number of BUSCO genes complete.

> <hands-on-title>Run Quast</hands-on-title>
>
> Run `Quast` as before with the new `Unicycler: Assembly` data
> 
> 1. Select the history item for our initial `Quast` job, then click the rerun {% icon dataset-rerun %} button. This will load the settings used for the previous `Quast` job.
> 
> 2. Change the *"Contigs/scaffolds file"* input to `Unicycler: Assembly`.
> 
> 3. Click the `Run Tool` button to submit the job.
>   
> 4. ***Rename*** the `tabular report` output to `Quast on Unicycler: Assembly`.
> 
> 5. At time of writing, these were the `Quast` results:
> 
> ![output quast report for unicycler assembly](../../images/denovo_assembly/quast_unicycler_assembly.png)
> 
{: .hands_on}

> <hands-on-title>Run Busco</hands-on-title>
>
> Run `Busco` as before with the new `Unicycler: Assembly` data
> 
> 1. Select the history item for our initial `Busco` job, then click the rerun {% icon dataset-rerun %} button. This will load the settings used for the previous `Busco` job.
> 
> 2. Change the *"Sequences to analyse"* input to `Unicycler: Assembly`.
> 
> 3. Click the `Run Tool` button to submit the job.
>
> 4. ***Rename*** the `short summary` output to `Busco on Unicycler: Assembly`.
> 
> 5. At time of writing, these were the `Busco` results:
> 
> ![output busco report for unicycler assembly](../../images/denovo_assembly/busco_unicycler_assembly.png)
> 
{: .hands_on}

We can now create an aggregate summary report from all of our `Quast` and `Busco` outputs using `MultiQC`, which enables us to easily compare the quality of our genome assemblies.

> <hands-on-title>Run MultiQC</hands-on-title>
>
> 1. {% tool [MultiQC](toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.27+galaxy4) %}:
>    - *"1: Results"*:
>      - *"Which tool was used generate logs?"*: `BUSCO`
>      - *"Output of BUSCO"*: `Busco on Reference`, `Busco on Fly: Assembly`, `Busco on Fly: Polished assembly`, `Busco on Unicycler: Assembly`, 
>    - Click `+ Insert Results` to add a second results field
>    - *"2: Results"*:
>      - *"Which tool was used generate logs?"*: `QUAST`
>      - *"Output of Quast"*: `Quast on Fly: Assembly`, `Quast on Fly: Polished assembly`, `Quast on Unicycler: Assembly`
>
> 2. View the `Webpage` output from `MultiQC`:
>    - At time of writing, these were the results.
>      - The number of BUSCOs identified in each assembly by `Busco`:
>        ![buscos identified in each assembly](../../images/denovo_assembly/busco_plot_bacillales_odb10.png)
>      - A table of `Quast` statistics for each assembly:
>        ![table of quast statistics](../../images/denovo_assembly/quast_table_statistics.png)
>      - The number of contigs identified in each assembly by `Quast`:
>        ![number of contigs in each assembly](../../images/denovo_assembly/quast_num_contigs.png) 
> 
{: .hands_on}


> <question-title></question-title>
> 
> How does the `Unicycler` assembly compare with the previous assembly produced using `Flye` + `Pilon`?
> 
> > <solution-title></solution-title>
> > 
> > It seems that the `Unicycler` assembly is much better than the `Flye` + `Pilon` polishing assembly, and produces the same report as the reference genome. Awesome! Looks like this is pretty good.
> > 
> > The `Unicycler` assembly has:
> > 
> > - A very high Genome fraction (coverage of the reference genome)
> > - Very few mismatches and indels per 100 kbp
> > - 15 contigs (compared with the reference genome only having 1 contig)
> > 
> {: .solution}
> 
{: .question}


These results indicate that our `Unicycler` assembly is high quality, and the organism is extremely similar to the organism of the reference genome.

To improve our assembly further to be publishable "complete" quality, we would want to reduce the number of contigs (ideally to a single contig). If the organism has plasmids, we would expect a handfull of contigs, likely below 10.

To do this, the best course of action would be to generate more long-read data.


### Section Questions

> <question-title></question-title>
> 
> Why did we select `Paired` for our Illumina reads in the `Unicycler` tool?
> 
> > <solution-title></solution-title>
> > 
> > Our short read set was `paired-end`. Short read technology can only sequence a few hundred base-pairs in a single read. To provide better structural information, paired-end sequencing was created, where longer fragments (fixed length) are used. A few hundred bp is sequenced at both ends of the fragment, leaving the middle section unsequenced. The reads produced (the mate-pair) from a single fragment are separated by a fixed length, so we know they are nearby in the genome.
> > 
> {: .solution}
> 
{: .question}

> <question-title></question-title>
> 
> Does `Unicycler` begin by using the long or short reads?
> 
> > <solution-title></solution-title>
> > 
> > Unicycler uses short reads first. It creates an assembly graph from short reads, then uses the long reads to provide better structural information of the genome.
> > 
> {: .solution}
> 
{: .question}

> <question-title></question-title>
> 
> How does `Unicycler` use long reads to improve its assembly graph?
> 
> > <solution-title></solution-title>
> > 
> > The assembly graph produced by short reads has tangled regions. When we don't know how sections of the genome are arranged, tangled regions appear in the graph. `Unicycler` uses Nanopore reads which overlap these tangled regions to resolve the proper structure of the genome.
> > 
> {: .solution}
> 
{: .question}


# Conclusion

We have covered two methods for hybrid de novo assembly. The combination of long- and short-read technology is clearly powerful, represented by our ability to create a good assembly with only 25x coverage (100Mb) of Nanopore and 50x coverage (200Mb) of Illumina reads.

To further improve our assembly, extra Nanopore read data may provide the most benefit. At 50x coverage (200Mb), we may achieve a single, or few contig assembly with high per-base accuracy.

In this tutorial, we had a reference genome for a closely related organism, which we could use when assessing the quality of our draft genome assemblies. Without access to a reference genome, `Quast` will produce fewer metrics for assessing the assembled genome quality. In this situation, you can look at the number of contigs in the assembled genome, the minimum number of contigs that produce 50\% and 90\% of the bases in the assembly and the GC\% content. For good quality assemblies, the total and minimum number (50\% and 90\%) of contigs should be small and the GC\% content should be near 50\%.

The development of new purpose-built tools for hybrid de novo assembly like `Unicycler` have improved the quality of assemblies we can produce. These tools are of great importance and while they already produce great results, they will continue to improve over time.

-------------------------------
## Additional reading
Links to additional recommended reading and suggestions for related tutorials.

- Flye: [https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#algorithm](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#algorithm)
- Pilon: [https://github.com/broadinstitute/pilon/wiki/Methods-of-Operation](https://github.com/broadinstitute/pilon/wiki/Methods-of-Operation)
- Unicycler: [https://github.com/rrwick/Unicycler](https://github.com/rrwick/Unicycler)
- Quast: [https://academic.oup.com/bioinformatics/article/29/8/1072/228832](https://academic.oup.com/bioinformatics/article/29/8/1072/228832)
- BUSCO analysis: [https://academic.oup.com/bioinformatics/article/31/19/3210/211866](https://academic.oup.com/bioinformatics/article/31/19/3210/211866)

