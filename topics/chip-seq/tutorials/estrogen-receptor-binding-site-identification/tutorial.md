---
layout: tutorial_hands_on
topic_name: chip-seq
tutorial_name: estrogen-receptor-binding-site-identification
---

# Introduction
{:.no_toc}

This exercise uses the dataset from the Nature publication by [Ross-Inness et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/22217937). The goal of this article was to identify the binding sites of the Estrogen receptor, a transcription factor known to be associated with different types of breast cancer.

To this end, ChIP-seq was performed in breast cancer cells from 4 patients of different outcomes (good and poor). For each ChIP-seq experiment there is a matching technical control, *i.e.* there are 8 samples in total:

Patient | Outcome | Treatment
--- | --- | ---
Patient 1 | Good | ChIP ER
Patient 1 | Good | input (no immunoprecipitation step)
Patient 2 | Good | ChIP ER
Patient 2 | Good | input (no immunoprecipitation step)
Patient 3 | Poor | ChIP ER
Patient 3 | Poor | input (no immunoprecipitation step)
Patient 4 | Poor | ChIP ER
Patient 4 | Poor | input (no immunoprecipitation step)

Half of which are the so-called 'input' samples for which the same treatment as the ChIP-seq samples was done except for the immunoprecipitation step. The input files are used to identify sequencing bias like open chromatin or GC bias.

Because of the long processing time for the large original files, we have downsample the original data for practice and provide already processed data for subsequent steps.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Quality control and treatment of the sequences

The first step of any ChIP-Seq data analysis is control of the data quality.

> ### {% icon hands_on %} Hands-on: Quality control
>
> 1. Create and name a new history for this tutorial
> 2. Import `patient1_input_good_outcome` from [Zenodo]() or from the data library into the history
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start** and **Close** the window
>    > * Click on the `pencil` icon once the file is imported
>    > * Click on **Datatype** in the middle panel
>    > * Select `fastqsanger` as **New Type**
>    {: .tip}
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > * Go into "Shared data" (top panel) then "Data libraries"
>    > * Click on "Training data" and then "Analyses of ChIP-Seq data"
>    > * Select interesting file
>    > * Click on "Import selected datasets into history"
>    > * Import in the history
>    {: .tip}
>
>    As default, Galaxy takes the link as name, so rename them. 
>
> 3. Inspect the file by clicking on the `eye` icon
> 
>    > ### {% icon question %} Questions
>    >
>    > 1. How are the DNA sequences stored?
>    > 2. What are the other entries?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The DNA sequences are stored in the second line of every 4-line groups</li>
>    >    <li>This file is a FastQ file. This type of file stores sequence information. Each sequence is represented by a group of 4 lines with the 1st line being the sequence id, the second the sequence of nucleotides, the third a transition line and the last one a sequence of quality score for each nucleotide</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 4. **FastQC** {% icon tool %} with
>    - "Short read data from your current history" to the imported file
>
>    Inspect the generated files
{: .hands_on}

It is often necessary to trim sequenced read, for example, to get rid of bases that were sequenced with high uncertainty (= 'low quality bases').

> ### {% icon hands_on %} Hands-on: Quality control
>
> 1. **Trim Galore!** {% icon tool %} with
>    - "Is this library paired- or single-end?" to `Single-end`
>    - "Reads in FASTQ format" to the imported file
>    - "Trim Galore! advanced settings" to `Full parameter list`
>    - "Trim low-quality ends from reads" to `15`
>    - "Overlap with adapter sequence required to trim a sequence" to `3`
>
>    > ### {% icon tip %} Tip: Importing data from a data library
>    >
>    > If your FASTQ files cannot be selected, you might check whether their format is FASTQ with Sanger-scaled quality values (`fastqsanger`). You can edit the data type by clicking on the `pencil` symbol.
>    {: .tip}
>
{: .hands_on}

# Mapping of the reads

In order to figure where the sequenced DNA fragments originated from in the genome, the short reads must be aligned to the reference genome. This is equivalent to solving a jigsaw puzzles, but unfortunately, not all pieces are unique. In principle, you could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way may, however, take a couple of weeks.

Nowadays, there are many read alignment programs for shotgun sequenced DNA, Bowtie2 being one of them.

> ### {% icon hands_on %} Hands-on: Mapping
>
> 1. **Bowtie2** {% icon tool %} with
>    - "Is this single or paired library" to `Single-end`
>    - "FASTA/Q file" to the Trim Galore! output with the trimmed reads
>    - "Will you select a reference genome from your history or use a built-in index?" to `Use a built-in genome index`
>    - "Select reference genome" to `Human (Homo sapiens): hg18`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. How many reads where mapped?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>This information can be accessed by clicking on the resulting history entry. You can see some basic mapping statistics once the alignment is completed. 16676 (66.96%) were mapped exactly 1 time and 7919 (31.80%) aligned >1 times. The overall alignment rate is then 98.76%.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

The read alignment step with bowtie2 resulted in a compressed, binary file (BAM) that is not human-readable. It's like the zipped version of a text file. 

We will show you two ways to inspect the file:

1. Visualization using a Genome Browser
2. Converting the binary format into its text file equivalent


## Visualization using a Genome Browser

> ### {% icon hands_on %} Hands-on: Visualization of the reads in IGV
>
> 1. Click on the `display with IGV local` to load the reads into the IGV browser
> 2. Zoom at the start of chromosome 11 (or `chr11:1,562,200-1,591,483`)
{: .hands_on}

The reads have a direction: they are mapped to the forward or reverse strand, respectively. When hovering over a read, extra information is displayed

> ### {% icon question %} Questions
>
> 1. Some reads have lines over them. Why?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>Try to zoom in in one of those lines to identify the reason for these</li>
>    </ol>
>    </details>
{: .question}

> ### {% icon comment %} Comments
> Because the number of reads over a region can be quite large, the IGV browser by default only allows to see the reads that fall into a small window. This behaviour can, in principle, be changed in the preferences panel.
{: .comment}


## Inspection of the SAM format

As mentioned above, you can convert the binary BAM file into a simple (but large!) text file, which is called a SAM (Sequence Alignment Map) file.

> ### {% icon hands_on %} Hands-on: Conversion into a SAM file
>
> 1. **BAM-to-SAM** {% icon tool %} with
>    - "BAM File to Convert" to the file generated by Bowtie2
>    - "Header options" to `Include header in SAM output`
>
> 2. Inspect the file by clicking on `eye` icon
>
{: .hands_on}

A SAM file is a file with

- A header with the chromosome names and lengths
- A file content as a tabular file with the location and other information of each read found in the FASTQ file and the mapping information

**Add a question box??**

# Control of the quality of the ChIP-seq

We checked the quality of the sequencing in the first step. We know would like to test the quality of the ChIP-seq preparation of the samples.

## Correlation between samples

To assess the similarity between the replicates of the ChIP-seq and the input, respectively, it is a common technique to calculate the correlation of read counts on different regions for the different samples. We expect that the replicates of the ChIP-seq experiments should be clustered more closely to each other than the replicates of the input sample.

**Not really clear for me... Add a bit more details? why such step? how such step is done?**

To do that, we need to run the previous steps (quality control and mapping) on each sample. To save time, we already did that and we can now work directly on the BAM files of the 8 samples

> ### {% icon hands_on %} Hands-on: Correlation between samples
>
> 1. Create a new history
> 2. Import the 8 BAM files from [Zenodo]() or from the data library into the history
> 3. **multiBamSummary** {% icon tool %} with
>    - "Sample order matters" to `No`
>    - "Bam file" to the 8 imported BAM files
>    - "Choose computation mode" to `Bins`
>    - "Bin size in bp" to `100`
>       
>       This corresponds to the length of the fragments that were sequenced; it is not the read length!
>
>    - "Distance between bins" to `500000` (to reduce the computation time for the tutorial)
>    - "Region of the genome to limit the operation to" to `chr1` (to reduce the computation time for the tutorial)
>  
>    This tool splits the reference genome into bins of equal size (*e.g.* 10kb) and counts the number of overlapping reads from each sample.
>
> 4. **plotCorrelation** {% icon tool %} with
>    - "Matrix file from the multiBamSummary tool" to the generated multiBamSummary output
>    
>    Feel free to try different parameters
{: .hands_on}

> ### {% icon question %} Questions
> 
> ![Output for plotCorrelation with the correlation scores between the 8 samples](../../images/estrogen-receptor-binding-site-identification/plotCorrelation_output.png "Correlation scores between the 8 samples")
>
> 1. How are clustered the samples? Does that correspond to the expectations?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li></li>
>    </ol>
>    </details>
{: .question}

> ### {% icon comment %} Comments
> More information on these two tools can be found at the [deepTools documentation page](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).
{: .comment}

## GC bias assessment

A common problem of PCR-based protocols is the observation that GC-rich regions tend to be amplified more readily than GC-poor regions. We need to check that our samples do not have more reads from regions of the genome with high GC. 

For practical reasons, we will focus here only on one of the BAM files, an input file.

> ### {% icon question %} Questions
> 
> 1. Can you guess why it makes more sense to check the input file?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>Only the bias induced by the PCR-based protocols, nothing added with the immunoprecipation step (because none in the input file) ????</li>
>    </ol>
>    </details>
{: .question}

> ### {% icon hands_on %} Hands-on: GC bias assessment
>
> 1. **computeGCbias** {% icon tool %} with
>    - "BAM file" to `patient1_input_good_outcome`
>    - "Reference genome" to `locally cached`
>    - "Using reference genome" to `Human (Homo sapiens): hg18`
>    - "Effective genome size" to `hg19 (2451960000)`
>    - "Fragment length used for the sequencing" to `300`
>    - "Region of the genome to limit the operation to" to `chr1` (to reduce the computation time for the tutorial)
>
>    > ### {% icon question %} Questions
>    >
>    > ![Output for computeGCbias with the GC bias estimation](../../images/estrogen-receptor-binding-site-identification/computeGCbias_output.png "Estimation of the GC bias for the input sample for the Patient 1")
>    >
>    > 1. Does this dataset have a GC bias?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>There is no significantly more reads in the GC-rich regions.</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}


## IP strength estimation

To evaluate the quality of the immuno-precipitation step, we can compute the IP strength. It determines how well the signal in the ChIP-seq sample can be differentiated from the background distribution of reads in the control sample. To do that we take the data for one Patient and compare the input sample and the ChIP-seq sample.

> ### {% icon comment %} Comments
> For more information on how the IP strength is estimated, you can check the [deepTools documentation page](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).
{: .comment}

> ### {% icon hands_on %} Hands-on: IP strength estimation
>
> 1. **plotFingerprint** {% icon tool %} with
>    - "Sample order matters" to `No`
>    - "Bam file" to `patient1_input_good_outcome` and `patient1_ChIP_ER_good_outcome`
>    - "Region of the genome to limit the operation to" to `chr1`
{: .hands_on}

The plotFingerprint tool generates a fingerprint plot. You need to intepret it to know the IP strength. The [deepTools documentation](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html) explains it clearly:

![A guide to interpret a fingerprint plot](../../images/estrogen-receptor-binding-site-identification/QC_fingerprint.png "How to interpret a fingerprint plot? Image extracted from the deepTools documentation")

> ### {% icon question %} Questions
>
> ![Output for plotFingerprint with the fingerprint plot to estimate the IP strength](../../images/estrogen-receptor-binding-site-identification/plotFingerprint_output.png "Fingerprint plot for the Patient 1 to estimate the IP strength")
>
> 1. What do you think about the quality of the IP for this experiment?
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>The difference between input and ChIP signal is not totally clear. 20% of the entire chromosome do not have any read (curve rise start)</li>
>    </ol>
>    </details>
{: .question}

> ### {% icon hands_on %} (Optional) Hands-on: GC bias assessment (other samples)
>
> 1. Run the same analysis on data of the 3 other patients
{: .hands_on}


# Information extraction (rename!!!)

We would like to know where are the binding site of the estrogen. We need then to extract which parts of the genome have been enriched (more reads mapped on them) within the samples that underwent immunoprecipitation. 

To extract such information, there is 3 possible processes:

1. Normalization by sequencing depth
2. Normalization by the coverage file
3. Enriched region calling

Here we will only do this extraction only on the data for the Patient 1. But the same steps can be done for all the patients.

## Generation of coverage files normalized by sequencing depth

We first need to make the samples comparable. Indeed, the different samples have usually a different sequencing depth, *i.e.* a different number of reads. These differences can bias the interpretation of the number of reads mapped on a genome portion.

> ### {% icon hands_on %} Hands-on: Coverage file normalization
>
> 1. **IdxStats** {% icon tool %} with
>    - "BAM file" to "Multiple datasets": `patient1_input_good_outcome` and `patient1_ChIP_ER_good_outcome`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What is the output of this tool?
>    > 2. How many reads has been mapped on chr2 for the input and for the ChIP-seq samples?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>This tool estimates how many reads mapped to which chromosome. Furthermore, it tells the chromosome lengths and naming convention (with or without 'chr' in the beginning)</li>
>    >    <li>1,089,370 for ChIP-seq samples and 1,467,480 for the input</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 2. **bamCoverage** {% icon tool %} with
>    - "BAM file" to "Multiple datasets": `patient1_input_good_outcome` and `patient1_ChIP_ER_good_outcome`
>    - "Bin size in bases" to `25`
>    - "Scaling/Normalization method" to `Normalize coverage to 1x`
>    - "Effective genome size" to `hg19 (2451960000)`
>    - "Coverage file format" to `bedgraph`
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What are the different columns 
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>...</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 3. **bamCoverage** {% icon tool %} with the same parameters but to generate a `bigWig` output file
> 4. **IGV** {% icon tool %} to inspect both signal coverages (input and ChIP samples) in IGV
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Add a question
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>And the answer</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
{: .hands_on}

## Generation of input-normalized coverage files and their visualization

To extract only the information induced by the immunoprecipitation, we normalize for each patient the coverage file for the sample that underwent immunoprecipitation by the coverage file for the input sample. Here we use the tool bamCompare which compare 2 BAM files while caring for sequencing depth normalization.

> ### {% icon hands_on %} Hands-on: Generation of input-normalized coverage files
>
> 1.  **bamCompare** {% icon tool %} with
>    - "First BAM file (e.g. treated sample)" to `patient1_ChIP_ER_good_outcome`
>    - "Second BAM file (e.g. control sample)" to `patient1_input_good_outcome`
>    - "Bin size in bases" to `50`
>    - "How to compare the two files" to `Compute log2 of the number of reads ratio`
>    - "Coverage file format" to `bedgraph`
>    - "Region of the genome to limit the operation to" to `chr11` (to reduce the computation time for the tutorial)
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What does mean a positive or a negative value in the 4th column?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>The 4th column contains the log2 of the number of reads ratio between the ChIP-seq sample and the input sample. A positive value means that the coverage on the portion is more important in the ChIP-seq sample than in the input sample</li>
>    >    </ol>
>    >    </details>
>    {: .question}
>
> 2. **bamCompare** {% icon tool %} with the same parameters but to generate a `bigWig` output file
> 3. **IGV** {% icon tool %} to inspect the log2 ratio
>
>    Remember that the bigWig file contains only the signal on chromosome 11!
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Add a question
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>And the answer</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

We would like now to visualize our coverage as an heatmap.

> ### {% icon hands_on %} Hands-on: Visualization of the coverage
>
> 1. **UCSC Main** {% icon tool %} with
>    - "assembly" to `hg18`
>    - "track" to `RefSeq genes`
>    - "region" to `position` with `chr11`
>    - "output format" to `BED`
>    - "Send output to" to `Galaxy`
>
> 2. **computeMatrix** {% icon tool %} with
>    - "Regions to plot" to the imported UCSC file
>    - "Score file" to the bigwig file generated by bamCompare
>    - "computeMatrix has two main output options" to `scale-regions`
>
>       This option stretches or shrinks all regions in the BED file (here: genes) to the same length (bp) as indicated by the user
>
>    - "Show advanced options" to `yes`
>    - "Convert missing values to 0?" to `Yes`
>
>    This tool prepares a file with scores per genomic region, which is required as input for the next tool.
>
> 3. **plotHeatmap** {% icon tool %} with
>    - "Matrix file from the computeMatrix tool" to the generated matrix
>    - "Show advanced options" to `yes`
>    - "Did you compute the matrix with more than one groups of regions?" to the correct setting
>
> 4. **computeMatrix** {% icon tool %} with the same parameters but `reference-point` as output option
>
>    With this option, it considers only those genomic positions before (downstream) and/or after (upstream) a reference point (*e.g.* TSS, which corresponds to the annotated gene start in our case)
{: .hands_on}

> ### {% icon question %} Questions
> 
> 1. Add a question to help heatmap interpretation
> 1. Add a question to relate the heatmap to the biological question...
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>Add the answer</li>
>    </ol>
>    </details>
{: .question}

## Enriched regions (peaks) calling

We can also call the enriched regions, or peaks, found in the ChIP-seq samples. 

**add more details on how it is done**

> ### {% icon hands_on %} Hands-on: Peak calling
>
> 1. **MACS2 callpeak** {% icon tool %} with
>    - "ChIP-Seq Treatment File" to `patient1_ChIP_ER_good_outcome`
>    - "ChIP-Seq Control File" to `patient1_input_good_outcome`
>    - "Effective genome size" to `H. sapiens (2,451,960,000)`
>    - "Outputs" to `Summary page (html)`
>
>    > ### {% icon comment %} Comments
>    > The advanced options may be adjusted, depending of the samples.
>    > If your ChIP-seq experiment targets regions of broad enrichment, *e.g.* non-punctuate histone modifications, select calling of broad regions.
>    > If your sample has a low duplication rate (*e.g.* below 10%), you might keep all duplicate reads (tags). Otherwise, you might use the 'auto' option to estimate the maximal allowed number of duplicated reads per genomic location.
>    {: .comment}
>
> 2. **IGV** {% icon tool %} to inspect with the signal coverage and log2 ratio tracks
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Do the called peaks match your expectation based on the signal coverage and log2 ratio tracks?
>    >
>    >    <details>
>    >    <summary>Click to view answers</summary>
>    >    <ol type="1">
>    >    <li>And the answer</li>
>    >    </ol>
>    >    </details>
>    {: .question}
{: .hands_on}

The called peak regions can be filtered by, *e.g.* fold change, FDR and region length for further downstream analysis.

> ### {% icon question %} Questions
> 
> 1. Add a question to relate to the biological question...
>
>    <details>
>    <summary>Click to view answers</summary>
>    <ol type="1">
>    <li>Add the answer</li>
>    </ol>
>    </details>
{: .question}

# Conclusion
{:.no_toc}

![Summary of the different steps of the tutorial and the generated files](../../images/estrogen-receptor-binding-site-identification/tutorial-scheme.png "Different steps of the tutorials with the generated files")

