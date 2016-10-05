#ChIP-seq Galaxy Workshop

## Content

* [Hands-on example](#example)
* [Useful literature](#literature)
  - [ChIP-seq in general](#chipseq)
  - [Peak calling](#peakcalling)
  - [Motif analysis](#motifs)

## Slides from the hands-on workshops at the University of Freiburg

The slides for part 1 of this session can be downloaded from here
[ChIP-seq1-galaxy_course_2015.pdf](https://drive.google.com/file/d/0B9urRnOAUUI8UmwzbTVpdmZucWM/view?usp=sharing)


The slides for part 2 can be downloaded from here
[ChIP-seq2-galaxy_course_2015.pdf](https://drive.google.com/file/d/0B9urRnOAUUI8cHpzYVBscjNKWEE/view?usp=sharing).  

<a name="example"/></a>
## Hands on example  

This exercise uses the dataset from the Nature publication by [Ross-Inness et al., 2012](http://www.ncbi.nlm.nih.gov/pubmed/22217937).
The goal was to identify the binding sites of the Estrogen receptor, a transcription factor known to be associated with different types of breast cancer.
To this end, ChIP-seq was performed in breast cancer cells from 4 patients of different outcomes (good and poor).
For each ChIP-seq experiment there is a matching technical control, i.e., there are 8 samples in total, half of which are the so-called 'input' samples for which the same treatment as the ChIP-seq samples was done except for the immunoprecipitation step.
The input files are used to identify sequencing bias like open chromatin or GC bias.

Because of the long processing time for the large original files, we have selected small samples for practice and provide already processed data for subsequent steps.

### Step 1: Quality control

Create a new history for this exercise.

- Import the [FASTQ file patient1_input_good_outcome_chr11_25k_reads.fastq](https://github.com/bgruening/training_data/raw/master/ChIPseq/Galaxy1-%5Bpatient1_input_good_outcome_chr11_25k_reads.fastq%5D.fastqsanger) into the history.

- Have a look at the file by clicking on the 'eye' icon. There is a lot of text, but can you spot where the DNA sequence is stored? Can you guess what the other entries mean?

- Run the tool `FastQC` on one of the two FASTQ files to control the quality of the reads. An explanation of the results can be found on the [FastQC web page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

**Step 2: Trimming and clipping of reads**

It is often necessary to trim sequenced read, for example, to get rid of bases that were sequenced with high uncertainty (= 'low quality bases').

- Apply the tool `TrimGalore` to your FASTQ file. Have a look at the full parameter list and then tell `TrimGalore` to:
    - __trim low quality bases__ with a cut-off of 15
    - __clip a possibly trailing adapter sequence__ (The default sequence that appears is the generic part of Illumina adaptors and should be used.) Instruct `TrimGalore` to trim adapter sequences with an overlap of at least 3.

-> If your FASTQ files cannot be selected, you might check whether their format is FASTQ with Sanger-scaled quality values (_fastqsanger_). You can edit the data type by clicking on the 'pencil' symbol.


**Step 3: Mapping of the reads**

In order to figure where the sequenced DNA fragments originated from in the genome, the short reads must be aligned to the reference genome. This is equivalent to solving a jigsaw puzzles, but unfortunately, not all pieces are unique. In principle, you could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way may, however, take a couple of weeks.
Nowadays, there are many read alignment programs for shot-gun sequenced DNA, `bowtie2` being one of them.

- Run `bowtie2` to map the single-end reads to the human genome (version GRCh38).

- By clicking on the resulting history entry, you can see some basic mapping statistics once the alignment is completed. How many reads where mapped?


**Step 4: Visualization of the reads in IGV**

The read alignment step with bowtie2 resulted in a compressed, binary file (BAM) that is not human-readable (it's like the zipped version of a text file). We will show you two ways to inspect the file, (1) by visualization using a Genome Browser and (2) by converting the binary format into its text file equivalent.

- Load the reads into the IGV browser by clicking the option 'display with IGV local'. To see the reads in IGV you will need to zoom in the start of chromosome 11. Try this region if you don't see any reads: chr11:1,562,200-1,591,483

- Notice that the reads have a _direction_ (i.e., they are mapped to the forward or reverse strand, respectively). When hovering over a read, extra information is displayed. Some reads have lines over them. Try to zoom in in one of those lines to identify the reason for these.

**Note**: because the number of reads over a region can be quite large, the IGV browser by default only allows to see the reads that fall into a small window. This behaviour can, in principle, be changed in the preferences panel.


**Step 5: Inspection of the SAM format**

As mentioned above, you can convert the binary BAM file into a simple (but large!) text file, which is called a SAM (Sequence Alignment Map) file.

- Go back to Galaxy and run the tool `BAM-to-SAM` using the BAM file that was created in step 2. Click 'include header in output'.

- Click on the view icon ('eye'). The first part of the file, the header, contains the chromosome names and lengths. After the header, the location and other information of each read found in the FASTQ file is given.


**Step 6: Correlation between samples**

To assess the similarity between the replicates of the ChIP-seq and the input, respectively, it is a common technique to calculate the correlation of read counts for the different samples.

We expect that the replicates of the ChIP-seq experiments should be clustered more closely to each other than the replicates of the input samples.

- For this step, we need to load all files with the aligned reads for each sample (.bam files) into the shared libraries. For this, please import into the current history all BAM files found in 'Galaxy courses' -> 'Hands-on'.

*bam files examples* -> link

- Next, run the tool `multiBamSummary` from the deepTools package. This tool will split the reference genome into bins of equal size (e.g. 10kb) and will count the number of overlapping reads from each sample.
    - Select only the 8 BAM files imported into the history.
    - Use as fragment length: 100 (this corresponds to the length of the fragments that were sequenced; it is not the _read_ length!)
    - To reduce the computation time, set the distance between bins to 500,000 and choose only one chromosome, for example 'chr1'.

- After the computation is done, run `plotCorrelation` from the deepTools package to visualize the results. Feel free to try different parameters.

More information on these two tools can be found at the [deepTools documentation page](http://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).

**Step 7: GC bias assessment**

A common problem of PCR-based protocols is the observation that GC-rich regions tend to be amplified more readily than GC-poor regions.

- Use the tool `computeGCbias` to check if the samples have more reads from regions of the genome with high GC.
    - for practical reasons select only one of the BAM files available, preferably an input file (can you guess why it makes more sense to check the input file?)
    - For _fragment size_, select 300.
    - limit the operation to only one chromosome (again, this is purely to speed up the analysis)

Does this dataset have a GC bias?

**Step 8: IP strength**

- To evaluate the quality of the immuno-precipitation step use the tool `plotFingerprint`.
    - Select one of the ChIP-seq samples and the matching input
    - Set as fragment size 100.
    - Limit the operation to only one chromosome.

What do you think about the quality of the IP for this experiment? If you are not sure how to interpret the resulting plots, please read the information [here](http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html#background)

**Step 9: Generate coverage files normalized by sequencing depth**

- If you have not done already, please import to the history the BAM files patient4_ChIP_ER_poor_outcome.bam and patient4_input_poor_outcome.bam. These files are provided in the data library under Galaxy courses -> ChIP-seq -> bam_files.

- The first file contains the reads from a ChIP-seq experiment for the transcription factor ER using human cells from a breast cancer patient. The second file contains the reads from the matching input sample. The sequencing reads were mapped to the human genome assembly hg18.

- Use the SAM Tools `IdxStats` tool to find out how many reads mapped to which chromosome. Furthermore, this will also tell you the chromosome lengths and naming convention (with or without 'chr' in the beginning).

- Run the tool `bamCoverage` to generate a signal coverage file for the ER ChIP sample normalized by sequencing depth. Set the fragment size to 100 and the bin size to 25. Normalize to 1x genomic coverage. The output file should be in human-readable format bedGraph. To speed up computation, limit the operation to chromosome 'chr11'.

Generally, you should adjust the effective genome size according to the used genome assembly. In our case, you however have to specify the size of chromosome chr11 only when limiting the computation to this region. We obtained this value just in the last step :-)

- Inspect the bedGraph output file.

- Re-run the tool and generate a *bigWig* output file. Inspect the signal coverage in IGV. Remember that the bigWig file contains only the signal on chromosome 11!


**Step 10: Generate input-normalized coverage files**

- Run the tool `bamCompare` to normalize the ChIP signal BAM file patient4_ChIP_ER_poor_outcome.bam by the input control provided by patient4_input_poor_outcome.bam.

- Set the fragment size to 100 again and the bin size to 50. Compute the log2 ratio of the read counts of ER ChIP vs. input sample. The output file should be in human-readable format bedGraph. To speed up computation, limit the operation to chromosome 'chr11'.

- Inspect the bedGraph output file.

- Re-run the tool and generate a bigWig output file. Inspect the log2 ratio in IGV. Remember that the bigWig file contains only the signal on chromosome 11!


**Step 11: Call enriched regions (peaks)**

To speed up computation, we will now restrict our analysis to chromosome 11.

- To this end, please import to the history the BAM files patient4_ChIP_ER_poor_outcome_chr11.bam and patient4_input_poor_outcome_chr11.bam. These files are provided in the data library under Galaxy courses -> ChIP-seq -> bam_files.

- Use the `MACS2 callpeak` tool to call peak regions from the alignment results in these two BAM files. Adjust the effective genome size to *Homo sapiens* genome assembly hg18, chromosome 11. Output the peaks as BED file and the HTML summary page.

You might adjust the advanced options depending on your sample:

-> If your ChIP-seq experiment targets regions of broad enrichment, e.g., non-punctuate histone modifications, select calling of broad regions.
-> If your sample has a low duplication rate (e.g. below 10%), you might keep all duplicate reads (tags). Otherwise, you might use the 'auto' option to estimate the maximal allowed number of duplicated reads per genomic location.

- Inspect the result files. When visualized in IGV, do the called peaks match your expectation based on the signal coverage and log2 ratio tracks?

The called peak regions can be filtered by, e.g., fold change, FDR and region length for further downstream analysis.


**Step 12: Visualize your ChIP-seq results as heatmap**

For the last exercise, we will use the data resulting from step 10.

- If you have not done so, please generate a bigWig file with the log2 ratio of ER ChIP over input signal restricted to chromosome 'chr11' as described in step 10.

- Furthermore, please download the following human gene annotation from **UCSC Table Browser**:
    - Human genome assembly hg18
    - RefSeq Genes track, limited to region chr11. Download the data a BED file and send it to Galaxy. This will import the regions of all genes on chromosome 11 that are annotated in RefSeq to your Galaxy history.

- Now, run the `computeMatrix` tool to prepare a file with scores per genomic region, which is required as input for the next tool.

Use the regions provided by the gene annotation file downloaded from UCSC and your log2 ratio score bigWig file obtained by `bamCompare` in step 10 as input. If you have trouble finding the correct input data: better rename datasets with self-explanatory names the next time ;-)

- You can run `computeMatrix` in two alternative modes:
    - **scale-regions**: stretches or shrinks all regions in the BED file (here: genes) to the same length (bp) as indicated by the user
    - **reference-point**: considers only those genomic positions before (downstream) and/or after (upstream) a reference point (e.g. TSS, which corresponds to the annotated gene start in our case).

- Run the tool in either mode.
  - Ensure that the advanced option 'Convert missing values to 0?' is activated.

- Generate a heatmap with `plotHeatmap` based on your output of `computeMatrix`.
  - Appropriately, set the advanced option 'Did you compute the matrix with more than one groups of regions?'

- How can you increase the "resolution" of the data in the resulting heatmap?

-------------------------------------------

<a name="literature"/></a>
##Useful literature

<a name="chipseq"/></a>
###ChIP-seq in general:

**Landt et al. (2012):** [ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia](http://genome.cshlp.org/content/22/9/1813.long), (doi:10.1101/gr.136184.111) - This is a very useful "encyclopedic" paper with many details about the tools the (mod)ENCODE consortia use. It also contains a long section about antibody validation etc.. It does not explain much of the reasoning behind the bioinformatics tools, though.

**Zentner and Henikoff (2012):** [Surveying the epigenomic landscape, one base at a time](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-250), (doi:10.1186/gb-2012-13-10-250) - Overview of popular *-seq techniques; very nice description of DNase-seq, MNase-seq, FAIRE-seq etc.

**Kidder et al. (2011):** [Technical considerations to obtaining high-quality data](http://www.nature.com/ni/journal/v12/n10/abs/ni.2117.html), (doi:10.1038/ni.2117) - Nice, readable introduction into all aspects of ChIP-seq experiments (from antibodies to cell numbers to replicates to data analysis)

**Leleu et al. (2010):** [Processing and analyzing ChIP-seq data](http://www.ncbi.nlm.nih.gov/pubmed/20861161), (doi: 10.1093/bfgp/elq022) - Fairly detailed review of key concepts of ChIP-seq data processing (less detailed on analysis)

**Peter Park (2009):** [ChIP-seq: Advantages and challenges of a maturing technology](http://www.nature.com/nrg/journal/v10/n10/full/nrg2641.html), (doi:10.1038/nrg2641)

**Kharchenko et al. (2008):** [Design and analysis of ChIP-seq experiments for DNA-binding proteins](http://www.ncbi.nlm.nih.gov/pubmed/19029915), (doi:10.1038/nbt.1508)

**Liu et al. (2010):** [Q&A: ChIP-seq technologies and the study of gene regulation](http://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-56), (doi:10.1186/1741-7007-8-56) - Short overview of several (typical) issues of ChIP-seq analysis

**Carroll et al. (2014):**  [Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data](http://journal.frontiersin.org/article/10.3389/fgene.2014.00075/full),(doi:10.3389/fgene.2014.00075)  

<a name="peakcalling"/></a>
###Peak Calling Methods (ChIP-seq)

**Pepke et al. (2009):** [Computation for ChIP-seq and RNA-seq studies](http://www.ncbi.nlm.nih.gov/pubmed/19844228), (doi: 10.1038/nmeth.1371) - First comparison of peak callers, focuses on the explanation of basic principles of ChIP-seq data processing and general workflows of peak calling algorithms

**Wilbanks et al. (2010):** [Evaluation of Algorithm Performance in ChIP-Seq Peak Detection](http://www.ncbi.nlm.nih.gov/pubmed/20628599), (doi: 10.1371/journal.pone.0011471) - Another comparison of peak callers - focuses more on the evaluation of the peak callers performances than Pepke et al. (2009)

**Micsinai et al. (2012):** [Picking ChIP-seq peak detectors for analyzing chromatin modification experiments](http://www.ncbi.nlm.nih.gov/pubmed/22307239), (doi: 10.1093/nar/gks048) - How to choose the best peak caller for your data set - their finding: default parameters, surprisingly, yield the most reproducible results regardless of the data set type

#### MACS

**Fen et al. (2012):** [Identifying ChIP-seq enrichment using MACS.](http://www.ncbi.nlm.nih.gov/pubmed/22936215), (doi:10.1038/nprot.2012.101) - How to use MACS - Nature Protocols

**Zhang et al. (2008):** [Model-based Analysis of ChIP-Seq (MACS)](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137), (doi:10.1186/gb-2008-9-9-r137) - The original publication of MACS

<a name="motifs"/></a>
### DNA motif analysis

**Das et al. (2007):** [A survey of DNA motif finding algorithms](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-S7-S21), (doi:10.1186/1471-2105-8-S7-S21) - Review of motif analysis tools

#### MEME (suite)

**Machanick and Bailey (2011):** [MEME-ChIP: motif analysis of large DNA datasets](http://www.ncbi.nlm.nih.gov/pubmed/21486936), (doi: 10.1093/bioinformatics/btr189) - MEME-ChIP-paper

**Bailey and Machanick (2012):** [Inferring direct DNA binding from ChIP-seq](http://www.ncbi.nlm.nih.gov/pubmed/22610855), (doi:10.1093/nar/gks433) - Centrimo: position-specific motif analysis, especially useful for ChIP-seq data

[TomTom](http://meme-suite.org/tools/tomtom) - Meme Suite Motif comparison tool: tool for the comparison of motifs from databases (not in Galaxy yet): [Manual](http://meme-suite.org/doc/tomtom.html?man_type=web)

#### TRAP

**Thomas-Chollier et al. (2012):** [Transcription factor binding predictions using TRAP for the analysis of ChIP-seq data and regulatory SNPs](http://www.ncbi.nlm.nih.gov/pubmed/22051799), (doi:10.1038/nprot.2011.409) - How to use TRAP - Nature Protocols

**Roider et al. (2006):** [Predicting transcription factor affinities to DNA from a biophysical model.](http://www.ncbi.nlm.nih.gov/pubmed/17098775), (doi:10.1093/bioinformatics/btl565) - Theoretical background of TRAP
