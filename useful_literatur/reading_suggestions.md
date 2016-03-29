#Sequencing biases

**Meacham et al. (2011):** [Identification and correction of systematic error in high-throughput sequence data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451), (doi:10.1186/1471-2105-12-451) - Correcting for systematic base pair errors in deep sequencing; important paper if you want to do look at any allele-specifically or if you're interested in SNPs

**Benjamini & Speed (2012):** GC bias of deeply sequenced samples - very good paper that systematically assesses many possible sourced of GC bias for deeply sequenced samples and eventually pinpoints it to the DNA Polymerase

A collection of papers on quality controls for various NGS applications: Frontiers in Genetics (2014)
---------------------------------------------------------------------------------------
Deep sequencing

Overview of popular *-seq techniques (2012) - by Zentner and Henikoff, very nice description of DNase-seq, MNase-seq, FAIRE-seq etc.

Paper on multiplexing - describes the individual steps of the Illumina deep sequencing protocols quite detailed

Illumina's tech report- focuses on Illumina's sequencing technique; nice educative figures
---------------------------------------------------------------------------------------
Mapping of short NGS reads

Intro to various aspects of NGS read mapping: informative slides with some nice reality checks

Fonseca et al. (2012): Tools for mapping high-throughput sequencing data - an excellent starting point despite its "old" age, you will learn a lot about the different philosophies behind the read alignment tools!

Hatem et al. (2013): Benchmarking short sequence mapping tools (spoiler alert: bowtie wins)

Engstrom et al. (2013): Systematic evaluation of spliced alignment programs for RNA-seq data
Genome Mappability

Lee and Schatz (2012): The reliability of short read mapping - very detailed paper about genome mappability issues that presents a new suite of tools for taking the mappability into account

mappability maps can be downloaded here
------------------------------------------------------------------------
ChIP-seq in general

Landt et al. (2012): ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia   - this is a very useful "encyclopedic" paper with many details about the tools the (mod)ENCODE consortia use. It also contains a long section about antibody validation etc., it does not explain much of the reasoning behind the bioinformatics tools, though.

Zentner et al. (2012): Surveying the epigenomic landscape, one base at a time - short review about MNase-Seq, DNase-Seq, ChIP-exo

Kidder et al. (2011): Technical considerations to obtaining high-quality data - nice, readable introduction into all aspects of ChIP-seq experiments (from antibodies to cell numbers to replicates to data analysis)

Leleu et al. (2010): Processing and analyzing ChIP-seq data - fairly detailed review of key concepts of ChIP-seq data processing (less detailed on analysis)

Peter Park (2009): Advantages and challenges of a maturing technology 

Kharchenko et al. (2008): Design and analysis of ChIP-seq experiments for DNA-binding proteins

Liu et al. (2010): Q&A: ChIP-seq technologies and the study of gene regulation - short overview of several (typical) issues of ChIP-seq analysis

Carroll et al. (2014):  Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data
---------------------------------------------------------------------------------------
RNA-seq

Nice graphical overview of the RNA-seq processing and analysis steps: http://figshare.com/articles/RNA_seq_Workflows_and_Tools/662782

Insights into properly planning your RNA-seq run! Read BEFORE your RNA-seq experiment! Auer and Doerge (2010)

A refreshingly honest view on the non-trivial aspects of RNA-seq analysis: Ian Korf, Nat Methods (2013)

Dillies et al. (2012): Systematic comparison of seven representative normalization methods for the differential analysis of RNA-seq data (Total Count, Upper Quartile, Median (Med), DESeq, edgeR, Quantile and Reads Per Kilobase per Million mapped reads (RPKM) normalization)

Rapaport et al. (2013): Evaluation of methods for differential gene expression analysis

Soneson et al. (2013): And another paper about different methods for differential gene expression analysis

Roberts et al. (2011): Fragment bias correction

Garber et al. (2011): Classical paper about the computational aspects of RNA-seq data analysis
---------------------------------------------------------------------------------------
Peak Calling Methods (ChIP-seq)

Pepke et al. (2009): First comparison of peak callers - focuses on the explanation of basic principles of ChIP-seq data processing and the general workflows of peak calling algorithms

Wilbanks et al. (2010): Another comparison of peak callers - focuses more on the evaluation of the peak callers performances than Pepke et al.

Micsinai et al. (2012): How to choose the best peak caller for your data set - their finding: default parameters, surprisingly, yield the most reproducible results regardless of the data set type
MACS

Fen et al. (2012): How to use MACS - Nature Protocols

Zhang et al. (2008): The original publication of MACS: Model-bases analysis of ChIP-seq data.
---------------------------------------------------------------------------------------
DNA motif analysis

Das et al. (2007): Review of Motif Analysis Tools
MEME (suite)

Machanick and Bailey (2011): MEME-ChIP - paper

Bailey and Machanick (2012): Centrimo - position-specific motif analysis, especially useful for ChIP-seq data

TomTom - tool for the comparison of motifs from databases (not in Galaxy yet): tutorial
TRAP

Thomas-Chollier et al. (2012): How to use TRAP - Nature Protocols

Roider et al. (2006): Theoretical background of TRAP
--------------------------------------------------------------------------------------
NGS data formats

UCSC has a very good overview with brief descriptions of BED, bedGraph, bigWig etc.: https://genome.ucsc.edu/FAQ/FAQformat.html

VCF format (encoding SNPs, indels etc.): Very readable, albeit not exhausting description

transcriptomes are often saved in GFF3 format (this is what TopHat needs, for example), but just to make things more complicated, GTF is another format used for transcriptome information, too (here are more information on GTF)
--------------------------------------------------------------------------------------
Bioinformatic Tools (Linux, R, BEDTools etc.) - Manuals, Courses, original papers

Why and how is bioinformatics software special? (Highly recommended read!)

A Field Guide to Genomics Research- very readable introduction about the different caveats of genomics research (with cute cartoons!)
Linux Command Line

Linux & Perl Primer for Biologists - very entertaining introduction to command line commands and perl scripts with a focus on bioinformatic application, i.e. handling of DNA sequences

Linux Tutorial for Beginners -  thorough, but concise online tutorial introducing the very basics of handling the Linux command line

Writing Linux shell scripts - useful for slightly more advanced Linux command line users
R

Hands on R course - for beginners - R is probably the most widely used open-source statistical software; through our epicenter website you can also access RStudio which provides are very nice interface to working and plotting with R. In fact, most of the plots generated within Galaxy are generated through R scripts, so if you're not happy with the default formats of the Galaxy graphs, definitely have a look at R yourself. The learning curve is steep, but it is worth it.
BEDTools

BEDTools Manual - when working with genomic intervals (e.g. genes, peaks, enriched regions...), BEDTools are invaluable! The manual is a very good read and we refer to it almost daily.
