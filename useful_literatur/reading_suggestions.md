#Sequencing biases

**Meacham et al. (2011):** [Identification and correction of systematic error in high-throughput sequence data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451), (doi:10.1186/1471-2105-12-451) - Correcting for systematic base pair errors in deep sequencing; important paper if you want to do look at any allele-specifically or if you're interested in SNPs

**Benjamini & Speed (2012):** [Summarizing and correcting the GC content bias in high-throughput sequencing](http://nar.oxfordjournals.org/content/40/10/e72.long), (doi: 10.1093/nar/gks001) - GC bias of deeply sequenced samples; very good paper that systematically assesses many possible sourced of GC bias for deeply sequenced samples and eventually pinpoints it to the DNA Polymerase

**A collection of papers on quality controls for various NGS applications:** [Frontiers in Genetics (2014)](http://journal.frontiersin.org/researchtopic/1683/quality-assessment-and-control-of-high-throughput-sequencing-data)


#Deep sequencing

**Zentner and Henikoff (2012):** [Surveying the epigenomic landscape, one base at a time](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-250), (doi:10.1186/gb-2012-13-10-250) - Overview of popular *-seq techniques; very nice description of DNase-seq, MNase-seq, FAIRE-seq etc.

**Son and Taylor (2011):** [Preparing DNA Libraries for Multiplexed Paired-End Deep Sequencing for Illumina GA Sequencers](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076644/), (doi:10.1002/9780471729259.mc01e04s20) - Paper on multiplexing; describes the individual steps of the Illumina deep sequencing protocols quite detailed

**Illumina's tech report** - focuses on [Illumina's sequencing technology](http://www.illumina.com/technology.html); nice educative figures


#Mapping of short NGS reads

**informative slides** [Mapping of sequencing reads](http://people.binf.ku.dk/krogh/tmp/Mapping_Krogh_Monday.pdf) - Intro to various aspects of NGS read mapping

**Fonseca et al. (2012):** [Tools for mapping high-throughput sequencing data](http://bioinformatics.oxfordjournals.org/content/28/24/3169.full), (doi:10.1093/bioinformatics/bts605) - an excellent starting point despite its "old" age, you will learn a lot about the different philosophies behind the read alignment tools!

**Hatem et al. (2013):**[Benchmarking short sequence mapping tools](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-184), (doi:10.1186/1471-2105-14-184) (spoiler alert: bowtie wins)

**Engstrom et al. (2013):** [Systematic evaluation of spliced alignment programs for RNA-seq data](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2722.html), (doi:10.1038/nmeth.2722) 

**Genome Mappability**

**Lee and Schatz (2012):** [The reliability of short read mapping](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413383/?report=reader), (doi:10.1093/bioinformatics/bts330)  - very detailed paper about genome mappability issues that presents a new suite of tools for taking the mappability into account

*mappability maps* can be downloaded [here](http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map/)


#ChIP-seq in general

**Landt et al. (2012):** [ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia](http://genome.cshlp.org/content/22/9/1813.long), (doi:10.1101/gr.136184.111) - this is a very useful "encyclopedic" paper with many details about the tools the (mod)ENCODE consortia use. It also contains a long section about antibody validation etc., it does not explain much of the reasoning behind the bioinformatics tools, though.

**Zentner and Henikoff (2012):** [Surveying the epigenomic landscape, one base at a time](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-250), (doi:10.1186/gb-2012-13-10-250) - Overview of popular *-seq techniques; very nice description of DNase-seq, MNase-seq, FAIRE-seq etc.

**Kidder et al. (2011):** [Technical considerations to obtaining high-quality data](http://www.nature.com/ni/journal/v12/n10/abs/ni.2117.html), (doi:10.1038/ni.2117) - nice, readable introduction into all aspects of ChIP-seq experiments (from antibodies to cell numbers to replicates to data analysis)

**Leleu et al. (2010):** [Processing and analyzing ChIP-seq data](http://www.ncbi.nlm.nih.gov/pubmed/20861161), (doi: 10.1093/bfgp/elq022) - fairly detailed review of key concepts of ChIP-seq data processing (less detailed on analysis)

**Peter Park (2009):** [ChIP-seq: Advantages and challenges of a maturing technology](http://www.nature.com/nrg/journal/v10/n10/full/nrg2641.html), (doi:10.1038/nrg2641)

**Kharchenko et al. (2008):** [Design and analysis of ChIP-seq experiments for DNA-binding proteins](http://www.ncbi.nlm.nih.gov/pubmed/19029915), (doi:10.1038/nbt.1508)

**Liu et al. (2010):** [Q&A: ChIP-seq technologies and the study of gene regulation](http://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-56), (doi:10.1186/1741-7007-8-56) - short overview of several (typical) issues of ChIP-seq analysis

**Carroll et al. (2014):**  [Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data](http://journal.frontiersin.org/article/10.3389/fgene.2014.00075/full),(doi:10.3389/fgene.2014.00075)


#RNA-seq

**Nice graphical overview of the RNA-seq processing and analysis steps:** http://figshare.com/articles/RNA_seq_Workflows_and_Tools/662782

**Auer and Doerge (2010):** [Statistical Design and Analysis of RNA Sequencing Data](http://www.genetics.org/content/185/2/405), (doi:10.1534/genetics.110.114983); Insights into properly planning your RNA-seq run! Read BEFORE your RNA-seq experiment! 

**Ian Korf (2013):**[Genomics: the state of the art in RNA-seq analysis](http://www.ncbi.nlm.nih.gov/pubmed/24296473), (doi:10.1038/nmeth.2735); A refreshingly honest view on the non-trivial aspects of RNA-seq analysis

**Dillies et al. (2012):** [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](http://bib.oxfordjournals.org/content/early/2012/09/15/bib.bbs046.long), (doi:10.1093/bib/bbs046) - Systematic comparison of seven representative normalization methods for the differential analysis of RNA-seq data (Total Count, Upper Quartile, Median (Med), DESeq, edgeR, Quantile and Reads Per Kilobase per Million mapped reads (RPKM) normalization)

**Rapaport et al. (2013):** [Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data](http://www.genomebiology.com/2013/14/9/R95/abstract), (doi:10.1186/gb-2013-14-9-r95) - Evaluation of methods for differential gene expression analysis

**Soneson et al. (2013):** [A comparison of methods for differential expression analysis of RNA-seq data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-91), (doi:10.1186/1471-2105-14-91)

**Roberts et al. (2011):** [Improving RNA-Seq expression estimates by correcting for fragment bias Fragment bias correction](http://www.ncbi.nlm.nih.gov/pubmed/21410973), (doi: 10.1186/gb-2011-12-3-r22)

**Garber et al. (2011):** [Computational methods for transcriptome annotation and quantification using RNA-seq](http://www.ncbi.nlm.nih.gov/pubmed/21623353), (doi: 10.1038/nmeth.1613) - Classical paper about the computational aspects of RNA-seq data analysis

#Peak Calling Methods (ChIP-seq)

**Pepke et al. (2009):** [Computation for ChIP-seq and RNA-seq studies](http://www.ncbi.nlm.nih.gov/pubmed/19844228), (doi: 10.1038/nmeth.1371) - First comparison of peak callers - focuses on the explanation of basic principles of ChIP-seq data processing and the general workflows of peak calling algorithms

**Wilbanks et al. (2010):** [Evaluation of Algorithm Performance in ChIP-Seq Peak Detection](http://www.ncbi.nlm.nih.gov/pubmed/20628599), (doi: 10.1371/journal.pone.0011471) - Another comparison of peak callers - focuses more on the evaluation of the peak callers performances than Pepke et al.

**Micsinai et al. (2012):** [Picking ChIP-seq peak detectors for analyzing chromatin modification experiments](http://www.ncbi.nlm.nih.gov/pubmed/22307239), (doi: 10.1093/nar/gks048) - How to choose the best peak caller for your data set - their finding: default parameters, surprisingly, yield the most reproducible results regardless of the data set type

**MACS**

**Fen et al. (2012):** [Identifying ChIP-seq enrichment using MACS.](http://www.ncbi.nlm.nih.gov/pubmed/22936215), (doi:10.1038/nprot.2012.101) - How to use MACS - Nature Protocols

**Zhang et al. (2008):** [Model-based Analysis of ChIP-Seq (MACS)](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137), (doi:10.1186/gb-2008-9-9-r137) - The original publication of MACS: Model-bases analysis of ChIP-seq data.

##DNA motif analysis

**Das et al. (2007):** [A survey of DNA motif finding algorithms](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-S7-S21), (doi:10.1186/1471-2105-8-S7-S21) - Review of Motif Analysis Tools 

**MEME (suite)**

**Machanick and Bailey (2011):** [MEME-ChIP: motif analysis of large DNA datasets](http://www.ncbi.nlm.nih.gov/pubmed/21486936), (doi: 10.1093/bioinformatics/btr189) - MEME-ChIP-paper

**Bailey and Machanick (2012):** [Inferring direct DNA binding from ChIP-seq](http://www.ncbi.nlm.nih.gov/pubmed/22610855), (doi:10.1093/nar/gks433) - Centrimo: position-specific motif analysis, especially useful for ChIP-seq data

[TomTom](http://meme-suite.org/tools/tomtom) - Meme Suite Motif comparison tool - tool for the comparison of motifs from databases (not in Galaxy yet): [Manual](http://meme-suite.org/doc/tomtom.html?man_type=web)

**TRAP**

**Thomas-Chollier et al. (2012):** [Transcription factor binding predictions using TRAP for the analysis of ChIP-seq data and regulatory SNPs](http://www.ncbi.nlm.nih.gov/pubmed/22051799), (doi:10.1038/nprot.2011.409) -How to use TRAP - Nature Protocols

**Roider et al. (2006):** [Predicting transcription factor affinities to DNA from a biophysical model.](http://www.ncbi.nlm.nih.gov/pubmed/17098775), (doi:10.1093/bioinformatics/btl565) - Theoretical background of TRAP

##NGS data formats

- UCSC has a very good overview with brief descriptions of BED, bedGraph, bigWig etc.: https://genome.ucsc.edu/FAQ/FAQformat.html

- [VCF format](http://gatkforums.broadinstitute.org/gatk/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk) (encoding SNPs, indels etc.): Very readable, albeit not exhausting description
- Transcriptomes are often saved in [GFF3 format](http://www.sequenceontology.org/gff3.shtml) (this is what TopHat needs, for example), but just to make things more complicated, GTF is another format used for transcriptome information, too ([here](http://gmod.org/wiki/GFF2) are more information on GTF)


#Bioinformatic Tools (Linux, R, BEDTools etc.) - Manuals, Courses, original papers

- Why and how is bioinformatics software special? **Altschul et a. (2013)** [The anatomy of successful computational biology software](http://www.ncbi.nlm.nih.gov/pubmed/24104757), (doi:10.1038/nbt.2721) **(Highly recommended to read!)**
- **Bild et al. (2014)** [A Field Guide to Genomics Research](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001744), (doi:10.1371/journal.pbio.1001744) - Very readable introduction about the different caveats of genomics research (with cute cartoons!)
 
**Linux Command Line**

- [Linux & Perl Primer for Biologists](http://korflab.ucdavis.edu/Unix_and_Perl/unix_and_perl_v3.1.1.html) - Very entertaining introduction to command line commands and perl scripts with a focus on bioinformatic application, i.e. handling of DNA sequences
- [Linux Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/) - Thorough, but concise online tutorial introducing the very basics of handling the Linux command line
- [Writing Linux shell scripts](http://www.freeos.com/guides/lsst/index.html) - Useful for slightly more advanced Linux command line users

**R**

- [Hands on R course](http://www.uwyo.edu/mdillon/hor.html) - For beginners - R is probably the most widely used open-source statistical software; through our epicenter website you can also access RStudio which provides are very nice interface to working and plotting with R. In fact, most of the plots generated within Galaxy are generated through R scripts, so if you're not happy with the default formats of the Galaxy graphs, definitely have a look at R yourself. The learning curve is steep, but it is worth it.

**BEDTools** 

- [BEDTools Manual](https://github.com/bgruening/training_data/raw/master/Literature/BEDTools-User-Manual.pdf) - When working with genomic intervals (e.g. genes, peaks, enriched regions...), BEDTools are invaluable! The manual is a very good read and we refer to it almost daily.
