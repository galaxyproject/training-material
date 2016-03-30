Training material
=================

- [Genome Annotation](https://github.com/bgruening/training-material/blob/master/genome-annotation/general-introduction/README.md)
- [RNA-seq](https://github.com/bgruening/training-material/blob/master/rna-seq/rna-seq.md)
- [ChIP-seq](https://github.com/bgruening/training-material/blob/master/ChIPseq/ChIPseq.md)
- [Exome-seq](https://github.com/bgruening/training-material/blob/master/Exome-Seq/Exome-Seq.md)
- [MethyC-seq](https://github.com/bgruening/training-material/blob/master/Methylation-Seq/Methylation-Seq.md)
- [Galaxy Introduction](https://github.com/bgruening/training-material/blob/master/Galaxy_Introduction/Galaxy_Introduction.md )
- [Data Sources](https://github.com/bgruening/training-material/blob/master/Data_Sources/Data_Sources.md)
  

# Literature

##Sequencing biases

**Meacham et al. (2011):** [Identification and correction of systematic error in high-throughput sequence data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451), (doi:10.1186/1471-2105-12-451) - Correcting for systematic base pair errors in deep sequencing; important paper if you want to look at any allele-specificy or if you're interested in SNPs

**Benjamini & Speed (2012):** [Summarizing and correcting the GC content bias in high-throughput sequencing](http://nar.oxfordjournals.org/content/40/10/e72.long), (doi: 10.1093/nar/gks001) - GC bias of deeply sequenced samples; very good paper that systematically assesses many possible sourced of GC bias for deeply sequenced samples and eventually pinpoints it to the DNA polymerase

**A collection of papers on quality controls for various NGS applications:** [Frontiers in Genetics (2014)](http://journal.frontiersin.org/researchtopic/1683/quality-assessment-and-control-of-high-throughput-sequencing-data)


##Deep sequencing

**Zentner and Henikoff (2012):** [Surveying the epigenomic landscape, one base at a time](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-250), (doi:10.1186/gb-2012-13-10-250) - Overview of popular *-seq techniques; very nice description of DNase-seq, MNase-seq, FAIRE-seq etc.

**Son and Taylor (2011):** [Preparing DNA Libraries for Multiplexed Paired-End Deep Sequencing for Illumina GA Sequencers](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076644/), (doi:10.1002/9780471729259.mc01e04s20) - Paper on multiplexing; describes the individual steps of the Illumina deep sequencing protocols quite in detail

**Illumina's technical report** - focuses on [Illumina's sequencing technology](http://www.illumina.com/technology.html); nice educative figures


##Mapping of short NGS reads

**informative slides** [Mapping of sequencing reads](http://people.binf.ku.dk/krogh/tmp/Mapping_Krogh_Monday.pdf) - Introduction to various aspects of NGS read mapping

**Fonseca et al. (2012):** [Tools for mapping high-throughput sequencing data](http://bioinformatics.oxfordjournals.org/content/28/24/3169.full), (doi:10.1093/bioinformatics/bts605) - An excellent starting point despite its "old" age, you will learn a lot about the different philosophies behind the read alignment tools!

**Hatem et al. (2013):** [Benchmarking short sequence mapping tools](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-184), (doi:10.1186/1471-2105-14-184) (spoiler alert: bowtie wins)

**Engstrom et al. (2013):** [Systematic evaluation of spliced alignment programs for RNA-seq data](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2722.html), (doi:10.1038/nmeth.2722) 

**Genome Mappability**

**Lee and Schatz (2012):** [The reliability of short read mapping](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413383/?report=reader), (doi:10.1093/bioinformatics/bts330)  - Very detailed paper about genome mappability issues that presents a new suite of tools for taking the mappability into account

*mappability maps* can be downloaded [here](http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map/)  

##NGS data formats

- UCSC has a very good overview with brief descriptions of BED, bedGraph, bigWig etc.: https://genome.ucsc.edu/FAQ/FAQformat.html

- [VCF format](http://gatkforums.broadinstitute.org/gatk/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk) (encoding SNPs, indels etc.): Very readable, albeit not exhausting description

- Transcriptomes are often saved in [GFF3 format](http://www.sequenceontology.org/gff3.shtml) (this is what TopHat needs, for example), but just to make things more complicated, GTF is another format used for transcriptome information, too ([here](http://gmod.org/wiki/GFF2) are more information on GTF)


##Bioinformatic Tools (Linux, R, BEDTools etc.) - Manuals, courses, original papers

- Why and how is bioinformatics software special? **Altschul et a. (2013)** [The anatomy of successful computational biology software](http://www.ncbi.nlm.nih.gov/pubmed/24104757), (doi:10.1038/nbt.2721) **(Highly recommended to read!)**
- **Bild et al. (2014)** [A Field Guide to Genomics Research](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001744), (doi:10.1371/journal.pbio.1001744) - Very readable introduction about the different caveats of genomics research (with cute cartoons!)
 
**Linux Command Line**

- [Linux & Perl Primer for Biologists](http://korflab.ucdavis.edu/Unix_and_Perl/unix_and_perl_v3.1.1.html) - Very entertaining introduction to command line commands and perl scripts with a focus on bioinformatic application, i.e. handling of DNA sequences
- [Linux Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/) - Thorough, but concise online tutorial introducing the very basics of handling the Linux command line
- [Writing Linux shell scripts](http://www.freeos.com/guides/lsst/index.html) - Useful for slightly more advanced Linux command line users

**R**

- [Hands on R course](http://www.uwyo.edu/mdillon/hor.html) - For beginners - R is probably the most widely used open-source statistical software; through our epicenter website you can also access RStudio which provides are very nice interface to working and plotting with R. In fact, most of the plots generated within Galaxy are generated through R scripts, so if you're not happy with the default formats of the Galaxy graphs, definitely have a look at R yourself. The learning curve is steep, but it is worth it.

**BEDTools** 

- [BEDTools Manual](http://bedtools.readthedocs.org) - When working with genomic intervals (e.g. genes, peaks, enriched regions...), BEDTools are invaluable! The manual is a very good read and we refer to it almost daily.



# License

This work is licensed under the [Creative Commons Attribution 3.0 Unported License](http://creativecommons.org/licenses/by/3.0/).
