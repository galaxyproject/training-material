Training material
=================

- [Genome Annotation](https://github.com/bgruening/training-material/blob/master/genome-annotation/general-introduction/README.md)
- [RNA-seq](https://github.com/bgruening/training-material/blob/master/rna-seq/rna-seq.md)
- [ChIP-seq](https://github.com/bgruening/training-material/blob/master/ChIPseq/ChIPseq.md)
- [Exome-seq](https://github.com/bgruening/training-material/blob/master/Exome-Seq/Exome-Seq.md)
- [MethyC-seq](https://github.com/bgruening/training-material/blob/master/Methylation-Seq/Methylation-Seq.md)
- [Galaxy Introduction](https://github.com/bgruening/training-material/blob/master/Galaxy_Introduction/Galaxy_Introduction.md )
- [Data Sources](https://github.com/bgruening/training-material/blob/master/Data_Sources/Data_Sources.md)
  

## Literature

#Sequencing biases

**Meacham et al. (2011):** [Identification and correction of systematic error in high-throughput sequence data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451), (doi:10.1186/1471-2105-12-451) - Correcting for systematic base pair errors in deep sequencing; important paper if you want to look at any allele-specificy or if you're interested in SNPs

**Benjamini & Speed (2012):** [Summarizing and correcting the GC content bias in high-throughput sequencing](http://nar.oxfordjournals.org/content/40/10/e72.long), (doi: 10.1093/nar/gks001) - GC bias of deeply sequenced samples; very good paper that systematically assesses many possible sourced of GC bias for deeply sequenced samples and eventually pinpoints it to the DNA polymerase

**A collection of papers on quality controls for various NGS applications:** [Frontiers in Genetics (2014)](http://journal.frontiersin.org/researchtopic/1683/quality-assessment-and-control-of-high-throughput-sequencing-data)


#Deep sequencing

**Zentner and Henikoff (2012):** [Surveying the epigenomic landscape, one base at a time](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-250), (doi:10.1186/gb-2012-13-10-250) - Overview of popular *-seq techniques; very nice description of DNase-seq, MNase-seq, FAIRE-seq etc.

**Son and Taylor (2011):** [Preparing DNA Libraries for Multiplexed Paired-End Deep Sequencing for Illumina GA Sequencers](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076644/), (doi:10.1002/9780471729259.mc01e04s20) - Paper on multiplexing; describes the individual steps of the Illumina deep sequencing protocols quite in detail

**Illumina's technical report** - focuses on [Illumina's sequencing technology](http://www.illumina.com/technology.html); nice educative figures


#Mapping of short NGS reads

**informative slides** [Mapping of sequencing reads](http://people.binf.ku.dk/krogh/tmp/Mapping_Krogh_Monday.pdf) - Introduction to various aspects of NGS read mapping

**Fonseca et al. (2012):** [Tools for mapping high-throughput sequencing data](http://bioinformatics.oxfordjournals.org/content/28/24/3169.full), (doi:10.1093/bioinformatics/bts605) - An excellent starting point despite its "old" age, you will learn a lot about the different philosophies behind the read alignment tools!

**Hatem et al. (2013):** [Benchmarking short sequence mapping tools](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-184), (doi:10.1186/1471-2105-14-184) (spoiler alert: bowtie wins)

**Engstrom et al. (2013):** [Systematic evaluation of spliced alignment programs for RNA-seq data](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2722.html), (doi:10.1038/nmeth.2722) 

**Genome Mappability**

**Lee and Schatz (2012):** [The reliability of short read mapping](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413383/?report=reader), (doi:10.1093/bioinformatics/bts330)  - Very detailed paper about genome mappability issues that presents a new suite of tools for taking the mappability into account

*mappability maps* can be downloaded [here](http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map/)  


### License

This work is licensed under the [Creative Commons Attribution 3.0 Unported License](http://creativecommons.org/licenses/by/3.0/).
