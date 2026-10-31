# Literature

##Deep sequencing

**Zentner and Henikoff (2012):** [Surveying the epigenomic landscape, one base at a time](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-250), (doi:10.1186/gb-2012-13-10-250) - Overview of popular *-seq techniques; very nice description of DNase-seq, MNase-seq, FAIRE-seq etc.

**Son and Taylor (2011):** [Preparing DNA Libraries for Multiplexed Paired-End Deep Sequencing for Illumina GA Sequencers](https://www.ncbi.nlm.nih.gov/pubmed/21400673), (doi:10.1002/9780471729259.mc01e04s20) - Paper on multiplexing; describes the individual steps of the Illumina deep sequencing protocols quite in detail

**Illumina's technical report** - focuses on [Illumina's sequencing technology](https://www.illumina.com/technology.html); nice educative figures

##NGS data formats

- UCSC has a very good overview with brief descriptions of BED, bedGraph, bigWig etc.: https://genome.ucsc.edu/FAQ/FAQformat.html

- [VCF format](https://gatkforums.broadinstitute.org/gatk/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk) (encoding SNPs, indels etc.): Very readable, albeit not exhausting description

- Transcriptomes are often saved in [GFF3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) (this is what TopHat needs, for example), but just to make things more complicated, GTF is another format used for transcriptome information, too


##Bioinformatic Tools (Linux, R, BEDTools etc.) - Manuals, courses, original papers

- Why and how is bioinformatics software special? **Altschul et a. (2013)** [The anatomy of successful computational biology software](https://www.ncbi.nlm.nih.gov/pubmed/24104757), (doi:10.1038/nbt.2721) **(Highly recommended to read!)**
- **Bild et al. (2014)** [A Field Guide to Genomics Research](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001744), (doi:10.1371/journal.pbio.1001744) - Very readable introduction about the different caveats of genomics research (with cute cartoons!)

**Linux Command Line**

- [Linux & Perl Primer for Biologists](http://korflab.ucdavis.edu/Unix_and_Perl/unix_and_perl_v3.1.1.html) - Very entertaining introduction to command line commands and perl scripts with a focus on bioinformatic application, i.e. handling of DNA sequences
- [Linux Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/) - Thorough, but concise online tutorial introducing the very basics of handling the Linux command line
- [Writing Linux shell scripts](http://www.freeos.com/guides/lsst/index.html) - Useful for slightly more advanced Linux command line users

**R**

- [Hands on R course](https://www.uwyo.edu/mdillon/hor.html) - For beginners - R is probably the most widely used open-source statistical software; through our epicenter website you can also access RStudio which provides are very nice interface to working and plotting with R. In fact, most of the plots generated within Galaxy are generated through R scripts, so if you're not happy with the default formats of the Galaxy graphs, definitely have a look at R yourself. The learning curve is steep, but it is worth it.

**BEDTools**

- [BEDTools Manual](https://bedtools.readthedocs.org) - When working with genomic intervals (e.g. genes, peaks, enriched regions...), BEDTools are invaluable! The manual is a very good read and we refer to it almost daily.
