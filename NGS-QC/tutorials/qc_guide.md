#NGS technologies most used
![NGS applications ranked by number of publications in PubMed](Pubmed_Sequencing_2016_06_04.svg)

#Collection of papers about NGS-QC and short summary

##A quality control system for profiles obtained by ChIP sequencing

[A quality control system for profiles obtained by ChIP sequencing][Mendoza-Parra2013]

[Mendoza-Parra2013]: http://nar.oxfordjournals.org/lookup/doi/10.1093/nar/gkt829 "A quality control system for profiles obtained by ChIP sequencing"

  * ChIP-Seq
  * RNA-Seq
  * Gro-Seq
  * footprinting
  * Epigen-Seq
  * IP-based 
      * CL-efficiency
      * Digestion
      * Shearing
      * Affinity and Accuracy of Antibody
  * Current method
      * Visual comparison
      * Peak finder
  * ENCODE
      * FRiP
      * IDR
  * Sequencing depth depends on
      * Peak profile
      * binding pattern
      * type of study
  * What depth saturated binding?
      * ChIP
          * sampling -> compare peaks
              *  does \# of peaks change of only peak height?
                  *  \# of peaks bad, peak heigth ok, as less signal was used
      *  RNA-Seq
          *  most stable as no manipulation like crosslinking or IP neccessary
      *  IP-based methods
          *  stability of peaks strongly depending on target

  *  How to get a hint if quality of samples is ok?
      *  Quality Control Generator
          *  [NGS Quality Control Generator][NGSQC]
			 [NGSQC]: http://www.ngs-qc.org/ "NGS Quality Control Generator"
          *  Allows comparison of sample with a database of samples scored by quality 
  

##A survey of best practices for RNA-seq data analysis

[A survey of best practices for RNA-seq data analysis][Conesa2016]
[Conesa2016]: http://genomebiology.com/2016/17/1/13 "A survey of best practices for RNA-seq data analysis"

  * Quality Control Checkpoints
      * Raw Reads
      * Read alignment
      * Quantification
      * Reproducibility
      * Transcript identification
      * Alignment
      * Transcript discovery
      * De novo transcript reconstruction
      * Transcript quantification
      * Differential gene expression analysis
      * Alternative splicing analysis
      * Visualization
      * Gene fusion discovery
      * Small RNAs
      * Functional profiling
      * Integration with other data types
          * DNA sequencing
          * DNA methylation
          * Chromatin features
          * MicroRNAs
          * Proteomics and metabolomics
          * Integration and visualization of multiple data types
      * Single cell RNA-Seq
      * Long-read sequencing
      * 

[RNA-Seq analysis roadmap][Conesa2016]
[Conesa2016]: https://static-content.springer.com/image/art%3A10.1186%2Fs13059-016-0881-8/MediaObjects/13059_2016_881_Fig1_HTML.gif "A survey of best practices for RNA-seq data analysis, Figure 1"
	
[Read mapping and transcript identification strategies][Conesa2016]
[Conesa2016]: https://static-content.springer.com/image/art%3A10.1186%2Fs13059-016-0881-8/MediaObjects/13059_2016_881_Fig2_HTML.gif "A survey of best practices for RNA-seq data analysis, Figure 2"

##Sequencing depth and coverage: key considerations in genomic analyses

[Sequencing depth and coverage: key considerations in genomic analyses][Sims2014]
[Sims2014]: http://www.nature.com/doifinder/10.1038/nrg3642 "Sequencing depth and coverage: key considerations in genomic analyses"


##Heat*seq: an interactive web tool for high-throughput sequencing experiment comparison with public data
[HeatSeq][HeatSeq]

[Devailly2016]: http://bioinformatics.oxfordjournals.org/lookup/doi/10.1093/bioinformatics/btw407 "HeatSeq"

##BaRC Standard operating procedures
[BaRC-SOP][BaRC-SOP]

[BaRC-SOP]: http://barcwiki.wi.mit.edu/wiki/SOPs "BaRC-SOP"

##BaRC Standard operating procedures
[RseQC][RseQC]

[RseQC]: http://rseqc.sourceforge.net/ "RseQC"


