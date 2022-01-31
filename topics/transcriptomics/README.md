RNA-Seq data analysis
======================

RNA-sequencing is a method used to reveal the presence and quantity of RNA in a biological sample at a given moment in time. Here, we propose some materials to learn how to analyze RNA-seq data.

# Slides

A deck of slides is available for this topic:

- [General introduction about RNA seq data analysis]({{site.url}}/topics/transcriptomics/slides/)

# Tutorials

A tutorial with hands-on is available for this topic:

- [Reference-based RNA-seq data analysis]({{site.url}}/topics/transcriptomics/tutorials/ref-based/tutorial.html)

## Input datasets

For de novo tutorial, data is available at [`Zenodo`](https://zenodo.org/record/254485).

For ref-based tutorial, the original data is available at NCBI Gene Expression Omnibus (GEO) under accession number [GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508). We will look at the 7 first samples (3 treated samples with Pasilla (PS) gene depletion: [GSM461179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461179), [GSM461180](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461180), [GSM461181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4611810) and 4 untreated samples: [GSM461176](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461176), [GSM461177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461177), [GSM461178](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461178), [GSM461182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461182)).

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/Dockerfile):

```
docker run -d -p 8080:80 bgruening/galaxy-rna-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080 ](http://localhost:8080).

# References

## Papers

**Shirley Pepke et al:** [Computation for ChIP-seq and RNA-seq studies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4121056/)


**Paul L. Auer & R. W. Doerge:** [Statistical Design and Analysis of RNA Sequencing Data](https://www.stat.purdue.edu/~doerge/BIOINFORM.D/SPRING10/auer_doerge_genetics_2010.pdf) DOI: 10.1534/genetics.110.114983

> Insights into proper planning of your RNA-seq run! To read before any RNA-seq experiment!

**Ian Korf:**[Genomics: the state of the art in RNA-seq analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4461013/)

> A refreshingly honest view on the non-trivial aspects of RNA-seq analysis

**Marie-Agnès Dillies et al:** [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](https://bib.oxfordjournals.org/content/14/6/671)

> Systematic comparison of seven representative normalization methods for the differential analysis of RNA-seq data (Total Count, Upper Quartile, Median (Med), DESeq, edgeR, Quantile and Reads Per Kilobase per Million mapped reads (RPKM) normalization)

**Franck  Rapaport et al:** [Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r95)

> Evaluation of methods for differential gene expression analysis

**Charlotte Soneson & Mauro Delorenzi:** [A comparison of methods for differential expression analysis of RNA-seq data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-91)

**Adam Roberts et al:** [Improving RNA-Seq expression estimates by correcting for fragment bias](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-3-r22)

**Manuel Garber et al:** [Computational methods for transcriptome annotation and quantification using RNA-seq](https://www.nature.com/nmeth/journal/v8/n6/abs/nmeth.1613.html)

> Classical paper about the computational aspects of RNA-seq data analysis

## Websites

**Stephen Turner:** [RNA-seq Workflows and Tools](https://figshare.com/articles/RNA_seq_Workflows_and_Tools/662782)

> Nice graphical overview of the RNA-seq processing and analysis step

# Contributors

This material is maintained by:

- Maintainer 1
- Maintainer 2

For any question related to this topic and the content, you can contact them.

The following individuals have contributed to this training material:

- Name 1
- Name 2
