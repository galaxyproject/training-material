RNA-Seq data analysis
======================

RNA-sequencing is a method used to reveal the presence and quantity of RNA in a biological sample at a given moment in time. Here, we propose some materials to learn how to analyze RNA-seq data.

# Slides

Several deck of slides are available for this topic:

- [General introduction about RNA seq data analysis](http://bgruening.github.io/training-material/RNA-Seq/slides/)

# Tutorials

A tutorial with hands-on is available for this topic:

- [Reference-based RNA-seq data analysis](tutorials/ref_based.md)

## Input datasets

The input datasets for the tutorials will be soon available on Zenodo.

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-rna-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080 ](http://localhost:8080).

# References

## Papers

**Shirley Pepke et al:** [Computation for ChIP-seq and RNA-seq studies](http://www.nature.com/nmeth/journal/v6/n11s/full/nmeth.1371.html)


**Paul L. Auer & R. W. Doerge:** [Statistical Design and Analysis of RNA Sequencing Data](http://www.genetics.org/content/185/2/405)

> Insights into proper planning of your RNA-seq run! To read before any RNA-seq experiment!

**Ian Korf:**[Genomics: the state of the art in RNA-seq analysis](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2735.html)

> A refreshingly honest view on the non-trivial aspects of RNA-seq analysis

**Marie-Agnès Dillies et al:** [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis](http://bib.oxfordjournals.org/content/14/6/671)

> Systematic comparison of seven representative normalization methods for the differential analysis of RNA-seq data (Total Count, Upper Quartile, Median (Med), DESeq, edgeR, Quantile and Reads Per Kilobase per Million mapped reads (RPKM) normalization)

**Franck  Rapaport et al:** [Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r95)

> Evaluation of methods for differential gene expression analysis

**Charlotte Soneson & Mauro Delorenzi:** [A comparison of methods for differential expression analysis of RNA-seq data](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-91)

**Adam Roberts et al:** [Improving RNA-Seq expression estimates by correcting for fragment bias](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-3-r22)

**Manuel Garber et al:** [Computational methods for transcriptome annotation and quantification using RNA-seq](http://www.nature.com/nmeth/journal/v8/n6/abs/nmeth.1613.html)

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
