Quality control of NGS data
===========================

Quality control is an important preprocessing step in NGS data analyses. Here, we propose some materials to learn how to control quality of your NGS data using Galaxy instance.

# Slides

Several deck of slides are available for this topic:

- [General introduction about quality control](http://bgruening.github.io/training-material/NGS-QC/slides/)
- Slide deck related to the tutorials:
    - [Dive into quality control](http://bgruening.github.io/training-material/Dev-Corner/slides/dive_into_qc.html)

# Tutorials

A tutorial with hands-on is available for this topic:

- [Dive into quality control](tutorials/dive_into_qc.md)

## Input datasets

The input datasets for the tutorial will be available soon on Zenodo

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-ngs-qc-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).

# References

## Papers

**Marco-Antonio Mendoza-Parra et al:** [*A quality control system for profiles obtained by ChIP sequencing*](http://nar.oxfordjournals.org/content/41/21/e196.short)

> Presents an approach that associates global and local QC indicators to ChIP-seq data sets as well as to a variety of enrichment-based studies by NGS.

**Ana Conesa et al:** [*A survey of best practices for RNA-seq data analysis*](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8)

> Highlights the challenges associated with each step of RNA-seq data analysis

**Yuval Benjamini1 et al:** [*Summarizing and correcting the GC content bias in high-throughput sequencing*](http://nar.oxfordjournals.org/content/40/10/e72.long)

> Summarizes the many possible sourced of GC bias for deeply sequenced samples

**David Sims et al:** [*Sequencing depth and coverage: key considerations in genomic analyses*](http://www.nature.com/nrg/journal/v15/n2/abs/nrg3642.html)

> Discuss the issue of sequencing depth in the design of next-generation sequencing experiments

**Frontiers in Genetics:** [*Quality assessment and control of high-throughput sequencing data*](http://journal.frontiersin.org/researchtopic/1683/quality-assessment-and-control-of-high-throughput-sequencing-data)

> Collection of papers on quality controls for various NGS applications

**Frazer Meacham et al:** [*Identification and correction of systematic error in high-throughput sequence data*](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451)

> Reviews and proposes correction for systematic base pair errors in deep sequencing

**Guillaume Devailly et al:** [*Heat seq: an interactive web tool for high-throughput sequencing experiment comparison with public data*](http://biorxiv.org/content/early/2016/04/18/049254.abstract)

> Presents a web-tool that allows genome scale comparison of high throughput experiments

## Websites

**Bioinformatics and Research Computing:** [*Quality Control and preprocessing of short reads*](http://barcwiki.wi.mit.edu/wiki/SOPs/qc_shortReads)

# Contributors

This material is maintained by:

- Maintainer 1
- Maintainer 2

For any question related to this topic and the content, you can contact them.

The following individuals have contributed to this training material:

- Name 1
- Name 2
