Mapping of NGS data
===================

After [quality control](../NGS-QC/README.md), mapping sequences on reference genomes is the second step of most NGS data analysis process. 

# Slides

Several deck of slides are available for this topic:

- [General introduction about mapping](http://bgruening.github.io/training-material/NGS-mapping/slides/)
- Slide deck related to the tutorials:
    - [Dive into mapping](http://bgruening.github.io/training-material/Dev-Corner/slides/dive_into_mapping.html)

# Tutorials

A tutorial with hands-on is available for this topic:

- [Dive into mapping](tutorials/dive_into_mapping.md)

## Input datasets

The input datasets for the tutorial will be available soon on Zenodo

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-ngs-mapping-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).

# References

## Papers

**Nuno A. Fonseca et al:** [*Tools for mapping high-throughput sequencing data*](http://bioinformatics.oxfordjournals.org/content/28/24/3169.full)

> Excellent starting point despite its "old" age to learn a lot about the different philosophies behind the read alignment tools

**Ayat Hatem et al:** [*Benchmarking short sequence mapping tools*](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-184)

> Compares several mapping tools

**Pär G Engström et al:** [*Systematic evaluation of spliced alignment programs for RNA-seq data*](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2722.html)

> Compared 26 mapping protocols based on 11 programs

**Hayan Lee and Michael C. Schatz:** [Genomic dark matter: the reliability of short read mapping illustrated by the genome mappability score](http://bioinformatics.oxfordjournals.org/content/28/16/2097.short)

> Very detailed paper about genome mappability issues that presents a new suite of tools for taking the mappability into account

# Contributors

This material is maintained by:

- Maintainer 1
- Maintainer 2

For any question related to this topic and the content, you can contact them.

The following individuals have contributed to this training material:

- Name 1
- Name 2
