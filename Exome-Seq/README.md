[![Docker Automated buil](https://img.shields.io/docker/automated/bgruening/galaxy-training-exome-seq.svg?maxAge=2592000)](https://hub.docker.com/r/bgruening/galaxy-training-exome-seq/)
[![Docker Pulls](https://img.shields.io/docker/pulls/bgruening/galaxy-training-exome-seq.svg?maxAge=2592000)](https://hub.docker.com/r/bgruening/galaxy-training-exome-seq/)
[![Docker Stars](https://img.shields.io/docker/stars/bgruening/galaxy-training-exome-seq.svg?maxAge=2592000)](https://hub.docker.com/r/bgruening/galaxy-training-exome-seq/)

Training course and material for Exome sequencing
====

Exome sequencing means that all protein-coding genes in a genome are sequenced.

# Slides

Slides can be found here: http://bgruening.github.io/training-material/Exome-Seq/slides/index.html

# Tutorials

Two tutorials are available for training on exome sequencing data analysis.

In the [first tutorial](tutorials/Exome-Seq.md), you will work with a family
where the parents are healthy but a child has a yet unknown disease. The goal is
to identify genetic variation that is responsible for the disease using the exome
sequencing data from both parents and the child.

The [second tutorial](tutorials/Diploid-variant-calling.md) is more detailed. It
follows a similar pipeline to first tutorial on the bottle dataset, but with
more details particularly on the theory behind.

## Input datasets

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.60520.svg)](http://dx.doi.org/10.5281/zenodo.60520)

The input datasets for both tutorials are available on
[Zenodo with a dedicated DOI](http://dx.doi.org/10.5281/zenodo.60520).

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-exome-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080 ](http://localhost:8080).
