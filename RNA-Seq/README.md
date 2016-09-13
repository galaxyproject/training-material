Training course and material for RNA-Seq analysis
====

Quality control is an important preprocessing step in NGS data analyses.

# Slides

A deck of slides will be available soon here: http://bgruening.github.io/training-material/RNA-Seq/slides/index.html

# Tutorials

Two tutorials are available:

In the [first tutorial](tutorials/RNA-Seq-hands-on.md), you will work from the study by Brooks et al. 2011, in which the pasilla gene in Drosophila melanogaster was depleted by RNAi and the effects on splicing events were analysed by RNA-seq.

The [second tutorial](tutorials/Reference-based-RNA-seq-detailed.md) is more detailed. It follows a similar pipeline to first tutorial on the bottle dataset, but with more details particularly on the theory behind.

## Input datasets

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.61771.svg)](http://dx.doi.org/10.5281/zenodo.61771)

The input datasets for first tutorial are available on
[Zenodo with a dedicated DOI](http://dx.doi.org/10.5281/zenodo.61771).

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-rna-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).
