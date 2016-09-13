Training course and material for RNA-Seq analysis
====

Quality control is an important preprocessing step in NGS data analyses.

# Slides

A deck of slides will be available soon here: http://bgruening.github.io/training-material/RNA-Seq/slides/index.html

# Tutorials

A [tutorial](tutorial/qc_guide.md) is the starting point for an exhaustive NGS-QC SOP collection.

## Input datasets

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.61771.svg)](http://dx.doi.org/10.5281/zenodo.61771)

The input datasets for both tutorials are available on
[Zenodo with a dedicated DOI](http://dx.doi.org/10.5281/zenodo.61771).

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-rna-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).
