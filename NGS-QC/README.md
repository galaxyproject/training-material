Training course and material for NGS quality control
====

Quality control is an important preprocessing step in NGS data analyses.

# Slides

A deck of slides will be available soon here: http://bgruening.github.io/training-material/NGS-QC/slides/index.html

# Tutorials

A [tutorial](tutorial/qc_guide.md) is the starting point for an exhaustive NGS-QC SOP collection.

## Input datasets

The input datasets for the tutorial will be available soon on Zenodo

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-ngs-qc-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).
