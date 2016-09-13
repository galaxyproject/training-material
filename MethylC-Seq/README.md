Training course and material for DNA Methylation data analysis
====

DNA methylation is an epigenetic mechanism used by higher eukaryotes and involved in eg. gene expression, X-Chromosome inactivating, imprinting, and gene silencing of germ-line specific gene and repetitive elements.

# Slides

A deck of slides will be available soon here: http://bgruening.github.io/training-material/MethylC-Seq/slides/index.html

# Tutorials

A [tutorial](tutorial/Methylation-Seq.md) uses a dataset from "Tissue-specific methylomes reveal epigenetic memory in adult mouse tissue" and call Bismarck BS Mapper to extract DNA methylation data.

## Input datasets

The input datasets for the tutorial will be available soon on Zenodo

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-methylc-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).
