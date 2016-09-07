Training course and material for ChIP-seq analysis
====

ChIP-sequencing is a method used to analyze protein interactions with DNA.
Here, you will material to training for ChIP-seq analysis.

# Slides

Two deck of slides are available:

- [High throughput sequencing analysis](https://drive.google.com/file/d/0B9urRnOAUUI8UmwzbTVpdmZucWM/view?usp=sharing)

    This slide deck is an introduction to high throughput sequencing analysis and ChIP-seq.

- [From aligned reads to read densities and peaks](https://drive.google.com/file/d/0B9urRnOAUUI8cHpzYVBscjNKWEE/view?usp=sharing)

    This slide deck describe how to analyze ChIP-seq data.

# Tutorials

The [tutorial](tutorial/ChIPseq.md) uses a dataset from a Nature publication, to identify the binding sites of the Estrogen receptor, a transcription factor known to be associated with different types of breast cancer.

## Input datasets

The input datasets for the tutorial will be available soon on Zenodo

## Galaxy instance

For this tutorial, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-chip-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).
