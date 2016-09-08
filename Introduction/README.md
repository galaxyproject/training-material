Training course and material for Galaxy introduction
====

Here, you will find the material for introduction to Galaxy.

# Slides

Two deck of slides are available:

- [Introduction to Galaxy](./slides/Introduction_to_Galaxy_Uni.pdf)
- [Galaxy course](./slides/Manke_2015.09.21a.pdf)

# Tutorials

Several tutorials are available:

- [Galaxy Introduction Exercise: From Peaks to Genes](./tutorials/Introduction.md)
- [Galaxy 101-1: The first thing you should try](./tutorials/Galaxy101-1.md)
- [Galaxy 101-2: Getting to know workflows](./tutorials/Galaxy101-2.md)
- [Overview of sequencing technologies](./tutorials/NGS-technologies.md)
- [Processing many samples at once](./tutorials/Processing-many-samples-at-once.md)
- [Using the Integrative Genomics Viewer](./tutorials/IGV_Introduction.md)


## Input datasets

The input datasets for the tutorials will be soon available on Zenodo

## Galaxy instance

For this tutorial, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-introduction-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).
