Galaxy introduction
===================

Galaxy is a scientific workflow, data integration, and data and analysis persistence and publishing platform that aims to make computational biology accessible to research scientists that do not have computer programming experience.

Here, you will find some material to learn how to use Galaxy.

# Slides

Several deck of slides are available for this topic:

- [General introduction about Galaxy](http://bgruening.github.io/training-material/Introduction/slides/)
- Slide deck related to the tutorials:
    - [Galaxy Introduction Exercise: From Peaks to Genes](http://bgruening.github.io/training-material/Dev-Corner/slides/introduction.html)
    - [Getting to know workflows](http://bgruening.github.io/training-material/Dev-Corner/slides/workflows.html)
    - [Processing many samples at once](http://bgruening.github.io/training-material/Dev-Corner/slides/processing_many_samples.html)
    - [Using the Integrative Genomics Viewer](http://bgruening.github.io/training-material/Dev-Corner/slides/igv.html)

# Tutorials

Several tutorials with hands-on are available for this topic:

- [Galaxy Introduction Exercise: From Peaks to Genes](tutorial/introduction.md)
- [Getting to know workflows](tutorial/workflows.md)
- [Processing many samples at once](tutorial/processing_many_samples.md)
- [Using the Integrative Genomics Viewer](./tutorials/igv.md)

## Input datasets

The input datasets for the tutorials will be soon available on Zenodo

## Galaxy instance

For this tutorial, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-introduction-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).

# Contributors

This material is maintained by:

- Maintainer 1
- Maintainer 2

For any question related to this topic and the content, you can contact them.

The following individuals have contributed to this training material:

- Name 1
- Name 2
