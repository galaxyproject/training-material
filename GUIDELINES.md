Guidelines
===

First of all, **thank you** for contributing.


By contributing, you agree that we may redistribute your work under [our license](LICENSE.md), and to follow the Guidelines in this document. In exchange, we will address your issues and/or assess your change proposal as promptly as we can, and help you become a member of our community.

# How to contribute

The easiest way to start contributing is to file an issue to tell us about a spelling mistake or a factual error. You can then introduce yourself and meet some of our community members.

Once, you are feeling more comfortable, you can propose changes to this training material. We strongly recommend you to read and follow [Software Carpentry's recommendations on lesson writing](http://swcarpentry.github.io/instructor-training/19-lessons/) if you planned to add or change some training material.

To manage changes, we use [GitHub flow](https://guides.github.com/introduction/flow/):

1. Create a fork of this repository on GitHub
2. Clone your fork of this repository to create a local copy on your computer
3. Create a new branch in your local copy for each significant change
4. Commit the changes in that branch
5. Push that branch to your fork on GitHub
6. Submit a pull request from that branch to the [master repository](https://github.com/bgruening/training-material)
7. If you receive feedback, make changes in your local clone and push to your branch on GitHub: the pull request will update automatically

For beginners, the GitHub interface will help you in the process in edition of one file. It will automatically create a fork of this repository where you can safely work and then submit the changes as a pull request.

If you have any question, you can reach us using the [Gitter chat]().

# Structure of the training material

Each each training material is related to a topic. All training material (slides, tutorials, ...) related to a topic is found in a dedicated directory (*e.g.* `Exome-seq` directory contains the material related to exome sequencing analysis). These repository have to have the following structure (as in `Exome-seq` directory):

```
├── docker
│   ├── Dockerfile
│   ├── README.md
│   ├── tools.yaml
├── images
├── slides
│   ├── index.html
├── tutorials
├── README.md
```

## `images` directory

The `images` directory collect all images/pictures needed for the training materials related to the topic, *i.e* pictures for the slides or the tutorials.

Images shared between several topics are in the `shared/images` directory at the root.

All images for the slides must be in `images` directory. The images must be in good quality. The sources (`svg` or other) of the images must also be added to the `images` directory. We encourage you to use [yEd](http://www.yworks.com/products/yed) to easily generate diagrams and [Inkscape](https://inkscape.org/en/) for any other image.

## `tutorials` directory

This directory collect the tutorials related to the topic. The tutorials must
be written in `markdown` and follow the [`shared/templates/tutorials/tutorial.md` template](shared/templates/tutorials/tutorial.md). The tutorials must be well structured (to help the trainees) and we recommend to have at least one scheme or diagram to sum up the pipeline used (for example).

The input data required for the tutorials must be upload on [Zenodo](https://zenodo.org/) to obtain a dedicated DOI (in the [Galaxy training network community](https://zenodo.org/communities/galaxy-training/?page=1&size=20)).

## `slides` directory

One slide deck must be available for each tutorial and it follows the same structure (but with a more detailed introduction). The slide decks are in the `slides` directory. The slides
must be in [`reveal.js`](https://github.com/hakimel/reveal.js/) format, in a
`html` file. A template for the `html` file can be found in
[`shared/templates/slides/slides.html` file](shared/templates/slides/slides.html).
You can also use the [`reveal.js` editor](https://slides.com/?ref=github)
to help you.

Once the slides are on the `master` branch, they will be available at `http://bgruening.github.io/training-material/<topic>/slides/<tutorial_name>.html`

## `docker` directory

For each topic, a flavored Docker image must integrate the needed tools for
the tutorials. The corresponding image must be based on official Galaxy Docker
images. We recommend you to use the content of [`shared/templates/docker`](shared/templates/docker) as a template.

The `docker` image must also integrate a Galaxy tour from the [`galaxy-tours` repository](https://github.com/galaxyproject/galaxy-tours)

## `README.md` file

The `README.md` introduces rapidly the topic, summarizes all this information about the training materials (where to find the slides, the tutorials and their
topics, the Docker image and the input dataset) and add the needed references. We recommend you to use the [`shared/templates/README.md`](shared/templates/README.md) as a template.
