Guidelines
===

First of all, **thank you** for contributing.

Here are a few rules to follow in order to add and maintain the training material.

Each each training material is related to a topic. All training material (slides, tutorials, ...) related
to a topic is found in a dedicated directory (*e.g.* `Exome-seq` directory contains
the material related to exome sequencing analysis). These repository have to have the
following structure (as in `Exome-seq` directory):

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

- `images` directory

  The `images` directory collect all images/pictures needed for the training materials
  related to the topic, *i.e* pictures for the slides or the tutorials.

- `tutorials` directory

  This directory collect the tutorials related to the topic. The tutorials must
  be written in `markdown`. The input data required for the tutorials must be
  upload on [Zenodo](https://zenodo.org/) to obtain a dedicated DOI.

- `slides` directory

  The slides related to the topic must be in the `slides` directory. The slides
  must be in [`reveal.js`](https://github.com/hakimel/reveal.js/) format, in a
  `index.html` file. A template for `index.html` file can be found in
  [`shared/shared/slide_template/index.html` file](shared/shared/slide_template/index.html).
  You can also use the [`reveal.js` editor](https://slides.com/?ref=github)
  to help you. The slides will be available at "http://bgruening.github.io/training-material/<topic>/slides/index.html"

- `docker` directory

  For each topic, a flavored Docker image must integrate the needed tools for
  the tutorials. The corresponding image must be based on official Galaxy Docker
  images. Check [`Exome-seq` Docker directory for an example](Exome-Seq/docker/).
  The `docker` image must also integrate a Galaxy tour from the [`galaxy-tours` repository](https://github.com/galaxyproject/galaxy-tours)

- `README.md` file

  The `README.md` introduces rapidly the topic and summarizes all this information
  with the links (where to find the slides, the tutorials and their
  topics, the Docker image and the input dataset). It also collect the interesting
  badges (Docker, DOI, ...). You can check [`README.md` of `Exome-seq`](Exome-Seq/README.md)
  for an example.
