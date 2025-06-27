---
layout: tutorial_hands_on

title: Featured Jupyter customizations
subtopic: practises
draft: true
time_estimation: 30m
questions:
  - Which customizations are battle-tested and recommended for teaching?
objectives:
  - An overview about customizations used by university of giessen.
key_points:
  - Don't reinvent the wheel.
contributors:
  - mittler-works
---

# Introduction

We present adaptions that the JLU Physics Department has been successfully using in its courses for over a year and consider to be a great addition to the default Jupyter notebook server.

# nbgitpuller

https://github.com/jupyterhub/nbgitpuller

Allows a teacher to host his teaching material in a git repository.

Teacher may update the repo while course is running in order to provide new exercises or to provide solutions for older exercises.

nbgitpuller will will sync if either

- the git puller url is opened
- the notebook is started and git puller is configured as default url

## Install

```dockerfile
FROM quay.io/jupyter/minimal-notebook:2025-06-23

RUN pip install nbgitpuller
```

## Configuration

Link generator: https://nbgitpuller.readthedocs.io/en/latest/link.html

No special config needed, you just generate a link and call it

Hint: when using nbgitpuller with a standalone JupyterLab and not JupyterHub, you need to remove the `hub/user-redirect/` part from the resulting URL!

So, instead of `http://localhost:8888/hub/user-redirect/git-pull?repo=...` you should use `http://localhost:8888/git-pull?repo=...`

## Use it

## Exercise: create a git puller url for the igv-notebook repo

From the earlier tutorial, expand the JupyterLab with igv-notebook installed: install nbgitpuller inside this JupyterLab.

Then goto link generator and create a git puller link for https://github.com/igvteam/igv-notebook

Configure it so the `BamFiles_Lab.ipynb` example openes automatically in the JupyterLab.

solution: configure the link generator like this:

| --- | --- |
| JupyterHub URL | http://localhost:8888 |
| Git Repository URL | https://github.com/igvteam/igv-notebook.git |
| branch | main |
| File to open | examples/BamFiles_Lab.ipynb |
| Application to Open | JupyterLab |
| Named Server to open | (leave it blank) |


Strip off the `hub/user-redirect/` part.

The link will be:

http://localhost:8888/git-pull?repo=https%3A%2F%2Fgithub.com%2Figvteam%2Figv-notebook.git&urlpath=tree%2Figv-notebook.git%2Fexamples%2FBamFiles_Lab.ipynb&branch=main

It should pull the repo and opens the example automatically. You should be able to run the example if you extended the igv-notebook tutorial.

# nbgrader

https://github.com/jupyter/nbgrader

## Install

## Configuration

## Usage

# Config customizations

## Dos and Don'ts

# Setup Scripts

## Where to place them what to do with them
