---
layout: tutorial_hands_on

title: Basics of Jupyter customization
subtopic: practises
draft: true
time_estimation: 1h
questions:
  - How can I make basic customizations to Jupyter?
objectives:
  - An introduction is provided on how to make basic Jupyter customizations, e.g. installing additional plugins.
key_points:
  - Basic customization costs litte effort and makes it easier for your course participants to get started.
contributors:
  - mittler-works
---

# Prerequisites

* We'll use containerized installation

* Examples with docker, as it is easy to install on all major platforms

# Start your first local Notebook Server with docker

* Show docker stacks

* Select an image

* WIP

```bash
sudo docker run --rm -it -p 127.0.0.1:8888:8888 quay.io/jupyter/minimal-notebook:2025-06-23
```

* Will print out a link to `http://localhost:8888/lab?token=<token>`

## Access your local Notebook Server

* Windows and Mac users will need verification work...

* WIP

## Missing packages

* Install a simple python package

* WIP

* You are now in a predefined environment. All Packages from mentioned in the docker stacks documentation are there, so you can easily open the `Python 3` Notebook presented on the start page and start working. As you are in a minimal notebook you just do have the python default libraries, but you could run this to get the current weekday:

```python
import datetime
now = datetime.datetime.now()
print(now.strftime("%A"))
```

* For scientific work you will likely need more than just the default python libraries. Reviewing the docker stacks is helpful and maybe an image like the `scipy-notebook` will have most of the packages you need. But of course these default collections cannot have suitable package selections for the very specific needs of all courses out in the world.

For example, in this minimal setup you might like to have the option to print out text in a fancier way. Luckily, there is a python package for this: `cowsay`. So let's go ahead, try running

```python
from cowsay import cow
cow('Muuuuh')
```

As you will see, this fails as cowsay is not available in minimal-notebook.


# Build your first custom Image

In order to customize this Container Image, we need a so-called `Dockerfile`. Go ahead, create a new directory as your playground. Then, create a file called `Dockerfile` in this new directory.

Fill it with content:

```dockerfile
# First, we say from which container image we want to start
# This is called the base image
FROM quay.io/jupyter/minimal-notebook:2025-06-23

# Now, we can run arbitrary commands in order to extend the base image
RUN pip install --quiet --no-cache-dir cowsay
```

* Build it, mind the period at the end of the command

```bash
sudo docker build -t my-custom-jupyterlab .
```


# Use your first custom Image

```bash
sudo docker run --rm -p 127.0.0.1:8888:8888 my-custom-jupyterlab
```

Again, you'll provided with a link to open you jupyterlab. Start a Pyhton 3 Notebook by clicking the tile.

Again, run:

```python
from cowsay import cow
cow('Muuuuh')
```

As you will see, the command now succeeds and prints a nice graphical text output.

# Exercise: build your own custom Image

Now, let's go ahead and try out some real usecase.

> Create a custom notebook from the `scipy-notebook` base image and install the `igv-browser` extension.

You succeeded, if you can run following snippet successfully and are provided with a genome browser with loaded track:

> <hands-on-title>Make this script succeed in your custom JupyterLab</hands-on-title>
>
>
> ```python
> # This example snippet is taken from https://github.com/igvteam/igv-notebook/blob/v3.1.4/README.md
> # License: MIT License
> 
> import igv_notebook
> igv_notebook.init()
> igv_browser= igv_notebook.Browser(
>     {
>         "genome": "hg19",
>         "locus": "chr22:24,376,166-24,376,456",
>         "tracks": [{
>             "name": "BAM",
>             "url": "https://s3.amazonaws.com/igv.org.demo/gstt1_sample.bam",
>             "indexURL": "https://s3.amazonaws.com/igv.org.demo/gstt1_sample.bam.bai",
>             "format": "bam",
>             "type": "alignment"
>         }],
>         "roi": [
>             {
>                 "name": "ROI set 1",
>                 "url": "https://s3.amazonaws.com/igv.org.test/data/roi/roi_bed_1.bed",
>                 "indexed": False,
>                 "color": "rgba(94,255,1,0.25)"
>             }
>         ]
>     }
> )
> ```
{: .question}

> <solution-title>Proposed approach</solution-title>
>
> Keep in mind, there are many ways to solve this exercise. This proposed one is only one way to get the job done.
> 
> 1. Create a new `Dockerfile`:
> ```dockerfile
> FROM quay.io/jupyter/scipy-notebook:2025-06-23
> RUN pip install --quiet --no-cache-dir igv-notebook
> ```
> 2. Build the new Container Image:
> ```bash
> sudo docker build -t my-custom-jupyterlab .
> ```
> 3. Run a container from your image
> ```bash
> sudo docker run --rm -p 127.0.0.1:8888:8888 my-custom-jupyterlab
> ```
> 4. Access the JupyterLab via URL provided by the container output
> 5. Open a Python 3 Jupyter Notebook, paste in the above code and run it. It should work now.
{: .solution}

# Docker Compose

While running all these docker commands into the command line directly works fine, it has significant disadvantages: it's verbose and therefore error-prone, it's difficult to reproduce and it's difficult to share and version.

A better approach is to use a `docker-compose` setup. All relevant config is written into a file called `compose.yml`, which is sharable, reproducable and versionizable.

```yaml
services:
  jupyterlab:
    image: my-custom-jupyterlab
    build:
      context: .
      dockerfile: Dockerfile
    ports:
      - 127.0.0.1:8888:8888
```

You may now build the image simply by executing:

```
sudo docker compose build
```

An you may run it easily by executing:

```
sudo docker compose up
```
