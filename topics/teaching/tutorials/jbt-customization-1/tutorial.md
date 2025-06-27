---
layout: tutorial_hands_on

title: Basics of JupyterLab customization
subtopic: practises
draft: true
time_estimation: 1h
questions:
  - How can I make basic customizations to JupyterLab?
objectives:
  - An introduction is provided on how to make basic JupyterLab customizations, e.g. installing additional plugins.
key_points:
  - Basic customization costs litte effort and makes it easier for your course participants to get started.
contributors:
  - mittler-works
---

# Prerequisites

One common way to serve JupyterLab is the container-based approach, especially when using a JupyterHub infrastructure provider.

This is why this tutorial also follows the container-based approach. As a teacher who wants to customize a container-based JupyterLab you need to install a container runtime as well as some related tools for building container images and more.

One of the easier ways to install all the required dependencies is the installation of docker. Docker is available on all major platforms, is very easy to use and has also the advantage, that your produced images will be likely working with any infrastructure provider.

All examples in this tutorial are therefore developed with docker in mind, but they should of course also work with other container platforms.

All commands that you should run in a terminal are given with linux in mind. On windows and mac, you might adapt them, e.g. it might not be neccessary (or possible) to use sudo on windows.

## Install docker

In order to install docker, please follow the official installation steps. Please read through the whole installation documentation first before starting the installation.

- for Linux, I recommend to install docker engine only\
[https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)

- for Windows, you need docker dekstop - Please use the WSL 2 backend to ensure proper function\
[https://docs.docker.com/desktop/setup/install/windows-install/](https://docs.docker.com/desktop/setup/install/windows-install/)

- for Mac, you also need docker desktop\
[https://docs.docker.com/desktop/setup/install/mac-install/](https://docs.docker.com/desktop/setup/install/mac-install/)

# Start your first local JupyterLab with docker

When you have installed docker successfully, you may now start your first containerized JupyterLab.

In order to run a container you need a container image to start from. The container image defines all the content of the running container. For JupyterLab, there are some prebuild container images for many usecases:

[https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html)

You may select an image that fits closest to your usecase. For this tutorial, we'll use the `jupyter/minimal-notebook` image.

To run this image, open up a terminal and execute:

```bash
sudo docker run --rm -p 127.0.0.1:8888:8888 quay.io/jupyter/minimal-notebook:2025-06-23
```

Let's tear this command apart:

|---|---|
| sudo | Executes the command with elevated privileges. This is neccassary as docker server is installed as privileged daemon by default. |
| docker | The docker binary. Allows to talk to the docker daemon. |
| run | Tells docker to create and start a new container |
| --rm | Tells docker to remove the container as soon as it's stopped |
| -p | Tells docker to setup a port forwarding |
| 127.0.0.1:8888:8888 | The argument for `-p`. A port forwarding on the local interface 127.0.0.1 from port 8888 to the container port 8888 |
| quay.io/jupyter/minimal-notebook:2025-06-23 | The image to create the container from |

## Access your local Notebook Server

The above command will print out a link to `http://localhost:8888/lab?token=<token>`. Copy this link into your browser or just right-click-open the link. You will then see your JupyterLab loading.

When the interface is booted, click on the `Python 3` tile in the Notebook tab to open up a new Notebook. Verify, that everything is working by printing out today's weekday:

```python
# Hint: press shift + enter to execute code in Jupyter Notebook
import datetime
now = datetime.datetime.now()
print(now.strftime("%A"))
```

## Missing packages

You are now in a predefined environment. The packages that are installed are defined by the container image you started, i.e. the `minimal-notebook`. All Packages from this image are available, which are the python default libraries. so you can easily succeed following above steps to print out the weekday, as the used packages are default packages.

For scientific work however, you will likely need more than just the default python libraries. Reviewing the docker stacks is helpful and maybe an image like the `scipy-notebook` will have most of the packages you need. But of course these default collections cannot have perfect fit package selections for the very specific needs of all courses that you could imagine of.

For example, in this minimal setup you might like to have the option to print out the weekday in a fancier way. Luckily, there is a python package for this: `cowsay`. So let's go ahead, try running:

```python
# Hint: press shift + enter to execute code in Jupyter Notebook
import datetime
from cowsay import cow

now = datetime.datetime.now()
cow(now.strftime("%A"))
```

As you will see, this fails with a `ModuleNotFoundError` as cowsay is not available in minimal-notebook.

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
