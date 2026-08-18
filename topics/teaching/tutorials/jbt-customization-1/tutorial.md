---
layout: tutorial_hands_on

title: Basic JupyterLab customizations
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

> <comment-title>Docker CLI</comment-title>
> The complete docker cli reference can be found at [https://docs.docker.com/reference/cli/docker/](https://docs.docker.com/reference/cli/docker/)
{: .comment}

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

Of course, you could run `pip install cowsay` in order to have the package available in your notebook. For a simple package installation like this, adding an install step inside the Notebook document itself might be feasible. But in general the installation of tooling should be done on beforehand. Not only is it more convinient, it also is also helps reproducing.

So we'll add the `cowsay` package as a simple example for customization. In order to customize the `minimal-notebook` Container Image, we need a so-called `Dockerfile`.

> <hands-on-title></hands-on-title>
> 
> Create a new directory as your playground. Then, create a file called `Dockerfile` in this new directory.
> 
> Fill the Dockerfile with this content:
> 
> ```dockerfile
> # First, we say from which container image we want to start
> # This is called the base image
> FROM quay.io/jupyter/minimal-notebook:2025-06-23
> 
> # Now, we can run arbitrary commands in order to extend the base image
> RUN pip install --no-cache-dir cowsay
> ```
{: .hands_on}

> <comment-title>Let's tear this command apart</comment-title>
> 
> |---|---|
> | FROM | This tells docker from which image to start building |
> | RUN | Runs an arbitrary command on top of the image and thus creates a new image layer |
{: .comment}

> <comment-title>Let's tear this command apart</comment-title>
> 
> The complete Dockerfile reference can be found at [https://docs.docker.com/reference/dockerfile/](https://docs.docker.com/reference/dockerfile/).
{: .comment}

> <hands-on-title></hands-on-title>
> 
> In order to build the customized image, run the following command. Please mind the period dot at the end of the command.
> 
> ```bash
> sudo docker build -t my-custom-jupyterlab .
> ```
{: .hands_on}

> <comment-title>docker build</comment-title>
> 
> The `-t` option tells docker to tag the built image with provided tag (`my-custom-jupyterlab`). The `docker build` command needs a build context as the last parameter. We tell docker to use the current directory (`.`) as build context.
{: .comment}



# Use your first custom Image

You may use the newly build image just the same way you have used the `minimal-notebook` image beforehand, just use the tag you have provided with your build command:

```bash
sudo docker run --rm -p 127.0.0.1:8888:8888 my-custom-jupyterlab
```

Again, you'll provided with a link to open you jupyterlab. Start a Pyhton 3 Notebook by clicking the tile.

Again, run:

```python
# Hint: press shift + enter to execute code in Jupyter Notebook
import datetime
from cowsay import cow

now = datetime.datetime.now()
cow(now.strftime("%A"))
```

As you will see, the command now succeeds and prints a ascii-cow with the current weekday.

# Build your own customized JupyterLab

You now know how to make basic customizations to your JupyterLab. Let's go ahead and try out some real usecase:

> <hands-on-title>Make this script succeed in your custom JupyterLab</hands-on-title>
> 
> Create a custom notebook from the `scipy-notebook` base image and install the `igv-browser` extension.
> 
> You find all relevant information about this extension at [https://github.com/igvteam/igv-notebook](https://github.com/igvteam/igv-notebook).
> 
> You have succeeded, if you can run following snippet in your JupyterLab and are provided with a genome browser with loaded track:
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
{: .hands_on}

> <solution-title>Possible solution</solution-title>
>
> Keep in mind, there are many ways to solve this exercise. The proposed solution is only one way to succeed.
> 
> 1. Create a new `Dockerfile`:
> ```dockerfile
> FROM quay.io/jupyter/scipy-notebook:2025-06-23
> RUN pip install --no-cache-dir igv-notebook
> ```
> 2. Build the new Container Image:
> ```bash
> sudo docker build -t my-igv-jupyterlab .
> ```
> 3. Run a container from your image
> ```bash
> sudo docker run --rm -p 127.0.0.1:8888:8888 my-igv-jupyterlab
> ```
> 4. Access the JupyterLab via URL provided by the container output
> 5. Open a Python 3 Jupyter Notebook, paste in the above snippet and run it
> 6. You should be provided with a genome browser with loaded track
{: .solution}

# Docker Compose

Even though running all these docker commands directly from the command line works fine, it has significant disadvantages:
- it's verbose and therefore error-prone
- it's difficult to reproduce
- it's difficult to share and version

A better approach is to use a container composition, e.g. with dockers `compose` plugin. All relevant config is written into a file called `compose.yml`, which is reviewable, sharable, reproducable and versionizable.

> <hands-on-title></hands-on-title>
> 
> Create a file called `compose.yml` inside the directory created above.
> 
> Fill the compose file with following content:
> 
> ```yaml
> services:
>   jupyterlab:
>     image: my-custom-jupyterlab
>     build:
>       context: .
>       dockerfile: Dockerfile
>     ports:
>       - 127.0.0.1:8888:8888
> ```
{: .hands_on}

> <comment-title>Docker Compose</comment-title>
> 
> A complete reference to the compose format can be found at [https://docs.docker.com/reference/compose-file/](https://docs.docker.com/reference/compose-file/)
{: .comment}

> <hands-on-title></hands-on-title>
> 
> Now, build the image by simply executing:
> 
> ```
> sudo docker compose build
> ```
> 
> An then run it by simply executing:
> 
> ```
> sudo docker compose up
> ```
{: .hands_on}

Using a container composition is the preferred way of using containers on your local machine. However, following tutorials will guide you without compose files. Feel free to practise composition and solve the exercises using container compositions.
