---
layout: tutorial_hands_on

title: Advanced JupyterLab customizations
subtopic: practises
draft: true
time_estimation: 90m
questions:
  - How can I tailor JupyterLab to my specific teaching needs?
objectives:
  - A guide to advanced JupyterLab customization is provided, e.g. how to install kernels and proxy applications.
key_points:
  - Advanced customitzations may be time-consuming, but offer many possibilities and reusability value.
contributors:
  - mittler-works
---

# Prerequisites

It is assumed that you either completed the previous tutorial "Basic JupyterLab customizations" or that you already have basic knowledge of containerization and JupyterLab customizing.

With this knowledge, we'll add some more advanced customizations to JupyterLab and look at some pitfalls.

# Installing an additional Kernel

One common szenario for customization is the need for an additional kernel, e.g. for Julia or R. This tutorial guides you through installing a rust kernel on top of `minimal-notebook` as an example.

> <comment-title>Jupyter Kernels</comment-title>
> 
> A list of availible Jupyter Kernels can be found at [https://github.com/jupyter/jupyter/wiki/Jupyter-kernels](https://github.com/jupyter/jupyter/wiki/Jupyter-kernels)
{: .comment}

If you look for `Rust` you'll find the `Evcxr` Jupyter Kernel which is the one we are going to add.

> <comment-title>Evcxr Kernel</comment-title>
> 
> All relevant information for the Evcxr Kernel can be found at [https://github.com/evcxr/evcxr/tree/main/evcxr_jupyter](https://github.com/evcxr/evcxr/tree/main/evcxr_jupyter)
{: .comment}

> <hands-on-title></hands-on-title>
> 
> Create a new directory. Inside, create a `Dockerfile` with the following content:
> 
> ```
> FROM quay.io/jupyter/minimal-notebook:2025-06-23
> 
> USER root
> RUN apt update \
>     && apt install -y build-essential \
>     && apt clean
> 
> USER ${NB_UID}
> ADD --chmod=755 https://sh.rustup.rs /tmp/rustup.sh
> RUN /tmp/rustup.sh -v -y
> 
> ENV PATH="${PATH}:${HOME}/.cargo/bin"
> RUN cargo install --locked evcxr_jupyter
> RUN evcxr_jupyter --install
> ```
{: .hands_on}

> <comment-title>Breakdown</comment-title>
> 
> Again, we are using `minimal-notebook` as base image. But this time we need to do more than just installing a python package:
> 1. we need to install system packages to meet requirements for Rust
> 2. we need to install Rust itself
> 3. we need to install and register the Evcxr kernel
> 
> In order to install system packages, we need a privileged user inside the docker build, which is why we are changing to the root user using the `USER` keyword. When installing packages, keep in mind that docker builds are non-interactive, thus you need to use the `-y` flag to bypass the confirmation prompt.
> 
> After installing system files, we can return to the JupyterLab default user, who is usually called `jovyan`. The variable `NB_UID` is inherited from the base image.
> 
> Next, the rustup utility will be downloaded and executed. Note, with the `--chmod` portion you can tell `ADD` in which filemod to save the downloaded file, thus you do not need to make the file executable in a separate step.
> The rustup script installs the binaries for rust in `${HOME}/.cargo/bin` which is why we need to add it to the path. We can do that simply by using the `ENV` keyword.
> 
> At least, we install and register the Evcxr as explained on the projects github page.
{: .question}

> <hands-on-title>Build your JupyterLab</hands-on-title>
> 
> ```bash
> sudo docker build -t my-rust-jupyterlab .
> ```
{: .hands_on}

> <hands-on-title>Run your JupyterLab</hands-on-title>
> 
> ```bash
> sudo docker run --rm -p 127.0.0.1:8888:8888 my-rust-jupyterlab .
> ```
{: .hands_on}

> <hands-on-title>Try your JupyterLab</hands-on-title>
> 
> Again, you'll be provided with a link to open your JupyterLab.
> 
> Open it, click on the Rust kernel and run following snippet:
> ```rust
> println!("Hello World!");
> ```
{: .hands_on}

Congratulations, you have successfully installed a custom kernel to your JupyterLab!

# Configuration Changes

There are a few files for configuring JupyterLab. One important one is the `jupyter_lab_config.py` file, where you can configure various different settings, e.g. the log level of the server application and much more.

> <comment-title></comment-title>
> 
> All relevant configuration options for `jupyter_lab_config.py` can be found at [https://jupyter-server.readthedocs.io/en/latest/other/full-config.html](https://jupyter-server.readthedocs.io/en/latest/other/full-config.html)
{: .comment}

You might include this file either as a user configuration in the user's home directory like this: `${HOME}/.jupyter/jupyter_lab_config.py`. But you better include it as a system-wide configuration by including it as `/etc/jupyter/jupyter_lab_config.py`. You can find out why this is the better solution in the [persistent data](#persistent-data) section.

> <hands-on-title>Include your custom configuration</hands-on-title>
> 
> First, create a snippet of custom config, e.g. let's edit the default name for Notebooks from "Untitled" to "HelloWorld". Save it as `jupyter_lab_config.py`:
>
> ```python
> c.ContentsManager.untitled_notebook = 'HelloWorld'
> ```
> 
> Then, add it as system wide config to your Dockerfile. Again, we're using the `minimal-notebook` as our base image:
> 
> ```dockerfile
> FROM quay.io/jupyter/minimal-notebook:2025-06-23
> 
> COPY --chown=root:root jupyter_lab_config.py /etc/jupyter/jupyter_lab_config.py
> ```
> 
> You may now build it as usual...
> 
> ```bash
> sudo docker build -t my-configured-jupyterlab .
> ```
>
> ...run it as usual...
> 
> ```bash
> sudo docker run --rm -p 127.0.0.1:8888:8888 my-configured-jupyterlab .
> ```
> 
> ...and access it with the printed URL. Now, open a new Notebook Document and you'll see it is now called `HelloWorld.ipynb` instead of `Untitled.ipynb`.
{: .hands_on}

Congratulations, you have successfully configured your JupyterLab with a custom configuration! While this example handles a very low impact customization, please review the available configuration options! You may also configure plugins that way, as you will learn later in the [Application Proxies](#application-proxies) section.

# Persistent Data

Even though persistent data is not strictly a topic of customization, it is a relevant topic that you have to keep in mind when customizing your JupyterLab.

## Volumes

Running JupyterLab in a container, means that your data may be lost when the container is removed. To avoid data loss, volumes can be used.

Volumes are only loosly coupled with containers at runtime and thus they are not tailored to a containers lifecycle. Volumes can be mounted into a container on an arbitrary path and multiple Volumes may by mounted.

Please note that data in the directory to which a volume is being mounted will be overlaid by the data in the volume in will wherefore not be easy accessible. That means, if you install user packages or configuration in your home directory when building the image and if you are mounting a volume as your home directory for persistence, all prebuild data in this home directory will not be available in the running container.

This is especially an issue if you are using a container-based JupyterHub Provider, as it is very likely that an empty volume will be mounted as users homedirectory by default. On your local machine you have the free choice on which directory a volume should be mounted. You could use `/mnt/data` or `/home/jovyan/data` in order not to override your complete home directory.

For your local machine, you have two basic options for using volumes with your containers: Named Volumes and Bind Mounts.

### Named Volumes

Named Volumes are managed by docker and will be automatically created on demand as soon as it is requested. The created Volume will be just empty. You may reference a named volume with it's name or id.

```bash
# This mounts the named volume my-data-volume to /data. If it does not exist, it will be created in first place.
sudo docker run --rm -v "my-data-volume:/data" -p 127.0.0.1:8888:8888 my-configured-jupyterlab
```

### Bind Mounts

Instead of having a volume managed by docker, you may also bind mount an existing directory into the container. This is useful when you have data to be shared on your local drive, e.g. for live coding or for analysis.

Please note, bind mounting data into the container may lead to problems with uid and gid numbers for the mounted directory, as the user inside docker not neccessarly has the same uid and gid numbers as your local machine user.

```bash
# This bind mounts the current directory to /data
sudo docker run --rm -v "$(pwd):/data" -p 127.0.0.1:8888:8888 my-configured-jupyterlab
```

## Large Data

Large Data must not be added to a container image, but should be mounted from a (shared) volume into the container.

Having big images wastes time and storage. It increases bulid, pull and push times and replicates the data on each host using this image.

On your local machine, you may use a bind mount for sharing. On a JupyterHub Provider your better use shared volumes or s3 in order to provide your JupyterLabs with data.

## Sensitive Data

As well as large data, sensitive data must not be added to a container image.

When it comes to sensitive data, it is not about resource consumption, but about data protection rules. It is very hard to use sensitive data inside a container image without violating GDPR, and especially in JupyterHub Provider szenarios it is hardly possible at all.

# Application Proxies

Sometimes it might be useful to start a server software with a web interface from within your JupyterLab. Great examples are RStudio or a Code Server. To access this web interface in a containerized setup, a proxy can and should be used. The proxy runs in your JupyterLab and hands all relevant traffic to your server application.

> <comment-title>Jupyter Server Proxy</comment-title>
> 
> All relevant info for Jupyter Server Proxy can be found at [https://github.com/jupyterhub/jupyter-server-proxy](https://github.com/jupyterhub/jupyter-server-proxy).
> 
{: .comment}

In this tutorial we'll install `code-server` and register it as a proxy application to be accessed directly from your JupyterLab's home screen. The roadmap for implementation is straight forward:
1. Install `code-server`
2. Install `jupyter_server_proxy`
3. Add Application Proxy to your JupyterLab configuration

> <comment-title>code server</comment-title>
> 
> All relevant info for code-server can be found at [https://github.com/coder/code-server](https://github.com/coder/code-server).
> 
{: .comment}

> <hands-on-title>Installing code server with application proxy</hands-on-title>
> 
> Create a new directory, inside create a `Dockerfile` with following content:
> 
> ```dockerfile
> FROM quay.io/jupyter/minimal-notebook:2025-06-23
> 
> USER root
> ARG CS_URL=https://github.com/coder/code-server/releases/download/v4.101.2/code-server-4.101.2-linux-amd64.tar.gz
> ARG CS_PATH=/opt/code-server
> RUN mkdir "${CS_PATH}"
> RUN wget -qO- "${CS_URL}" | tar xzvf - -C "${CS_PATH}" --strip-components 1
> COPY --chown=root:root jupyter_lab_config.py /etc/jupyter/jupyter_lab_config.py
> 
> USER ${NB_UID}
> RUN pip install jupyter_server_proxy
> ```
> 
> Next to the Dockerfile, create a configuration file `jupyter_lab_config.py` with following content:
> 
> ```python
> c.ServerProxy.servers = {
>     "code-server": {
>         "command": [
>           "/opt/code-server/bin/code-server",
>           "--auth=none",
>           "--socket={unix_socket}",
>           "--disable-telemetry",
>           "--disable-update-check"
>         ],
>         "unix_socket": True,
>         "timeout": 30,
>         "absolute_url": False,
>         "raw_socket_proxy": False,
>         "launcher_entry": {
>           "enabled": True,
>           "title": "Code",
>           "icon_path": "/opt/code-server/src/browser/media/favicon.svg"
>         }
>     }
> }
> ```
> 
> Now, build it...
> 
> ```bash
> sudo docker build -t my-code-jupyterlab .
> ```
>
> ...and run it...
> 
> ```bash
> sudo docker run --rm -p 127.0.0.1:8888:8888 my-code-jupyterlab .
> ```
> 
> ...and access it with the printed URL. You can now see a new Launcher "Code". Try it out and you'll be provided with a web-based IDE.
{: .hands_on}

Congratulations, you have successfully installed a `code-server` applications proxy into your JupyterLab. You may install any other arbitrary server applications this way too, e.g. Shiny Server or RStudio.

# Exercise 1: Persistent Data

> <question-title></question-title>
> 
> Why might be the rust installation above, in section 1 of this tutorial, problematic?
> 
> > <solution-title></solution-title>
> > 
> > With the installation above, all rust related tools will be installed and configured in the users homedir and thus might be overriden when using volumes. Try it out:
> > 
> > ```bash
> > sudo docker run --rm -v "jovyan-home:/home/jovyan" -p 127.0.0.1:8888:8888 my-rust-jupyterlab .
> > ```
> > 
> > Access the notebook and see that your kernel and all rust related packages are missing.
> {: .solution}
> 
{: .question}

> <question-title></question-title>
>
> Build a customized JupyterLab container that solves this problem.
> 
> > <solution-title></solution-title>
> > 
> > Keep in mind, that the proposed solution is only one possible way of solving the problem. You might find other solutions as well.
> > 
> > First, create a Dockerfile with following content:
> > 
> > ```dockerfile
> > FROM quay.io/jupyter/minimal-notebook:2025-06-23
> > 
> > USER root
> > RUN apt update \
> >     && apt install -y build-essential \
> >     && apt clean
> > 
> > ARG RUSTUP_URL=https://sh.rustup.rs
> > ARG RUSTUP_INIT=/tmp/rustup.sh
> > 
> > ENV RUSTUP_HOME=/opt/rustup
> > ENV CARGO_HOME=/opt/cargo
> > ENV JUPYTER_PATH=/usr/local/share/jupyter
> > ENV PATH="${PATH}:${CARGO_HOME}/bin"
> > 
> > ADD --chmod=755 "${RUSTUP_URL}" "${RUSTUP_INIT}"
> > RUN "${RUSTUP_INIT}" -v -y
> > RUN cargo install --locked evcxr_jupyter
> > RUN evcxr_jupyter --install
> > 
> > RUN rm "${RUSTUP_INIT}"
> > 
> > USER ${NB_UID}
> > ```
> >
> > Build it, run it with volume mounted on homedir...
> > 
> > ```bash
> > sudo docker build -t my-rust-jupyterlab .
> > sudo docker run --rm -v "jovyan-home:/home/jovyan" -p 127.0.0.1:8888:8888 my-rust-jupyterlab .
> > ```
> > 
> > Access it and verify that it's working.
> {: .solution}
> 
{: .question}

# Exercise 2: Proxy Application

In this exercise, you will configure your own Proxy Application in order to access a static web page served from within your JupyterLab.

> <hands-on-title></hands-on-title>
> 
> Create a simple html page like this:
> 
> ```html
> <!DOCTYPE html>
> <html>
> <head>
>   <title>Hello World</title>
>   <meta charset="UTF-8" />
> </head>
> <body>
>   <h1>Hello World!</h1>
>   <p>I am served by your python server application!</p>
> </body>
> </html>
> ```
{: .hands_on}

> <question-title></question-title>
> 
> Include the above html snippet as `index.html` in a custom folder in your container image.
> 
> Create a proxy application serving this index.html file with the simple python builtin `http.server` package.
> 
> > <comment-title>http.server</comment-title>
> > 
> > All relevant information about the http-server package can be found at [https://docs.python.org/3/library/http.server.html](https://docs.python.org/3/library/http.server.html).
> {: .comment}
> 
> > <tip-title></tip-title>
> > 
> > ```bash
> > python3 -m http.server -d {path_to_dir_containing_the_html_file} {PORT}
> > ```
> {: .tip}
> 
{: .question}

> <solution-title></solution-title>
> 
> Create an index.html file according to example above.
> 
> Create a JupyterLab config `jupyter_lab_config.py` like this:
> 
> ```python
> c.ServerProxy.servers = {
>     "hello-world-server": {
>         "command": [
>           "python3",
>           "-m", "http.server",
>           "-d", "/srv/html",
>           "{port}"
>         ],
>         "launcher_entry": {
>           "enabled": True,
>           "title": "Hello World Server"
>         }
>     }
> }
> ```
> 
> Create a Dockerfile like this:
> 
> ```dockerfile
> FROM quay.io/jupyter/minimal-notebook:2025-06-23
> 
> COPY --chown=root:root jupyter_lab_config.py /etc/jupyter/jupyter_lab_config.py
> COPY --chown=root:root index.html /srv/html/index.html
> 
> RUN pip install jupyter_server_proxy
> ```
> 
> Build and run it:
> 
> ```bash
> sudo docker build -t my-rust-jupyterlab .
> sudo docker run --rm -v "jovyan-home:/home/jovyan" -p 127.0.0.1:8888:8888 my-rust-jupyterlab .
> ```
> 
> Access it and verify that it's working.
{: .solution}
