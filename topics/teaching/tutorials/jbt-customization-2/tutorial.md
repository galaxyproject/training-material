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

Running JupyterLab in a container, means that your data may be lost when the container is removed. To avoid data loss, volumes can be used.

--> On local machine, mount it wherever you want

--> Bind mount vs actual volume

--> With JupyterHub it's common to overlay homedir

## Large Data

* Don't include it into your Image!

* S3

* WIP

## Sensitive Data (?)

* Just avoid at all cost

* Maybe WIP, but maybe put that into `jbt-intro`

# Application Proxies

Install the pip `jupyter_server_proxy`

Download code-server from Github https://github.com/coder/code-server/releases

```
wget -qO- 'https://github.com/coder/code-server/releases/download/v4.101.2/code-server-4.101.2-linux-amd64.tar.gz' | tar xzvf - -C code-server --strip-components 1
```

Configure the proxy server application in `/etc/jupyter/jupyter_lab_config.py`:

```python
c.ServerProxy.servers = {
    "code-server": {
        "command": [
          "/opt/code-server/bin/code-server",
          "--auth=none",
          "--socket={unix_socket}",
          "--disable-telemetry",
          "--disable-update-check"
        ],
        "unix_socket": True,
        "timeout": 30,
        "absolute_url": False,
        "raw_socket_proxy": False,
        "launcher_entry": {
          "enabled": True,
          "title": "Code",
          "icon_path": "/opt/code-server/src/browser/media/favicon.svg"
        }
    }
}
```

Extend your Dockerfile:

```dockerfile
FROM quay.io/jupyter/minimal-notebook:2025-06-23

USER root
RUN apt update \
    && apt install -y build-essential \
    && apt clean

USER ${NB_USER}
ADD --chmod=755 https://sh.rustup.rs /tmp/rustup.sh
RUN /tmp/rustup.sh -v -y
RUN . "$HOME/.cargo/env" && cargo install --locked evcxr_jupyter
RUN . "$HOME/.cargo/env" && evcxr_jupyter --install

USER root
ARG CS_URL=https://github.com/coder/code-server/releases/download/v4.101.2/code-server-4.101.2-linux-amd64.tar.gz
ARG CS_PATH=/opt/code-server
RUN mkdir "${CS_PATH}"
RUN wget -qO- "${CS_URL}" | tar xzvf - -C "${CS_PATH}" --strip-components 1
COPY --chown=root:root jupyter_lab_config.py /etc/jupyter/jupyter_lab_config.py

USER ${NB_UID}
RUN pip install jupyter_server_proxy
```

You may install any other arbitrary server applications this way too, e.g. RStudio.

# Exercise 1: Persistent Data

Why is the rust installation in step 1 of this tutorial problematic in case this notebook is served via JupyterHub?

--> it installs rust in homedir


HINT:

This method installs all the kernel related resources inside the default homedir.

Therefore this will not work inside deployments where an empty volume is mounted as homedir, e.g. in default z2jh JupyterHub Deployments.


Build a Customized JupyterLab with an available rust kernel even in JupyterHub szenario.

```dockerfile
FROM quay.io/jupyter/minimal-notebook:2025-06-23

USER root
RUN apt update \
    && apt install -y build-essential \
    && apt clean

ARG RUSTUP_URL=https://sh.rustup.rs
ARG RUSTUP_INIT=/tmp/rustup.sh

ENV RUSTUP_HOME=/opt/rustup
ENV CARGO_HOME=/opt/cargo
ENV JUPYTER_PATH=/usr/local/share/jupyter
ENV PATH="${PATH}:${CARGO_HOME}/bin"

ADD --chmod=755 "${RUSTUP_URL}" "${RUSTUP_INIT}"
RUN "${RUSTUP_INIT}" -v -y
RUN cargo install --locked evcxr_jupyter
RUN evcxr_jupyter --install

RUN rm "${RUSTUP_INIT}"

USER ${NB_UID}
```

# Exercise 2: Proxy Application

Create a simple html page like this:

```html
<!DOCTYPE html>
<html>
<head>
  <title>Hello World</title>
  <meta charset="UTF-8" />
</head>
<body>
  <h1>Hello World!</h1>
  <p>I am served by your python server application!</p>
</body>
</html>
```

Include it as `index.html` in a separate folder in your container image.

Create a proxy application serving this index.html file with the simple python builtin `http.server` package.

Tip:

```bash
python3 -m http.server -d {path_to_dir_containing_the_html_file} {PORT}
```



Serve it via simple python 

Solution:

Create the index.html as explained above.

Create a JupyterLab config `jupyter_lab_config.py` like this:

```python
c.ServerProxy.servers = {
    "hello-world-server": {
        "command": [
          "python3",
          "-m", "http.server",
          "-d", "/srv/html",
          "{port}"
        ],
        "launcher_entry": {
          "enabled": True,
          "title": "Hello World Server"
        }
    }
}
```

Create a Dockerfile like this:

```dockerfile
FROM quay.io/jupyter/minimal-notebook:2025-06-23

USER root

COPY --chown=root:root jupyter_lab_config.py /etc/jupyter/jupyter_lab_config.py
COPY --chown=root:root index.html /srv/html/index.html

USER ${NB_UID}
RUN pip install jupyter_server_proxy
```
