---
layout: tutorial_hands_on

title: Advanced Jupyter customization
subtopic: practises
draft: true
time_estimation: 90m
questions:
  - How can I tailor Jupyter to my specific teaching needs?
objectives:
  - A guide to advanced Jupyter customization is provided, e.g. how to install proxy applications.
key_points:
  - Advanced customitzations may be time-consuming, but offer great potential and reusability value.
contributors:
  - mittler-works
---

# Introduction

Now we know how to include changes, we will build something more complex and look at some pitfalls.

# Install an additional Kernel

https://github.com/jupyter/jupyter/wiki/Jupyter-kernels

Rust: https://github.com/evcxr/evcxr/tree/main/evcxr_jupyter


```
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
```

Build It.

Run it.

Try it out:

```rust
println!("Hello World!");
```

HINT:

This method installs all the kernel related resources inside the default homedir.

Therefore this will not work inside deployments where an empty volume is mounted as homedir, e.g. in default z2jh JupyterHub Deployments.

# Config Changes

* Where to put them
You may put it in `~/.jupyter/jupyter_lab_config.py`

But you should put it in `/etc/jupyter/jupyter_lab_config.py` because homedir will be overriden in some infrastructures.


* What to put there
https://docs.jupyter.org/en/latest/use/config.html

You find all configuration options here: https://jupyter-server.readthedocs.io/en/latest/other/full-config.html

* How to put there
Create a new local file `jupyter_lab_config.py` and fill it with content.

Then, in your Dockerfile add it either as system wide config (/etc) or user config (~/)

```dockerfile
FROM quay.io/jupyter/minimal-notebook:2025-06-23

COPY --chown=root:root jupyter_lab_config.py /etc/jupyter/jupyter_lab_config.py

```

# Persistent Data

* HomeDir will be overlayed with a mounted Volume

* WIP

## External persistent Data

* Depends on your environment

* WIP

## Large Data

* Don't include it into your Image!

* S3

* WIP

## Sensitive Data (?)

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
