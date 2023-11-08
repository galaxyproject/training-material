---
layout: tutorial_hands_on

title: "Running the GTN website locally"
questions:
  - "How to setup the infrastructure to build training webpages?"
objectives:
  - "Installing packages needed for rendering the webpage"
  - "Running the GTN material website locally"
  - "Tracking changes to the content live in the webbrowser"
time_estimation: "15m"
subtopic: writing
key_points:
  - "Checking the generated website can be done locally"
contributors:
  - bebatut
  - bgruening
  - shiltemann
---

# Introduction


If you want to run the entire GTN material website locally or test your new training material you can do this!

Currently, the website is generated from the metadata and the tutorials using Jekyll, a simple static site builder.
We can use Jekyll to run a server to check if the tutorial is correctly added and rendered.

> <agenda-title></agenda-title>
>
> In this tutorial, you will learn how to run a local instance of the GTN website:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Installation of the requirements

The first step is to install the needed tools inside a conda environment. A conda environment is a directory that contains a specific collection of packages. For example here to run the website, we need ruby, pandas, requests, etc. By creating a conda environment and installing the needed tools there, we do not affect your main installation.

This step has to be done once.

We also need to make sure that a couple of other utilities and build requirements are present: **git**, **curl** & **make**. The easiest way to install these is with your package manager of choice - i.e. yum, apt, brew etc.

> <hands-on-title>Install the requirements</hands-on-title>
>
> 1. Open a Terminal
> 2. Use your package manager to install `git`, `curl` and `make`
>
>    - For Debian/Ubuntu: `sudo apt update && sudo apt install git curl make`
>    - For Fedora/CentOs/RedHat: `sudo yum install git curl make`
>
> 3. (If not done yet) Clone the training material GitHub repository: `git clone https://github.com/galaxyproject/training-material.git`
> 4. Navigate to the `training-material/` folder with `cd`
> 5. Set up the conda environment
>
>     It will install some needed tools (ruby, nodejs, etc) in a protected environment, without interfering with the existing tools or versions.
>
>     1. Install conda (if not already installed): `make install-conda`
>     2. (You may need to exit the terminal and re-open for conda to be recognised. Navigate back to the same place.)
>     3. Create the `galaxy_training_material` conda environment: `make create-env`
>
> 6. Install Jekyll and related modules into the conda environment: `make install`
{: .hands_on}

> <details-title>Troubleshooting `libxml2` errors</details-title>
> If you encounter an error about libxml2 on Linux, please try to install `libxml2-dev` (executing `sudo apt install libxml2-dev`) if on Debian/Ubuntu or `libxml2-devel` (executing `sudo yum install libxml2-devel`) if on Fedora/RedHat/CentOS, and re-run `make install` .
{: .details}


# Checking the website generation

Once Jekyll and its modules are installed in our conda environment, we can check the generation of the website locally:

> <hands-on-title>Checking the website generation locally</hands-on-title>
>
> 1. Run a local Jekyll server with `make serve-quick`
> 2. Visualize at [http://localhost:4000/training-material/ ](http://localhost:4000/training-material/)
> 3. Edit one of the tutorials:
>    - For example, open `topics/introduction/tutorials/galaxy-intro-peaks2genes/tutorial.md` in a text editor of your choice.
>    - Make some changes to the *Introduction* paragraph, and save the file.
>    - Refresh the tutorial page in your browser until you can see the changes you made.
>        - this may take a little bit of time; in the terminal you can monitor when the regeneration is complete
>
{: .hands_on}

With `make serve-quick`, a local Jekyll server will run in background. It will check the changes and regenerate the website accordingly. You may need to reload the page to see the changes (and sometimes to wait 1-2 minutes).


> <tip-title>Tips</tip-title>
>
> 1. Use `make serve` instead of `make serve-quick` to get all plugins, but also configure the post, host and pass additional flags. This however can be quite slow.
>
> 2. Need to speed up the cloning step? You coud fetch only the latest commit of the main branch:
>
>    ```
>    $ git clone {{ site.github_repository }}.git --depth 1 --branch main
>    ```
>
> 3. Running on a VM or remote machine?
>
>    If you are not running this on your local machine, but e.g. on a VM, you may need to configure a webserver to serve the website.
>
>    Below is an example NGINX configuration (e.g. in `/etc/nginx/sites-enabled/default`)
>
>    ```
>    location /training-material/ {
>      root /home/ubuntu/training-material/_site/;
>    }
>    ```
>
>    Change the `root` path above to wherever you cloned the training material folder
>
{: .tip}


# Stopping the server

Once you are done, you can stop the server and clean your repository.

> <hands-on-title>Stopping the server</hands-on-title>
>
> 1. Stop the server with <kbd>CTRL</kbd>-<kbd>C</kbd>
> 2. Clean the repository: `make clean`
>
{: .hands_on}

# Conclusion

