---
layout: tutorial_hands_on

title: "Galaxy Tool Management"
zenodo_link: ""
questions:
  - How to install, update, and maintain Galaxy tools?
  - How to extract a list of tools from a workflow or Galaxy instance?
objectives:
  - Learn about Ephemeris
  - Extract a list of tools from a workflow
  - Install these tools into a given Galaxy
time_estimation: "45m"
key_points:
  - Ephemeris and automation help with the tool management on Galaxy
  - There are tool management best practices you can learn from
  - Do not manage your Galaxy's tools manually
contributors:
  - martenson
  - erasche
---

# Overview
{:.no_toc}

This tutorial will introduce you to one of the Galaxy's projects - the [Ephemeris](https://github.com/galaxyproject/ephemeris). It aims to limit the manual interactions admins have to do when maintaining a Galaxy instance.

This is a hands-on training that will walk you through a real-life-ish user story.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# The situation

You are an administrator of your lab's Galaxy. A colleague approached you with a request to run a specific [Galaxy workflow](https://raw.githubusercontent.com/galaxyproject/training-material/master/topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga)  on the lab's data. You want to find all tools used in the workflow, review and install the missing ones on your instance.

# What is Ephemeris?

Ephemeris is a small Python library and set of scripts for managing the bootstrapping of Galaxy plugins - tools, index data, and workflows.

# Hands on

> ### {% icon hands_on %} Hands-on: Setting up our workspace
> For Python library a virtualenv is the preferred way to set up a single-purpose environment.
>
> 1. Install virtualenv `$ pip install virtualenv`.
>
> 2. Create an empty directory and `cd` into it.
>
> 3. Create virtualenv with `$ virtualenv .v_ephemeris`.
>
> 4. Activate virtualenv with `$ source .v_ephemeris/bin/activate`.
>
> 5. Install ephemeris `$ pip install ephemeris`


> ### {% icon hands_on %} Hands-on: Extracting tool list from a workflow
>
> 1. Download an exported Galaxy workflow using `$ wget <url>` command
>
> 2. Use proper Ephemeris command to extract tool list from this workflow into a file named `workflow_tools.yml`. Consult the [docs](https://ephemeris.readthedocs.io)
>
>    > ### {% icon solution %} Solution
>    > Example command:
>    > ```console
>    > $ workflow-to-tools -w mapping.ga -o workflow_tools.yml -l "mapping tools"
>    > ```
>    {: .solution }
>
> 3. Inspect the tool list file.


> ### {% icon hands_on %} Hands-on: Installing tools from a tool list
>
> 1. Identify url and port your Galaxy is running on.
>
> 2. (optional) Use `tail` command to watch the installation log.
>
> 3. Use proper Ephemeris command to install all tools from the `workflow_tools.yml` file.
>
>    > ### {% icon solution %} Solution
>    > Example command:
>    > ```console
>    > $ shed-tools install -t workflow_tools.yml -g "http://127.0.0.1:8080" -a "5cfd0d5f88c8addd5700b6a522a6a983"
>    > ```
>    {: .solution }
>
> 4. Load the Galaxy interface and check that the tools have been loaded.
>
> 5. Using the UI import the workflow from file `mapping.ga`.

> ### {% icon hands_on %} Hands-on: Test the installed tools
>
> 1. Use proper Ephemeris command to *test* all tools from the `workflow_tools.yml` file.
>
>    > ### {% icon solution %} Solution
>    > Example command:
>    > ```console
>    > $ shed-tools install -t workflow_tools.yml -g "http://127.0.0.1:8080" -a "5cfd0d5f88c8addd5700b6a522a6a983" --test
>    > ```
>    {: .solution }

> ### {% icon hands_on %} Hands-on: Get Galaxy's full tool list
>
> 1. Use proper Ephemeris command to obtain a tool list of all tools installed into your Galaxy.
>
>    {: .solution }>
>    > ### {% icon solution %} Solution
>    > Example command:
>    > ```console
>    > $ get-tool-list -g "http://127.0.0.1:8080" -o "galaxy_tool_list.yml"
>    > ```
>    {: .solution }

# Production (best) practices

Servers in `usegalaxy.*` network use Ephemeris extensively to manage its large tool sets.

You can explore the usegalaxy.eu [approach](https://github.com/usegalaxy-eu/usegalaxy-eu-tools) and the usegalaxy.org.au [way](https://github.com/usegalaxy-au/usegalaxy-au-tools/tree/current) of doing things.

There is also an existing Ansible [role](https://github.com/galaxyproject/ansible-galaxy-tools) and a sample [playbook](https://github.com/afgane/galaxy-tools-playbook) that can help automate some tasks even more.

