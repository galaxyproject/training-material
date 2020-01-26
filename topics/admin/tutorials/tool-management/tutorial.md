---
layout: tutorial_hands_on

title: "Ephemeris for Galaxy Tool Management"
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
  - hexylena
subtopic: features
tags:
  - tools
---

# Overview
{:.no_toc}

This tutorial will introduce you to one of Galaxy's associated projects - [Ephemeris](https://github.com/galaxyproject/ephemeris). Ephemeris is a small Python library and set of scripts for managing the bootstrapping of Galaxy plugins - tools, index data, and workflows. It aims to help automate, and limit the quantity of manual actions admins must do in order to maintain a Galaxy instance.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background

You are an administrator of your lab's Galaxy. A colleague has approached you with a request to run a specific [Galaxy workflow]({% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}) on the lab's data. In order to run this workflow you have to accomplish several substeps first. You will need to:

- identify what tools are required for the workflow
- and install these on your Galaxy instance


# Requirements

To run this tutorial, first you will need to [install Ephemeris](https://ephemeris.readthedocs.io/en/latest/installation.html). We recommend installing it in a virtualenv.


# Extracting Tools

Galaxy workflow files are complex json documents, and the process of mapping the tool IDs to a ToolShed repository and revision is not trivial. Workflow files contain tool IDs which look like:

```
toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.71
toolshed.g2.bx.psu.edu/repos/bgruening/trim_galore/trim_galore/0.4.3.1
toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.6
toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.3.4.2
toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.1
toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.4.1
toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.1
```

However in order to actually install these tools, we need to convert the ID into a `revision` number. Ephemeris takes care of this process for us, and identifies the exact revision which corresponds to the version listed in your workflow. For example, FastQC version 0.71 corresponds to revision `ff9530579d1f` in the ToolShed.

```yaml
- name: fastqc
  owner: devteam
  revisions:
  - ff9530579d1f
  tool_panel_section_label: Tools from workflows
  tool_shed_url: https://toolshed.g2.bx.psu.edu
```

Let's try running this on a real worfklow.

> ### {% icon hands_on %} Hands-on: Extracting a list of tools from a workflow
>
> 1. Download the mapping workflow:
>    ```console
>    $ wget {{ site.url }}{% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}
>    ```
>
> 2. Use the Ephemeris command [`workflow-to-tools`](https://ephemeris.readthedocs.io/en/latest/commands/workflow-to-tools.html) to extract the tool list from this workflow, into a file named `workflow_tools.yml`.
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > ```console
>    > > $ workflow-to-tools -w mapping.ga -o workflow_tools.yml -l "mapping tools"
>    > > ```
>    > {: .solution }
>    {: .question}
>
> 3. Inspect the tool list file.
>
{: .hands_on }

# Installing Tools

Now that you have extracted a list of tools, let's install these to Galaxy. In order to accomplish this section you will need:

- The URL of your Galaxy instance
- To be an admin of this Galaxy
- The API key for your account

There are additionally two ways to install tools, we'll show both:

> ### {% icon hands_on %} Hands-on: Installing tools from a specific tool name
>
> 1. Use the Ephemeris command [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) to install the tool `bwa`, owned by user `devteam` into a section named `Mapping`
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > Note that your API key and URL will probably be different than in the example command below:
>    > >
>    > > ```console
>    > > $ shed-tools install -g https://your-galaxy -a <api-key> --name bwa --owner devteam --section_label Mapping
>    > > ```
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

This provides an easy way to do one-off installation of tools, but is less convenient if you want to install many tools. For that, you can install from a yaml file:

> ### {% icon hands_on %} Hands-on: Installing tools from a tool list
>
> 1. (optional) Use the `tail` command to watch the installation proceed
>
> 2. Use the Ephemeris command [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) to install all of the tools from the `workflow_tools.yml` file to your Galaxy.
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > Note that your API key and URL will probably be different than in the example command below:
>    > >
>    > > ```console
>    > > $ shed-tools install -t workflow_tools.yml -g "http://127.0.0.1:8080" -a "5cfd0d5f88c8addd5700b6a522a6a983"
>    > > ```
>    > {: .solution}
>    {: .question}
>
> 4. Open your Galaxy's admin interface and check that the tools have been installed.
>
> 5. Using the UI import the workflow file that you used, `mapping.ga`.
{: .hands_on}

Occasionally this will fail due to network issues, if it does just re-run the `shed-tools` installation process until it succeeds. This is a known issue the developers are working on.


# Tool Testing

Having the tools installed is a good first step, but usually your users will expect that they actually work as well. You can use Ephemeris to automatically test all of the installed tools

> ### {% icon hands_on %} Hands-on: Test the installed tools
>
> 1. Use the Ephemeris command [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) to test all of the tools from the `workflow_tools.yml` file on your Galaxy.
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > Note that your API key and URL will probably be different than in the example command below:
>    > >
>    > > ```console
>    > > $ shed-tools install -t workflow_tools.yml -g "http://127.0.0.1:8080" -a "5cfd0d5f88c8addd5700b6a522a6a983" --test
>    > > ```
>    > {: .solution}
>    {: .question}
{: .hands_on}

This can give you some more confidence that things are working correctly. Oftentimes, users provide workflows for biological domains that we are not familiar with, so knowing how we can test these tools is impossible for us as admins. Leveraging the built-in tool test cases can give you reassurance that things are functional before you inform your users of the new tools.

# Obtaining a Tool List

Sometimes your users might request that you install all of the same tools as they were previously using in a domain specific server. Ephemeris offers functionality to obtain a `tool_list.yaml` for all of the tools installed on an instance.

> ### {% icon hands_on %} Hands-on: Obtain UseGalaxy.eu's tool list
>
> 1. Use the Ephemeris command [`get-tool-list`](https://ephemeris.readthedocs.io/en/latest/commands/get-tool-list.html) to obtain the full set of tools installed to UseGalaxy.eu
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > This command does not require authentication and can be used to obtain the tool list from any public Galaxy server:
>    > >
>    > > ```console
>    > > $ get-tool-list -g "https://usegalaxy.eu" -o "eu_tool_list.yaml"
>    > > ```
>    > {: .solution}
>    {: .question}
{: .hands_on}

We will not install the tools from that server as the EU Galaxy server has more tools than most other known Galaxies, but it is useful as an example of how you can use Ephemeris to help you mirror another Galaxy instance to meet user needs.

# Production Best Practices

Servers in `usegalaxy.*` network use Ephemeris extensively to manage their large tool sets.

UseGalaxy.eu and UseGalaxy.org.au have different approaches:
- AU maintains [separate yaml files](https://github.com/usegalaxy-au/usegalaxy-au-tools/tree/current) per toolbox category, which allows easily identifying where tools should be added or found to add new revisions. They follow a cycle of adding the tool to their yaml file, triggering installation, and then updating the yaml files from the current server status.
- EU maintains [yaml files roughly per domain](https://github.com/usegalaxy-eu/usegalaxy-eu-tools), it is not as clear of an ordering. They maintain a yaml file where humans add the tools which should be installed in a given category, and lock files are automatically generated from these with the latest revision if it is missing. They follow a cycle of updating the lock files with the latest available revisions of tools, and then installing from these lock files any missing revisions. They use a [Jenkins server](https://build.galaxyproject.eu/job/usegalaxy-eu/job/install-tools/) to automatically run tool installation weekly.
- Together, `usegalaxy.*` are working on a collaborative approach at [galaxyproject/usegalaxy-tools](https://github.com/galaxyproject/usegalaxy-tools) but this is not consumption ready yet.

If running ephemeris at the command line is not your preference, there is also an Ansible [role](https://github.com/galaxyproject/ansible-galaxy-tools) and a sample [playbook](https://github.com/afgane/galaxy-tools-playbook) that can help automate some tasks.
