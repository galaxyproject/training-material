---
layout: tutorial_hands_on

title: "Ephemeris for Galaxy Tool Management"
questions:
  - How to install, update, and maintain Galaxy tools?
  - How to extract a list of tools from a workflow or Galaxy instance?
objectives:
  - Learn about Ephemeris
  - Extract a list of tools from a workflow
  - Install these tools on a given Galaxy
time_estimation: "45m"
key_points:
  - Ephemeris and automation help with the tool management on Galaxy
  - There are tool management best practices you can learn from
  - Do not manage your Galaxy tools manually
contributors:
  - martenson
  - hexylena
  - nsoranzo
subtopic: features
tags:
  - tools
---

# Overview
{:.no_toc}

This tutorial will introduce you to one of Galaxy's associated projects - [Ephemeris](https://ephemeris.readthedocs.io/). Ephemeris is a small Python library and set of scripts for managing the bootstrapping of Galaxy plugins - tools, index data, and workflows. It aims to help automate, and limit the quantity of manual actions admins have to do in order to maintain a Galaxy instance.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background

You are an administrator of a Galaxy server. A colleague has approached you with a request to run a specific [Galaxy workflow]({% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}) on their data. In order to enable this workflow for your users, you will have to:

- identify what tools are required for the workflow
- install these tools and their dependencies on your Galaxy instance.


# Requirements

To run this tutorial, you will need to [install Ephemeris](https://ephemeris.readthedocs.io/en/latest/installation.html). You would normally install it on your workstation, but during training courses we recommend to install it on the same virtual machine used for the Galaxy server.

> ### {% icon tip %} Installing Ephemeris in a Python virtual environment
>
> 1. Install `virtualenv` if it is not already available. On Ubuntu this can be done with `sudo apt install virtualenv`
> 2. Create a virtual environment just for ephemeris, activate it and install ephemeris inside it:
>    ```console
>    virtualenv -p python3 ephemeris_venv
>    . ephemeris_venv/bin/activate
>    pip install ephemeris
>    ```
{: .tip}


# Extracting Tools

Galaxy workflow files are complex JSON documents, and the process of mapping the tool IDs to a ToolShed repository and revision is not trivial. Workflow files contain tool IDs which look like:

```
toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.71
toolshed.g2.bx.psu.edu/repos/bgruening/trim_galore/trim_galore/0.4.3.1
toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.6
toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.3.4.2
toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.1
toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.4.1
toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.1
```

However, in order to actually install these tools, we need to convert each tool ID into a ToolShed repository name and revision. For example, FastQC version 0.71 corresponds to revision `ff9530579d1f` in the ToolShed.

```yaml
- name: fastqc
  owner: devteam
  revisions:
  - ff9530579d1f
  tool_panel_section_label: Tools from workflows
  tool_shed_url: https://toolshed.g2.bx.psu.edu
```

Ephemeris can take care of this process. Let's practice this on a real worfklow.

> ### {% icon hands_on %} Hands-on: Extracting a list of tools from a workflow
>
> 1. Download the mapping workflow:
>    ```console
>    $ wget {{ site.url }}{% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}
>    ```
>
> 2. Use the Ephemeris [`workflow-to-tools`](https://ephemeris.readthedocs.io/en/latest/commands/workflow-to-tools.html) command to extract the tool list from this workflow into a file named `workflow_tools.yml`.
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
> 3. Inspect the `workflow_tools.yml` file, which contains a tool list in YAML format.
>
{: .hands_on }

# Installing Tools

Now that you have extracted a list of tools, let's install these on your Galaxy instance. In order to accomplish this, you will need:

- The URL of your Galaxy server
- The API key for your account, which must be an admin

> ### {% icon tip %} Get the API key of an admin account
>
> Galaxy admin accounts are specified as a comma-separated email list in the `admin_users` directive of `galaxy.yml` . If you have set up your Galaxy server using the [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, this is set to `admin@example.org` .
> 1. In your browser, open your Galaxy homepage
> 2. Log in using the admin email, or register a new account with it if it is the first time you use it
> 3. Go to `User -> Preferences` in the top menu bar, then click on `Manage API key`
> 4. If there is no current API key available, click on `Create a new key` to generate it
> 5. Copy your API key to somewhere convenient, you will need it throughout this tutorial
{: .tip}

There are two ways to install tools, depending on how you specify the tools to install:

> ### {% icon hands_on %} Hands-on: Installing a single tool
>
> 1. Use the Ephemeris [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) command to install the tool `bwa`, owned by user `devteam` into a section named `Mapping`
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > Use your Galaxy URL and API key in the example command below:
>    > >
>    > > ```console
>    > > $ shed-tools install -g https://your-galaxy -a <api-key> --name bwa --owner devteam --section_label Mapping
>    > > ```
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

> ### {% icon tip %} Certificate issues
>
> If your Galaxy instance is served via the HTTPS protocol (as it should be!), ephemeris will use the [requests](https://requests.readthedocs.io) Python library to encrypt the communication with Galaxy. Therefore, if your Galaxy uses a self-signed SSL certificate, `shed-tools` may fail with a `CERTIFICATE_VERIFY_FAILED` error.
>
> Under Ubuntu, you can allow the use of the unrecognized certificate as follows:
> 1. Get hold of the Certificate Authority (CA) certificate used to sign your Galaxy SSL certificate. For a [Galaxy Admin Training](https://github.com/galaxyproject/dagobah-training) course, this is usually the [Fake LE Root X1 certificate](https://letsencrypt.org/certs/fakelerootx1.pem).
> 2. Copy the CA certificate file into `/usr/local/share/ca-certificates/` with a `.crt` extension.
> 3. Run `update-ca-certificates` as root.
> 4. Execute `export REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt`, as explained in [requests docs](https://requests.readthedocs.io/en/master/user/advanced/#ssl-cert-verification).
>
> Now you should be able to execute successfully the `shed-tools` commands.
{: .tip}


This provides an easy way to do a one-off installation of a tool, but is not very convenient if you want to install many tools.
For that, you can install from a YAML file:

> ### {% icon hands_on %} Hands-on: Installing tools from a tool list
>
> 1. (optional) Watch the installation proceed by running `journalctl -f -u galaxy` in a separate remote shell.
>
> 2. Use the Ephemeris [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) command to install all of the tools from the `workflow_tools.yml` file on your Galaxy.
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > Use your Galaxy URL and API key in the example command below:
>    > >
>    > > ```console
>    > > $ shed-tools install -g https://your-galaxy -a <api-key> -t workflow_tools.yml
>    > > ```
>    > {: .solution}
>    {: .question}
>
> 4. Open your Galaxy's admin interface and check that the tools have been installed.
>
> 5. Using the UI import the workflow file that you used, [mapping.ga]({% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}) .
{: .hands_on}

Occasionally the tool installation may fail due to network issues; if it does, just re-run the `shed-tools` installation process until it succeeds. This is a known issue the developers are working on.


# Tool Testing

Having the tools installed is a good first step, but your users will expect that they actually work as well. You can use Ephemeris to automatically test all of the installed tools.

> ### {% icon hands_on %} Hands-on: Test the installed tools
>
> 1. Use the Ephemeris [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) command to test the `bamtools_filter` tool on your Galaxy.
>
>    > ### {% icon question %} Question
>    > What did your command look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > Use your Galaxy URL and API key in the example command below:
>    > >
>    > > ```console
>    > > $ shed-tools test -g https://your-galaxy -a <api-key> --name bamtools_filter --owner devteam
>    > > ```
>    > {: .solution}
>    {: .question}
{: .hands_on}

This can give you some more confidence that things are working correctly. Oftentimes, users provide workflows for biological domains that we are not familiar with, so knowing how we can test these tools is impossible for us as admins. Leveraging the built-in tool test cases can give you reassurance that things are functional before you inform your users of the new tools.

# Obtaining a Tool List

Sometimes a user might ask you to install all the tools they were previously using on another Galaxy instance. Ephemeris can produce a `tool_list.yaml` file for all the tools installed on a server.

> ### {% icon hands_on %} Hands-on: Obtain UseGalaxy.eu's tool list
>
> 1. Use the Ephemeris [`get-tool-list`](https://ephemeris.readthedocs.io/en/latest/commands/get-tool-list.html) command to obtain the full set of tools installed on UseGalaxy.eu
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

We will not install all the tools from the EU Galaxy server as that server likely has more tools than any other Galaxy instance, but it is useful as an example of how you can use Ephemeris to facilitate the mirroring of another Galaxy instance.

# Production Best Practices

The servers which are part of the `usegalaxy.*` network use Ephemeris extensively to manage their large tool sets.

Interestingly, UseGalaxy.eu and UseGalaxy.org.au have different approaches:
- AU maintains [a separate YAML file](https://github.com/usegalaxy-au/usegalaxy-au-tools/tree/current) per each toolbox category, which allows to easily identify where tools should be added or found to add new revisions. They follow a cycle of adding the tool to their YAML file, triggering installation, and then updating the YAML files from the current server status.
- EU maintains [YAML files roughly per domain](https://github.com/usegalaxy-eu/usegalaxy-eu-tools), it is not as clear of an ordering. They maintain a YAML file where humans add the tools which should be installed in a given category, and lock files are automatically generated from these with the latest revision if it is missing. They follow a cycle of updating the lock files with the latest available revisions of tools, and then installing from these lock files any missing revisions. They use a [Jenkins server](https://build.galaxyproject.eu/job/usegalaxy-eu/job/install-tools/) to automatically run tool installation weekly.
- Together, `usegalaxy.*` are working on a collaborative approach at [galaxyproject/usegalaxy-tools](https://github.com/galaxyproject/usegalaxy-tools) but this is not consumption ready yet.

If running ephemeris directly is not your preference, there is also an Ansible [role](https://github.com/galaxyproject/ansible-galaxy-tools) and a sample [playbook](https://github.com/afgane/galaxy-tools-playbook) that can help automate some tasks.
