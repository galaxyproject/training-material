---
layout: tutorial_hands_on

title: "Galaxy Tool Management with Ephemeris"
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
  - git-gat
---

This tutorial will introduce you to one of Galaxy's associated projects - [Ephemeris](https://ephemeris.readthedocs.io/). Ephemeris is a small Python library and set of scripts for managing the bootstrapping of Galaxy plugins - tools, index data, and workflows. It aims to help automate, and limit the quantity of manual actions admins have to do in order to maintain a Galaxy instance.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="tool-management" %}

# Background

You are an administrator of a Galaxy server. A colleague has approached you with a request to run a specific [Galaxy workflow]({% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}) on their data. In order to enable this workflow for your users, you will have to:

- identify what tools are required for the workflow
- install these tools and their dependencies on your Galaxy instance.


# Requirements

To run this tutorial, you will need to [install Ephemeris](https://ephemeris.readthedocs.io/en/latest/installation.html). You would normally install it on your workstation, but during training courses we recommend to install it on the same virtual machine used for the Galaxy server.

> <tip-title>Installing Ephemeris in a Python virtual environment</tip-title>
>
> 1. Install `virtualenv` if it is not already available. On Ubuntu this can be done with `sudo apt install virtualenv`
> 2. Create a virtual environment just for ephemeris, activate it and install ephemeris inside it:
>    ```console
>    virtualenv -p python3 ~/ephemeris_venv
>    . ~/ephemeris_venv/bin/activate
>    pip install ephemeris
>    ```
{: .tip}


# Extracting Tools

A common request you will experience as an administrator is "I want to run this workflow". Since this is such a common workflow, Galaxy has a built in way to accomplish it. We can use Ephemeris to extract a list of tools from a Galaxy workflow document, and then use Ephemeris to install these tools and specific versions into your Galaxy.

However, Galaxy workflow files are complex JSON documents, and the process of mapping the tool IDs to a ToolShed repository and revision is not trivial. Workflow files contain tool IDs which look like:

```
toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.71
toolshed.g2.bx.psu.edu/repos/bgruening/trim_galore/trim_galore/0.4.3.1
toolshed.g2.bx.psu.edu/repos/iuc/multiqc/multiqc/1.6
toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.3.4.2
toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.1
toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.4.1
toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.1
```

In order to actually install these tools, we need to convert each tool ID into a ToolShed repository name and revision. For example, FastQC version 0.71 corresponds to revision `ff9530579d1f` in the ToolShed.

```yaml
- name: fastqc
  owner: devteam
  revisions:
  - ff9530579d1f
  tool_panel_section_label: Tools from workflows
  tool_shed_url: https://toolshed.g2.bx.psu.edu
```

Ephemeris can take care of this process. Let's practice this on a real workflow.

> <hands-on-title>Extracting a list of tools from a workflow</hands-on-title>
>
> 1. Download the mapping workflow:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > wget {{ site.url }}{% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 2. Use the Ephemeris [`workflow-to-tools`](https://ephemeris.readthedocs.io/en/latest/commands/workflow-to-tools.html) command to extract the tool list from this workflow into a file named `workflow_tools.yml`.
>
>    > <question-title></question-title>
>    > What did your command look like?
>    >
>    > > <solution-title></solution-title>
>    > > > <code-in-title>Bash</code-in-title>
>    > > > ```bash
>    > > > workflow-to-tools -w mapping.ga -o workflow_tools.yml -l Mapping
>    > > > ```
>    > > > {: data-cmd="true"}
>    > > {: .code-in}
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

{% snippet faqs/galaxy/preferences_admin_api_key.md admin=true %}

There are two ways to install tools, depending on how you specify the tools to install:

> <hands-on-title>Installing a single tool</hands-on-title>
>
> 1. Use the Ephemeris [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) command to install the tool `bwa`, owned by `devteam` into a section named `Mapping`
>
>    > <question-title></question-title>
>    > What did your command look like?
>    >
>    > > <solution-title></solution-title>
>    > > Use your Galaxy URL and API key in the example command below:
>    > >
>    > > ```bash
>    > > shed-tools install -g https://your-galaxy -a <api-key> --name bwa --owner devteam --section_label Mapping
>    > > ```
>    > > {: data-cmd="true"}
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

> <tip-title>Certificate issues</tip-title>
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

> <hands-on-title>Installing tools from a tool list</hands-on-title>
>
> 1. (optional) Watch the installation proceed by running `journalctl -f -u galaxy` in a separate remote shell.
>
> 2. Use the Ephemeris [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html) command to install all of the tools from the `workflow_tools.yml` file on your Galaxy.
>
>    > <question-title></question-title>
>    > What did your command look like?
>    >
>    > > <solution-title></solution-title>
>    > > Use your Galaxy URL and API key in the example command below:
>    > >
>    > > ```bash
>    > > shed-tools install -g https://your-galaxy -a <api-key> -t workflow_tools.yml
>    > > ```
>    > > {: data-cmd="true"}
>    > {: .solution}
>    {: .question}
>
> 4. Open your Galaxy's admin interface and check that the tools have been installed.
>
> 5. Using the UI import the workflow file that you used, [mapping.ga]({% link topics/sequence-analysis/tutorials/mapping/workflows/mapping.ga %}).
>    1. Right-click or Ctrl-click on the link above and copy the link address
>    2. On your Galaxy instance click on `Workflow`, then `Import`.  Paste the URL into the `Import Archived URL` field.
{: .hands_on}

Occasionally the tool installation may fail due to network issues; if it does, just re-run the `shed-tools` installation process until it succeeds. This is a known issue the developers are working on.

> <tip-title>Opening a split screen in byobu</tip-title>
> <kbd>Shift-F2</kbd>: Create a horizontal split
>
> <kbd>Shift-Left/Right/Up/Down</kbd>: Move focus among splits
>
> <kbd>Ctrl-F6</kbd>:  Close split in focus
>
> <kbd>Ctrl-D</kbd>:  (Linux, Mac users) Close split in focus
>
> There are more byobu commands described in this [gist](https://gist.github.com/devhero/7b9a7281db0ac4ba683f)
{: .tip}

> <tip-title>Can I install tools without restarting?</tip-title>
> Yes. The default tool config (`config/tool_conf.xml.sample`, copy to `config/tool_conf.xml` to modify) has an option, `monitor="true"` set in the root `<toolbox>` tag. This instructs Galaxy to watch the tool files referenced in that config and load or reload them as necessary. It will also add any tools you have added.
{: .tip}

> <tip-title>Can I install tools without a ToolShed?</tip-title>
> Yes. The `galaxy_local_tools` option for the `galaxyproject.galaxy` Ansible role can be used to install local tools, or you can manage them in another way that fits your workflow better. UseGalaxy.eu, for example, maintains a repository of tools that are not installed from the ToolShed to aid their local developers. This is deployed to the server using the `git` module, rather than the Galaxy Ansible role.
{: .tip}


# Tool Testing

Having the tools installed is a good first step, but your users will expect that they actually work as well. You can use Ephemeris to automatically test all of the installed tools.

> <hands-on-title>Test the installed tools</hands-on-title>
>
> 1. Use the Ephemeris [`shed-tools`](https://ephemeris.readthedocs.io/en/latest/commands/shed-tools.html#test) command to test the `bamtools_filter` tool on your Galaxy.
>
>    > <question-title></question-title>
>    > What did your command look like?
>    >
>    > > <solution-title></solution-title>
>    > > Use your Galaxy URL and API key in the example command below:
>    > >
>    > > ```bash
>    > > shed-tools test -g https://your-galaxy -a <api-key> --name bamtools_filter --owner devteam
>    > > ```
>    > > {: data-cmd="true"}
>    > {: .solution}
>    {: .question}
> 2. Shed-tools test outputs a file with details of `tool_test_output.json` with details of jobs that have run.  Have a look at this file.
{: .hands_on}

This can give you some more confidence that things are working correctly. Oftentimes, users provide workflows for biological domains that we are not familiar with, so knowing how we can test these tools is impossible for us as admins. Leveraging the built-in tool test cases can give you reassurance that things are functional before you inform your users of the new tools.

> <tip-title>Tool test reports</tip-title>
> The ephemeris `shed-tools test` command produces an output file `tool_test_output.json` with information about the test jobs that have run. Galaxyproject's [Planemo](https://github.com/galaxyproject/planemo) can be used to generate [test reports](https://planemo.readthedocs.io/en/latest/commands/test_reports.html?highlight=test_reports#test-reports-command) from tool_test_output.json in HTML and other formats.
{: .tip}


# Obtaining a Tool List

Sometimes a user might ask you to install all the tools they were previously using on another Galaxy instance. Ephemeris can produce a `tool_list.yaml` file for all the tools installed on a server.

> <hands-on-title>Obtain UseGalaxy.eu's tool list</hands-on-title>
>
> 1. Use the Ephemeris [`get-tool-list`](https://ephemeris.readthedocs.io/en/latest/commands/get-tool-list.html) command to obtain the full set of tools installed on UseGalaxy.eu
>
>    > <question-title></question-title>
>    > What did your command look like?
>    >
>    > > <solution-title></solution-title>
>    > > This command does not require authentication and can be used to obtain the tool list from any public Galaxy server:
>    > >
>    > > > <code-in-title>Bash</code-in-title>
>    > > > ```bash
>    > > > get-tool-list -g "https://usegalaxy.eu" -o "eu_tool_list.yaml"
>    > > > ```
>    > > > {: data-cmd="true"}
>    > > {: .code-in}
>    > {: .solution}
>    {: .question}
> 2. Inpect the first few lines of tool list: Run `head -n 20 eu_tool_list.yaml`.
{: .hands_on}

We will not install all the tools from the EU Galaxy server as that server likely has more tools than any other Galaxy instance, but it is useful as an example of how you can use Ephemeris to facilitate the mirroring of another Galaxy instance.

> <tip-title>Non-shed tools</tip-title>
> The output of `get-tool-list` only includes ToolShed tools, not local non-TS tools.
{: .tip}

> <tip-title>But how many tools is it really?</tip-title>
> If you've seen the [European Galaxy tools view](https://usegalaxy.eu/tools/view) (this is available on any galaxy! Just access `/tools/view`)
> you'll notice they report somewhere over 2700 tools, however the `get-tool-list` output lists significantly fewer tools.
>
> This is for a couple reasons:
> - That's the number of repositories that are installed, and some repositories include multiple tools (e.g. [`circos`](https://toolshed.g2.bx.psu.edu/view/iuc/circos/df7356989ac1) has quite a few)
> - This only lists ToolShed tools, while EU also includes a number of non-TS tools
> - There are a number of tools built into Galaxy (e.g. the Collection Operation tools)
{: .tip}

# Production Best Practices

The servers which are part of the `usegalaxy.*` network use Ephemeris extensively to manage their large tool sets.

Interestingly, UseGalaxy.eu and UseGalaxy.org.au have different approaches:
- EU maintains [YAML files roughly per domain](https://github.com/usegalaxy-eu/usegalaxy-eu-tools), it is not as clear of an ordering. They maintain a YAML file where humans add the tools which should be installed in a given category, and lock files are automatically generated from these with the latest revision if it is missing. They follow a cycle of updating the lock files with the latest available revisions of tools, and then installing from these lock files any missing revisions. They use a [Jenkins server](https://build.galaxyproject.eu/job/usegalaxy-eu/job/install-tools/) to automatically run tool installation weekly.
- AU maintains [a separate YAML file per tool panel section](https://github.com/usegalaxy-au/usegalaxy-au-tools/tree/master/usegalaxy.org.au) as a record of all installed tools on the server.  They accept tool requests as pull requests and a Jenkins server is notified when a pull request is merged so that the tool can be automatically installed and tested. Like EU, they run automatic updates of installed tools once a week.
- Together, `usegalaxy.*` are working on a collaborative approach at [galaxyproject/usegalaxy-tools](https://github.com/galaxyproject/usegalaxy-tools) but this is not consumption ready yet.

If running ephemeris directly is not your preference, there is an Ansible [role](https://github.com/galaxyproject/ansible-galaxy-tools) and a sample [playbook](https://github.com/afgane/galaxy-tools-playbook) that can help automate some tasks.

> <tip-title>What if environments are not working</tip-title>
> It sometimes happens in Galaxy, that one environment isn't working anymore. It mostly happens from the start when it does happen. You can remove the environment on disk, or use the "Manage Dependencies" interface, select the environment, and delete it. Then re-install the dependency through the same Manage Dependencies interface.
{: .tip}

> <tip-title>Can you install multiple tools simultaneously?</tip-title>
> Previous experience with this is not good, there was a lot of unsafe code that would do things simultaneously that would destroy config files. With conda this is not recommended either, as this can corrupt conda environments. A solution for this is how [UseGalaxy.eu does it](https://github.com/usegalaxy-eu/usegalaxy-eu-tools), where they keep the full list of tools they want to install, and then a CI server (Jenkins) installs them. In this way they can enforce that only a single install process is running at any time.
{: .tip}

> <tip-title>Certificate Issues (GAT Only)</tip-title>
> If a student is running `shed-tools` on the VM, then it should work without certificate issues, because we installed the Fake LE X1 CA certificate, meaning that to your VM, the certificate chain is valid.
> We cannot (and would not recommend) setting that certificate on your local machine. That is the first way that comes to mind, that running an ephemeris command could generate that error.
{: .tip}

> <tip-title>Uninstalling tools</tip-title>
> While there is a function to accomplish this in [BioBlend](https://github.com/galaxyproject/bioblend), it has not been [included in Ephemeris yet.](https://github.com/galaxyproject/ephemeris/issues/83) If you're looking for a way to contribute to Galaxy, this would be great :)
{: .tip}

> <tip-title>"Recommended way" to install tools</tip-title>
> There is no one recommended way. Different people like to do different things. Some of us install things through the GUI still, some of us use ephemeris just for the automation.
>
> EU and others try and force all tool installation to go through ephemeris to make sure that rebuilding a server in case of disaster is easy. Anything not managed through their files will be lost.
>
> Others do a more mixed approach:
>
> > typically use Ephemeris, when I have to re-deploy a Galaxy instance (and thus install all tools) or during a Galaxy maintenance, when I want to update all installed toolshed tools in one go. Otherwise, when I need to install a specific version of one tool I use the GUI.
> {: .quote}
>
> There is no one right answer.
{: .tip}

> <tip-title>Tools not showing up? Refresh the toolbox!</tip-title>
> Sometimes the toolbox will fail to reload. You can correct for this by manually triggering the toolbox reload with a query:
>
> ```console
> curl -X PUT https://<your-galaxy>/api/configuration/toolbox -H "x-api-key: $GALAXY_API_KEY"
> ```
>
> This will request the toolbox to reload and you can check after if it's discovered your newly installed tools.
{: .tip}
