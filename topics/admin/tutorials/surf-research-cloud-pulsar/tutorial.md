---
layout: tutorial_hands_on

title: "Pulsar usage on SURF Research Cloud"
zenodo_link: ""
questions:
  - How do I start a Pulsar instance on SURF Research Cloud?
  - How do I connect to Pulsar?
objectives:
  - Be able to attach a Pulsar node to Galaxy
  - Send jobs to Pulsar
key_points:
  - With SRC you can start your own Pulsar on-demand instance in a secure environment.
  - The Pulsar node is publicly accessible for a Galaxy with the credentials to use.
  - You can send jobs to a Pulsar node from either a Galaxy instance running on SRC, or even inside your own network.
requirements:
  - type: "none"
    title: Access to the SURF Research Cloud
  - type: "none"
    title: SRAM authentication
  - type: "none"
    title: An SSH key connected to your account
  - type: "none"
    title: Have a basic understanding of how Galaxy works
time_estimation: "30m"
level: Introductory
contributions:
  authorship:
  - hexylena
  - mirelaminkova
  funding:
  - surf

priority: 11

edam_ontology:
- topic_0605 # Informatics
- topic_3071 # Data Management

abbreviations:
  SRC: SURF Research Cloud
  GDPR: General Data Protection Regulations
  SRAM: SURF Research Access Management
  CO: Collaborative Organisation

follow_up_training:
  - type: "internal"
    topic_name: admin
    tutorials:
      - surf-research-cloud-galaxy
subtopic: cloud
tags:
 - deploying
---


Using Pulsar via the {SRC} allows researchers to start Pulsar instances on-demand to expand their computational resources and even access GPUs to help and analyze their data in a secure environment following the {GDPR}.

There are two main use cases we envision this role being useful for:

Saving costs on SRC
:  Maybe you're already running Galaxy in SRC, but you don't want to run a GPU node because it is very expensive. By using the SRC Pulsar Catalog Item, you can launch a node to do computations and then shut it down when you're done, saving money. Pulsar instances can be started and stopped on demand, depending on personal cases and requirements, giving you a lot of freedom!

Accessing a GPU from a local (in UMC/hospital Galaxy)
:  If you do not have a GPU easily available within your institute, it may be attractive to send jobs securely to SRC, by launching a Pulsar node in SRC and attaching it to your institute's Galaxy instance.


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prerequisites

The instance provides secure authentication, where users must have a SURF Research account prior to this tutorial, have set the {SRAM} authentication method, and connect an SSH key to their accounts. In case you are not familiar with {SRC} and need help in setting up your accounts, please follow the instructions on the [SURF Knowledge Base](https://servicedesk.surf.nl/wiki/display/WIKI/SURF+Research+Cloud)

Inside the SRC members should have access to all publicly available catalog items. If you are not able to create a catalog item, please [contact SURF servicedesk](mailto:servicedesk@surf.nl).

This tutorial assumes you are member of a {CO} in {SRAM} that has access to {SRC} and a wallet with budget in SRC with enough sources to create Galaxy and Pulsar catalog items. (For more information please refer to the [SURF Knowledge Base](https://servicedesk.surf.nl/wiki/display/WIKI/Budgets%2C+wallets%2C+contracts).

You should have previous experience working with data inside Galaxy.

> <tip-title>SSH Access Required</tip-title>
> SSH access is required to reconfigure your Galaxy instance. Please make sure you set an SSH key in your SRAM account if you are planning to use Galaxy in SRC for this tutorial.
{: .tip}

> <hands-on-title>Log In to SRC</hands-on-title>
> 1. Log in to [SURF Research Cloud](https://portal.live.surfresearchcloud.nl/)
>
>    ![screenshot of SRC login portal with a bright yellow login button](images/surf1.png)
{: .hands_on}

# Launching Pulsar in SRC

> <hands-on-title>Launching Pulsar</hands-on-title>
> 1. In the **Workspaces Tab** on the bottom half of the screen, you'll find a **Plus Button** at right to add a new workspace
>
>    ![A small plus button is hovered over which says Add. Below galaxy-test is shown running with a yellow Access button](images/surf3.png)
>
> 2. Clicking that will let you choose any of the *Catalog Items* from SRC. They've got a wide selection but we're only interested in the two Pulsar Catalog Items
>
>    ![the two pulsar components in SRC: Galaxy Pulsar GPU Node (CUDA) and Galaxy Pulsar Node are shown. The second is expanded showing an author, Helena Rasche, and a description: let's you run galaxy jobs on another node](images/surf4.png)
>
>    > <warning-title>GPUs are Expensive</warning-title>
>    > The GPU nodes are *expensive*. In fact it
>    > was the motivating reason for building this catalog item: to enable you to
>    > launch a node, run some computations, and shut it down, saving you money.
>    {: .warning}
>
> 3. Creating a "workspace" from a catalog item (a template) is easy, most of the options are fixed for you, you just need to choose the size of the item. Pick an appropriate size for whatever computations you need to do.
>
>    ![workspace creation screen, choose the cloud provider is locked to SURF HPC, flavour is locked to 22.04, and the size options are available from 1GB/8CPU to 60C/750GB Hi-mem](images/surf5.png)
>
> 6. Pick a name, it can be anything, it does not matter. Check the expiration date to ensure it is just enough time for your computation and no more.
>
>    > <tip-title>Expiration Date</tip-title>
>    > The standard life-time of the VM is 5 days. If you need it for longer, this option can be changed once the machine is running.
>    > Note, that once the machine is expired and deleted it cannot be restored! Plan accordingly and migrate your data in time to prevent data loss!
>    >
>    > This is an incredibly useful feature as it saves you from forgetting to destroy a VM. Especially for GPU nodes it can help you ensure that they disappear after your computation is complete.
>    {:.tip}
>
>    ![almost there! some final details reads this page, asking for name and description. both have been filled out with 'pulsar'. a yellow submit button waits at the bottom](images/surf6.png)
>
> 7. Click submit when you are happy with the configuration.
>
{: .hands_on}

Once done, the workspace will be created for you. You'll need to wait ~5 minutes usually. Go for a beverage ☕️

![workspace list showing a workspace named pulsar being created.](images/surf7.png)

## Using Pulsar on SRC

> <hands-on-title>Access the Information Page</hands-on-title>
>
> 1. Once the workspace is up, you'll see an **Access** link:
>
>    ![A small plus button is hovered over which says Add. Below galaxy-test is shown running with a yellow Access button](images/surf3.png)
>
> 2. Click that will show you a Pulsar information page. This page is running on your pulsar node itself, and is restricted to ensure only authorised members can access the contents. It includes some configuration you will need to copy to your Galaxy node in order to make use of the Pulsar node.
>
>    ![pulsar configuration information page showing an about with admins and metadata like workspace fqdn. Configuration for galaxy is shown below including XML changes](images/surf8.png)
{: .hands_on}

This information page should have more than enough information to connect this Pulsar instance to your Galaxy server. You will need to reference information from this page in the following steps:

> <hands-on-title>Configuring Galaxy to use SRC Pulsar</hands-on-title>
> 1. Collect the requirements for accessing the Galaxy machine. You will need:
>
>    - your username from the first step
>    - your SSH key that is associated with your SRAM account
>
> 2. SSH into your *Galaxy* machine (not pulsar!).
>
>    ```
>    ssh -i path/to/your/sram-ssh-key username@galaxy.src-winter-school.src.surf-hosted.nl
>    ```
>
> 3. You will need to `sudo su` to do anything useful. Do that now.
> 4. `cd /srv/galaxy/` to move into the directory Galaxy configuration is stored.
> 5. The configuration is discussed fully in the Pulsar information, but it will be briefly covered here as well. Generally there are a few steps that must be followed:
>
>    - A runner must be registered
>    - A destination/environment must be added with the pulsar details
>    - Some tools should be redirected to this Pulsar
>
>    Here is an example of what those changes *might* look like in your Galaxy node.  In this example our pulsar node was called `p20` but that will be different for you.
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
>
>    ```diff
>     runners:
>       local:
>         load: galaxy.jobs.runners.local:LocalJobRunner
>         workers: 4
>       condor:
>         load: galaxy.jobs.runners.condor:CondorJobRunner
>    +  pulsar:
>    +    load: galaxy.jobs.runners.pulsar:PulsarRESTJobRunner
>
>     execution:
>       default: docker_dispatch
>       environments:
>         local_destination:
>           runner: local
>         # ... probably some more environments here.
>
>    +    remote_p20:
>    +       runner: pulsar
>    +       url: https://p20.src-sensitive-i.src.surf-hosted.nl
>    +       private_token: ySgfM1rnGIsiVN8Xlfk
>    +       dependency_resolution: remote
>    +       manager: _default_
>    +       # Uncomment the following to enable interactive tools:
>    +       docker_enabled: true
>    +       docker_set_user: null
>    +       docker_memory: "8G"
>    +       singularity_enabled: false
>    +       tmp_dir: true
>    +       outputs_to_working_directory: false
>    +       container_resolvers:
>    +       - type: explicit
>    +       require_container: True
>    +       container_monitor_command: /mnt/pulsar/venv/bin/galaxy-container-monitor
>    +       container_monitor_result: callback
>    +       container_monitor_get_ip_method: command:echo p20.src-sensitive-i.src.surf-hosted.nl
>
>
>     tools:
>     - class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>       environment: local_env
>     - id: Cut1
>       environment: condor_1x1
>    +- id: interactive_tool_jupyter_notebook
>    +  environment: remote_p20
>    +- id: interactive_tool_rstudio
>    +  environment: remote_p20
>    ```
{: .hands_on}

With the above, the minimal configuration is done, the Galaxy server will be aware of the Pulsar node, and two tools will be sent there: the Jupyter and RStudio Interactive Tools.

While you will simply copy-paste the runner and environment, you will need to identify yourself which tools should go to this Pulsar node.
You can find the tool ID from the dropdown at the top right, just to the left of "Execute" at the top:

![url bar and tool interface for the Cut1 tool](./images/id2.png)

> <tip-title>An Easy Configuration Option: Send Every Job to Pulsar</tip-title>
>
> If you are running jobs for a limited period of time, you might consider making this pulsar node the default destination. Remember to use the `remote_...` name of your pulsar node, based on what you copied. Not `remote_p20`.
>
> ```diff
>  execution:
> -  default: docker_dispatch
> +  default: remote_p20
>    environments:
>      local_destination:
>        runner: local
> ```
{: .tip}

With that, you're done, and for the length of time your node is running, your chosen tools (or everything) will be executed on that Pulsar node with more memory and CPU than the Galaxy host, and maybe a GPU as well!

> <hands-on-title>Launch a Job on Pulsar</hands-on-title>
> 1. Login to your Galaxy
> 2. Run one of the tools you have decided to send to Pulsar
> 3. On the pulsar machine, you can check that it runs by following the logs:
>
>    ```bash
>    sudo journalctl -fu pulsar
>    ```
>
{: .hands_on}

Congratulations on launching Pulsar in SRC! 🌌
