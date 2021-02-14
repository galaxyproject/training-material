---
layout: tutorial_hands_on

title: "Reference Data with CVMFS"
zenodo_link: ""
questions:
objectives:
  - Have an understanding of what CVMFS is and how it works
  - Install and configure the CVMFS client on a linux machine and mount the Galaxy reference data repository
  - Configure your Galaxy to use these reference genomes and indices
  - Use an Ansible playbook for all of the above.
time_estimation: "1h"
key_points:
contributors:
  - slugger70
  - hexylena
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

# Overview
{:.no_toc}

The CernVM-FS is a distributed filesystem perfectly designed for sharing readonly data across the globe. We use it in the [Galaxy Project](https://galaxyproject.org) for sharing things that a lot of Galaxy servers need. Namely:
* **Reference Data**
    * Genome sequences for hundreds of useful species.
    * Indices for the genome sequences
    * Various bioinformatic tool indices for the available genomes
* **Tool containers**
    * [Singularity](https://www.sylabs.io/) containers of everything stored in [Biocontainers](https://biocontainers.pro/) (A bioinformatic tool container repository.) You get these for free every time you build a [Bioconda](https://bioconda.github.io/) recipe/package for a tool.
* Others too..

From the Cern website:

> The CernVM File System provides a scalable, reliable and low-maintenance software distribution service. It was developed to assist High Energy Physics (HEP) collaborations to deploy software on the worldwide-distributed computing infrastructure used to run data processing applications. CernVM-FS is implemented as a POSIX read-only file system in user space (a FUSE module). Files and directories are hosted on standard web servers and mounted in the universal namespace /cvmfs."
>
> -- [https://cernvm.cern.ch/portal/filesystem](https://cernvm.cern.ch/portal/filesystem)
{: .quote}

A slideshow presentation on this subject can be found [here](slides.html). More details on the usegalaxy.org (Galaxy Main's) reference data setup and CVMFS system can be found [here](https://galaxyproject.org/admin/reference-data-repo/#usegalaxyorg-reference-data)

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Ansible-CVMFS and Galaxy

The Galaxy project supports a few CVMFS repositories.


| Repository                 | Repository Address              | Contents                                                                         |
| ----------                 | ------------------              | --------                                                                         |
| Reference Data and Indices | `data.galaxyproject.org`        | Genome sequences and their tool indices, Galaxy `.loc` files for them as well    |
| Singularity Containers     | `singularity.galaxyproject.org` | Singularity containers for everything in Biocontainers for use in Galaxy systems |
| Galaxy Main Configuration  | `main.galaxyproject.org`        | The configuration files etc for Galaxy Main (usegalaxy.org)                      |

You can browse the contents of `data.galaxyproject.org` at the [datacache](http://datacache.galaxyproject.org/).

## Installing and configuring Galaxy's CVMFS reference data with Ansible

Luckily for us, the Galaxy Project has a lot of experience with using and configuring CVMFS and we are going to leverage off that. To get CVMFS working on our Galaxy server, we will use the Ansible role for CVMFS written by the Galaxy Project. Firstly, we need to install the role and then write a playbook for using it.

If the terms "Ansible", "role" and "playbook" mean nothing to you, please checkout [the Ansible introduction slides]({% link topics/admin/tutorials/ansible/slides.html %}) and [the Ansible introduction tutorial]({% link topics/admin/tutorials/ansible/tutorial.md %})

{% include snippets/ansible_local.md %}

> ### {% icon hands_on %} Hands-on: Installing CVMFS with Ansible
>
> 1. In your working directory, add the CVMFS role to your `requirements.yml`
>
>    ```yaml
>    - src: galaxyproject.cvmfs
>      version: 0.2.13
>    ```
>
> 2. Install the role with:
>
>    ```console
>    ansible-galaxy role install -p roles -r requirements.yml
>    ```
>
> 3. The variables available in this role are:
>
>    | Variable             | Type          | Description                                                                                                                                                                    |
>    | ----------           | -------       | -------------                                                                                                                                                                  |
>    | `cvmfs_role`         | string        | Type of CVMFS host: `client`, `stratum0`, `stratum1`, or `localproxy`. Controls what packages are installed and what configuration is performed.                               |
>    | `cvmfs_keys`         | list of dicts | Keys to install on hosts of all types.                                                                                                                                         |
>    | `cvmfs_server_urls`  | list of dicts | CVMFS server URLs, the value of `CVMFS_SERVER_URL` in `/etc/cvmfs/domain.d/<domain>.conf`.                                                                                     |
>    | `cvmfs_repositories` | list of dicts | CVMFS repository configurations, `CVMFS_REPOSITORIES` in `/etc/cvmfs/default.local` plus additional settings in `/etc/cvmfs/repositories.d/<repository>/{client,server}.conf`. |
>    | `cvmfs_quota_limit`  | integer in MB | Size of CVMFS client cache. Default is `4000`.                                                                                                                                 |
>
>    But, luckily for us, the Galaxy Project CVMFS role has a lot of defaults for these variables which we can use by just setting `galaxy_cvmfs_repos_enabled` to `config-repo`. We will also set the `cvmfs_quota_limit` to something sensible (500MB) as we have relatively small disks on our instances. In a production setup, you should size this appropriately for the client.
>
>    Add the following lines to your `group_vars/all.yml` file, creating it if it doesn't exist:
>
>    ```yaml
>    # CVMFS vars
>    cvmfs_role: client
>    galaxy_cvmfs_repos_enabled: config-repo
>    cvmfs_quota_limit: 500
>    ```
>
>    > ### {% icon tip %} Why all.yml?
>    > We've integrated the cvmfs and pulsar tutorials better, such that CVMFS will be used for Pulsar as well, this configuration will be needed on all of our machines. This mirrors real life where you want CVMFS on every node that does computation.
>    {: .tip}
>
> 4. Add the new role to the list of roles under the `roles` key in your playbook, `galaxy.yml`:
>
>    ```yaml
>    - hosts: galaxyservers
>      become: true
>      roles:
>        # ... existing roles ...
>        - galaxyproject.cvmfs
>    ```
>
> 5. Run the playbook
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
{: .hands_on}

Congratulations, you've set up CVMFS.

## Exploring the CVMFS Installation

> ### {% icon hands_on %} Hands-on: Exploring CVMFS
>
> 1. SSH into your machine
>
> 2. Change directory into `/cvmfs/` and list the files in that folder
>
>    > ### {% icon question %} Question
>    >
>    > What do you see?
>    >
>    > > ### {% icon solution %} Solution
>    > > You should see nothing, as CVMFS uses `autofs` in order to mount paths only upon request. Once you `cd` into the directory, autofs will automatically mount the repository and files will be listed.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>
> 3. Change directory into `/cvmfs/data.galaxyproject.org/`. Have a browse through the contents. You'll see `.loc` files, genomes and indices.
>
>    And just like that we all have access to all the reference genomes and associated tool indices thanks to the Galaxy Project's and mostly Nate's hard work!
>
{: .hands_on}

## Configuring Galaxy to use the CVMFS references.

Now that we have mounted the CVMFS repository we need to tell Galaxy how to find it and use it.

There are two primary directories in the reference data repository:

| Directory   | Contents                                                                                                                                                                  |
| ----------- | ----------                                                                                                                                                                |
| `/managed`  | Data generated with Galaxy Data Managers, organized by data table (index format), then by genome build.                                                                   |
| `/byhand`   | Data generated prior to the existence/use of Data Managers, manually curated. (For legacy reasons, this directory is shared as `/indexes` on the HTTP and rsync servers.) |

These directories have somewhat different structures:

* `/managed` is organized by index type, then by genome build (Galaxy dbkey)
* `/byhand` is organzied by genome build, then by index type

Both directories contain a location subdirectory, and each of these contain a `tool_data_table_conf.xml` file:

* `/managed/location/tool_data_table_conf.xml`
* `/byhand/location/tool_data_table_conf.xml`

Galaxy consumes these `tool_data_table_conf.xml` files and the `.loc` "location" files they reference. The paths contained in these files are valid if the data is mounted via CVMFS.

Examples of data include:

* twoBit (`.2bit`) and FASTA (`.fa`) sequence files
* Bowtie 2 and BWA indexes
* Mutation Annotation Format (`.maf`) files
* SAMTools FASTA indexes (`.fai`)

Now all we need to do is tell Galaxy how to find it! This tutorial assumes that you have run the tutorial in the requirements, [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}). The hands-on below will use the Galaxy Project Ansible role to configure everything.

> ### {% icon hands_on %} Hands-on: Configuring Galaxy to use CVMFS
>
> 1. Edit the `group_vars/galaxyservers.yml` file and add a `tool_data_table_config_path` entry under the `galaxy` key of the `galaxy_config` section in the `group_vars/galaxyservers.yml` file. This new entry should be a list containing the paths to both `tool_data_table_conf.xml` files referenced above.
>
>    > ### {% icon question %} Question
>    >
>    > How does your final configuration look?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > ```yaml
>    > > galaxy_config:
>    > >   galaxy:
>    > >     # ... existing configuration options in the `galaxy` section ...
>    > >     tool_data_table_config_path: /cvmfs/data.galaxyproject.org/byhand/location/tool_data_table_conf.xml,/cvmfs/data.galaxyproject.org/managed/location/tool_data_table_conf.xml
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>
> 2. Re-run the playbook (`ansible-playbook galaxy.yml`)
>
> 3. Install the BWA-MEM tool, if needed.
>
>    {% include snippets/install_tool.md query="bwa" name="Map with BWA-MEM" section="Mapping" %}
>
> 4. In your Galaxy server, open the **Map with BWA-MEM** {% icon tool %} tool. Now check that there are a lot more reference genomes available for use!
>
>    ![available_genomes.png](../../images/available_genomes.png)
>
> 5. Login to Galaxy as the admin user, and go to **Admin → Data Tables → bwa_mem indexes**
>
>    ![bwa_mem_indices.png](../../images/bwa_mem_indices.png)
>
{: .hands_on}

# Other Aspects

## Development

If you are developing a new tool, and want to add a reference genome, we recommend you [talk to us on Gitter](https://gitter.im/galaxy-iuc/iuc). You can also look at one of the tools that uses reference data, and try and copy from that. If you’re developing the location files completely new, you need to write the data manager.

## Automation

You can automate the process of installing and setting up data managers and data with ephemeris. We're working in the [IDC](https://github.com/galaxyproject/idc) to democratise this CVMFS repository, and make this a community-controlled resource. You can also look here for ideas on automating your data management.

## Access Control

It is not easily possible to filter access to reference data depending on the user's role or group.

You could set up a tool per user/group, [secure access to running this tool](https://galaxyproject.org/admin/config/access-control/), and then allow this private tool to access a private tool data table. But you will not get tool updates, you will have to copy and edit this tool every time it gets updated. Or write more advanced job control rules to reject specific jobs which use specific datasets.

## Proxying Recap

The client talks directly to the stratum 1 (or to a proxy), and manages the data, and exposes it to the user. The proxy stores an opaque cache, that can't really be used, except as a proxy to a stratum 1.

## Plant Data

If you are working with plants, you can find separate reference data here: [frederikcoppens/galaxy_data_management](https://github.com/frederikcoppens/galaxy_data_management)
