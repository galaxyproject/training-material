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

{% snippet topics/admin/faqs/ansible_local.md %}

> ### {% icon hands_on %} Hands-on: Installing CVMFS with Ansible
>
> 1. In your working directory, add the CVMFS role to your `requirements.yml`
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -18,3 +18,5 @@
>       version: 048c4f178077d05c1e67ae8d9893809aac9ab3b7
>     - src: gantsign.golang
>       version: 2.6.3
>    +- src: galaxyproject.cvmfs
>    +  version: 0.2.13
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement"}
>
> 2. Install the role with:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    {: .code-in data-cmd="true"}
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
>    > ### {% icon tip %} What is a good size for this?
>    > In production UseGalaxy.org.au uses 100GB, different sites have different needs and you can make your cache smaller depending on your usage. E.g. if your users only use one dataset from the reference data (e.g. just hg38) then perhaps you don't need such a large cache.
>    {: .tip}
>
>    Add the following lines to your `group_vars/all.yml` file, creating it if it doesn't exist:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/all.yml
>    @@ -0,0 +1,4 @@
>    +# CVMFS vars
>    +cvmfs_role: client
>    +galaxy_cvmfs_repos_enabled: config-repo
>    +cvmfs_quota_limit: 500
>    {% endraw %}
>    ```
>    {: data-commit="Configure CVMFS variables"}
>
>    > ### {% icon tip %} Why all.yml?
>    > We've integrated the cvmfs and pulsar tutorials better, such that CVMFS will be used for Pulsar as well, this configuration will be needed on all of our machines. This mirrors real life where you want CVMFS on every node that does computation.
>    {: .tip}
>
> 4. Add the new role to the list of roles under the `roles` key in your playbook, `galaxy.yml`:
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -27,3 +27,4 @@
>           become_user: "{{ galaxy_user.name }}"
>         - usegalaxy_eu.galaxy_systemd
>         - galaxyproject.nginx
>    +    - galaxyproject.cvmfs
>    {% endraw %}
>    ```
>    {: data-commit="Add role to playbook"}
>
> 5. Run the playbook
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > ansible-playbook galaxy.yml
>    > ```
>    {: .code-in data-cmd="true"}
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
>    > > You should see nothing, as CVMFS uses `autofs` in order to mount paths only upon request.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>
> 3. Change directory into `/cvmfs/data.galaxyproject.org/`.
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > cd /cvmfs/data.galaxyproject.org/
>    > ls
>    > ls byhand
>    > ls managed
>    > ```
>    {: .code-in}
>
>    > ### {% icon question %} Question
>    >
>    > What do you see now?
>    >
>    > > ### {% icon solution %} Solution
>    > >  You'll see `.loc` files, genomes and indices.
>    > > AutoFS only mounts the files when they're accessed, so it appears like there is no folder there.
>    > {: .solution }
>    >
>    {: .question}
>
>    And just like that we all have access to all the reference genomes and associated tool indices thanks to the Galaxy Project, IDC, and Nate's hard work!
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
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -29,6 +29,7 @@ miniconda_manage_dependencies: false
>
>     galaxy_config:
>       galaxy:
>    +    tool_data_table_config_path: /cvmfs/data.galaxyproject.org/byhand/location/tool_data_table_conf.xml,/cvmfs/data.galaxyproject.org/managed/location/tool_data_table_conf.xml
>         dependency_resolvers_config_file: "{{ galaxy_config_dir }}/dependency_resolvers_conf.xml"
>         containers_resolvers_config_file: "{{ galaxy_config_dir }}/container_resolvers_conf.xml"
>         brand: "🧬🔬🚀"
>    {% endraw %}
>    ```
>    {: data-commit="Add tool_data_table_config_path to group variables"}
>
>
> 2. Re-run the playbook
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > ansible-playbook galaxy.yml
>    > ```
>    {: .code-in data-cmd="true"}
>
> 3. Install the BWA-MEM tool, if needed.
>
>    {% snippet topics/admin/faqs/install_tool.md query="bwa" name="Map with BWA-MEM" section="Mapping" %}
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

# Common Production Questions

> ### {% icon question %} For the most used datasets (for ex. hg38) could we have a local copy, or would that be irrelevant?
> This would be irrelevant, the most used datasets will stay in the cache. CVMFS uses a Least Recently Used (LRU) cache (see their [docs](https://cvmfs.readthedocs.io/en/latest/cpt-details.html#disk-cache)), so whenever it runs out of space, it will remove the least recently used file. If you have a file that is very commonly used, it will remain in the cache.
{: .question}

> ### {% icon question %} Could you explain how to calculate a good cache space?
> Here are two approaches, there are others:
> 1. Allocate some cache, see how it is, make it larger if it is fully used + users complain of speed.
> 2. Enable reference data, and collect a week or two of data, analyse which reference datasets are being used, and allocate enough space for all of them.
>
> Essentially you just need data on how your users will behave and what reference data they want, combined with "when will they accept a wait period" to determine how much space you must allocate.
{: .question}

> ### {% icon question %} If I use a cluster, will I need to configure this FS in each node (given that the folder is at / directly)?
> Yes. Often admins with a cluster keep a smaller cache local to each compute node, and then setup a Squid proxy to hold the most commonly accessed data on a machine with more storage. E.g. each compute node could have 10-50GB of CVMFS storage while you might setup a Squid proxy with 200-300 GB of storage that will store everything your site uses.
{: .question}

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
