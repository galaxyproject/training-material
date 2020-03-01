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

There are two sections to this exercise. The first shows you how to use Ansible to setup and configure CVMFS for Galaxy. The second shows you how to do everything manually. It is recommended that you use the Ansible method. The manual method is included here mainly for a more in depth understanding of what is happening.

If you really want to perform all these tasks manually, go [here](#cvmfs-and-galaxy-without-ansible), otherwise just follow along.

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
>    ---
>    - src: galaxyproject.cvmfs
>      version: 0.2.8
>    ```
>
> 2. Install the requirements with `ansible-galaxy`:
>
>    ```console
>    ansible-galaxy role install -p roles -r requirements.yml
>    ```
>
> 3. Edit the group variables file, `group_vars/galaxyservers.yml`:
>
>    <br/>
>
>    The variables available in this role are:
>
>    | Variable             | Type          | Description                                                                                                                                                                    |
>    | ----------           | -------       | -------------                                                                                                                                                                  |
>    | `cvmfs_role`         | string        | Type of CVMFS host: `client`, `stratum0`, `stratum1`, or `localproxy`. Controls what packages are installed and what configuration is performed.                               |
>    | `cvmfs_keys`         | list of dicts | Keys to install on hosts of all types.                                                                                                                                         |
>    | `cvmfs_server_urls`  | list of dicts | CVMFS server URLs, the value of `CVMFS_SERVER_URL` in `/etc/cvmfs/domain.d/<domain>.conf`.                                                                                     |
>    | `cvmfs_repositories` | list of dicts | CVMFS repository configurations, `CVMFS_REPOSITORIES` in `/etc/cvmfs/default.local` plus additional settings in `/etc/cvmfs/repositories.d/<repository>/{client,server}.conf`. |
>    | `cvmfs_quota_limit`  | integer in MB | Size of CVMFS client cache. Default is `4000`.                                                                                                                                 |
>
>    <br/>
>
>    But, luckily for us, the Galaxy Project CVMFS role has a lot of defaults for these variables which we can use by just setting `galaxy_cvmfs_repos_enabled` to `config-repo`. We will also set the `cvmfs_quota_limit` to something sensible (500MB) as we have relatively small disks on our instances. In a production setup, you should size this appropriately for the client.
>
>    <br/>
>
>    Add the following lines to your `group_vars/galaxyservers.yml` file:
>
>    ```yaml
>    # CVMFS vars
>    cvmfs_role: client
>    galaxy_cvmfs_repos_enabled: config-repo
>    cvmfs_quota_limit: 500
>    ```
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
>    > >     tool_data_table_config_path:
>    > >       - /cvmfs/data.galaxyproject.org/byhand/location/tool_data_table_conf.xml
>    > >       - /cvmfs/data.galaxyproject.org/managed/location/tool_data_table_conf.xml
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>
> 2. Re-run the playbook (`ansible-playbook galaxy.yml`)
>
> 3. In your Galaxy interface, open the **BWA** {% icon tool %}, **BWA-MEM** {% icon tool %} or **Bowtie2** {% icon tool %} tool interface (whichever you may have installed). Now check that there are a lot more Genomes available for use!
>
>    ![available_genomes.png](../../images/available_genomes.png)
>
> 4. Login to Galaxy as the admin user, and go to **Admin → Data Tables → bwa_mem indexes**
>
>    ![bwa_mem_indices.png](../../images/bwa_mem_indices.png)
>
{: .hands_on}

You've now finished the tutorial, and you can [jump to the end](#feedback-google) or read on to learn about configuring CVMFS without Ansible.

# CVMFS and Galaxy without Ansible

> ### {% icon comment %} Manual version of Ansible Commands
> If you wish to perform the same thing that we've just done, but by building the ansible script manually, follow these instructions. Otherwise, you have already done everything below and do not need to re-do it.
{: .comment}

We are going to setup a CVMFS mount to the Galaxy reference data repository on our machines. To do this we have to install and configure the CVMFS client and then mount the appropriate CVMFS repository using the publicly available keys.

> ### {% icon hands_on %} Hands-on: Installing the CVMFS Client
>
> 1. On your remote machine, we need to first install the Cern software apt repo and then the CVMFS client and config utility:
>
>    ```bash
>    sudo apt install lsb-release
>    wget https://ecsft.cern.ch/dist/cvmfs/cvmfs-release/cvmfs-release-latest_all.deb
>    sudo dpkg -i cvmfs-release-latest_all.deb
>    rm -f cvmfs-release-latest_all.deb
>    sudo apt-get update
>
>    sudo apt install cvmfs cvmfs-config
>    ```
>
> 2. Now we need to run the CVMFS setup script.
>
>    ```bash
>    sudo cvmfs_config setup
>    ```
>
{: .hands_on}


## Configuring CVMFS

The configuration is not complex for CVMFS:

> ### {% icon hands_on %} Hands-on: Configuring CVMFS
>
> 1. Create a `/etc/cvmfs/default.local` file with the following contents:
>
>    ```
>    CVMFS_REPOSITORIES="data.galaxyproject.org"
>    CVMFS_HTTP_PROXY="DIRECT"
>    CVMFS_QUOTA_LIMIT="500"
>    CVMFS_CACHE_BASE="/srv/cvmfs/cache"
>    CVMFS_USE_GEOAPI=yes
>    ```
>
>    This tells CVMFS to mount the Galaxy reference data repository and use a specific location for the cache which is limited to 500MB in size and to use the instance's geo-location to choose the best CVMFS repo server to connect to.
>
> 2. Create a `/etc/cvmfs/domain.d/galaxyproject.org.conf` file with the following contents:
>
>    ```
>    CVMFS_SERVER_URL="http://cvmfs1-tacc0.galaxyproject.org/cvmfs/@fqrn@;http://cvmfs1-iu0.galaxyproject.org/cvmfs/@fqrn@;http://cvmfs1-psu0.galaxyproject.org/cvmfs/@fqrn@;http://galaxy.jrc.ec.europa.eu:8008/cvmfs/@fqrn@;http://cvmfs1-mel0.gvl.org.au/cvmfs/@fqrn@;http://cvmfs1-ufr0.galaxyproject.eu/cvmfs/@fqrn@"
>    ```
>
>    This is a list of the available stratum 1 servers that have this repo.
>
> 3. Create a `/etc/cvmfs/keys/data.galaxyproject.org.pub` file with the following contents:
>
>    ```
>    -----BEGIN PUBLIC KEY-----
>    MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEA5LHQuKWzcX5iBbCGsXGt
>    6CRi9+a9cKZG4UlX/lJukEJ+3dSxVDWJs88PSdLk+E25494oU56hB8YeVq+W8AQE
>    3LWx2K2ruRjEAI2o8sRgs/IbafjZ7cBuERzqj3Tn5qUIBFoKUMWMSIiWTQe2Sfnj
>    GzfDoswr5TTk7aH/FIXUjLnLGGCOzPtUC244IhHARzu86bWYxQJUw0/kZl5wVGcH
>    maSgr39h1xPst0Vx1keJ95AH0wqxPbCcyBGtF1L6HQlLidmoIDqcCQpLsGJJEoOs
>    NVNhhcb66OJHah5ppI1N3cZehdaKyr1XcF9eedwLFTvuiwTn6qMmttT/tHX7rcxT
>    owIDAQAB
>    -----END PUBLIC KEY-----
>    ```
>
> 4. Make a directory for the cache files
>
>    ```
>    sudo mkdir /srv/cvmfs
>    ```
{: .hands_on}


## Testing it out

Probe the connection.

> ### {% icon hands_on %} Hands-on: Testing it out
>
> 1. Run `sudo cvmfs_config probe data.galaxyproject.org`
>
>    > ### {% icon question %} Question
>    >
>    > What does it output?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > ```
>    > > OK
>    > > ```
>    > >
>    > > If this doesn't return `OK` then you may need to restart autofs: `sudo systemctl restart autofs`
>    > >
>    > {: .solution }
>    >
>    {: .question}
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

#### Step 4: Look at the repository

Now to configure Galaxy to use the CVMFS references we have just installed, see [here](#configuring-galaxy-to-use-the-cvmfs-references)
