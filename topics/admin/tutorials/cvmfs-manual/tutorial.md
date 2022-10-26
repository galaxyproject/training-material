---
layout: tutorial_hands_on

title: "Reference Data with CVMFS without Ansible"
zenodo_link: ""
questions:
objectives:
  - Have an understanding of what CVMFS is and how it works
  - Install and configure the CVMFS client on a linux machine and mount the Galaxy reference data repository
  - Configure your Galaxy to use these reference genomes and indices
time_estimation: "1h"
key_points:
contributors:
  - slugger70
  - hexylena
subtopic: data
---

# Overview


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

A slideshow presentation on this subject can be found [here]({% link topics/admin/tutorials/cvmfs/slides.html %}). More details on the usegalaxy.org (Galaxy Main's) reference data setup and CVMFS system can be found [here](https://galaxyproject.org/admin/reference-data-repo/#usegalaxyorg-reference-data)

There are two sections to this exercise. The first shows you how to use Ansible to setup and configure CVMFS for Galaxy. The second shows you how to do everything manually. It is recommended that you use the Ansible method. The manual method is included here mainly for a more in depth understanding of what is happening.

If you really want to perform all these tasks manually, go [here](#cvmfs-and-galaxy-without-ansible), otherwise just follow along.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}



# CVMFS and Galaxy without Ansible

> <comment-title>Manual version of Ansible Commands</comment-title>
> If you wish to perform the same thing that we've just done, but by building the ansible script manually, follow these instructions. Otherwise, you have already done everything below and do not need to re-do it.
{: .comment}

We are going to setup a CVMFS mount to the Galaxy reference data repository on our machines. To do this we have to install and configure the CVMFS client and then mount the appropriate CVMFS repository using the publicly available keys.

> <hands-on-title>Installing the CVMFS Client</hands-on-title>
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

> <hands-on-title>Configuring CVMFS</hands-on-title>
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
>    This tells CVMFS to mount the Galaxy reference data repository and use a specific location for the cache which is limited to 500MB in size and to use the instance's geo-location to choose the best CVMFS repo server to connect to. You can use the `cvmfs_quota_limit` role variable to control this setting.
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

> <hands-on-title>Testing it out</hands-on-title>
>
> 1. Run `sudo cvmfs_config probe data.galaxyproject.org`
>
>    > <question-title></question-title>
>    >
>    > What does it output?
>    >
>    > > <solution-title></solution-title>
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
>    > <question-title></question-title>
>    >
>    > What do you see?
>    >
>    > > <solution-title></solution-title>
>    > > You should see nothing, as CVMFS uses `autofs` in order to mount paths only upon request.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>
> 3. Change directory into `/cvmfs/data.galaxyproject.org/`.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > cd /cvmfs/data.galaxyproject.org/
>    > ls
>    > ls byhand
>    > ls managed
>    > ```
>    {: .code-in}
>
>    > <question-title></question-title>
>    >
>    > What do you see now?
>    >
>    > > <solution-title></solution-title>
>    > >  You'll see `.loc` files, genomes and indices.
>    > > AutoFS only mounts the files when they're accessed, so it appears like there is no folder there.
>    > {: .solution }
>    >
>    {: .question}
>
>    And just like that we all have access to all the reference genomes and associated tool indices thanks to the Galaxy Project, IDC, and Nate's hard work!
>
>    > <tip-title>Contributing Reference Genomes</tip-title>
>    > If you are developing a new tool, and want to add a reference genome, we recommend you [talk to us on Gitter](https://gitter.im/galaxy-iuc/iuc). You can also look at one of the tools that uses reference data, and try and copy from that. If you’re developing the location files completely new, you need to write the data manager.
>    {: .tip}
>
{: .hands_on}

## Look at the repository

Now to configure Galaxy to use the CVMFS references we have just installed, see [the Ansible tutorial.]({% link topics/admin/tutorials/cvmfs/tutorial.md %}#configuring-galaxy-to-use-the-cvmfs-references)

