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
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
voice:
  id: Olivia
  lang: en-AU
  neural: true
subtopic: data
tags:
  - ansible
  - git-gat
---

> These words come from a transcript of Simon Gladman teaching this course. He
> is a bioinformatician at the University of Melbourne in Australia and also
> one of the administrators of Galaxy Australia.
>
> Hello everybody and welcome back to the Galaxy
> administrators course. in this session we're going to be talking about
> reference data in Galaxy using CVMFS.
>
> Hopefully by now you've all seen the video of
> the slides - if not I suggest you go and do
> look at those first as that explains how reference
> data works in Galaxy and also what CVMFS does to
> help us with that reference data.
>
> Some of the requirements from that we hope you've
> already completed: Galaxy server installation
> with Ansible and hopefully you
> understand what Ansible is and how
> it's used.
{: .spoken data-visual="gtn" data-target="#top-navbar" }

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
{: .quote id="cvmfs-quote"}

> All right, just a quick recap. CVMFS or
> Cern-VMFS is a distributed file system perfectly
> designed for sharing read-only data across the
> globe and we use it extensively in the Galaxy
> project for sharing things that
> a lot of Galaxy servers need.
> Namely all the reference data; so the genome
> sequences for all the different genomes that
> we need to think about in in Galaxy and
> bioinformatics. Things like a human genome,
> mouse genome, etc., etc. And we need a lot of
> indices for the genome sequences and lots of
> tool indices for all those genomes. And we also
> have a repository that contains tool containers.
> So singularity containers for all of the
> bioinformatics tools that we might want to use.
> There's a tutorial idea already done this week
> or coming up soon this week that explains how to
> use those singularity containers within Galaxy as
> well. So we're just going to get on with things.
{: .spoken data-visual="gtn" data-target="#cvmfs-quote"}

A slideshow presentation on this subject can be found [here](slides.html). More details on the usegalaxy.org (Galaxy Main's) reference data setup and CVMFS system can be found [here](https://galaxyproject.org/admin/reference-data-repo/#usegalaxyorg-reference-data)

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="cvmfs" %}

> The agenda we're going to follow today is: We're going to install and
> configure Galaxy CVMFS reference data using ansible. We're going to explore
> the CVMFS installation and then we're going to configure Galaxy to use it.
{: .spoken data-visual="gtn" data-target="#agenda"}

# Ansible-CVMFS and Galaxy

> There are a few different repositories that um Galaxy project has created and
> shared with everybody. The first one is the reference data and indices and
> this is in data.galaxyproject.org. And then we have another one called
> singularity.galaxyproject.org.
>
> If you want to have a look at what's in there you can click on the "data
> cache", and this link here will show you all the things that are inside this
> data.galaxyproject.org
{: .spoken data-visual="gtn" data-target="#ansible-cvmfs-and-galaxy" }

The Galaxy project supports a few CVMFS repositories.


| Repository                 | Repository Address              | Contents                                                                         |
| ----------                 | ------------------              | --------                                                                         |
| Reference Data and Indices | `data.galaxyproject.org`        | Genome sequences and their tool indices, Galaxy `.loc` files for them as well    |
| Singularity Containers     | `singularity.galaxyproject.org` | Singularity containers for everything in Biocontainers for use in Galaxy systems |
| Galaxy Main Configuration  | `main.galaxyproject.org`        | The configuration files etc for Galaxy Main (usegalaxy.org)                      |

You can browse the contents of `data.galaxyproject.org` at the [datacache](http://datacache.galaxyproject.org/).

## Installing and Configuring

Luckily for us, the Galaxy Project has a lot of experience with using and configuring CVMFS and we are going to leverage off that. To get CVMFS working on our Galaxy server, we will use the Ansible role for CVMFS written by the Galaxy Project. Firstly, we need to install the role and then write a playbook for using it.

If the terms "Ansible", "role" and "playbook" mean nothing to you, please checkout [the Ansible introduction slides]({% link topics/admin/tutorials/ansible/slides.html %}) and [the Ansible introduction tutorial]({% link topics/admin/tutorials/ansible/tutorial.md %})

{% snippet topics/admin/faqs/ansible_local.md %}

> Okay, so now we're going to move on to installing and configuring Galaxy's
> CVMFS reference data with Ansible. We are going to do some Ansible here and
> we're going to install CVMFS onto our Galaxy server. Hopefully you all have
> access to a Galaxy server.
{: .spoken data-visual="gtn" data-target="#installing-and-configuring"}

> Here is mine. And on it I am an admin user and I have access to the admin
> page. This was done as part of the installation of Galaxy um and hopefully
> you've installed a tool.
{: .spoken data-visual="galaxy" data-target="/" data-action="goto"}

> I have bwa and bwa-mem installed under Mapping.
{: .spoken data-visual="galaxy" data-target=".search-input input" data-action="fill" data-value="bwa-mem"}

> So if I click on that you can see
> here that I have bwa-mem but there are no options available for reference
> genomes, so we want to fix that. We want to connect to all of Galaxy's
> pre-built references and so we're going to use Galaxy's CVMFS system to let
> our own Galaxies connect and get access to all the pre-built caches and
> everything we already have.
{: .spoken data-visual="galaxy" data-target="/?tool_id=bwa_mem" data-action="goto"}

> Okay, so let's get started. If we go back to our
> tutorial here, it says that we need to install a CVMFS role into our
> requirements.yml and then add it to our Ansible.
{: .spoken data-visual="gtn" data-target="#hands-on-installing-cvmfs-with-ansible"}

> <hands-on-title>Installing CVMFS with Ansible</hands-on-title>
>
> 1. In your working directory, add the CVMFS role to your `requirements.yml`
>
>    > Well the first thing I need to do is log into my Galaxy machine in the
>    > terminal.
>    {: .spoken data-visual="terminal" data-cmd="whoami; hostname -f"}
>
>    > Have a look at he contents of this directory and I'll go into galaxy and here we have all of the Ansible scripts that hopefully everybody already has.
>    {: .spoken data-visual="terminal" data-cmd="ls -al"}
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -14,3 +14,5 @@
>       version: 0.1.5
>     - name: galaxyproject.tusd
>       version: 0.0.1
>    +- src: galaxyproject.cvmfs
>    +  version: 0.2.13
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement" data-ref="add-req"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
>    > Okay, so the first thing I'm going to do is I'm going to add the CVMFS
>    > role to the requirements.yml.
>    > Edit requirements.yml and we need to add this to the bottom of that file. Copy. Paste. And save it.
>    {: .spoken data-visual="terminal" data-ref="add-req"}
>
> 2. Install the role with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true" data-ref="req-install"}
>    {: .code-in}
>
>    > And now install the role into our local Ansible scripts using the
>    > ansible-galaxy command. And as you can see, it's downloading the CVMFS
>    > role.
>    {: .spoken data-visual="terminal" data-ref="req-install"}
>
>    > And if we look into roles now you can see that we have galaxyproject.cmfs
>    {: .spoken data-visual="terminal" data-cmd="ls roles/"}
>
>    > Right, clear the screen.
>    {: .spoken data-visual="terminal" data-cmd="clear"}
>
> 3. The variables available in this role are:
>
>    > Okay, now what we need to do is we need to run this role and that will
>    > install the CVMFS client onto our Galaxy server. So the first thing we
>    > need to do is edit our group files galaxyservers file. There are a bunch
>    > of different variables that we can set. We can set the cvmfs role to be
>    > a client or stratum zero, stratum one server. URLs - so where we're
>    > getting all this data from, and if you remember in the slideshow, the
>    > data can come from Europe or America or Australia. Which repositories
>    > we want to have installed. And this is an important one - the quota
>    > limit. This is basically saying that CVMFS will cache some data on your
>    > local machine and the quota limit is the maximum size of that cache.
>    {: .spoken  data-visual="gtn" data-target="#variables-table"}
>
>    | Variable             | Type          | Description                                                                                                                                                                    |
>    | ----------           | -------       | -------------                                                                                                                                                                  |
>    | `cvmfs_role`         | string        | Type of CVMFS host: `client`, `stratum0`, `stratum1`, or `localproxy`. Controls what packages are installed and what configuration is performed.                               |
>    | `cvmfs_keys`         | list of dicts | Keys to install on hosts of all types.                                                                                                                                         |
>    | `cvmfs_server_urls`  | list of dicts | CVMFS server URLs, the value of `CVMFS_SERVER_URL` in `/etc/cvmfs/domain.d/<domain>.conf`.                                                                                     |
>    | `cvmfs_repositories` | list of dicts | CVMFS repository configurations, `CVMFS_REPOSITORIES` in `/etc/cvmfs/default.local` plus additional settings in `/etc/cvmfs/repositories.d/<repository>/{client,server}.conf`. |
>    | `cvmfs_quota_limit`  | integer in MB | Size of CVMFS client cache. Default is `4000`.                                                                                                                                 |
>    {: id="variables-table"}
>
>    But, luckily for us, the Galaxy Project CVMFS role has a lot of defaults for these variables which we can use by just setting `galaxy_cvmfs_repos_enabled` to `config-repo`. We will also set the `cvmfs_quota_limit` to something sensible (500MB) as we have relatively small disks on our instances. In a production setup, you should size this appropriately for the client.
>
>    > <tip-title>What is a good size for this?</tip-title>
>    > In production UseGalaxy.org.au uses 100GB, different sites have different needs and you can make your cache smaller depending on your usage. E.g. if your users only use one dataset from the reference data (e.g. just hg38) then perhaps you don't need such a large cache.
>    {: .tip}
>
>    > Okay but instead of just modifying galaxyservers.yml and adding in some
>    > of these variables - instead what we're going to do here is - we're
>    > going to create a new group file called all.yml.
>    >
>    > Because one of the things that we may want to do in the future - is we
>    > may want to create other machines. We might want to have worker nodes
>    > for our Galaxy cluster; or we may want to have other machines that we
>    > want to be able to create using these Ansible scripts that also have the
>    > CVMFS role installed. And instead of reproducing these variables in each
>    > of the group var files for those particular machines, we can create a
>    > special group vars file called all.yml. And whatever we put in there
>    > will be automatically available to all machines that we create with
>    > Ansible from this directory. Hopefully that makes a bit of sense and
>    > we're going to use that a bit later on if you come along to the pulsar
>    > tutorial where we will also be installing CVMFS on another machine - on
>    > a remote machine to run Pulsar.
>    {: .spoken data-visual="terminal" data-cmd="ls group_vars/"}
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
>    {: data-commit="Configure CVMFS variables" data-ref="vars-all"}
>
>    > So what we're going to do, is we're going to create a new file called
>    > groupvars/all.yml and we're going to put some of these CVMFS variables
>    > inside it. So i'll just copy that. Okay, groupvars all dot yaml. And
>    > I'll paste this in. So basically the CVMFS role we want this machine to
>    > have is client. Which means that we just want it to be able to access
>    > all of our CVMFS reference data. And then we want we're going to set
>    > this one here to say that, "yes we want to set this up for Galaxy." And
>    > the config repo is the one that tells CVMFS how to set everything else
>    > up. And then this is the other important one - the CVMFS quota limit.
>    > We're setting to 500 megabytes and that's just so we don't fill the root
>    > disk of these machines. So I'll save that.
>    {: .spoken data-visual="terminal" data-ref="vars-all"}
>
>    > <tip-title>Why all.yml?</tip-title>
>    > We've integrated the cvmfs and pulsar tutorials better, such that CVMFS will be used for Pulsar as well, this configuration will be needed on all of our machines. This mirrors real life where you want CVMFS on every node that does computation.
>    {: .tip}
>
> 4. Add the new role to the list of roles under the `roles` key in your playbook, `galaxy.yml`:
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -20,3 +20,4 @@
>           become_user: "{{ galaxy_user.name }}"
>         - galaxyproject.nginx
>         - galaxyproject.tusd
>    +    - galaxyproject.cvmfs
>    {% endraw %}
>    ```
>    {: data-commit="Add role to playbook" data-ref="pb"}
>
>    > And now we need to add the role - we need to add the role to our Galaxy
>    > playbook. So we edit galaxy.yml which is our playbook and we just need to
>    > add the galaxyproject.cvmfs to the bottom of this. galaxyproject.cvmfs.
>    > And that's pretty much it.
>    {: .spoken data-visual="terminal" data-ref="pb"}
>
> 5. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true" data-ref="run-pb1"}
>    {: .code-in }
{: .hands_on}

> Okay now we just run the playbook. ansible-playbook. We want to run the
> Galaxy one. So. All right, we're about to get to the CVMFS repo here. We are
> now it's installing the um the apt package of CVMFS. It's going to get
> that out of this special um Cern apt repository. Okay, it's installing
> it. Hopefully it won't take too long. Okay, now setting up the
> repositories. And it's done.
{: .spoken data-visual="terminal" data-ref="run-pb1"}

Congratulations, you've set up CVMFS.

## Exploring the CVMFS Installation


> Okay, that was completed. We've now installed CVMFS client onto our machine
> and we've told it to go looking for certain repositories. Now to get access
> to them. We'll see what's in them. uh They'll be located at slash cvmfs. So
> under the cvmfs directory in your root directory. So, we can go to that. cd
> /
{: .spoken data-visual="terminal" data-cmd="cd /"}

> Do an ll. You can see here there's a directory here called cvmfs.
{: .spoken data-visual="terminal" data-cmd="ls -al"}

>  So we'll go in there and have a look and see what's in there.
{: .spoken data-visual="terminal" data-cmd="cd /cvmfs"}

> So let's have a look. Oh there's nothing in there!
> Well actually, what what's going to happen is; as soon as we go looking for
> something in this directory, so autofs will automatically mount the
> particular thing we're looking for. And so what we're going to do here is I'm
> going to go: cd data.galaxyproject.org, because I know that's one of the one
> of the repositories that should have been installed.
{: .spoken data-visual="terminal" data-cmd="ls -al"}

> And when I do that, autofs is automatically going to mount it for me on the
> fly. Like that.
{: .spoken data-visual="terminal" data-cmd="cd data.galaxyproject.org"}

> And now I've cd'd into it and if I do an ll, you can see here I've got some
> things in here now. I've got byhand and managed.
{: .spoken data-visual="terminal" data-cmd="ls -al"}

>  If I go into byhand you can see here that I have quite a lot of different
>  genomes and their tool indices.
{: .spoken data-visual="terminal" data-cmd="cd byhand; ls"}

> So if I'm going to, say, sacCer2, I can see in here there are bowtie index,
> the bwa index, the uh the picard index, the sam index; a whole bunch of
> other different things - including the original sequences etc. So yeah,
> quite a lot of data and we just have access to that on the fly. And then as
> soon as we try and look at any of these files, what will happen is CVMFS
> will automatically cache it to the local disk within within that 500
> megabyte cache that we uh we set up earlier. This is really cool. And we can
> tell Galaxy to look at all of this data and use it as its reference data.
{: .spoken data-visual="terminal" data-cmd="ls -al sacCer2/"}

> <hands-on-title>Exploring CVMFS</hands-on-title>
>
> 1. SSH into your machine
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
* Multiple Alignment Format (`.maf`) files
* SAMTools FASTA indexes (`.fai`)

Now all we need to do is tell Galaxy how to find it! This tutorial assumes that you have run the tutorial in the requirements, [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}). The hands-on below will use the Galaxy Project Ansible role to configure everything.
{: id="spoken-7"}

> So that's what we're going to do now. Okay. So now we're going to try and
> configure Galaxy to use this CVMFS data. And and to have it so that we can run
> things like bwa and bwa mem and run them against the human genome or the bee
> genome or the mouse genome and take advantage of the fact that a lot of other
> people in the Galaxy community have done a lot of work for reference data for
> us already.
>
> So the way to do this is we're going to edit the groupvars galaxyservers
> file and we're going to add a variable called tool_data_table_config_path.
> And then we're going to point it to the two files that are in - there's one
> in byhand and one in managed.
{: .spoken  data-visual="gtn" data-target="#spoken-7"}


> If we go into byhand you can see here. And then we're going to location. As
> you can see in here, are all the lock files, but you'll also see there's an
> xml file here called tool_data_table_conf xml and we're going to point Galaxy
> at this file and there's another one in the same position in managed.
{: .spoken data-visual="terminal" data-cmd="cd /cvmfs/data.galaxyproject.org/byhand;  ls location;"}

> And you can see here there's another one in managed there. And so we're going to
> add both of these files to our Galaxy configuration and then Galaxy will be
> able to use all of the data contained within this repository.
{: .spoken data-visual="terminal" data-cmd="cd ../; ls managed/location"}

> Okay, so we'll go back to our Ansible directory.
{: .spoken data-visual="terminal" data-cmd="cd ~/galaxy/"}

> <hands-on-title>Configuring Galaxy to use CVMFS</hands-on-title>
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
>         brand: "🧬🔬🚀"
>         admin_users: admin@example.org
>         database_connection: "postgresql:///galaxy?host=/var/run/postgresql"
>    {% endraw %}
>    ```
>    {: data-commit="Add tool_data_table_config_path to group variables" data-ref="gvconf"}
>
>    > This time I'm going to edit the groupvars galaxyservers.yml file. And in
>    > our Galaxy section. Which is here. Which is here. At the bottom of that
>    > I'm going to add a variable called tool_data_table_config_path. And I'm
>    > going to point it to the locations that we um I showed you before. um I
>    > can't remember what they are off the top of my head but luckily they're
>    > inside this solution box and so I will just copy them. And paste. And as
>    > you can see pointing to /cvmfs/data.galaxyproject.org/byhand/location
>    > and then that tool data table conf xml file. And then we have a list
>    > here and we separate it by commas and then we point it to the second
>    > one. Right, so we save this file.
>    {: .spoken data-visual="terminal" data-ref="gvconf"}
>
> 2. Re-run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true" data-ref="pb-run2"}
>    {: .code-in }
>
>    > So this time what we're going to do, all we're doing is making a minor change
>    > to the Galaxy yaml file in the Galaxy config to add that one line and then
>    > we're going to restart Galaxy. And yeah, Galaxy will suddenly automatically
>    > have access to all of that data. So you can see here we've changed the Galaxy
>    > configuration file. And then Galaxy is now restarting and it's done.
>    {: .spoken data-visual="terminal" data-ref="pb-run2"}
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

> Okay let's go and have a look at our Galaxy server and see if bwa can
> suddenly see all of those - that stuff. All right so we're back on our Galaxy
> server. I'll click on Analyze Data to just to reload the page.
{: .spoken data-visual="galaxy" data-target="#analysis a" data-action="click"}

> I'll go back to Mapping
{: .spoken data-visual="galaxy" data-target=".search-input input" data-action="click"}

>  and load bwa-mem.
{: .spoken data-visual="galaxy" data-target=".search-input input" data-action="fill" data-value="bwa-mem"}

> And then suddenly, instead of
> having no options available, you can see here we've got the Bee genome.
{: .spoken data-visual="galaxy" data-target="/?tool_id=bwa_mem" data-action="goto"}

>  Now click on that. Oh look at that, there are lots and lots and lots of
>  available genomes now including: lots of human, mouse, rat, yeast, all sorts
>  of things. And in fact if you want to see the list of all the different
>  available genomes now, that we have available to us.
>
>  If you go to admin. We go to data tables over here. You can see here that we
>  have um a couple of data tables for managed and for all fasta. So if we
>  click on that one, you can see that we have a lot of genomes available now
>  in the all fasta data table that Galaxy can get access to. If we go back to
>  the data tables again, and go down to bwa indexes or bwa mem indexes here.
>  You can see we have access to a lot of pre-built indexes for bwa for all of
>  these different genomes. That is pretty powerful.
>
>  So what did that take us? Maybe 30 minutes? And uh suddenly our Galaxy
>  server has access to all the uh the data the reference data and the tool
>  indices that the community have built over a number of years and it's super
>  simple.
{: .spoken data-visual="galaxy" data-target="s2id_field-uid-12_select a" data-action="click"}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}

{% snippet topics/admin/faqs/missed-something.md step=5 %}

# Common Production Questions

> <question-title>For the most used datasets (for ex. hg38) could we have a local copy, or would that be irrelevant?</question-title>
> This would be irrelevant, the most used datasets will stay in the cache. CVMFS uses a Least Recently Used (LRU) cache (see their [docs](https://cvmfs.readthedocs.io/en/latest/cpt-details.html#disk-cache)), so whenever it runs out of space, it will remove the least recently used file. If you have a file that is very commonly used, it will remain in the cache.
{: .question}

> <question-title>Could you explain how to calculate a good cache space?</question-title>
> Here are two approaches, there are others:
> 1. Allocate some cache, see how it is, make it larger if it is fully used + users complain of speed.
> 2. Enable reference data, and collect a week or two of data, analyse which reference datasets are being used, and allocate enough space for all of them.
>
> Essentially you just need data on how your users will behave and what reference data they want, combined with "when will they accept a wait period" to determine how much space you must allocate.
{: .question}

> <question-title>If I use a cluster, will I need to configure this FS in each node (given that the folder is at / directly)?</question-title>
> Yes. Often admins with a cluster keep a smaller cache local to each compute node, and then setup a Squid proxy to hold the most commonly accessed data on a machine with more storage. E.g. each compute node could have 10-50GB of CVMFS storage while you might setup a Squid proxy with 200-300 GB of storage that will store everything your site uses.
{: .question}

> <tip-title>Debugging failed mounting</tip-title>
>
> Are you having issues mounting your CVMFS mount? Is it giving strange errors like "Endpoint not connected"
> Try running this command as root:
>
> > > <code-in-title>Bash</code-in-title>
> > > ```console
> > > /usr/bin/cvmfs2 -d -o rw,system_mount,fsname=cvmfs2,allow_other,grab_mountpoint singularity.galaxyproject.org /mnt
> > > ```
> > {: .code-in}
> >
> > > <code-out-title>Consolue</code-out-title>
> > > ```console
> > > Debug: using library /usr/lib/libcvmfs_fuse3_stub.so
> > > CernVM-FS: running in debug mode
> > > CernVM-FS: loading Fuse module... (cvmfs) Parsing config file /etc/cvmfs/default.conf    [07-21-2022 11:11:20 UTC]
> > > (cvmfs) execve'd /bin/sh (PID: 280373)    [07-21-2022 11:11:20 UTC]
> > > (cvmfs) Parsing config file /etc/cvmfs/default.d/50-cern-debian.conf    [07-21-2022 11:11:20 UTC]
> > > (cvmfs) execve'd /bin/sh (PID: 280375)    [07-21-2022 11:11:20 UTC]
> > > (cvmfs) Parsing config file /etc/cvmfs/default.d/80-ansible-galaxyproject-cvmfs.conf    [07-21-2022 11:11:20 UTC]
> > > (cvmfs) execve'd /bin/sh (PID: 280378)    [07-21-2022 11:11:20 UTC]
> > > [...]
> > > (dns) empty hostname    [07-21-2022 11:11:20 UTC]
> > > (download) installed 1 proxies in 1 load-balance groups    [07-21-2022 11:11:20 UTC]
> > > (cvmfs) DNS roaming is disabled for this repository.    [07-21-2022 11:11:20 UTC]
> > > (catalog) constructing client catalog manager    [07-21-2022 11:11:20 UTC]
> > > (catalog) Initialize catalog    [07-21-2022 11:11:20 UTC]
> > > (cache) unable to read local checksum    [07-21-2022 11:11:20 UTC]
> > > (download) escaped http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/.cvmfspublished to http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/.cvmfspublished    [07-21-2022 11:11:20 UTC]
> > > (download) Verify downloaded url /.cvmfspublished, proxy DIRECT (curl error 0)    [07-21-2022 11:11:20 UTC]
> > > (cache) miss ./e2/ab48b0984729d99951cb62c4312f501b3ddc6b (-2)    [07-21-2022 11:11:20 UTC]
> > > (download) escaped http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/data/e2/ab48b0984729d99951cb62c4312f501b3ddc6bX to http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/data/e2/ab48b0984729d99951cb62c4312f501b3ddc6bX    [07-21-2022 11:11:20 UTC]
> > > (download) Verify downloaded url /data/e2/ab48b0984729d99951cb62c4312f501b3ddc6bX, proxy DIRECT (curl error 0)    [07-21-2022 11:11:20 UTC]
> > > (download) escaped http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/.cvmfswhitelist to http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/.cvmfswhitelist    [07-21-2022 11:11:20 UTC]
> > > [...]
> > > (download) escaped http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/data/c7/f1555f421b1868b979291dc23f34a83132eadbC to http://cvmfs1-psu0.galaxyproject.org/cvmfs/singularity.galaxyproject.org/data/c7/f1555f421b1868b979291dc23f34a83132eadbC    [07-21-2022 11:11:20 UTC]
> > > (download) Verify downloaded url /data/c7/f1555f421b1868b979291dc23f34a83132eadbC, proxy DIRECT (curl error 0)    [07-21-2022 11:11:25 UTC]
> > > (cache) finished downloading of /data/c7/f1555f421b1868b979291dc23f34a83132eadbC    [07-21-2022 11:11:25 UTC]
> > > (cache) commit ./c7/f1555f421b1868b979291dc23f34a83132eadb ./txn/fetchJWcwtt    [07-21-2022 11:11:25 UTC]
> > > (quota) pin into lru c7f1555f421b1868b979291dc23f34a83132eadb, path file catalog at singularity.galaxyproject.org:/ (c7f1555f421b1868b979291dc23f34a83132eadb)    [07-21-2022 11:11:25 UTC]
> > > (cache) commit failed: cannot pin c7f1555f421b1868b979291dc23f34a83132eadb    [07-21-2022 11:11:25 UTC]
> > > (catalog) failed to load catalog '' (2 - not enough space to load catalog)    [07-21-2022 11:11:25 UTC]
> > > (catalog) failed to initialize root catalog    [07-21-2022 11:11:25 UTC]
> > > Failed to initialize root file catalog (16 - file catalog failure)
> > > (cache) unpinning / unloading all catalogs    [07-21-2022 11:11:25 UTC]
> > > ```
> > {: .code-out}
> {: .code-2col}
{: .tip}

# Other Aspects

> Right, we'll go back to our tutorial. um Yeah. Just finally, just before we
> finish up. um If we are developing a new tool and you want to add a reference
> genome or a different index just give us - drop us a line on Gitter and we'll
> be able to add it into our - into the reference data for the community. um
> We're looking at automating the process of building all of this material
> using data managers and ephemeris. And we're working with a group of people
> called the IDC, which is the Intergalactic Data Commission, which is a funny
> name for everyone in Galaxy - in the Galaxy community who likes reference
> data. And we're looking at making a community controlled resource that will
> be semi-automatic. One of the other things that you can do is have automatic
> fallbacks. So, if say, you're in Australia and you're hooked up to the
> Australian mirror of the CVMFS repository and the Australian mirror dies, the
> CVMFS client is smart enough to automatically go to the next closest one and
> so you won't lose anything. If you're interested in looking at plant data
> there's a link here for that.
{: .spoken data-visual="gtn" data-target="#other-aspects" }

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

> And finally if you could please click on this link here and give us some
> feedback on how you think the tutorial went, whether it was useful, if you
> enjoyed it or um if you have any criticisms, could you please put them in
> here as well. And if you end up using this to build a Galaxy server and you
> publish that Galaxy server um could you cite the tutorial for us please. That
> would be make a big difference to us. All right, thank you very much and I
> hope you enjoyed it and hopefully I'll get to meet some of you in person one
> day soon at a Galaxy conference. Thank you and goodbye.
{: .spoken data-visual="gtn" data-target="#feedback"}
