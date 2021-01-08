---
layout: tutorial_hands_on

title: "Use Singularity containers for running Galaxy jobs"
zenodo_link: ""
questions:
objectives:
  - Configure your Galaxy to use Singularity and BioContainers for running jobs
  - Use an Ansible playbook for all of the above
time_estimation: "1h"
key_points:
contributors:
  - torfinnnome
  - mvdbeek
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - cvmfs
---

# Overview
{:.no_toc}

In this tutorial you will learn how to configure Galaxy to run jobs using [Singularity](https://sylabs.io/singularity/) containers
provided by the [BioContainers](https://biocontainers.pro/) community.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# BioContainers

From the [BioContainers website](https://biocontainers.pro/):
> BioContainers is a community-driven project that provides the infrastructure and basic guidelines to create, manage and distribute bioinformatics packages (e.g conda) and containers (e.g docker, singularity). BioContainers is based on the popular frameworks Conda, Docker and Singularity.
>
> -- [https://biocontainers-edu.readthedocs.io/en/latest/what_is_biocontainers.html](https://biocontainers-edu.readthedocs.io/en/latest/what_is_biocontainers.html)
{: .quote}

# Singularity and Galaxy

From the Sylabs website:

> Singularity is a container platform. It allows you to create and run containers that package up pieces of software in a way that is portable and reproducible.
>
> -- [https://sylabs.io/guides/3.7/user-guide/introduction.html](https://sylabs.io/guides/3.7/user-guide/introduction.html)
{: .quote}

## Installing Singularity

First, we will install Singularity using Ansible.

> ### {% icon hands_on %} Hands-on: Installing Singularity with Ansible
>
>    > ### {% icon tip %} CentOS7
>    > If you are using CentOS7, you can skip this hands-on section and instead install the `epel-release` and `singularity` system packages.
>    {: .tip}
>
> 1. In your working directory, add the Singularity role to your `requirements.yml` file:
>
>    ```yaml
>    - src: cyverse-ansible.singularity
>    ```
>
> 2. Install the requirements with `ansible-galaxy`:
>
>    ```console
>    ansible-galaxy role install -p roles -r requirements.yml
>    ```
>
> 4. Specify which version of Singularity you want to install, in `group_vars/galaxyservers.yml`:
>
>    ```yaml
>    # Singularity target version
>    singularity_version: "3.7.0"
>    ```
> 4. Add the new role to the list of roles under the `roles` key in your playbook, `galaxy.yml`:
>
>    ```yaml
>    - hosts: galaxyservers
>      become: true
>      roles:
>        # ... existing roles ...
>        - cyverse-ansible.singularity
>    ```
>
> 5. Add the Go compiler to be installed using the Ansible `galaxy.yml` playbook:
>
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -4,7 +4,7 @@
>       pre_tasks:
>         - name: Install Dependencies
>           package:
>    -        name: ['git', 'make', 'python3-psycopg2', 'virtualenv', 'tar', 'bzip2']
>    +        name: ['git', 'make', 'python3-psycopg2', 'virtualenv', 'tar', 'bzip2', 'golang-go']
>       handlers:
>         - name: Restart Galaxy
>    ```
>
> 5. Run the playbook
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
>
> 6. Singularity should now be installed on your Galaxy server. You can test this by connecting
> to your server and run the following command:
>
>    ```console
>    $ singularity exec docker://hello-world
>    ```
{: .hands_on}

## Configure Galaxy to use Singularity

Now, we will configure Galaxy to run tools using Singularity containers, which will be
fetched from [Quay.io](https://quay.io/organization/biocontainers).

> ### {% icon hands_on %} Hands-on: Configure Galaxy to use Singularity
>
> 1. Edit the `group_vars/galaxyservers.yml` file and add a `dependency_resolvers_config_file` entry:
>{% raw %}
>    ```yaml
>    galaxy_config:
>      galaxy:
>        dependency_resolvers_config_file: "{{ galaxy_config_dir }}/dependency_resolvers_conf.xml"
>    ```
>{% endraw %}
>
> 2. Also in `group_vars/galaxyservers.yml`, add a `galaxy_config_files` entry:
> {% raw %}
>    ```yaml
>    galaxy_config_files:
>    - src: files/galaxy/config/dependency_resolvers_conf.xml
>      dest: "{{ galaxy_config_dir }}/dependency_resolvers_conf.xml"
>    ```
>{% endraw %}
>
> 3. Create the new file `files/galaxy/config/dependency_resolvers_conf.xml`:
>
>    ```xml
>    <dependency_resolvers>
>       <tool_shed_packages />
>       <galaxy_packages />
>       <galaxy_packages versionless="true" />
>    </dependency_resolvers>
>    ```
>
> 3. Now, we want to make Galaxy run jobs using Singularity. Modify the file `templates/galaxy/config/job_conf.xml.j2`, by adding the `singularity_enabled` parameter:
>
>    ```diff
>    --- a/templates/galaxy/config/job_conf.xml.j2
>    +++ b/templates/galaxy/config/job_conf.xml.j2
>    @@ -3,7 +3,9 @@
>         </plugins>
>         <destinations>
>             <destination id="local" runner="local"/>
>    -        <destination id="local_destination" runner="local_plugin"/>
>    +        <destination id="local_destination" runner="local_plugin">
>    +            <param id="singularity_enabled">true</param>
>    +        </destination>
>         </destinations>
>         <tools>
>         </tools>
>    ```
>
>    > ### {% icon tip %} Environment variables
>    > You might want to pass environment variables to Singularity, which can be configured like this:
>    >{% raw %}
>    >    ```xml
>    >    <destination id="local" runner="local">
>    >       <param id="singularity_enabled">true</param>
>    >       <!-- Some potentially useful environment modifications here: -->
>    >       <env id="LC_ALL">C</env>
>    >       <env id="SINGULARITY_CACHEDIR">/tmp/singularity</env>
>    >       <env id="SINGULARITY_TMPDIR">/tmp</env>
>    >    </destination>
>    >    ```
>    >{% endraw %}
>    {: .tip}
>
> 4. Re-run the playbook (`ansible-playbook galaxy.yml`)
>
> 5. In your Galaxy admin interface, install the minimap2 tool.
>
> 6. Upload a fasta file, and run the minimap2 tool using this fasta file as both reference and target file:
>
>    ```
>    >testing
>    GATTACAGATHISISJUSTATESTGATTACA
>    ```
>
> Your job should be executed using Singularity with a BioContainer!
>
{: .hands_on}

## Use Singularity containers from CVMFS

Galaxy can be configured to use pre-made Singularity containers available from /cvmfs/singularity.galaxyproject.org/.
In order to do so, you will first need to set up CVMFS by doing the [CVMFS]({{ site.baseurl }}/topics/admin/tutorials/cvmfs/tutorial.html) tutorial.
After finishing the CVMFS tutorial, come back, and do this hands-on.

> ### {% icon hands_on %} Optional: Hands-on: Configure Galaxy to use Singularity containers from CVMFS
>
> 1. Edit the `group_vars/galaxyservers.yml` file and add `containers_resolvers_config_file` and `galaxy_singularity_images_cvmfs_path`:
>{% raw %}
>    ```yaml
>    galaxy_singularity_images_cvmfs_path: "/cvmfs/singularity.galaxyproject.org/all/"
>    galaxy_config:
>      galaxy:
>        ...
>        containers_resolvers_config_file: "{{ galaxy_config_dir }}/container_resolvers_conf.xml"
>    ```
>{% endraw %}
>
> 2. Also in `group_vars/galaxyservers.yml`, add a `galaxy_config_templates` entry:
>{% raw %}
>    ```yaml
>    galaxy_config_templates:
>      - src: templates/galaxy/config/container_resolvers_conf.xml.j2
>        dest: "{{ galaxy_config_dir }}/container_resolvers_conf.xml"
>    ```
>{% endraw %}
>
> 3. Create the new file `templates/galaxy/config/container_resolvers_conf.xml.j2`:
>{% raw %}
>    ```xml
>    <containers_resolvers>
>        <explicit_singularity />
>        <cached_mulled_singularity cache_directory="{{ galaxy_singularity_images_cvmfs_path }}"/>
>    </containers_resolvers>
>    ```
>{% endraw %}
>
> 4. Re-run the playbook (`ansible-playbook galaxy.yml`)
>
> 5. Rerun your minimap2 job from the previous section. This time it should use a container served by CVMFS.
>
{: .hands_on}



