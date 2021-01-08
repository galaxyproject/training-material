---
layout: tutorial_hands_on

title: "Use Singularity containers for running Galaxy jobs"
zenodo_link: ""
questions:
objectives:
  - Configure your Galaxy to use Singularity and BioContainers for running jobs
time_estimation: "1h"
key_points:
contributors:
  - torfinnnome
  - mvdbeek
  - hexylena
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

In this tutorial you will learn how to configure Galaxy to run jobs using [Singularity](https://sylabs.io/singularity/) containers provided by the [BioContainers](https://biocontainers.pro/) community.

## Background

> BioContainers is a community-driven project that provides the infrastructure and basic guidelines to create, manage and distribute bioinformatics packages (e.g conda) and containers (e.g docker, singularity). BioContainers is based on the popular frameworks Conda, Docker and Singularity.
>
> -- [https://biocontainers-edu.readthedocs.io/en/latest/what_is_biocontainers.html](https://biocontainers-edu.readthedocs.io/en/latest/what_is_biocontainers.html)
{: .quote}

Singularity is an alternative to Docker that is much friendlier for HPCs

> Singularity is a container platform. It allows you to create and run containers that package up pieces of software in a way that is portable and reproducible.
>
> -- [https://sylabs.io/guides/3.7/user-guide/introduction.html](https://sylabs.io/guides/3.7/user-guide/introduction.html)
{: .quote}

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Installing Singularity

First, we will install Singularity using Ansible. On most operating systems there is no package for singularity yet, so we must use a role which will compile it from source. If you're on CentOS7/8, it is available through the EPEL repository.

> ### {% icon tip %} CentOS7
> If you are using CentOS7, you can skip this hands-on section and instead install the `epel-release` and `singularity` system packages in your `pre_tasks`.
{: .tip}

> ### {% icon hands_on %} Hands-on: Installing Singularity with Ansible
>
> 1. In your working directory, add the Singularity role to your `requirements.yml` file:
>
>    ```yaml
>    - src: cyverse-ansible.singularity
>      version: ad4de5e4b0bb3f8a43de0f3565757aa156853485
>    - src: gantsign.golang
>      version: 2.6.3
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
>    # Golang
>    golang_gopath: '/opt/workspace-go'
>    # Singularity target version
>    singularity_version: "3.7.0"
>    singularity_go_path: "{{ golang_install_dir }}"
>    ```
> 4. Add the new roles to your `galaxy.yml` playbook, before the Galaxy server itself. We'll do this bceause it's a dependency of Galaxy to run, so it needs to be there before Galaxy starts.
>
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -16,6 +16,8 @@
>           become: true
>           become_user: postgres
>         - geerlingguy.pip
>    +    - gantsign.golang
>    +    - cyverse-ansible.singularity
>         - galaxyproject.galaxy
>         - role: uchida.miniconda
>           become: true
>    ```
>
> 5. Run the playbook
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > ansible-playbook galaxy.yml
>    > ```
>    {: .code-in}
>
> 6. Singularity should now be installed on your Galaxy server. You can test this by connecting
> to your server and run the following command:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > singularity exec docker://hello-world /hello
>    > ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output: Bash
>    > ```
>    > INFO:    Converting OCI blobs to SIF format
>    > INFO:    Starting build...
>    > Getting image source signatures
>    > Copying blob 0e03bdcc26d7 done
>    > Copying config b23a8f6569 done
>    > Writing manifest to image destination
>    > Storing signatures
>    > 2021/01/08 11:25:12  info unpack layer: sha256:0e03bdcc26d7a9a57ef3b6f1bf1a210cff6239bff7c8cac72435984032851689
>    > INFO:    Creating SIF file...
>    > WARNING: passwd file doesn't exist in container, not updating
>    > WARNING: group file doesn't exist in container, not updating
>    >
>    > Hello from Docker!
>    > This message shows that your installation appears to be working correctly.
>    > ...
>    > ```
>    {: .code-out}
{: .hands_on}

## Configure Galaxy to use Singularity

Now, we will configure Galaxy to run tools using Singularity containers, which will be automatically fetched from [the BioContainers repository](https://quay.io/organization/biocontainers).

> ### {% icon hands_on %} Hands-on: Configure Galaxy to use Singularity
>
> 1. Edit the `group_vars/galaxyservers.yml` file and add a `dependency_resolvers_config_file` entry and a corresponding `galaxy_config_files` entry:
>
>    {% raw %}
>    ```yaml
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -28,6 +28,7 @@ miniconda_manage_dependencies: false
>
>     galaxy_config:
>       galaxy:
>    +    dependency_resolvers_config_file: "{{ galaxy_config_dir }}/dependency_resolvers_conf.xml"
>         brand: "🧬🔬🚀"
>         admin_users: admin@example.org
>         database_connection: "postgresql:///galaxy?host=/var/run/postgresql"
>    @@ -65,6 +66,11 @@ galaxy_config_templates:
>       - src: templates/galaxy/config/job_conf.xml.j2
>         dest: "{{ galaxy_config.galaxy.job_config_file }}"
>
>    +galaxy_config_files:
>    +- src: files/galaxy/config/dependency_resolvers_conf.xml
>    +  dest: "{{ galaxy_config_dir }}/dependency_resolvers_conf.xml"
>    ```
>    {% endraw %}
>
> 2. Create the `files/galaxy/config` directory if it doesn't exist:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > mkdir -p files/galaxy/config
>    > ```
>    {: .code-in}
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
>    -        <destination id="local_destination" runner="local_plugin"/>
>    +        <destination id="local_destination" runner="local_plugin">
>    +            <param id="singularity_enabled">true</param>
>    +            <!-- Ensuring a consistent collation environment is good for reproducibility. -->
>    +            <env id="LC_ALL">C</env>
>    +            <!-- The cache directory holds the docker containers that get converted. -->
>    +            <env id="SINGULARITY_CACHEDIR">/tmp/singularity</env>
>    +            <!-- Singularity uses a temporary directory to build the squashfs filesystem. -->
>    +            <env id="SINGULARITY_TMPDIR">/tmp</env>
>    +        </destination>
>         </destinations>
>         <tools>
>         </tools>
>    ```
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



