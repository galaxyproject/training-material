---
layout: tutorial_hands_on
redirect_from:
- /topics/admin/tutorials/singularity/tutorial

title: "Use Apptainer containers for running Galaxy jobs"
zenodo_link: ""
questions:
objectives:
  - Configure your Galaxy to use Apptainer and BioContainers for running jobs
time_estimation: "1h"
key_points:
contributors:
  - torfinnnome
  - mvdbeek
  - bernt-matthias
  - hexylena
  - mira-miracoli
subtopic: jobs
tags:
  - jobs
  - ansible
  - git-gat
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

In this tutorial you will learn how to configure Galaxy to run jobs using [Apptainer](https://apptainer.org) containers provided by the [BioContainers](https://biocontainers.pro/) community.
Make sure to read the documentation on [container in Galaxy](https://docs.galaxyproject.org/en/master/admin/special_topics/mulled_containers.html) and in particular [container resolvers in Galaxy](https://docs.galaxyproject.org/en/master/admin/container_resolvers.html).


## Background

> BioContainers is a community-driven project that provides the infrastructure and basic guidelines to create, manage and distribute bioinformatics packages (e.g Conda) and containers (e.g Docker, Apptainer). BioContainers is based on the popular frameworks Conda, Docker and Apptainer.
>
{: .quote cite="https://biocontainers-edu.readthedocs.io/en/latest/what_is_biocontainers.html"}

Apptainer is an alternative to Docker that is much friendlier for HPCs

> <comment-title>Apptainer, Singularity, SingularityCE?</comment-title>
>
> Name                  | Singularity                                          | SingularityCE                                                  | Apptainer
> --------------------- | ---------------------------------------------------- | ---------------------                                          | ----------------
> Origin                | Original name, used by the project until 2021        | Name of Sylabs' Fork (CE for <b>C</b>ommunity <b>E</b>dition) | Name change when the project joined the Linux Foundation
> Status                | renamed in 2021                                      | currently active                                               | currently active
> RPM Package available | discontinued                                         | ❌                                                             | ✅
> CLI name              | `singularity`                                        | `singularity`                                                  | `apptainer` or `singularity`
>
> Many people still know Apptainer under its former name, Singularity.
> Singularity was forked in 2021, with the non-commercial fork being renamed to Apptainer and joining the Linux Foundation.
> Sylabs maintains their own free and open source fork, 'SingularityCE' (for <u>C</u>ommunity <u>E</u>dition), as well as a commercial version, 'SingularityPRO'.
>
> We will use the Name Apptainer for our training material, because most rpm packages are now named Apptainer.
{: .comment}

> Apptainer is a container platform. It allows you to create and run containers that package up pieces of software in a way that is portable and reproducible.
{: .quote cite="https://sylabs.io/guides/3.7/user-guide/introduction.html"}

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="apptainer" %}

# Installing Apptainer

First, we will install Apptainer using Ansible. Since there is a package available for major Linux distros now, we could simply install it in the pre tasks. However, we would have to enable additional repos, so we decided to create a role for this:

> <hands-on-title>Installing Apptainer with Ansible</hands-on-title>
>
> 1. In your working directory, add the Apptainer role to your `requirements.yml` file:
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -20,3 +20,6 @@
>     # CVMFS Support
>     - src: galaxyproject.cvmfs
>       version: 0.2.21
>    +# Singularity/Apptainer
>    +- src: usegalaxy_eu.apptainer
>    +  version: 0.0.1
>    {% endraw %}
>    ```
>    {: data-commit="Add Apptainer ansible roles"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Install the requirements with `ansible-galaxy`:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 4. Add the new roles to your `galaxy.yml` playbook, before the Galaxy server itself. We'll do this because it's a dependency of Galaxy to run, so it needs to be there before Galaxy starts.
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -31,6 +31,7 @@
>           when: ansible_os_family == 'Debian'
>       roles:
>         - galaxyproject.tusd
>    +    - usegalaxy_eu.apptainer
>         - galaxyproject.galaxy
>         - role: galaxyproject.miniconda
>           become: true
>    {% endraw %}
>    ```
>    {: data-commit="Add the roles to the playbook"}
>
> 5. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 6. Apptainer should now be installed on your Galaxy server. You can test this by connecting
>    to your server and run the following command:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > apptainer run docker://hello-world
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
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

## Configure Galaxy to use Apptainer

Now, we will configure Galaxy to run tools using Apptainer containers, which will be automatically fetched from [the BioContainers repository](https://quay.io/organization/biocontainers).  

> <warning-title>Galaxy Calls It Singularity, not Apptainer</warning-title>
> Galaxy still uses `singularity` in most variables, they will be replaced successively.
{: .warning}

> <hands-on-title>Configure Galaxy to use Apptainer</hands-on-title>
>
> 1. Edit the `group_vars/galaxyservers.yml` file and add a `dependency_resolvers_config_file` entry and a corresponding `galaxy_config_templates` entry:
>
> 1. Edit the `group_vars/galaxyservers.yml` file and add a `container_resolvers_config_file` entry and a corresponding `galaxy_config_templates` entry:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -72,6 +72,9 @@ galaxy_config:
>         tus_upload_store: "{{ galaxy_tus_upload_store }}"
>         # CVMFS
>         tool_data_table_config_path: /cvmfs/data.galaxyproject.org/byhand/location/tool_data_table_conf.xml,/cvmfs/data.galaxyproject.org/managed/location/tool_data_table_conf.xml
>    +    # Tool Dependencies
>    +    dependency_resolvers_config_file: "{{ galaxy_config_dir }}/dependency_resolvers_conf.xml"
>    +    container_resolvers_config_file: "{{ galaxy_config_dir }}/container_resolvers_conf.yml"
>       gravity:
>         process_manager: systemd
>         galaxy_root: "{{ galaxy_root }}/server"
>    @@ -111,6 +114,12 @@ galaxy_config_files:
>       - src: files/galaxy/themes.yml
>         dest: "{{ galaxy_config.galaxy.themes_config_file }}"
>     
>    +galaxy_config_templates:
>    +  - src: templates/galaxy/config/container_resolvers_conf.yml.j2
>    +    dest: "{{ galaxy_config.galaxy.container_resolvers_config_file }}"
>    +  - src: templates/galaxy/config/dependency_resolvers_conf.xml
>    +    dest: "{{ galaxy_config.galaxy.dependency_resolvers_config_file }}"
>    +
>     galaxy_extra_dirs:
>       - /data
>     
>    {% endraw %}
>    ```
>    {: data-commit="Configure the container and dependency resolvers"}
>
> 2. Create the `templates/galaxy/config` directory if it doesn't exist:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > mkdir -p templates/galaxy/config
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Create the new file `templates/galaxy/config/dependency_resolvers_conf.xml`. This will not enable any dependency resolvers like the legacy toolshed packages or Galaxy packages, and instead everything will be resolved through Apptainer.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/galaxy/config/dependency_resolvers_conf.xml
>    @@ -0,0 +1,2 @@
>    +<dependency_resolvers>
>    +</dependency_resolvers>
>    {% endraw %}
>    ```
>    {: data-commit="Configure the dependency resolvers"}
>
> 3. Create the new file `templates/galaxy/config/container_resolvers_conf.yml.j2`, this specifies the order in which to attempt container resolution.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/galaxy/config/container_resolvers_conf.yml.j2
>    @@ -0,0 +1,11 @@
>    +
>    +- type: cached_explicit_singularity
>    +  cache_directory: "{{ galaxy_mutable_data_dir }}/cache/singularity/explicit/"
>    +- type: cached_mulled_singularity
>    +  cache_directory: "{{ galaxy_mutable_data_dir }}/cache/singularity/mulled/"
>    +- type: mulled_singularity
>    +  auto_install: False
>    +  cache_directory: "{{ galaxy_mutable_data_dir }}/cache/singularity/mulled/"
>    +- type: build_mulled_singularity
>    +  auto_install: False
>    +  cache_directory: "{{ galaxy_mutable_data_dir }}/cache/singularity/built/"
>    {% endraw %}
>    ```
>    {: data-commit="Configure the container resolver"}
>
> 3. Now, we want to make Galaxy run jobs using Apptainer. Modify the file `group_vars/galaxyservers.yml`, by adding the `singularity_enabled` parameter:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -21,11 +21,24 @@ galaxy_job_config:
>       handling:
>         assign: ['db-skip-locked']
>       execution:
>    -    default: local_env
>    +    default: singularity
>         environments:
>           local_env:
>             runner: local_runner
>             tmp_dir: true
>    +      singularity:
>    +        runner: local_runner
>    +        singularity_enabled: true
>    +        env:
>    +        # Ensuring a consistent collation environment is good for reproducibility.
>    +        - name: LC_ALL
>    +          value: C
>    +        # The cache directory holds the docker containers that get converted
>    +        - name: APPTAINER_CACHEDIR
>    +          value: /tmp/singularity
>    +        # Apptainer uses a temporary directory to build the squashfs filesystem
>    +        - name: APPTAINER_TMPDIR
>    +          value: /tmp
>       tools:
>         - class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>           environment: local_env
>    {% endraw %}
>    ```
>    {: data-commit="Update the job_conf with singularity destination"}
>
> 4. Re-run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 5. In your Galaxy admin interface, install the minimap2 tool.
>
>    - Login to Galaxy as the admin user
>    - Click the "admin" menu at the top
>    - Under "Tool Management" on the left select "Install and Uninstall"
>    - search for `minimap2` and install the latest version with the Target Section "Mapping"
>
>    ![Screenshot of the install interface, minimap2 is entered in the search box and the latest revision shows it is currently cloning](../../images/install-minimap2.png)
>
> 6. Upload the following fasta file
>
>    ```
>    >testing
>    GATTACAGATHISISJUSTATESTGATTACA
>    ```
>
> 2. **Map with minimap2** {% icon tool %} with the following parameters
>    - *"Will you select a reference genome from your history or use a built-in index"*: `Use a genome from history and build index`
>    - *"Use the following dataset as the reference sequence"*: The fasta file you uploaded
>    - *"Single or Paired-end reads"*: `Single`
>        - {% icon param-file %} *"Select fastq dataset"*: The fasta file you uploaded
>
>    Your job should be executed using Apptainer with a BioContainer! You can watch the logs of Galaxy to see this happening.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > journalctl -f
>    > ```
>    {: .code-in}
>
>    > <code-out-title></code-out-title>
>    > ```
>    > gunicorn[1190010]: galaxy.tool_util.deps.containers INFO 2021-01-08 13:37:30,342 [p:1190010,w:0,m:2] [LocalRunner.work_thread-1] Checking with container resolver [MulledSingularityContainerResolver[namespace=biocontainers]] found description [ContainerDescription[identifier=docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:e1ea28074233d7265a5dc2111d6e55130dff5653-0,type=singularity]]
>    > gunicorn[1190010]: galaxy.jobs.command_factory INFO 2021-01-08 13:37:30,418 [p:1190010,w:0,m:2] [LocalRunner.work_thread-1] Built script [/srv/galaxy/jobs/000/23/tool_script.sh] for tool command [minimap2 --version > /srv/galaxy/jobs/000/23/outputs/COMMAND_VERSION 2>&1; ln -f -s '/data/000/dataset_22.dat' reference.fa && minimap2           -t ${GALAXY_SLOTS:-4} reference.fa '/data/000/dataset_22.dat' -a | samtools sort -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O BAM -o '/data/000/dataset_23.dat' > '/data/000/dataset_23.dat']
>    > gunicorn[1190010]: galaxy.jobs.runners DEBUG 2021-01-08 13:37:30,441 [p:1190010,w:0,m:2] [LocalRunner.work_thread-1] (23) command is: mkdir -p working outputs configs
>    > gunicorn[1190010]: if [ -d _working ]; then
>    > gunicorn[1190010]:     rm -rf working/ outputs/ configs/; cp -R _working working; cp -R _outputs outputs; cp -R _configs configs
>    > gunicorn[1190010]: else
>    > gunicorn[1190010]:     cp -R working _working; cp -R outputs _outputs; cp -R configs _configs
>    > gunicorn[1190010]: fi
>    > gunicorn[1190010]: cd working; SINGULARITYENV_GALAXY_SLOTS=$GALAXY_SLOTS SINGULARITYENV_HOME=$HOME SINGULARITYENV__GALAXY_JOB_HOME_DIR=$_GALAXY_JOB_HOME_DIR SINGULARITYENV__GALAXY_JOB_TMP_DIR=$_GALAXY_JOB_TMP_DIR SINGULARITYENV_TMPDIR=$TMPDIR SINGULARITYENV_TMP=$TMP SINGULARITYENV_TEMP=$TEMP singularity -s exec -B /srv/galaxy/server:/srv/galaxy/server:ro -B /srv/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/iuc/minimap2/8c6cd2650d1f/minimap2:/srv/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/iuc/minimap2/8c6cd2650d1f/minimap2:ro -B /srv/galaxy/jobs/000/23:/srv/galaxy/jobs/000/23 -B /srv/galaxy/jobs/000/23/outputs:/srv/galaxy/jobs/000/23/outputs -B /srv/galaxy/jobs/000/23/configs:/srv/galaxy/jobs/000/23/configs -B /srv/galaxy/jobs/000/23/working:/srv/galaxy/jobs/000/23/working -B /data:/data -B /srv/galaxy/var/tool-data:/srv/galaxy/var/tool-data:ro -B /srv/galaxy/var/tool-data:/srv/galaxy/var/tool-data:ro --home $HOME:$HOME docker://quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:e1ea28074233d7265a5dc2111d6e55130dff5653-0 /bin/bash /srv/galaxy/jobs/000/23/tool_script.sh > ../outputs/tool_stdout 2> ../outputs/tool_stderr; return_code=$?; cd '/srv/galaxy/jobs/000/23';
>    > ```
>    {: .code-out.code-max-300}
>
{: .hands_on}

> ```bash
> 1-run-minimap2.sh
> ```
> {: data-test="true"}
{: .hidden}

> <comment-title>Manage dependencies menu</comment-title>
> You can manually pull one or many containers for tools in the admin menu. Go to the admin menu, click Manage Dependencies and select the Containers tab. This will list all tools, their dependencies and whether containers are already pulled or can be pulled on demand.
>
> When a container has been resolved through Apptainer, you'll see something like this:
> ![Image of a table entry with minimap2 having requirements minimap2+singularity, a resolved column with a green checkmark next to via singularity, the resolver is mulled_singularity, and a container column with a path to /srv/galaxy/var/cache/singularity/mulled and some long hash.](../../images/singularity-resolved.png)
{: .comment}

> <tip-title>Apptainer, Conda, something else?</tip-title>
> We often hear
>
> > What would be the best practice, use Conda or Apptainer?
> {: .quote}
>
> Many of us are moving towards Apptainer. Conda environments can resolve differently if they were installed at different times, which isn't great for reproducibility. Apptainer images are never updated after generation which makes them fantastic. Also the isolation that's there by default is an incredible improvement for less-trustworthy binaries.
{: .tip}

> <tip-title>Does Apptainer fix issues with Conda dependencies resolution?</tip-title>
> Yes and no. Apptainer images are built from Conda environments. Only now you are no longer responsible for solving the conda environment, or ensuring that all of the dependencies are installed. The Galaxy project uses a system called "mulling" to bring together multiple Conda dependencies together in a single environment, and Apptainer images are produced for these dependencies as well. That said, complex or unresolvable Conda environments are not solved by Apptainer, because Apptainer is really just packaging Conda's environment into a single binary file.
{: .tip}


> <tip-title>Gateway Time-out (504) in Dependencies view</tip-title>
> When you open "Admin -> Tool Management -> Manage Dependencies -> Containers", it sometimes shows "Gateway Time-out (504)"
>
> Resolving all dependencies for all tools can take a bit, you can increase your timeout with the `proxy_read_timeout` setting in `templates/nginx/galaxy.j2`
{:.tip}

> <tip-title>Resolution is "unresolved"</tip-title>
> In "Admin -> Tool Management -> Manage Dependencies -> Dependencies", the Resolution for minimap2 @ 2.24 (as well as samtools @1.14) is "unresolved". How can I resolve this issue?
>
> Because our training uses containers for resolution it is expected that the non-container dependencies show as "unresolved". There is not currently a view which indicates if the containers have been resolved.
{: .tip}





<!--
## Use Apptainer containers from CVMFS

Galaxy can be configured to use pre-made Apptainer containers available from /cvmfs/singularity.galaxyproject.org/.
In order to do so, you will first need to set up CVMFS by doing the [CVMFS]({{ site.baseurl }}/topics/admin/tutorials/cvmfs/tutorial.html) tutorial.
After finishing the CVMFS tutorial, come back, and do this hands-on.

> <hands-on-title>Optional: Configure Galaxy to use Apptainer containers from CVMFS</hands-on-title>
>
> 1. Edit the `group_vars/galaxyservers.yml` file and add `container_resolvers_config_file` and `galaxy_singularity_images_cvmfs_path`:
>    {% raw %}
>    ```yaml
>    galaxy_singularity_images_cvmfs_path: "/cvmfs/singularity.galaxyproject.org/all/"
>    galaxy_config:
>      galaxy:
>        ...
>        container_resolvers_config_file: "{{ galaxy_config_dir }}/container_resolvers_conf.yml"
>    ```
>    {% endraw %}
>
> 2. Also in `group_vars/galaxyservers.yml`, add a `galaxy_config_templates` entry:
>    {% raw %}
>    ```yaml
>    galaxy_config_templates:
>      - src: templates/galaxy/config/container_resolvers_conf.yml.j2
>        dest: "{{ galaxy_config_dir }}/container_resolvers_conf.yml"
>    ```
>    {% endraw %}
>
> 3. Create the new file `templates/galaxy/config/container_resolvers_conf.yml.j2`:
>    {% raw %}
>    ```xml
>    <containers_resolvers>
>        <explicit_singularity />
>        <cached_mulled_singularity cache_directory="{{ galaxy_singularity_images_cvmfs_path }}"/>
>    </containers_resolvers>
>    ```
>    {% endraw %}
>
> 4. Re-run the playbook (`ansible-playbook galaxy.yml`)
>
> 5. Rerun your minimap2 job from the previous section. This time it should use a container served by CVMFS.
>
{: .hands_on}

> ```bash
> 1-run-minimap2.sh
> ```
> {: data-test="true"}
{: .hidden}

-->

{% snippet topics/admin/faqs/git-commit.md page=page %}

{% snippet topics/admin/faqs/missed-something.md step=6 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="apptainer" %}
