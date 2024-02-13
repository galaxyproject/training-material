---
layout: tutorial_hands_on

title: "Mapping Jobs to Destinations using TPV"
questions:
  - How can I configure job dependent resources, like cores, memory for my DRM?
  - How can I map jobs to resources and destinations
objectives:
  - Know how to map tools to job destinations
  - Be able to use the dynamic job runner to make arbitrary destination mappings
  - Understand the job resource selector config and dynamic rule creation
  - The various ways in which tools can be mapped to destinations, both statically and dynamically
  - How to write a dynamic tool destination (DTD)
  - How to write a dynamic python function destination
  - How to use the job resource parameter selection feature
time_estimation: "2h"
key_points:
  - Dynamic Tool Destinations are a convenient way to map
  - Job resource parameters can allow you to give your users control over job resource requirements, if they are knowledgeable about the tools and compute resources available to them.
contributions:
  authorship:
    - natefoo
    - bgruening
    - nuwang
    - mira-miracoli
  editing: # And reviewing
    - hexylena
    - afgane
  funding: []
  testing:
    - cat-bro
    - edmontosaurus
tags:
  - jobs
  - git-gat
subtopic: jobs
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - connect-to-compute-cluster
abbreviations:
  TPV: Total Perspective Vortex
---


This tutorial heavily builds on the [Connecting Galaxy to a compute cluster]({% link topics/admin/tutorials/connect-to-compute-cluster/tutorial.md %}) and it's expected you have completed this tutorial first.

Now that you have a working scheduler, we will start configuring which jobs are sent to which destinations.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="job-destinations" %}

# Mapping jobs to destinations

In order to run jobs in Galaxy, you need to assign them to a resource manager that can handle the task. This involves specifying the appropriate amount of memory and CPU cores. For production installations, the jobs
must be routed to a resource manager like SLURM, HTCondor, or Pulsar. Some tools may need specific resources such as GPUs or multi-core machines to work efficiently.

Sometimes, your available resources are spread out across multiple locations and resource managers. In such cases, you need a way to route your jobs to the appropriate location. Galaxy offers several methods for
routing jobs, ranging from simple static mappings to custom Python functions via [dynamic job destinations](https://docs.galaxyproject.org/en/latest/admin/jobs.html#dynamic-destination-mapping).

Recently, the Galaxy project has introduced a library named [{TPV}](https://total-perspective-vortex.readthedocs.io/) to simplify this process. {TPV} provides a admin-friendly YAML configuration that works for most scenarios.
For more complex cases, TPV also allows you to embed Python code into the configuration YAML file and implement fine-grained control over jobs. 

Lastly, TPV offers a shared global database of default resource requirements (more below). By leveraging this database, admins don't have to figure out the requirements for each tool
separately.

## Writing a testing tool

To demonstrate a real-life scenario and {TPV}'s role in it, let's plan on setting up a configuration where the VM designated for training jobs doesn't run real jobs and hence doesn't get overloaded. To start, we'll create a "testing" tool that we'll use in our configuration. This testing tool can run quickly, and without overloading our small machines.

> <hands-on-title>Deploying a Tool</hands-on-title>
>
> 1. Create the directory `files/galaxy/tools/` if it doesn't exist and edit a new file in `files/galaxy/tools/testing.xml` with the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/tools/testing.xml
>    @@ -0,0 +1,11 @@
>    +<tool id="testing" name="Testing Tool">
>    +    <command>
>    +        <![CDATA[echo "Running with '\${GALAXY_SLOTS:-1}' threads" > "$output1"]]>
>    +    </command>
>    +    <inputs>
>    +        <param name="input1" type="data" format="txt" label="Input Dataset"/>
>    +    </inputs>
>    +    <outputs>
>    +        <data name="output1" format="txt" />
>    +    </outputs>
>    +</tool>
>    {% endraw %}
>    ```
>    {: data-commit="Add testing tool"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Add the tool to the Galaxy group variables under the new item `galaxy_local_tools` :
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -152,6 +152,9 @@ galaxy_config_templates:
>     galaxy_extra_dirs:
>       - /data
>     
>    +galaxy_local_tools:
>    +- testing.xml
>    +
>     # Certbot
>     certbot_auto_renew_hour: "{{ 23 |random(seed=inventory_hostname)  }}"
>     certbot_auto_renew_minute: "{{ 59 |random(seed=inventory_hostname)  }}"
>    {% endraw %}
>    ```
>    {: data-commit="Deploy testing tool"}
>
> 3. Run the Galaxy playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 4. Reload Galaxy in your browser and the new tool should now appear in the tool panel. If you have not already created a dataset in your history, upload a random text dataset. Once you have a dataset, click the tool's name in the tool panel, then click Execute.
>
>    > <question-title></question-title>
>    >
>    > What is the tool's output?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > ```
>    > > Running with '1' threads
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
{: .hands_on}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}

Of course, this tool doesn't actually *use* the allocated number of cores. In a real tool, you would call the tools's underlying command with whatever flag that tool provides to control the number of threads or processes it starts, such as `samtools sort -@ \${GALAXY_SLOTS:-1}`.

## Safeguard: TPV Linting

If we want to change something in production, it is always a good idea to have a safeguard in place. In our case, we would like to check the {TPV} configuration files for syntax errors, so nothing will break when we deploy broken yaml files or change them quickly on the server. TPV-lint-and-copy works with two separate locations:
- one where you can safely edit your files
- and the actual production config directory that galaxy reads.

Once you are done with your changes, you can run the script and it will automatically lint and copy over the files, if they are correct *and* mentioned in your job_conf.yml file, or in your `group_vars/galaxyservers.yml` inline `job_conf`.
And of course, Galaxy has an Ansible Role for that.

> <hands-on-title>Adding automated TPV-lind-and-copy-script</hands-on-title>
>
> 1. Add the role to your `requirements.yml`.
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -28,3 +28,6 @@
>       version: 0.0.3
>     - src: galaxyproject.slurm
>       version: 1.0.2
>    +# TPV Linting
>    +- name: usegalaxy_eu.tpv_auto_lint
>    +  version: 0.4.3
>    {% endraw %}
>    ```
>    {: data-commit="Add tpv-auto-lint to requirements"}
>
> 2. Install the missing role
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Change your `group_vars/galaxyservers.yml`. We need to create a new directory where the TPV configs will be stored after linting, and add that directory name as variable for the role. The default name is 'TPV_DO_NOT_TOUCH' for extra safety 😉. If you want a different name, you need to change the `tpv_config_dir_name` variable, too. We also need to create a directory, `tpv_mutable_dir` (a role default variable), where TPV configs are copied before linting.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -135,6 +135,8 @@ galaxy_config:
>               - job-handlers
>               - workflow-schedulers
>     
>    +galaxy_job_config_file: "{{ galaxy_config_dir }}/galaxy.yml"
>    +
>     galaxy_config_files_public:
>       - src: files/galaxy/welcome.html
>         dest: "{{ galaxy_mutable_config_dir }}/welcome.html"
>    @@ -151,6 +153,11 @@ galaxy_config_templates:
>     
>     galaxy_extra_dirs:
>       - /data
>    +  - "{{ galaxy_config_dir }}/{{ tpv_config_dir_name }}"
>    +
>    +galaxy_extra_privsep_dirs:
>    +  - "{{ tpv_mutable_dir }}"
>    +tpv_privsep: true
>     
>     galaxy_local_tools:
>     - testing.xml
>    {% endraw %}
>    ```
>    {: data-commit="Add TPV config dir"}
> 4. Add the role to your `galaxy.yml` playbook.
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -38,6 +38,7 @@
>         - galaxyproject.slurm
>         - usegalaxy_eu.apptainer
>         - galaxyproject.galaxy
>    +    - usegalaxy_eu.tpv_auto_lint
>         - role: galaxyproject.miniconda
>           become: true
>           become_user: "{{ galaxy_user_name }}"
>    {% endraw %}
>    ```
>    {: data-commit="Add tpv-auto-lint role"}
> 5. At this point we won't run the modified playbook just yet. Because TPV itself has not yet been installed,
>    the tpv_auto_lint role would fail at this point. So first, we'll have to install and configure TPV itself before the
>    linter can work.
>
{: .hands_on}
## Configuring TPV

We want our tool to run with more than one core. To do this, we need to instruct Slurm to allocate more cores for this job. First however, we need to configure Galaxy to use TPV.

> <hands-on-title>Adding TPV to your job configuration</hands-on-title>
>
> 1. Edit `group_vars/galaxyservers.yml` and add the following destination under your `job_config` section to route all jobs to TPV.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -24,34 +24,18 @@ galaxy_job_config:
>       handling:
>         assign: ['db-skip-locked']
>       execution:
>    -    default: slurm
>    +    default: tpv_dispatcher
>         environments:
>           local_env:
>             runner: local_runner
>             tmp_dir: true
>    -      slurm:
>    -        runner: slurm
>    -        singularity_enabled: true
>    -        env:
>    -        - name: LC_ALL
>    -          value: C
>    -        - name: APPTAINER_CACHEDIR
>    -          value: /tmp/singularity
>    -        - name: APPTAINER_TMPDIR
>    -          value: /tmp
>    -      singularity:
>    -        runner: local_runner
>    -        singularity_enabled: true
>    -        env:
>    -        # Ensuring a consistent collation environment is good for reproducibility.
>    -        - name: LC_ALL
>    -          value: C
>    -        # The cache directory holds the docker containers that get converted
>    -        - name: APPTAINER_CACHEDIR
>    -          value: /tmp/singularity
>    -        # Apptainer uses a temporary directory to build the squashfs filesystem
>    -        - name: APPTAINER_TMPDIR
>    -          value: /tmp
>    +      tpv_dispatcher:
>    +        runner: dynamic
>    +        type: python
>    +        function: map_tool_to_destination
>    +        rules_module: tpv.rules
>    +        tpv_config_files:
>    +          - "{{ tpv_config_dir }}/tpv_rules_local.yml"
>       tools:
>         - class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>           environment: local_env
>    @@ -144,6 +128,8 @@ galaxy_config_files_public:
>     galaxy_config_files:
>       - src: files/galaxy/themes.yml
>         dest: "{{ galaxy_config.galaxy.themes_config_file }}"
>    +  - src: files/galaxy/config/tpv_rules_local.yml
>    +    dest: "{{ tpv_mutable_dir }}/tpv_rules_local.yml"
>     
>     galaxy_config_templates:
>       - src: templates/galaxy/config/container_resolvers_conf.yml.j2
>    {% endraw %}
>    ```
>    {: data-commit="Add TPV to job config"}
>
>  Note that we set the default execution environment to the tpv_dispatcher, added the tpv_dispatcher itself as a dynamic runner, and removed all other destinations.
>  Adding TPV as a runner will cause Galaxy to automatically install the `total-perspective-vortex` package on startup as a conditional dependency.
>  Finally, we added a new config file named `tpv_rules_local.yml`, which we will create next.
>
> 2. Create a new file named `tpv_rules_local.yml` in the `files/galaxy/config/` folder of your ansible playbook, so that it is copied to the config folder on the target.
>    The file should contain the following content:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -0,0 +1,30 @@
>    +tools:
>    +  .*testing.*:
>    +    cores: 2
>    +    mem: cores * 4
>    +
>    +destinations:
>    +  local_env:
>    +    runner: local_runner
>    +    max_accepted_cores: 1
>    +    params:
>    +      tmp_dir: true
>    +  singularity:
>    +    runner: local_runner
>    +    max_accepted_cores: 1
>    +    params:
>    +      singularity_enabled: true
>    +    env:
>    +      # Ensuring a consistent collation environment is good for reproducibility.
>    +      LC_ALL: C
>    +      # The cache directory holds the docker containers that get converted
>    +      APPTAINER_CACHEDIR: /tmp/singularity
>    +      # Singularity uses a temporary directory to build the squashfs filesystem
>    +      APPTAINER_TMPDIR: /tmp
>    +  slurm:
>    +    inherits: singularity
>    +    runner: slurm
>    +    max_accepted_cores: 16
>    +    params:
>    +      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    +
>    {% endraw %}
>    ```
>    {: data-commit="Add TPV local rules"}
>
>    In this TPV config, we have specified that the testing tool should use `2` cores. Memory has been defined as an expression and should be 4 times as much as cores, which, in this case, is 8GB.
>    Note that the tool id is matched via a regular expression against the full tool id. For example, a full tool id for hisat may look like: `toolshed.g2.bx.psu.edu/repos/iuc/hisat2/hisat2/2.1.0+galaxy7`
>    This enables complex matching, including matching against specific versions of tools.
>
>    Destinations must also be defined in TPV itself. Importantly, note that any destinations defined in the job conf are ignored by TPV. Therefore, we have moved all destinations from the job conf to TPV. In addition, we have removed some
>    redundancy by using the "inherits" clause in the `slurm` destination. This means that slurm will inherit all of the settings defined for singularity, but selectively override some settings. We have additionally
>    defined the `native_specification` param for SLURM, which is what SLURM uses to allocate resources per job. Note the use of the `{cores}`
>    parameter within the native specification, which TPV will replace at runtime with the value of cores assigned to the tool.
>
>    Finally, we have also defined a new property named `max_accepted_cores`, which is the maximum amount of cores this destination will accept. Since the testing tool requests 2 cores, but only the `slurm`
>    destination is able to accept jobs greater than 1 core, TPV will automatically route the job to the best matching destination, in this case, slurm.
>
> 3. Run the Galaxy playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 4. Click the rerun button on the last history item, or click **Testing Tool** in the tool panel, and then click the tool's Run Tool button.
>
>    > <question-title></question-title>
>    >
>    > What is the tool's output?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > ```
>    > > Running with '2' threads
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}

> ```bash
> 2.sh
> ```
> {: data-test="true"}
{: .hidden}


### Configuring defaults

Now that we've configured the resource requirements for a single tool, let's see how we can configure defaults for all tools, and reuse those defaults to reduce repetition.


> <hands-on-title>Configuring defaults and inheritance</hands-on-title>
>
> 1. Edit your `files/galaxy/config/tpv_rules_local.yml` and add the following settings.
>
>    {% raw %}
>    ```diff
>    --- a/files/galaxy/config/tpv_rules_local.yml
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -1,4 +1,11 @@
>    +global:
>    +  default_inherits: default
>    +
>     tools:
>    +  default:
>    +    abstract: true
>    +    cores: 1
>    +    mem: cores * 4
>       .*testing.*:
>         cores: 2
>         mem: cores * 4
>    @@ -27,4 +34,3 @@ destinations:
>         max_accepted_cores: 16
>         params:
>           native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    -
>    {% endraw %}
>    ```
>    {: data-commit="Add TPV default inherits"}
>
>    We have defined a `global` section specifying that all tools and destinations should inherit from a specified `default`. We have then defined a tool named `default`, whose properties
>    are implicitly inherited by all tools at runtime. This means that our `testing` tool will also inherit from this default tool, but it explicitly overrides cores.
>    We can also explicitly specify an `inherits` clause if we wish to extend a specific tool or destination, as previously shown in the destinations section.
>
> 2. Run the Galaxy playbook. When the new `tpv_rules_local.yml` is copied, TPV will automatically pickup the changes without requiring a restart of Galaxy.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on}

### TPV reference documentation

Please see [TPV's dedicated documentation](https://total-perspective-vortex.readthedocs.io/en/latest/) for more information.


# Configuring the TPV shared database

The Galaxy Project maintains a [shared database of TPV rules](https://github.com/galaxyproject/tpv-shared-database) so that admins do not have to independently rediscover ideal resource allocations for specific tools. These rules are based
on settings that have worked well in the usegalaxy.* federation. The rule file can simply be imported directly, with local overrides applied on top.

> <hands-on-title>Adding the TPV shared database to job_conf</hands-on-title>
>
> 1. Edit `group_vars/galaxyservers.yml` and add the location of the TPV shared rule file to the `tpv_dispatcher` destination.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -35,6 +35,7 @@ galaxy_job_config:
>             function: map_tool_to_destination
>             rules_module: tpv.rules
>             tpv_config_files:
>    +          - https://raw.githubusercontent.com/galaxyproject/tpv-shared-database/main/tools.yml
>               - "{{ tpv_config_dir }}/tpv_rules_local.yml"
>       tools:
>         - class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>    {% endraw %}
>    ```
>    {: data-commit="Importing TPV shared database via job conf"}
>
>    Note how TPV allows the file to be imported directly via its http url. As many local and remote rule files as necessary can be combined, with rule files specified later overriding
>    any previously specified rule files. The TPV shared database does not define destinations, only cores and mem settings, as well as any required environment vars.
>    Take a look at the shared database of rules and note that some tools have very large recommended memory settings, which may or may not be available within your local cluster.
>    Nevertheless, you may still wish to execute these tools with memory adjusted to suit your cluster's capabilities. 
>
> 2. Edit your `files/galaxy/config/tpv_rules_local.yml` and make the following changes.
>
>    {% raw %}
>    ```diff
>    --- a/files/galaxy/config/tpv_rules_local.yml
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -31,6 +31,9 @@ destinations:
>       slurm:
>         inherits: singularity
>         runner: slurm
>    -    max_accepted_cores: 16
>    +    max_accepted_cores: 24
>    +    max_accepted_mem: 256
>    +    max_cores: 2
>    +    max_mem: 8
>         params:
>           native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    {% endraw %}
>    ```
>    {: data-commit="TPV clamp max cores and mem"}
>
>    These changes indicate that the destination will accept jobs that are up to `max_accepted_cores: 24` and `max_accepted_mem: 256`. If the tool requests resources that exceed these limits, the tool will be rejected
>    by the destination. However, once accepted, the resources will be forcibly clamped down to 2 and 8 at most because of the `max_cores` and `max_mem` clauses. (E.g. a tool requesting 24 cores would only be submitted with 16 cores at maximum.) Therefore, a trick that can be used here to support
>    job resource requirements in the shared database that are much larger than your destination can actually support, is to combine `max_accepted_cores/mem/gpus` with `max_cores/mem/gpus` to accept the job and then
>    clamp it down to a supported range. This allows even the largest resource requirement in the shared database to be accomodated.
>
>    > <comment-title>Clamping in practice</comment-title>
>    > For the purposes of this tutorial, we've clamped down from 16 cores to 2 cores, and mem from 256 to 8, which is unlikely to work in practice. In production, you will probably need to manually test
>    > any tools that exceed your cluster's capabilities, and decide whether you want those tools to run in the first place.
>    {: .comment}
>
> 3. Run the Galaxy playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on}

# Basic access controls

You may wish to apply some basic restrictions on which users are allowed to run specific tools. {TPV} accomodates user and role specific rules. In addition, TPV supports tagging of tools, users, roles and destinations. These tags
can be matched up so that only desired combinations are compatible with each other. While these mechanisms are detailed in the TPV documentation, we will choose a different problem that highlights some other capabilities
- restricting a tool so that only an admin can execute that tool.

> <hands-on-title>Using conditionals to restrict a tool to admins only</hands-on-title>
>
> 1. Edit your `files/galaxy/config/tpv_rules_local.yml` and add the following rule.
>
>    {% raw %}
>    ```diff
>    --- a/files/galaxy/config/tpv_rules_local.yml
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -9,6 +9,15 @@ tools:
>       .*testing.*:
>         cores: 2
>         mem: cores * 4
>    +    rules:
>    +      - id: admin_only_testing_tool
>    +        if: |
>    +          # Only allow the tool to be executed if the user is an admin
>    +          admin_users = app.config.admin_users
>    +          # last line in block must evaluate to a value - which determines whether the TPV if conditional matches or not
>    +          not user or user.email not in admin_users
>    +        fail: Unauthorized. Only admins can execute this tool.
>    +
>     
>     destinations:
>       local_env:
>    {% endraw %}
>    ```
>    {: data-commit="TPV admin only tool"}
>
> Note the use of the `if` rule, which allows for conditional actions to be taken in TPV. An `if` block is evaluated as a multi-line python block, and can execute arbitrary code, but the last line of the block must evaluate to a value. That value
> determines whether the `if` condtional is matched or not. If the conditional is matched, the `fail` clause is executed in this case, and the message specified in the `fail` clause is displayed to the user.
> It is similarly possible to conditionally add job parameters, modify cores/mem/gpus and take other complex actions.
>
> 2. Run the Galaxy playbook to update the TPV rules.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Try running the tool as both an admin user and a non-admin user, non-admins should not be able to run it. You can start a private browsing session to test as a non-admin, anonymous user. Anonymous users were enabled in your Galaxy configuration.
>
{: .hands_on}

> ```bash
> 3.sh
> ```
> {: data-test="true"}
{: .hidden}

# Job Resource Selectors

Certain tools can benefit from allowing users to select appropriate job resource parameters, instead of having admins decide resource allocations beforehand. For example, a user might know that a particular set of parameters and inputs to a certain tool needs a larger memory allocation than the standard amount given to that tool.
Galaxy provides functionality to have extra form elements in the tool execution form to specify these additional job resource parameters.
This of course assumes that your users are well behaved enough not to choose the maximum whenever available, although such concerns can be mitigated somewhat by the use of concurrency limits on larger memory destinations.

Such form elements can be added to tools without modifying each tool's configuration file through the use of the **job resource parameters configuration file**

> <hands-on-title>Configuring a Resource Selector</hands-on-title>
>
> 1. Create and open `templates/galaxy/config/job_resource_params_conf.xml.j2`
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/galaxy/config/job_resource_params_conf.xml.j2
>    @@ -0,0 +1,7 @@
>    +<parameters>
>    +    <param label="Cores" name="cores" type="select" help="Number of cores to run job on.">
>    +        <option value="1">1 (default)</option>
>    +        <option value="2">2</option>
>    +    </param>
>    +  <param label="Time" name="time" type="integer" size="3" min="1" max="24" value="1" help="Maximum job time in hours, 'walltime' value (1-24). Leave blank to use default value." />
>    +</parameters>
>    {% endraw %}
>    ```
>    {: data-commit="Add job resource params configuration"}
>
>    This defines two resource fields, a select box where users can choose between 1 and 2 cores, and a text entry field where users can input an integer value from 1-24 to set the walltime for a job.
>
> 2. As usual, we need to instruct Galaxy of where to find this file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -37,9 +37,17 @@ galaxy_job_config:
>             tpv_config_files:
>               - https://raw.githubusercontent.com/galaxyproject/tpv-shared-database/main/tools.yml
>               - "{{ tpv_config_dir }}/tpv_rules_local.yml"
>    +  resources:
>    +    default: default
>    +    groups:
>    +      default: []
>    +      testing: [cores, time]
>       tools:
>         - class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>           environment: local_env
>    +    - id: testing
>    +      environment: tpv_dispatcher
>    +      resources: testing
>     
>     galaxy_config:
>       galaxy:
>    @@ -56,6 +64,7 @@ galaxy_config:
>         object_store_store_by: uuid
>         id_secret: "{{ vault_id_secret }}"
>         job_config: "{{ galaxy_job_config }}" # Use the variable we defined above
>    +    job_resource_params_file: "{{ galaxy_config_dir }}/job_resource_params_conf.xml"
>         # SQL Performance
>         slow_query_log_threshold: 5
>         enable_per_request_sql_debugging: true
>    @@ -137,6 +146,8 @@ galaxy_config_templates:
>         dest: "{{ galaxy_config.galaxy.container_resolvers_config_file }}"
>       - src: templates/galaxy/config/dependency_resolvers_conf.xml
>         dest: "{{ galaxy_config.galaxy.dependency_resolvers_config_file }}"
>    +  - src: templates/galaxy/config/job_resource_params_conf.xml.j2
>    +    dest: "{{ galaxy_config.galaxy.job_resource_params_file }}"
>     
>     galaxy_extra_dirs:
>       - /data
>    {% endraw %}
>    ```
>    {: data-commit="Configure resources in job conf"}
>
>    We have added a resources section. The group ID will be used to map a tool to job resource parameters, and the text value of the `<group>` tag is a comma-separated list of `name`s from `job_resource_params_conf.xml` to include on the form of any tool that is mapped to the defined `<group>`.
>
>    We have also listed the `testing` tool as a tool which uses resource parameters. When listed here, the resource parameter selection form will be displayed.
>
>    Finally, we have specified that the `job_resource_params_conf.xml.j2` should be copied across.
>
{: .hands_on}

This will set everything up to use the function. We have:

- A set of "job resources" defined which will let the user select the number of cores and walltime.
- A job configuration which says:
    - that our testing tool should allow selection of the cores and time parameters
    - directs it to TPV's `tpv_dispatcher` destination

This is a lot but we're still missing the last piece for it to work:

## Configuring TPV to process resource parameters

Lastly, we need to write a rule in TPV that will read the value of the job resource parameter form fields and decide how to submit the job.

> <hands-on-title>Processing job resource parameters in TPV</hands-on-title>
>
> 1. Create and edit `files/galaxy/config/tpv_rules_local.yml`. Create it with the following contents:
>
>    {% raw %}
>    ```diff
>    --- a/files/galaxy/config/tpv_rules_local.yml
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -6,6 +6,8 @@ tools:
>         abstract: true
>         cores: 1
>         mem: cores * 4
>    +    params:
>    +      walltime: 8
>       .*testing.*:
>         cores: 2
>         mem: cores * 4
>    @@ -17,7 +19,13 @@ tools:
>               # last line in block must evaluate to a value - which determines whether the TPV if conditional matches or not
>               not user or user.email not in admin_users
>             fail: Unauthorized. Only admins can execute this tool.
>    -
>    +      - id: resource_params_defined
>    +        if: |
>    +          param_dict = job.get_param_values(app)
>    +          param_dict.get('__job_resource', {}).get('__job_resource__select') == 'yes'
>    +        cores: int(job.get_param_values(app)['__job_resource']['cores'])
>    +        params:
>    +           walltime: "{int(job.get_param_values(app)['__job_resource']['time'])}"
>     
>     destinations:
>       local_env:
>    @@ -45,4 +53,4 @@ destinations:
>         max_cores: 2
>         max_mem: 8
>         params:
>    -      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    +      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores} --time={params['walltime']}:00:00
>    {% endraw %}
>    ```
>    {: data-commit="process resource params in TPV"}
>
>    We define a conditional rule and check that the job_resource_params have in fact been defined by the user. If defined, we override the cores value with the one specified by the user. We also define a walltime param, setting it to 8 hours by default.
>    We also override walltime if the user has specified it.
>    It is important to note that **you are responsible for parameter validation, including the job resource selector**. This function only handles the job resource parameter fields, but it could do many other things - examine inputs, job queues, other tool parameters, etc.
>
>    Finally, we pass the walltime as part of the native specification.
>
>    > <comment-title>Rules for TPV code evaluation</comment-title>
>    > Notice how the `walltime` parameter is wrapped in braces, whereas the `cores` value isn't, yet they are both expressions. The reason for this difference is that when TPV evaluates an expression, all string fields (e.g. env, params) are evaluated as Python f-strings,
>    > while all non-string fields (integer fields like `cores` and `gpus`, float fields like `mem`, boolean fields like `if`) are evaluated as Python code-blocks. This rule enables the YAML file to be much more readable, but requires you to keep this simple rule in mind.
>    > Note also that Python code-blocks can be multi-line, and that the final line must evaluate to a value that can be assigned to the field.
>    {: .comment}
>
> 2. Run the Galaxy playbook to update the TPV rules.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Run the **Testing Tool** with various resource parameter selections
>
>    - Use default job resource parameters
>    - Specify job resource parameters:
>      - 1 core
>      - 2 cores
>      - Some value for walltime from 1-24
>
{: .hands_on}

The cores parameter can be verified from the output of the tool. The walltime can be verified with `scontrol`:

> <code-in-title>Bash</code-in-title>
> Your job number may be different.
> ```
> scontrol show job 24
> ```
{: .code-in}

> <code-out-title></code-out-title>
> Your output may look slightly different. Note that the `TimeLimit` for this job (which I gave a 12 hour time limit) was set to `12:00:00`.
> ```console
> JobId=24 JobName=g24_multi_anonymous_10_0_2_2
>    UserId=galaxy(999) GroupId=galaxy(999)
>    Priority=4294901747 Nice=0 Account=(null) QOS=(null)
>    JobState=COMPLETED Reason=None Dependency=(null)
>    Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
>    RunTime=00:00:05 TimeLimit=12:00:00 TimeMin=N/A
>    SubmitTime=2016-11-05T22:01:09 EligibleTime=2016-11-05T22:01:09
>    StartTime=2016-11-05T22:01:09 EndTime=2016-11-05T22:01:14
>    PreemptTime=None SuspendTime=None SecsPreSuspend=0
>    Partition=debug AllocNode:Sid=gat2016:1860
>    ReqNodeList=(null) ExcNodeList=(null)
>    NodeList=localhost
>    BatchHost=localhost
>    NumNodes=1 NumCPUs=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
>    TRES=cpu=1,node=1
>    Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
>    MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
>    Features=(null) Gres=(null) Reservation=(null)
>    Shared=OK Contiguous=0 Licenses=(null) Network=(null)
>    Command=(null)
>    WorkDir=/srv/galaxy/server/database/jobs/000/24
>    StdErr=/srv/galaxy/server/database/jobs/000/24/galaxy_24.e
>    StdIn=StdIn=/dev/null
>    StdOut=/srv/galaxy/server/database/jobs/000/24/galaxy_24.o
>    Power= SICP=0
> ```
{: .code-out}

{% snippet topics/admin/faqs/missed-something.md step=10 %}


# More on TPV

The goal of this tutorial is to provide a quick overview of some of the basic capabilities of TPV. However, there are numerous features that we have not convered such as:
a. Custom code blocks - execute arbitrary python code blocks, access additional context variables etc.
a. User and Role Handling - Add scheduling constraints based on the user's email or role
b. Metascheduling support - Perform advanced querying and filtering prior to choosing an appropriate destination
c. Job resubmissions - Resubmit a job if it fails for some reason
d. Linting, formatting and dry-run - Automatically format tpv rule files, catch potential syntax errors and perform a dry-run to check where a tool would get scheduled.

These features are covered in detail in the [TPV documentation](https://total-perspective-vortex.readthedocs.io/en/latest/).

## Further Reading

- The [sample dynamic tool destination config file](https://github.com/galaxyproject/galaxy/blob/dev/config/tool_destinations.yml.sample) fully describes the configuration language
- [Dynamic destination documentation](https://docs.galaxyproject.org/en/latest/admin/jobs.html)
- Job resource parameters are not as well documented as they could be, but the [sample configuration file](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/test/files/galaxy/config/job_resource_params_conf.xml) shows some of the possibilities.
- [usegalaxy.org's job_conf.yml](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/common/templates/galaxy/config/job_conf.yml.j2) is publicly available for reference.
- [usegalaxy.eu's job_conf.xml](https://github.com/usegalaxy-eu/infrastructure-playbook/search?l=YAML&q=galaxy_jobconf) is likewise (see the `group_vars/galaxy.yml` result)

{% snippet topics/admin/faqs/git-gat-path.md tutorial="job-destinations" %}
