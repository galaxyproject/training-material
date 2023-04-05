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
contributors:
  - natefoo
  - bgruening
  - hexylena
tags:
  - jobs
  - git-gat
subtopic: jobs
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - connect-to-compute-cluster
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
routing jobs, ranging from simple static mappings to custom Python functions called dynamic job destinations.

Recently, the Galaxy project has introduced a library named Total-Perspective-Vortex (TPV) to simplify this process. TPV provides a user-friendly YAML configuration that works for most scenarios.
For more complex cases, TPV allows you to embed Python code into the YAML file. Additionally, TPV shares a global database of resource requirements, so admins don't have to figure out the requirements for each tool
separately.

## Writing a testing tool

To check that resources are being allocated appropriately, We don't want to overload our training VMs trying to run real tools, so to demonstrate how to map a multicore tool to a multicore destination using TPV,
we'll create a fake tool.

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
>    @@ -126,6 +126,9 @@ galaxy_config_templates:
>       - src: templates/galaxy/config/dependency_resolvers_conf.xml
>         dest: "{{ galaxy_config.galaxy.dependency_resolvers_config_file }}"
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

## Configuring TPV

We want our tool to run with more than one core. To do this, we need to instruct Slurm to allocate more cores for this job. First however, we need to configure Galaxy to use TPV.

#TODO: Are we going to rely on tpv being installed as a conditional requirement? Or add it to the playbook as a package?


> <hands-on-title>Adding TPV to your job configuration</hands-on-title>
>
> 1. Edit your `templates/galaxy/config/job_conf.yml.j2` and add the following destination to route all jobs to TPV. # TODO: set default execution context so all jobs go to TPV
>
>    {% raw %}
>    ```diff
>    --- a/templates/galaxy/config/job_conf.yml.j2
>    +++ b/templates/galaxy/config/job_conf.yml.j2
>    @@ -20,4 +20,10 @@ execution:
>             value: /tmp/singularity
>           - name: SINGULARITY_TMPDIR
>             value: /tmp
>    +    tpv_dispatcher:
>    +      runner: dynamic
>    +      type: python
>    +      function: map_tool_to_destination
>    +      rules_module: tpv.rules
>    +      tpv_config_files:
>    +      - config/tpv_rules_local.yml
>         singularity:
>           runner: local_runner
>           singularity_enabled: true
>    {% endraw %}
>    ```
>    {: data-commit="Configure testing tool in job conf"}
>
> 2. Create a new file named `tpv_rules_local.yml` in the "files" folder of your ansible playbook, so that it is copied to the config folder on the target.
>    The file should contain the following content:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -0,0 +1,16 @@
>    +tools:
>    +  .*testing.*:
>    +    cores: 2
>    +    mem: cores * 4
>    +
>    +destinations:
>    +  slurm:
>    +    runner: slurm
>    +    params:
>    +      singularity_enabled: "true"
>    +      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    +    env:
>    +      LC_ALL: C
>    +      SINGULARITY_CACHEDIR: /tmp/singularity
>    +      SINGULARITY_TMPDIR: /tmp
>    +
>    {% endraw %}
>
>    In this TPV config, we have specified that the testing tool should use `2` cores, and that memory is an expression and should be 4 times as much as cores, or in this case, 8GB.
>    Note that the tool id is matched via a regular expression against the full tool id. For example, a full tool id for hisat may look like: toolshed.g2.bx.psu.edu/repos/iuc/hisat2/hisat2/2.1.0+galaxy7:
>    This enables complex matching, including matching against specific versions of tools.
>
>    Destinations must also be defined in TPV itself (Destinations need not be defined in job_conf.yml and are ignored by TPV). We have designated SLURM as an
>    available destination, and specified that the `native_specification` param, which is what SLURM uses to allocate resources per job. Note the use of the `{cores}`
>    parameter, which TPV will replace at runtime with the value of cores assigned to the tool.
>
>
> 3. Run the Galaxy playbook. Because we modified `job_conf.yml`, Galaxy will be restarted to reread its config files.
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

Now that we've configured the resource requirements for a single tool, let's see how we can configure defaults for all tools, and reuse those defaults to reduce repetition and tedium.


> <hands-on-title>Configuring defaults and inheritance</hands-on-title>
>
> 1. Edit your `files/galaxy/config/tpv_rules_local.yml` and add the following settings.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -0,0 +1,16 @@
>    +global:
>    +  default_inherits: default
>    +
>    +tools:
>    +  default:
>    +    abstract: true
>    +    cores: 1
>    +    mem: cores * 4
>    +  testing:
>    +    cores: 2
>    -    mem: cores * 4
>    +
>    +destinations:
>    +  default:
>    +    abstract: true
>    +    params:
>    +      singularity_enabled: "true"
>    +      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    +  slurm:
>    +    inherits: default
>    +    runner: slurm
>    -    params:
>    -      singularity_enabled: "true"
>    -      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    +    env:
>    +      LC_ALL: C
>    +      SINGULARITY_CACHEDIR: /tmp/singularity
>    +      SINGULARITY_TMPDIR: /tmp
>    +
>    {% endraw %}
>
> We have defined a `global` section specifying that all tools and destinations should inherit from a specified `default`. We have then defined a tool named `default`, whose properties
> are implicitly inherited by all tools at runtime. This means that our `testing` tool will also inherit from this default tool, but it explicitly overrides cores
> We can also explicitly specify an `inherits` clause if we wish to extend a specific tool or destination, as shown in the destinations section. 
>
>

### Selecting between multiple destinations

1. using tags
2. Using max_accepted_cores


### TPV reference documentation

TPV has dedicated documentation at: https://total-perspective-vortex.readthedocs.io/en/latest/


# Configuring the TPV shared database

The Galaxy Project maintains a [shared database of TPV rules](https://github.com/galaxyproject/tpv-shared-database) so that admins do not have to independently rediscover ideal resource allocations for specific tools. These rules are based
on settings that have worked well in the usegalaxy.* federation. The rule file can simply be imported directly, with local overrides applied on top.

> <hands-on-title>Configuring defaults and inheritance</hands-on-title>
>
> 1. Edit your `templates/galaxy/config/job_conf.yml.j2` and add the location of the TPV shared rule file.
>
>    {% raw %}
>    ```diff
>    --- a/templates/galaxy/config/job_conf.yml.j2
>    +++ b/templates/galaxy/config/job_conf.yml.j2
>    @@ -20,4 +20,10 @@ execution:
>             value: /tmp/singularity
>           - name: SINGULARITY_TMPDIR
>             value: /tmp
>         tpv_dispatcher:
>           runner: dynamic
>           type: python
>           function: map_tool_to_destination
>           rules_module: tpv.rules
>           tpv_config_files:
>    +      - https://raw.githubusercontent.com/galaxyproject/tpv-shared-database/main/tools.yml
>    +      - config/tpv_rules_local.yml
>         singularity:
>           runner: local_runner
>           singularity_enabled: true
>    {% endraw %}
>    ```
>    {: data-commit="Importing TPV shared database via job conf"}
>
> Note how TPV allows the file to be imported directly via its http url. As many local and remote rule files as necessary can be combined, with rule files specified later overriding
> any previously specified rule files. The TPV shared database does not define destinations, only cores and mem settings, as well as any required environment vars.
> Take a look at the shared database of rules and note that some tools have very large recommended memory settings, which may or may not be available within your local cluster.
> Nevertheless, you may still wish to execute these tools with memory adjusted to suit your cluster's capabilities. 
>
> 2. Edit your `files/galaxy/config/tpv_rules_local.yml` and make the following changes.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -0,0 +1,16 @@
>    +global:
>    +  default_inherits: default
>    +
>    +tools:
>    +  default:
>    +    abstract: true
>    +    cores: 1
>    +    mem: cores * 4
>    +  testing:
>    +    cores: 2
>    -    mem: cores * 4
>    +
>    +destinations:
>    +  default:
>    +    abstract: true
>    +    params:
>    +      singularity_enabled: "true"
>    +      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    +  slurm:
>    +    inherits: default
>    +    runner: slurm
>    +    max_accepted_cores: 24
>    +    max_accepted_mem: 256
>    +    max_cores: 16
>    +    max_mem: 128
>    -    params:
>    -      singularity_enabled: "true"
>    -      native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores}
>    +    env:
>    +      LC_ALL: C
>    +      SINGULARITY_CACHEDIR: /tmp/singularity
>    +      SINGULARITY_TMPDIR: /tmp
>    +
>    {% endraw %}
>
>
> These changes indicate that the destination will accept jobs that are up to `max_accepted_cores: 24` and `max_accepted_mem: 256`. However, once accepted, it will forcibly clamp it down to 16 and 128 at most using the `max_cores` and `max_mem` clauses.
> If the tool does not exceed the specified limit, it will run happily. If not, it will be clamped down to the desired range. This allows even the largest resource requirement in the shared database to be accomodated.
>

# Basic access controls

You may wish to apply some basic restrictions on which users are allowed to run specific tools. TPV accomodates user and role specific rules. Combined with tags, these allow fine-grained control over


# More reading

The tutorial provides a quick overview of some of the basic features in TPV. However, there are numerous features that we have not convered such as:
a. Custom code blocks - execute arbitrary python code blocks, access additional context variables etc.
a. User and Role Handling - Add scheduling constraints based on the user's email or role
b. Metascheduling support - Perform advanced querying and filtering prior to choosing an appropriate destination
c. Job resubmissions - Resubmit jobs on failure
d. Linting, formatting and dry-run - Automatically format tpv rule files, catch potential syntax errors and perform a dry-run to check where a tool would get scheduled.

These features are covered in detail in the [TPV documentation](https://total-perspective-vortex.readthedocs.io/en/lates/).


# Job Resource Selectors

You may find that certain tools can benefit from having form elements added to them to allow for controlling certain job parameters, so that users can select based on their own knowledge. For example, a user might know that a particular set of parameters and inputs to a certain tool needs a larger memory allocation than the standard amount for a given tool. This of course assumes that your users are well behaved enough not to choose the maximum whenever available, although such concerns can be mitigated somewhat by the use of concurrency limits on larger memory destinations.

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
>    @@ -148,6 +148,8 @@ galaxy_config_templates:
>         dest: "{{ galaxy_config.galaxy.dependency_resolvers_config_file }}"
>       - src: templates/galaxy/config/tool_destinations.yml
>         dest: "{{ galaxy_config.galaxy.tool_destinations_config_file }}"
>    +  - src: templates/galaxy/config/job_resource_params_conf.xml.j2
>    +    dest: "{{ galaxy_config.galaxy.job_resource_params_file }}"
>     
>     galaxy_local_tools:
>     - testing.xml
>    {% endraw %}
>    ```
>    {: data-commit="Deploy job resource params configuration"}
>
> 3. Next, we define a new section in `job_conf.yml`: `<resources>`. This groups together parameters that should appear together on a tool form. Add the following section to your `templates/galaxy/config/job_conf.yml.j2`:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -69,6 +69,11 @@ galaxy_job_config:
>           dtd:
>             runner: dynamic
>             type: dtd
>    +  resources:
>    +    default: default
>    +    groups:
>    +      default: []
>    +      testing: [cores, time]
>       tools:
>         - class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>           environment: local_env
>    {% endraw %}
>    ```
>    {: data-commit="Configure resources in job conf"}
>
>    The group ID will be used to map a tool to job resource parameters, and the text value of the `<group>` tag is a comma-separated list of `name`s from `job_resource_params_conf.xml` to include on the form of any tool that is mapped to the defined `<group>`.
>
>
> 4. Finally, in `job_conf.yml`, move the previous `<tool>` definition for the `testing` tool into the comment and define a new `<tool>` that defines the `resources` for the tool:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -78,7 +78,8 @@ galaxy_job_config:
>         - class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>           environment: local_env
>         - id: testing
>    -      environment: dtd
>    +      environment: dynamic_cores_time
>    +      resources: testing
>     
>     galaxy_config:
>       galaxy:
>    {% endraw %}
>    ```
>    {: data-commit="Configure resources in job conf"}
>
> 5. We have assigned the `testing` tool to a new destination: `dynamic_cores_time`, but this destination does not exist. We need to create it. Add the following destination in your job conf:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -69,6 +69,9 @@ galaxy_job_config:
>           dtd:
>             runner: dynamic
>             type: dtd
>    +      dynamic_cores_time:
>    +        runner: dynamic
>    +        function: dynamic_cores_time
>       resources:
>         default: default
>         groups:
>    {% endraw %}
>    ```
>    {: data-commit="Add dynamic_cores_time destination"}
>
>    This will be another dynamic destination. Galaxy will load all python files in the {% raw %}`{{ galaxy_dynamic_rule_dir }}`{% endraw %}, and all functions defined in those will be available `dynamic_cores_time` to be used in the `job_conf.yml`
>
{: .hands_on}

This will set everything up to use the function. We have:

- A set of "job resources" defined which will let the user select the number of cores and walltime.
- A job configuration which says:
    - that our testing tool should allow selection of the cores and time parameters
    - directs it to use a new, `dynamic_cores_time` destination
    - and a has a new destination, `dynamic_cores_time`, which is defined as a dynamic destination which will call a python function we will load.

This is a lot but we're still missing the last piece for it to work:

## A dynamic destination

Lastly, we need to write the rule that will read the value of the job resource parameter form fields and decide how to submit the job.

> <hands-on-title>Writing a dynamic destination</hands-on-title>
>
> 1. Create and edit `files/galaxy/dynamic_job_rules/map_resources.py`. Create it with the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/dynamic_job_rules/map_resources.py
>    @@ -0,0 +1,42 @@
>    +import logging
>    +from galaxy.jobs.mapper import JobMappingException
>    +
>    +log = logging.getLogger(__name__)
>    +
>    +DESTINATION_IDS = {
>    +    1 : 'slurm',
>    +    2 : 'slurm-2c'
>    +}
>    +FAILURE_MESSAGE = 'This tool could not be run because of a misconfiguration in the Galaxy job running system, please report this error'
>    +
>    +
>    +def dynamic_cores_time(app, tool, job, user_email):
>    +    destination = None
>    +    destination_id = 'slurm'
>    +
>    +    # build the param dictionary
>    +    param_dict = job.get_param_values(app)
>    +
>    +    if param_dict.get('__job_resource', {}).get('__job_resource__select') != 'yes':
>    +        log.info("Job resource parameters not seleted, returning default destination")
>    +        return destination_id
>    +
>    +    # handle job resource parameters
>    +    try:
>    +        # validate params
>    +        cores = int(param_dict['__job_resource']['cores'])
>    +        time = int(param_dict['__job_resource']['time'])
>    +        destination_id = DESTINATION_IDS[cores]
>    +        destination = app.job_config.get_destination(destination_id)
>    +        # set walltime
>    +        if 'nativeSpecification' not in destination.params:
>    +            destination.params['nativeSpecification'] = ''
>    +        destination.params['nativeSpecification'] += ' --time=%s:00:00' % time
>    +    except:
>    +        # resource param selector not sent with tool form, job_conf.yml misconfigured
>    +        log.warning('(%s) error, keys were: %s', job.id, param_dict.keys())
>    +        raise JobMappingException(FAILURE_MESSAGE)
>    +
>    +    log.info('returning destination: %s', destination_id)
>    +    log.info('native specification: %s', destination.params.get('nativeSpecification'))
>    +    return destination or destination_id
>    {% endraw %}
>    ```
>    {: data-commit="Add map_resources python"}
>
>    It is important to note that **you are responsible for parameter validation, including the job resource selector**. This function only handles the job resource parameter fields, but it could do many other things - examine inputs, job queues, other tool parameters, etc.
>
>
> 2. As usual, we need to instruct Galaxy of where to find this file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -127,6 +127,7 @@ galaxy_config:
>         user_library_import_dir: /libraries/user
>         # Tool Destination Configuration
>         tool_destinations_config_file: "{{ galaxy_config_dir }}/tool_destinations.yml"
>    +    job_resource_params_file: "{{ galaxy_config_dir }}/job_resource_params_conf.xml"
>       gravity:
>         process_manager: systemd
>         galaxy_root: "{{ galaxy_root }}/server"
>    @@ -164,6 +165,7 @@ galaxy_local_tools:
>     - testing.xml
>     galaxy_dynamic_job_rules:
>     - my_rules.py
>    +- map_resources.py
>     
>     # Certbot
>     certbot_auto_renew_hour: "{{ 23 |random(seed=inventory_hostname)  }}"
>    {% endraw %}
>    ```
>    {: data-commit="Deploy map_resources.py"}
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
> 4. Run the **Testing Tool** with various resource parameter selections
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

{% snippet topics/admin/faqs/git-commit.md page=page %}

{% snippet topics/admin/faqs/missed-something.md step=9 %}

## Further Reading

- The [sample dynamic tool destination config file](https://github.com/galaxyproject/galaxy/blob/dev/config/tool_destinations.yml.sample) fully describes the configuration language
- [Dynamic destination documentation](https://docs.galaxyproject.org/en/latest/admin/jobs.html)
- Job resource parameters are not as well documented as they could be, but the [sample configuration file](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/test/files/galaxy/config/job_resource_params_conf.xml) shows some of the possibilities.
- [usegalaxy.org's job_conf.yml](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/common/templates/galaxy/config/job_conf.yml.j2) is publicly available for reference.
- [usegalaxy.eu's job_conf.xml](https://github.com/usegalaxy-eu/infrastructure-playbook/search?l=YAML&q=galaxy_jobconf) is likewise (see the `group_vars/galaxy.yml` result)
