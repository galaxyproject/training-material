---
layout: tutorial_hands_on

title: "Connecting Galaxy to a compute cluster"
questions:
  - How to connect Galaxy to a compute cluster?
  - How can I configure job dependent resources, like cores, memory for my DRM?
objectives:
  - Be familiar with the basics of installing, configuring, and using Slurm
  - Understand all components of the Galaxy job running stack
  - Understand how the `job_conf.xml` file controls Galaxy's jobs subsystem
  - Have a strong understanding of Galaxy job destinations
  - Know how to map tools to job destinations
  - Be able to use the dynamic job runner to make arbitrary destination mappings
  - Understand the job resource selector config and dynamic rule creation
  - The various ways in which tools can be mapped to destinations, both statically and dynamically
  - How to write a dynamic tool destination (DTD)
  - How to write a dynamic python function destination
  - How to use the job resource parameter selection feature
time_estimation: "4h"
key_points:
  - Galaxy supports a variety of different DRMs.
  - Dynamic Tool Destinations are a convenient way to map
  - Job resource parameters can allow you to give your users control over job resource requirements, if they are knowledgeable about the tools and compute resources available to them.
contributors:
  - natefoo
  - bgruening
  - hexylena
tags:
  - jobs
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

# Running Galaxy Jobs with Slurm

{% include snippets/warning_results_may_vary.md %}

The tools that are added to Galaxy can have a wide variance in the compute resources that they require and work efficiently on.
To account for this, Galaxy's job configuration needs to be tuned to run these tools properly. In addition, site-specific variables must
be taken into consideration when choosing where to run jobs and what parameters to run them with.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Installing Slurm

> ### {% icon hands_on %} Hands-on: Installing Slurm
>
> 1. Create and edit a file in your working directory called `requirements.yml` and include the following contents:
>
>    ```yaml
>    - galaxyproject.repos
>    - galaxyproject.slurm
>    ```
>
>    The `galaxyproject.repos` role adds the [Galaxy Packages for Enterprise Linux (GPEL)](https://depot.galaxyproject.org/yum/) repository for RedHat/CentOS, which provides both Slurm and Slurm-DRMAA (neither are available in standard repositories or EPEL). For Ubuntu versions 18.04 or newer, it adds the [Slurm-DRMAA PPA](https://launchpad.net/~natefoo/+archive/ubuntu/slurm-drmaa) (Slurm-DRMAA was removed from Debian/Ubuntu in buster/bionic).
>
> 2. In the same directory, run `ansible-galaxy install -p roles -r requirements.yml`. This will install all of the required modules for this training into the `roles/` folder. We choose to install to a folder to give you easy access to look through the different roles when you have questions on their behaviour.
>
> 3. Create the hosts file if you have not done so, include a group for `[galaxyservers]` with the address of the host where you will install Slurm
>
> 4. Create a playbook, `slurm.yml` which looks like the following:
>
>    ```yaml
>    - hosts: galaxyservers
>      become: true
>      vars:
>        slurm_roles: ['controller', 'exec']
>        slurm_nodes:
>        - name: localhost
>          CPUs: 2                              # Here you would need to figure out how many cores your machine has. (Hint, `htop`)
>        slurm_config:
>          FastSchedule: 2                      # Ignore errors if the host actually has cores != 2
>          SelectType: select/cons_res
>          SelectTypeParameters: CR_CPU_Memory  # Allocate individual cores/memory instead of entire node
>      roles:
>        - galaxyproject.repos
>        - galaxyproject.slurm
>    ```
>
> 5. Run the playbook (`ansible-playbook -i hosts slurm.yml`)
>
{: .hands_on}


Note that the above Slurm config options are only those that are useful for this training exercise. In production, you would want to use a more appropriate configuration specific to your cluster (and setting `FastSchedule` to `2` is not recommended).

Installed with Slurm is MUNGE (MUNGE Uid 'N Gid Emporium...) which authenticates users between cluster hosts. You would normally need to ensure the same Munge key is distributed across all cluster hosts (in `/etc/munge/munge.key`) - A great task for Ansible. However, the installation of the munge package has created a random key for you, and you will not need to distribute this since you'll run jobs only on a single host.

You can now check that all of the daemons are running with the command `systemctl status munge slurmd slurmctld`

```console
$ sudo systemctl status munge slurmd slurmctld
● munge.service - MUNGE authentication service
   Loaded: loaded (/usr/lib/systemd/system/munge.service; enabled; vendor preset: disabled)
   Active: active (running) since Sa 2019-01-26 22:38:13 CET; 28min ago
     Docs: man:munged(8)
 Main PID: 22930 (munged)
    Tasks: 4
   Memory: 128.0K
   CGroup: /system.slice/munge.service
           └─22930 /usr/sbin/munged

Jan 26 22:38:13 helena-test.novalocal systemd[1]: Starting MUNGE authentication service...
Jan 26 22:38:13 helena-test.novalocal systemd[1]: Started MUNGE authentication service.

● slurmd.service - Slurm node daemon
   Loaded: loaded (/usr/lib/systemd/system/slurmd.service; enabled; vendor preset: disabled)
   Active: active (running) since Sa 2019-01-26 23:04:21 CET; 2min 25s ago
  Process: 15051 ExecStart=/usr/sbin/slurmd $SLURMD_OPTIONS (code=exited, status=0/SUCCESS)
 Main PID: 15054 (slurmd)
    Tasks: 1
   Memory: 628.0K
   CGroup: /system.slice/slurmd.service
           └─15054 /usr/sbin/slurmd

Jan 26 23:04:21 helena-test.novalocal systemd[1]: Starting Slurm node daemon...
Jan 26 23:04:21 helena-test.novalocal systemd[1]: PID file /var/run/slurmd.pid not readable (yet?) after start.
Jan 26 23:04:21 helena-test.novalocal systemd[1]: Started Slurm node daemon.

● slurmctld.service - Slurm controller daemon
   Loaded: loaded (/usr/lib/systemd/system/slurmctld.service; enabled; vendor preset: disabled)
   Active: active (running) since Sa 2019-01-26 23:04:20 CET; 2min 26s ago
  Process: 15040 ExecStart=/usr/sbin/slurmctld $SLURMCTLD_OPTIONS (code=exited, status=0/SUCCESS)
 Main PID: 15042 (slurmctld)
    Tasks: 7
   Memory: 1.1M
   CGroup: /system.slice/slurmctld.service
           └─15042 /usr/sbin/slurmctld

Jan 26 23:04:20 helena-test.novalocal systemd[1]: Starting Slurm controller daemon...
Jan 26 23:04:20 helena-test.novalocal systemd[1]: PID file /var/run/slurmctld.pid not readable (yet?) after start.
Jan 26 23:04:20 helena-test.novalocal systemd[1]: Started Slurm controller daemon.
```

Running the playbook, the Slurm configuration, `/etc/slurm/slurm.conf` (or `/etc/slurm-llnl/slurm.conf` on Debian-based distributions) was created for you automatically. All of the variables were set by default. If you need to override the configuration yourself, Slurm provides [an online tool](https://slurm.schedmd.com/configurator.html) which will help you configure it.

## Using Slurm

You should now be able to see that your Slurm cluster is operational with the `sinfo` command. This shows the state of nodes and partitions (synonymous with queues in other DRMs). The "node-oriented view" provided with the `-N` flag is particularly useful:

```console
$ sinfo
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
debug*       up   infinite      1   idle localhost
$ sinfo -Nel
Fri Nov  4 16:51:24 2016
NODELIST   NODES PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT FEATURES REASON
localhost      1    debug*        idle    1    2:1:1      1        0      1   (null) none
```

If your node state is not `idle`, something has gone wrong. If your node state ends with an asterisk \*, the Slurm controller is attempting to contact the Slurm execution daemon but has not yet been successful (the \* next to the partition name is normal, it indicates the default partition).

We want to ensure that Slurm is actually able to run jobs. There are two ways this can be done:

- `srun`: Run an interactive job (e.g. a shell, or a specific program with its stdin, stdout, and stderr all connected to your terminal.
- `sbatch`: Run a batch job, with stdin closed and stdout/stderr redirected to a file.

Galaxy runs `sbatch` jobs but we can use both `srun` and `sbatch` to test:


> ### {% icon hands_on %} Hands-on: Running commands with `srun`
>
> 1. Use [`srun`](https://slurm.schedmd.com/srun.html) to run the command `uname -a`
>
>
>    > ### {% icon question %} Question
>    >
>    > How did the output look?
>    >
>    > > ### {% icon solution %} Solution
>    > > Your output may look slightly different:
>    > > ```console
>    > > $ srun uname -a
>    > > Linux helena-test.novalocal 3.10.0-862.14.4.el7.x86_64 #1 SMP Wed Sep 26 15:12:11 UTC 2018 x86_64 x86_64 x86_64 GNU/Linux
>    > > ```
>    > {: .solution }
>    {: .question}
{: .hands_on}

Although it looks like this command ran as if I had not used `srun`, it was in fact routed through Slurm.

> ### {% icon hands_on %} Hands-on: Running commands with `sbatch`
>
> 1. Create a test job script somewhere, such as in `~/sbatch-test.sh`. It should be a batch script which runs `uname -a`, `uptime`, and sleeps for 30 seconds.
>
>    > ### {% icon question %} Question
>    >
>    > What does your shell script look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > ```bash
>    > > #!/bin/bash
>    > > uname -a
>    > > uptime
>    > > sleep 30
>    > > ```
>    > {: .solution }
>    {: .question}
>
> 2. Make the script executable, `chmod +x ~/sbatch-test.sh`
>
> 3. Use [`sbatch`](https://slurm.schedmd.com/sbatch.html) to submit the job script
>
>    > ### {% icon question %} Question
>    >
>    > What command did you run?
>    >
>    > > ### {% icon solution %} Solution
>    > > ```console
>    > > $ sbatch ~/sbatch-test.sh
>    > > ```
>    > {: .solution }
>    {: .question}
>
> 4. Use [`squeue`](https://slurm.schedmd.com/squeue.html) to check the queue
>
>    > ### {% icon question %} Question
>    >
>    > What did the output look like?
>    >
>    > > ### {% icon solution %} Solution
>    > > ```console
>    > >JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
>    > >    3     debug sbatch-t   ubuntu  R       0:22      1 localhost
>    > > ```
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

If you've made it this far, your Slurm installation is working!

## Slurm-DRMAA

Above Slurm in the stack is slurm-drmaa, a library that provides a translational interface from the Slurm API to the generalized DRMAA API in C.

> ### {% icon hands_on %} Hands-on: Installing Slurm-DRMAA
>
> 1. Add a `post_task` to your playbook to install `slurm-drmaa1` (Debian/Ubuntu) or `slurm-drmaa` (RedHat/CentOS), and additionally include the `galaxyproject.repos` role
>
>    ```yaml
>    - hosts: galaxyservers
>      become: true
>      vars:
>        slurm_roles: ['controller', 'exec']
>      roles:
>        - galaxyproject.repos
>        - galaxyproject.slurm
>      post_tasks:
>        - name: Install slurm-drmaa
>          package:
>            name: slurm-drmaa1
>    ```
>
> 2. Run the playbook (`ansible-playbook -i hosts slurm.yml`)
>
{: .hands_on}

Moving one level further up the stack, we find DRMAA Python. This is a Galaxy framework *conditional dependency*. Conditional dependencies are only installed if, during startup, a configuration option is set that requires that dependency. The `galaxyproject.galaxy` Ansible role will install these conditional dependencies, automatically.

# Galaxy and Slurm

At the top of the stack sits Galaxy. Galaxy must now be configured to use the cluster we've just set up. The DRMAA Python documentation (and Galaxy's own documentation) instruct that you should set the `$DRMAA_LIBRARY_PATH` environment variable so that DRMAA Python can find `libdrmaa.so` (aka slurm-drmaa). Because Galaxy runs under supervisor, the environment that Galaxy starts under is controlled by the `environment` option in `/etc/supervisor/conf.d/galaxy.conf`. The galaxy task should thus be updated to refer to the path to slurm-drmaa, which is `/usr/lib/slurm-drmaa/lib/libdrmaa.so.1`:


> ### {% icon hands_on %} Hands-on: Making Galaxy aware of DRMAA
>
> 1. Open your group variables and edit the supervisor task, update the environment variable:
>
>    ```yaml
>    supervisor_programs:
>      - name: galaxy
>        ...
>        configuration: |
>          ...
>          {% raw %}environment=HOME="{{ galaxy_mutable_data_dir }}",VIRTUAL_ENV="{{ galaxy_venv_dir }}",PATH="{{ galaxy_venv_dir }}/bin:%(ENV_PATH)s",DRMAA_LIBRARY_PATH="/usr/lib/slurm-drmaa/lib/libdrmaa.so.1"{% endraw %}
>   ```
>
> 2. We need to modify `job_conf.xml` to instruct Galaxy's job handlers to load the Slurm job runner plugin, and set the Slurm job submission parameters. A job runner plugin definition must have the `id`, `type`, and `load` attributes. The entire `<plugins>` tag group should look like:
>
>    If the folder does not exist, create `files/galaxy/config` next to your `playbook.yml` (`mkdir -p files/galaxy/config/`)
>
>    Create `files/galaxy/config/job_conf.xml` with the following contents:
>
>    ```xml
>    <job_conf>
>        <plugins workers="4">
>            <plugin id="local" type="runner" load="galaxy.jobs.runners.local:LocalJobRunner"/>
>            <plugin id="slurm" type="runner" load="galaxy.jobs.runners.slurm:SlurmJobRunner"/>
>        </plugins>
>        <destinations>
>            <destination id="local" runner="local"/>
>        </destinations>
>    </job_conf>
>    ```
> 3. Next, we need to add a new destination for the Slurm job runner. This is a basic destination with no parameters, Galaxy will do the equivalent of submitting a job as `sbatch /path/to/job_script.sh`. Note that we also need to set a default destination now that more than one destination is defined. In a `<destination>` tag, the `id` attribute is a unique identifier for that destination and the `runner` attribute must match the `id` of defined plugin:
>
>    ```xml
>    <destinations default="slurm">
>        <destination id="slurm" runner="slurm"/>
>        <destination id="local" runner="local"/>
>    </destinations>
>    ```
>
> 4. Inform `galaxyproject.galaxy` of where you would like the `job_conf.xml` to reside in your group variables:
>
>    ```yaml
>    galaxy_config:
>      galaxy:
>        job_config_file: {% raw %}"{{ galaxy_config_dir }}/job_conf.xml"{% endraw %}
>    ```
>
>    And then deploy the new config file using the `galaxy_config_files` var in your group vars
>
>    ```yaml
>    galaxy_config_files:
>      - src: files/galaxy/config/job_conf.xml
>        dest: {% raw %}"{{ galaxy_config['galaxy']['job_config_file'] }}"{% endraw %}
>    ```
>
>      The variable `galaxy_config_files` is an array of hashes, each with `src` and `dest`, the files from src will be copied to dest on the server. `galaxy_template_files` exist to template files out.
>
> 5. Run your *Galaxy* playbook (`ansible-playbook -i hosts playbook.yml`)
>
> 6. Follow the logs with `supervisorctl tail -f galaxy stderr`
>
> 7. Rerun the playbook. Because we updated the supervisor config, Galaxy will automatically be restarted.
>
{: .hands_on}

Two sections of the log output are of interest. First, when Galaxy parses `job_conf.xml`:

```
galaxy.jobs DEBUG 2016-11-05 14:07:12,649 Loading job configuration from /srv/galaxy/config/job_conf.xml
galaxy.jobs DEBUG 2016-11-05 14:07:12,650 Read definition for handler 'handler0'
galaxy.jobs DEBUG 2016-11-05 14:07:12,651 Read definition for handler 'handler1'
galaxy.jobs DEBUG 2016-11-05 14:07:12,652 <handlers> default set to child with id or tag 'handlers'
galaxy.jobs DEBUG 2016-11-05 14:07:12,652 <destinations> default set to child with id or tag 'slurm'
galaxy.jobs DEBUG 2016-11-05 14:07:12,653 Done loading job configuration
```

Second, when Galaxy loads job runner plugins:

```
galaxy.jobs.manager DEBUG 2016-11-05 14:07:22,341 Starting job handler
galaxy.jobs INFO 2016-11-05 14:07:22,347 Handler 'handler0' will load all configured runner plugins
galaxy.jobs.runners DEBUG 2016-11-05 14:07:22,355 Starting 4 LocalRunner workers
galaxy.jobs DEBUG 2016-11-05 14:07:22,367 Loaded job runner 'galaxy.jobs.runners.local:LocalJobRunner' as 'local'
pulsar.managers.util.drmaa DEBUG 2016-11-05 14:07:22,434 Initializing DRMAA session from thread MainThread
galaxy.jobs.runners DEBUG 2016-11-05 14:07:22,443 Starting 4 SlurmRunner workers
galaxy.jobs DEBUG 2016-11-05 14:07:22,455 Loaded job runner 'galaxy.jobs.runners.slurm:SlurmJobRunner' as 'slurm'
galaxy.jobs.handler DEBUG 2016-11-05 14:07:22,455 Loaded job runners plugins: slurm:local
```

## Running a Job

You should now be able to run a Galaxy job through Slurm. The simplest way to test is using the upload tool to upload some text.

> ### {% icon hands_on %} Hands-on: Testing a Slurm Job
>
> 1. If you're not still following the log files with `tail`, do so now.
> 2. Click the upload button at the top of the tool panel (on the left side of the Galaxy UI).
> 3. In the resulting modal dialog, click the "Paste/Fetch data" button.
> 4. Type some random characters into the text field that has just appeared.
> 5. Click "Start" and then "Close"
>
{: .hands_on}

In your `tail` terminal window you should see the following messages:

```
galaxy.jobs DEBUG 2016-11-05 14:07:22,862 (2) Persisting job destination (destination id: slurm)
galaxy.jobs.runners DEBUG 2016-11-05 14:07:22,958 Job [2] queued (328.180 ms)
galaxy.jobs.handler INFO 2016-11-05 14:07:22,996 (2) Job dispatched
galaxy.tools.deps DEBUG 2016-11-05 14:07:23,621 Building dependency shell command for dependency 'samtools'
  ...
galaxy.tools.deps WARNING 2016-11-05 14:07:23,631 Failed to resolve dependency on 'samtools', ignoring
galaxy.jobs.command_factory INFO 2016-11-05 14:07:23,674 Built script [/srv/galaxy/server/database/jobs/000/2/tool_script.sh] for tool command[python /srv/galaxy/server/tools/data_source/upload.py /srv/galaxy/server /srv/galaxy/server/database/tmp/tmpkiMZKd /srv/galaxy/server/database/tmp/tmpJuSMo5 2:/srv/galaxy/server/database/jobs/000/2/dataset_2_files:/srv/galaxy/server/database/datasets/000/dataset_2.dat]
galaxy.tools.deps DEBUG 2016-11-05 14:07:24,033 Building dependency shell command for dependency 'samtools'
  ...
galaxy.tools.deps WARNING 2016-11-05 14:07:24,038 Failed to resolve dependency on 'samtools', ignoring
galaxy.jobs.runners DEBUG 2016-11-05 14:07:24,052 (2) command is: mkdir -p working; cd working; /srv/galaxy/server/database/jobs/000/2/tool_script.sh; return_code=$?; cd '/srv/galaxy/server/database/jobs/000/2'; python "/srv/galaxy/server/database/jobs/000/2/set_metadata_CALKH0.py" "/srv/galaxy/server/database/tmp/tmpkiMZKd" "/srv/galaxy/server/database/jobs/000/2/working/galaxy.json" "/srv/galaxy/server/database/jobs/000/2/metadata_in_HistoryDatasetAssociation_2_nnti4M,/srv/galaxy/server/database/jobs/000/2/metadata_kwds_HistoryDatasetAssociation_2_sN3gVP,/srv/galaxy/server/database/jobs/000/2/metadata_out_HistoryDatasetAssociation_2_jIhXJJ,/srv/galaxy/server/database/jobs/000/2/metadata_results_HistoryDatasetAssociation_2_v4v_dv,/srv/galaxy/server/database/datasets/000/dataset_2.dat,/srv/galaxy/server/database/jobs/000/2/metadata_override_HistoryDatasetAssociation_2_OQwwTH" 5242880; sh -c "exit $return_code"
galaxy.jobs.runners.drmaa DEBUG 2016-11-05 14:07:24,125 (2) submitting file /srv/galaxy/server/database/jobs/000/2/galaxy_2.sh
galaxy.jobs.runners.drmaa INFO 2016-11-05 14:07:24,172 (2) queued as 7
galaxy.jobs DEBUG 2016-11-05 14:07:24,172 (2) Persisting job destination (destination id: slurm)
galaxy.jobs.runners.drmaa DEBUG 2016-11-05 14:07:24,539 (2/7) state change: job is queued and active
```

At this point the job has been accepted by Slurm and is awaiting scheduling on a node. Once it's been sent to a node and starts running, Galaxy logs this event:

```
galaxy.jobs.runners.drmaa DEBUG 2016-11-05 14:07:25,559 (2/7) state change: job is running
```

Finally, when the job is complete, Galaxy performs its job finalization process:

```
galaxy.jobs.runners.drmaa DEBUG 2016-11-05 14:07:30,883 (2/7) state change: job finished normally
galaxy.model.metadata DEBUG 2016-11-05 14:07:31,132 loading metadata from file for: HistoryDatasetAssociation 2
galaxy.jobs INFO 2016-11-05 14:07:31,336 Collecting metrics for Job 2
galaxy.jobs DEBUG 2016-11-05 14:07:31,370 job 2 ended (finish() executed in (411.821 ms))
galaxy.model.metadata DEBUG 2016-11-05 14:07:31,375 Cleaning up external metadata files
```

Note a few useful bits in the output:
- `Persisting job destination (destination id: slurm)`: Galaxy has selected the `slurm` destination we defined
- `submitting file /srv/galaxy/server/database/jobs/000/2/galaxy_2.sh`: This is the path to the script that is submitted to Slurm as it would be with `sbatch`
- `(2) queued as 7`: Galaxy job id "2" is Slurm job id "7".
- If `job <id> ended` is reached, the job should show as done in the UI

Slurm allows us to query the exit state of jobs for a time period of the value of Slurm's `MinJobAge` option, which defaults to 300 (seconds, == 5 minutes):

```console
$ scontrol show job 7
JobId=7 JobName=g2_upload1_anonymous_10_0_2_2
   UserId=galaxy(999) GroupId=galaxy(999)
   Priority=4294901754 Nice=0 Account=(null) QOS=(null)
   JobState=COMPLETED Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:06 TimeLimit=UNLIMITED TimeMin=N/A
   SubmitTime=2016-11-05T14:07:24 EligibleTime=2016-11-05T14:07:24
   StartTime=2016-11-05T14:07:24 EndTime=2016-11-05T14:07:30
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=debug AllocNode:Sid=gat2016:16025
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=localhost
   BatchHost=localhost
   NumNodes=1 NumCPUs=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
   Features=(null) Gres=(null) Reservation=(null)
   Shared=0 Contiguous=0 Licenses=(null) Network=(null)
   Command=(null)
   WorkDir=/srv/galaxy/server/database/jobs/000/2
   StdErr=/srv/galaxy/server/database/jobs/000/2/galaxy_2.e
   StdIn=StdIn=/dev/null
   StdOut=/srv/galaxy/server/database/jobs/000/2/galaxy_2.o
   Power= SICP=0
```

After the job has been purged from the active jobs database, a bit of information (but not as much as `scontrol` provides) can be retrieved from Slurm's logs. However, it's a good idea to set up Slurm's accounting database to keep old job information in a queryable format.

## Further Reading

- [Galaxy's cluster documentation](https://docs.galaxyproject.org/en/latest/admin/cluster.html) describes in detail alternative cluster configurations.
- [The job_conf.xml documentation](https://docs.galaxyproject.org/en/latest/admin/jobs.html) fully describes the syntax of the job configuration file.
- The [Distributed Resource Management Application API (DRMAA)](https://www.drmaa.org/) page contains the DRMAA specification as well as documentation for various implementations. It also includes a list of DRMs supporting DRMAA.
- The [Slurm documentation](http://slurm.schedmd.com/) is extensive and covers all the features and myriad of ways in which you can configure slurm.
- [PSNC slurm-drmaa](http://apps.man.poznan.pl/trac/slurm-drmaa)'s page includes documentation and the SVN repository, which has a few minor fixes since the last released version. PSNC also wrote the initial implementations of the DRMAA libraries for PBSPro and LSF, so all three are similar.
- [Our own fork of slurm-drmaa](http://github.com/natefoo/slurm-drmaa) includes support for Slurms `-M`/`--clusters` multi-cluster functionality.
- [Slurm Accounting documentation](http://slurm.schedmd.com/accounting.html) explains how to set up SlurmDBD.

# Galaxy and Slurm - Statically Mapping a Job

We don't want to overload our training VMs trying to run real tools, so to demonstrate how to map a multicore tool to a multicore destination, we'll create a fake tool.

## Writing a testing tool

> ### {% icon hands_on %} Hands-on: Deploying a Tool
>
> 1. Create the directory `files/galaxy/tools/` if it doesn't exist and edit a new file in `files/galaxy/tools/testing.xml` with the following contents:
>
>    ```xml
>    <tool id="testing" name="Multicore Tool">
>        <command>
>            <![CDATA[echo "Running with '\${GALAXY_SLOTS:-1}' threads" > "$output1"]]>
>        </command>
>        <inputs>
>            <param name="input1" type="data" format="txt" label="Input Dataset"/>
>        </inputs>
>        <outputs>
>            <data name="output1" format="txt" />
>        </outputs>
>    </tool>
>    ```
>    {: .question}
>
> 2. Add the tool to the Ansible automated tool conf, `galaxy_local_tools`
>
>    ```yaml
>    galaxy_local_tools:
>    - testing.xml
>    ```
>
> 3. Run the playbook
>
> 4. Reload Galaxy in your browser and the new tool should now appear in the tool panel. If you have not already created a dataset in your history, upload a random text dataset. Once you have a dataset, click the tool's name in the tool panel, then click Execute.
>
>    > ### {% icon question %} Question
>    >
>    > What is the tool's output?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > ```
>    > > Running with '1' threads
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
{: .hands_on}

Of course, this tool doesn't actually *use* the allocated number of cores. In a real tool, you would call the tools's underlying command with whatever flag that tool provides to control the number of threads or processes it starts, such as `samtools sort -@ \${GALAXY_SLOTS:-1}`.

## Running with more resources

We want our tool to run with more than one core. To do this, we need to instruct Slurm to allocate more cores for this job. This is done in the job configuration file.


> ### {% icon hands_on %} Hands-on: Allocating more resources
>
> 1. Edit your `files/galaxy/config/job_conf.xml` and add the following destination:
>
>    ```xml
>    <destination id="slurm-2c" runner="slurm">
>        <param id="nativeSpecification">--nodes=1 --ntasks=2</param>
>    </destination>
>    ```
> 2. Then, map the new tool to the new destination using the tool ID (`<tool id="testing">`) and destination id (`<destination id="slurm-2c">`) by adding a new section to the job config, `<tools>`:
>
>    ```xml
>        <tools>
>            <tool id="testing" destination="slurm-2c"/>
>        </tools>
>    ```
>
> 3. Run the playbook. Because we modified `job_conf.xml`, Galaxy will be restarted to reread its config files.
>
> 4. Click the rerun button on the last history item, or click **Testing Tool** in the tool panel, and then click the tool's Execute button.
>
>    > ### {% icon question %} Question
>    >
>    > What is the tool's output?
>    >
>    > > ### {% icon solution %} Solution
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


# Dynamic Job Destinations

Dynamic destinations allow you to write custom python code to dispatch jobs based on whatever rules you like. For example, UseGalaxy.eu at one point used a very complex custom dispatching configuration to handle sorting jobs between multiple clusters. Galaxy has [extensive documentation](https://docs.galaxyproject.org/en/latest/admin/jobs.html#dynamic-destination-mapping-python-method) on how to write these sort of destinations.

> ### {% icon hands_on %} Hands-on: Writing a dynamic job destination
>
> 1. Create and open `files/galaxy/dynamic_job_rules/my_rules.py`
>
>    ```python
>    from galaxy.jobs import JobDestination
>    from galaxy.jobs.mapper import JobMappingException
>    import os
>
>    def admin_only(app, user_email):
>        # Only allow the tool to be executed if the user is an admin
>        admin_users = app.config.get( "admin_users", "" ).split( "," )
>        if user_email not in admin_users:
>            raise JobMappingException("Unauthorized.")
>        return JobDestination(runner="slurm")
>    ```
>
>    This destination will check that the `user_email` is in the set of `admin_users` from your config file.
>
> 2. As usual, we need to instruct Galaxy of where to find this file:
>
>    - Edit your group variables file and add the following:
>
>      ```yml
>      galaxy_dynamic_job_rules:
>        - my_rules.py
>      ```
>
> 3. We next need to configure this plugin in our job configuration:
>
>    ```xml
>    <destination id="dynamic_admin_only" runner="dynamic">
>        <param id="type">python</param>
>        <param id="function">admin_only</param>
>    </destination>
>    ```
>
>    This is a **Python function dynamic destination**. Galaxy will load all python files in the {% raw %}`{{ galaxy_dynamic_rule_dir }}`{% endraw %}, and all functions defined in those will be available `my_rules.py` to be used in the `job_conf.xml`
>
> 4. Finally, in `job_conf.xml`, update the `<tool>` definition and point it to this destination:
>
>    ```xml
>    <tools>
>        <tool id="testing" destination="dynamic_admin_only" />
>    </tools>
>    ```
>
> 5. Run the playbook / restart Galaxy
>
{: .hands_on}


Try running the tool as both an admin user and a non-admin user, non-admins should not be able to run it. You can imagine extending this to complex logic for permissions, or for destination mapping depending on numerous factors. We did not cover it, but in the documentation you can add additional variables to your function signature, and they will be automatically supplied. Some useful variables are `tool`, `user`, `job`, and `app` if you need to load configuration information.

# Dynamically map a tool to a job destination

If you don't want to write dynamic destinations yourself, Dynamic Tool Destinations (DTDs) utilize the dynamic job runner to provide dynamic job mapping functionality without having to explicitly write code to perform the mapping. The mapping functionality is mostly limited to input sizes, but often input size is the most important factor in deciding what resources to allocate for a job.

## Writing a Dynamic Tool Destination

> ### {% icon hands_on %} Hands-on: Writing a DTD
>
> 1. Dynamic tool destinations are configured via a YAML file. As before, we'll use a fake example but this is extremely useful in real-life scenarios. Create the file `files/galaxy/config/tool_destinations.yml` with the following contents:
>
>    ```yaml
>    ---
>    tools:
>      testing:
>        rules:
>          - rule_type: file_size
>            lower_bound: 16
>            upper_bound: Infinity
>            destination: slurm-2c
>        default_destination: slurm
>    default_destination: local
>    verbose: True
>    ```
>
>    The rule says:
>    - If the tool has ID `testing`:
>      - If the input dataset is >=16 bytes, run on the destination `slurm-2c`
>      - If the input dataset is <16 bytes, run on the destination `slurm`
>    - Else, run on the destination `local`
>
> 2. We also need to inform Galaxy of the path to the file we've just created, which is done using the `tool_destinations_config_file` in `galaxy_config` > `galaxy`. Additionally we need to add a `galaxy_config_files` entry to ensure it is deployed.
>
>    ```yml
>    galaxy_config:
>      galaxy:
>        tool_destinations_config_file: {% raw %}"{{ galaxy_config_dir }}/tool_destinations.yml"{% endraw %}
>    ...
>    galaxy_config_files:
>        ...
>        - src: files/galaxy/config/tool_destinations.yml
>          dest: {% raw %}"{{ galaxy_config['galaxy']['tool_destinations_config_file'] }}"{% endraw %}
>    ```
>
> 3. We need to update Galaxy's job configuration to use this rule. Open `files/galaxy/config/job_conf.xml` and add a DTD destination:
>
>    ```xml
>    <destination id="dtd" runner="dynamic">
>        <param id="type">dtd</param>
>    </destination>
>    ```
>
>    Also, comment out the previous `<tool>` definition for the `testing` tool, and replace it with a mapping to the dtd destination like so:
>
>    ```xml
>    <tools>
>    <!--
>        <tool id="testing" destination="slurm-2c"/>
>        <tool id="testing" destination="dynamic_admin_only" />
>    -->
>        <tool id="testing" destination="dtd"/>
>    </tools>
>    ```
>
> 4. Run the playbook and restart Galaxy
>
{: .hands_on}

## Testing the DTD

Our rule specified that any invocation of the `testing` tool with an input dataset with size <16 bytes would run on the 1 core destination, whereas any with >= 16 bytes would run on the 2 core destination.

> ### {% icon hands_on %} Hands-on: Testing the DTD
>
> 1. Create a dataset using the upload paste tool with a few (<16) characters
>
> 2. Create a dataset using the upload paste tool with >16 characters
>
> 3. Run the `Testing Tool` on both datasets.
>
{: .hands_on}

You can imagine using this to run large blast jobs on compute hardware with more resources, or giving them more CPU cores. Some tools require more memory as job inputs increase, you can use this to run tools with a larger memory limit, if you know it will need it to process a certain size of inputs.

# Job Resource Selectors

You may find that certain tools can benefit from having form elements added to them to allow for controlling certain job parameters, so that users can select based on their own knowledge. For example, a user might know that a particular set of parameters and inputs to a certain tool needs a larger memory allocation than the standard amount for a given tool. This of course assumes that your users are well behaved enough not to choose the maximum whenever available, although such concerns can be mitigated somewhat by the use of concurrency limits on larger memory destinations.

Such form elements can be added to tools without modifying each tool's configuration file through the use of the **job resource parameters configuration file**

> ### {% icon hands_on %} Hands-on: Configuring a Resource Selector
>
> 1. Create and open `files/galaxy/config/job_resource_params_conf.xml`
>
>    ```xml
>    <parameters>
>        <param label="Cores" name="cores" type="select" help="Number of cores to run job on.">
>            <option value="1">1 (default)</option>
>            <option value="2">2</option>
>        </param>
>      <param label="Time" name="time" type="integer" size="3" min="1" max="24" value="1" help="Maximum job time in hours, 'walltime' value (1-24). Leave blank to use default value." />
>    </parameters>
>    ```
>
>    This defines two resource fields, a select box where users can choose between 1 and 2 cores, and a text entry field where users can input an integer value from 1-24 to set the walltime for a job.
>
> 2. As usual, we need to instruct Galaxy of where to find this file:
>
>    - Edit your group variables file and add the following:
>
>      ```yml
>      galaxy_config:
>        galaxy:
>          job_resource_params_file: {% raw %}"{{ galaxy_config_dir }}/job_resource_params_conf.xml"{% endraw %}
>      ...
>      galaxy_config_files:
>        ...
>        - src: files/galaxy/config/job_resource_params_conf.xml
>          dest: {% raw %}"{{ galaxy_config['galaxy']['job_resource_params_file'] }}"{% endraw %}
>      ```
>
> 3. Next, we define a new section in `job_conf.xml`: `<resources>`. This groups together parameters that should appear together on a tool form. Add the following section to your `files/galaxy/config/job_conf.xml`:
>
>    ```xml
>    <resources>
>        <group id="testing">cores,time</group>
>    </resources>
>    ```
>
>    The group ID will be used to map a tool to job resource parameters, and the text value of the `<group>` tag is a comma-separated list of `name`s from `job_resource_params_conf.xml` to include on the form of any tool that is mapped to the defined `<group>`.
>
>
> 4. Finally, in `job_conf.xml`, move the previous `<tool>` definition for the `testing` tool into the comment and define a new `<tool>` that defines the `resources` for the tool:
>
>    ```xml
>    <tools>
>        <!--
>        <tool id="testing" destination="slurm-2c"/>
>        <tool id="testing" destination="dtd"/>
>        -->
>        <tool id="testing" destination="dynamic_cores_time" resources="testing"/>
>    </tools>
>    ```
> 5. We have assigned the `testing` tool to a new destination: `dynamic_cores_time`, but this destination does not exist. We need to create it. Add the following destination in your job conf:
>
>    ```xml
>    <destination id="dynamic_cores_time" runner="dynamic">
>        <param id="type">python</param>
>        <param id="function">dynamic_cores_time</param>
>    </destination>
>    ```
>
>    This will be another dynamic destination. Galaxy will load all python files in the {% raw %}`{{ galaxy_dynamic_rule_dir }}`{% endraw %}, and all functions defined in those will be available `dynamic_cores_time` to be used in the `job_conf.xml`
>
{: .hands_on}

This will set everything up to use the function. We have:

- A set of "job resources" defined which will let the user select the number of cores and walltime.
- A job configuration which says:
    -  that our testing tool should allow selection of the cores and time parameters
    - directs it to use a new, `dynamic_cores_time` destination
    - and a has a new destination, `dynamic_cores_time`, which is defined as a dynamic destination which will call a python function we will load.

This is a lot but we're still missing the last piece for it to work:

## A dynamic destination

Lastly, we need to write the rule that will read the value of the job resource parameter form fields and decide how to submit the job.

> ### {% icon hands_on %} Hands-on: Writing a dynamic destination
>
> 1. Create and edit `files/galaxy/dynamic_job_rules/map_resources.py`. Create it with the following contents:
>
>    ```python
>    import logging
>    from galaxy.jobs.mapper import JobMappingException
>
>    log = logging.getLogger(__name__)
>
>    DESTINATION_IDS = {
>        1 : 'slurm',
>        2 : 'slurm-2c'
>    }
>    FAILURE_MESSAGE = 'This tool could not be run because of a misconfiguration in the Galaxy job running system, please report this error'
>
>
>    def dynamic_cores_time(app, tool, job, user_email):
>        destination = None
>        destination_id = 'slurm'
>
>        # build the param dictionary
>        param_dict = job.get_param_values(app)
>
>        # handle job resource parameters
>        try:
>            # validate params
>            cores = int(param_dict['__job_resource']['cores'])
>            time = int(param_dict['__job_resource']['time'])
>            destination_id = DESTINATION_IDS[cores]
>            destination = app.job_config.get_destination(destination_id)
>            # set walltime
>            if 'nativeSpecification' not in destination.params:
>                destination.params['nativeSpecification'] = ''
>            destination.params['nativeSpecification'] += ' --time=%s:00:00' % time
>        except:
>            # resource param selector not sent with tool form, job_conf.xml misconfigured
>            log.warning('(%s) error, keys were: %s', job.id, param_dict.keys())
>            raise JobMappingException(FAILURE_MESSAGE)
>
>        log.info('returning destination: %s', destination_id)
>        return destination or destination_id
>    ```
>
>    It is important to note that **you are responsible for parameter validation, including the job resource selector**. This function only handles the job resource parameter fields, but it could do many other things - examine inputs, job queues, other tool parameters, etc.
>
>
> 2. As usual, we need to instruct Galaxy of where to find this file:
>
>    - Edit your group variables file and add the following:
>
>      ```yml
>      galaxy_dynamic_job_rules:
>        - my_rules.py
>        - map_resources.py
>      ```
>
> 3. Run the playbook, restart Galaxy
>
> 4. Run the **Multicore Tool** with various resource parameter selections
>
>    - Use default job resource parameters
>    - Specify job resource parameters:
>      - 1 core
>      - 2 cores
>      - Some value for walltime from 1-24
>
{: .hands_on}

The cores parameter can be verified from the output of the tool. The walltime can be verified with `scontrol`:

```console
$ scontrol show job 24
JobId=24 JobName=g24_multi_anonymous_10_0_2_2
   UserId=galaxy(999) GroupId=galaxy(999)
   Priority=4294901747 Nice=0 Account=(null) QOS=(null)
   JobState=COMPLETED Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:05 TimeLimit=12:00:00 TimeMin=N/A
   SubmitTime=2016-11-05T22:01:09 EligibleTime=2016-11-05T22:01:09
   StartTime=2016-11-05T22:01:09 EndTime=2016-11-05T22:01:14
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=debug AllocNode:Sid=gat2016:1860
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=localhost
   BatchHost=localhost
   NumNodes=1 NumCPUs=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
   Features=(null) Gres=(null) Reservation=(null)
   Shared=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=(null)
   WorkDir=/srv/galaxy/server/database/jobs/000/24
   StdErr=/srv/galaxy/server/database/jobs/000/24/galaxy_24.e
   StdIn=StdIn=/dev/null
   StdOut=/srv/galaxy/server/database/jobs/000/24/galaxy_24.o
   Power= SICP=0
```

Note that the `TimeLimit` for this job (which I gave a 12 hour time limit) was set to `12:00:00`.

## Further Reading

- The [sample dynamic tool destination config file](https://github.com/galaxyproject/galaxy/blob/dev/config/tool_destinations.yml.sample) fully describes the configuration language
- [Dynamic destination documentation](https://docs.galaxyproject.org/en/latest/admin/jobs.html)
- Job resource parameters are not as well documented as they could be, but the [sample configuration file](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/test/files/galaxy/config/job_resource_params_conf.xml) shows some of the possibilities.
- [usegalaxy.org's job_conf.xml](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/main/templates/galaxy/config/job_conf.xml.j2) is publicly available for reference.
