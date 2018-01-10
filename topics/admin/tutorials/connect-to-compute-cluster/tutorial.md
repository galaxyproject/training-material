---
layout: tutorial_hands_on
topic_name: admin
tutorial_name: connect-to-compute-cluster
---

## Running Galaxy Jobs with Slurm

### Introduction

The tools that are added to Galaxy can have a wide variance in the compute resources that they require and work efficiently on.
To account for this, Galaxy's job configuration needs to be tuned to run these tools properly. In addition, site-specific variables must
be taken into consideration when choosing where to run jobs and what parameters to run them with.

### Section 1 - Install and configure Slurm

**Part 1 - Install Slurm**

Install Slurm with apt:

```console
$ sudo apt-get install -y slurm-wlm
Reading package lists... Done
Building dependency tree       
Reading state information... Done
The following additional packages will be installed:
  ...
$
```

Installed with Slurm is MUNGE (MUNGE Uid 'N Gid Emporium...) which authenticates users between cluster hosts. You would normally need to ensure the same Munge key is distributed across all cluster hosts (in `/etc/munge/munge.key`) - A great task for Ansible. However, the installation of the munge package has created a random key for you, and you will not need to distribute this since you'll run jobs locally.

Verify that MUNGE is now running with `systemctl status munge`:

```console
$ systemctl status munge
● munge.service - MUNGE authentication service
   Loaded: loaded (/etc/systemd/system/munge.service; enabled; vendor preset: enabled)
   Active: active (running) since Wed 2018-01-10 13:15:33 UTC; 3min 46s ago
     Docs: man:munged(8)
 Main PID: 4805 (munged)
   CGroup: /system.slice/munge.service
           └─4805 /usr/sbin/munged --syslog

Jan 10 13:15:33 galaxy-admin-ws-big-35 systemd[1]: Starting MUNGE authentication service...
Jan 10 13:15:33 galaxy-admin-ws-big-35 systemd[1]: Started MUNGE authentication service.
$
```

You can also see that Slurm's controller and execution daemon processes are configured to start automatically, and that they attempted to start, but failed:

```console
$ systemctl status slurmctld
● slurmctld.service - Slurm controller daemon
   Loaded: loaded (/lib/systemd/system/slurmctld.service; enabled; vendor preset: enabled)
   Active: inactive (dead)
Condition: start condition failed at Fri 2016-11-04 16:05:30 EDT; 39s ago
$ systemctl status slurmd
● slurmd.service - Slurm node daemon
   Loaded: loaded (/lib/systemd/system/slurmd.service; enabled; vendor preset: enabled)
   Active: inactive (dead)
Condition: start condition failed at Fri 2016-11-04 16:05:29 EDT; 43s ago
```

The start condition that failed was the missing slurm config file.

**Part 2 - Configure Slurm**

Under Ubuntu, Slurm configs are stored in `/etc/slurm-llnl`<sup>[1]</sup>. No config is created by default.

Slurm provides a tool to create a configuration file. This is available online for the latest version, but Ubuntu 16.04 ships with Slurm 15.08. There's a copy of the configurator in `/usr/share/doc/slurmctld/slurm-wlm-configurator.html`. I've copied that to the training repository:

[Slurm Version 15.08 Configuration Tool](./slurm-wlm-configurator.html)

Enter the following values into the configuration tool (leaving others at their defaults):
- ControlMachine: `localhost`
- NodeName: `localhost`
- CPUs: 2
- SelectType: Cons_res

Then click **Submit** at the bottom of the form. You should now see the contents of a `slurm.conf` which you can copy and paste into `/etc/slurm-llnl/slurm.conf`.

We need to make one change to the configuration that was generated: Uncomment `SelectTypeParameters` and set its value to `CR_CPU`.

Your VM should have 2 CPUs allocated to it, but if you're running this exercise on a different VM that only has one core, set `FastSchedule` to `2` or else Slurm will set the node state to `DRAINING` because it does not match the configuration.

**Part 3 - Start Slurm daemons**

It should now be possible to start Slurm's daemons. Begin by starting `slurmctld`, The Slurm controller daemon (only one host runs the controller, this orchestrates job scheduling, dispatching, and completion):

```console
$ sudo systemctl start slurmctld
$ systemctl status slurmctld
● slurmctld.service - Slurm controller daemon
   Loaded: loaded (/lib/systemd/system/slurmctld.service; enabled; vendor preset: enabled)
   Active: active (running) since Fri 2016-11-04 16:46:08 EDT; 4s ago
  Process: 5134 ExecStart=/usr/sbin/slurmctld $SLURMCTLD_OPTIONS (code=exited, status=0/SUCCESS)
 Main PID: 5138 (slurmctld)
    Tasks: 10
   Memory: 884.0K
      CPU: 11ms
   CGroup: /system.slice/slurmctld.service
           └─5138 /usr/sbin/slurmctld

Nov 04 16:46:08 gat2016 systemd[1]: Starting Slurm controller daemon...
Nov 04 16:46:08 gat2016 systemd[1]: Started Slurm controller daemon.
```

Next, start up `slurmd`, the Slurm execution daemon. Every host that will execute jobs runs slurmd, which manages the processes that the slurm controller dispatches to it:

```console
$ sudo systemctl start slurmd
$ systemctl status slurmd
● slurmd.service - Slurm node daemon
   Loaded: loaded (/lib/systemd/system/slurmd.service; enabled; vendor preset: enabled)
   Active: active (running) since Fri 2016-11-04 16:50:35 EDT; 2s ago
  Process: 5169 ExecStart=/usr/sbin/slurmd $SLURMD_OPTIONS (code=exited, status=0/SUCCESS)
 Main PID: 5173 (slurmd)
    Tasks: 1
   Memory: 2.0M
      CPU: 10ms
   CGroup: /system.slice/slurmd.service
           └─5173 /usr/sbin/slurmd

Nov 04 16:50:35 gat2016 systemd[1]: Starting Slurm node daemon...
Nov 04 16:50:35 gat2016 systemd[1]: slurmd.service: PID file /var/run/slurm-llnl/slurmd.pid not readable (yet?) after start: No such
Nov 04 16:50:35 gat2016 systemd[1]: Started Slurm node daemon.
```

You should now be able to see that your slurm cluster is operational with the `sinfo` command. This shows the state of nodes and partitions (synonymous with queues in other DRMs). The "node-oriented view" provided with the `-N` flag is particularly useful:

```console
$ sinfo
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
debug*       up   infinite      1   idle localhost
$ sinfo -Nel
Fri Nov  4 16:51:24 2016
NODELIST   NODES PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT FEATURES REASON
localhost      1    debug*        idle    1    2:1:1      1        0      1   (null) none
```

If your node state is not `idle`, something has gone wrong. If your node state ends with an asterisk `*`, the slurm controller is attempting to contact the slurm execution daemon but has not yet been successful.

### Section 2 - Get Slurm ready for Galaxy

**Part 1 - Test Slurm**

We want to ensure that Slurm is actually able to run jobs. There are two ways this can be done:

- `srun`: Run an interactive job (e.g. a shell, or a specific program with its stdin, stdout, and stderr all connected to your terminal.
- `sbatch`: Run a batch job, with stdin closed and stdout/stderr redirected to a file.

Galaxy runs `sbatch` jobs but we can use both `srun` and `sbatch` to test:

```console
$ srun uname -a
Linux gat2016 4.4.0-31-generic #50-Ubuntu SMP Wed Jul 13 00:07:12 UTC 2016 x86_64 x86_64 x86_64 GNU/Linux
```

Although it looks like this command ran as if I had not used `srun`, it was in fact routed through Slurm.

Next, create a test job script somewhere, such as in `~/sbatch-test.sh`. This should be a shell script and must include the shell "shebang" line:

```bash
#!/bin/sh

uname -a
uptime
cat /etc/issue
sleep 30
```

Submit it with `sbatch` and monitor it with `squeue`:

```console
$ sbatch sbatch-test.sh
Submitted batch job 3
$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                 3     debug sbatch-t galaxygu  R       0:03      1 localhost
$ cat slurm-3.out
Linux gat2016 4.4.0-31-generic #50-Ubuntu SMP Wed Jul 13 00:07:12 UTC 2016 x86_64 x86_64 x86_64 GNU/Linux
 17:09:18 up  1:28,  2 users,  load average: 0.00, 0.00, 0.00
Ubuntu 16.04.1 LTS \n \l

$
```

If you've made it this far, your Slurm installation is working!

**Part 2 - Install slurm-drmaa**

Above Slurm in the stack sits slurm-drmaa, a library that provides a translational interface from the Slurm API to the generalized DRMAA API in C. Thankfully, Ubuntu has a package for it as well:

```console
$ sudo apt-get install slurm-drmaa1
Reading package lists... Done
Building dependency tree
Reading state information... Done
The following additional packages will be installed:
  libslurm29
The following NEW packages will be installed:
  libslurm29 slurm-drmaa1
0 upgraded, 2 newly installed, 0 to remove and 91 not upgraded.
Need to get 574 kB of archives.
After this operation, 1,676 kB of additional disk space will be used.
Do you want to continue? [Y/n]
Get:1 http://us.archive.ubuntu.com/ubuntu xenial/universe amd64 libslurm29 amd64 15.08.7-1build1 [522 kB]
Get:2 http://us.archive.ubuntu.com/ubuntu xenial/universe amd64 slurm-drmaa1 amd64 1.0.7-1build3 [52.3 kB]
Fetched 574 kB in 0s (1,102 kB/s)
Selecting previously unselected package libslurm29.
(Reading database ... 60829 files and directories currently installed.)
Preparing to unpack .../libslurm29_15.08.7-1build1_amd64.deb ...
Unpacking libslurm29 (15.08.7-1build1) ...
Selecting previously unselected package slurm-drmaa1.
Preparing to unpack .../slurm-drmaa1_1.0.7-1build3_amd64.deb ...
Unpacking slurm-drmaa1 (1.0.7-1build3) ...
Processing triggers for libc-bin (2.23-0ubuntu3) ...
Setting up libslurm29 (15.08.7-1build1) ...
Setting up slurm-drmaa1 (1.0.7-1build3) ...
Processing triggers for libc-bin (2.23-0ubuntu3) ...
$
```

### Section 3 - Run Galaxy jobs through Slurm

**Part 1 - Install DRMAA Python**

Moving one level further up the stack, we find DRMAA Python. This is a Galaxy framework *conditional dependency*. Conditional dependencies are only installed if, during startup, a configuration option is set that requires that dependency. Galaxy will automatically install it into the virtualenv if we're using the `run.sh` (which calls `scripts/common_startup.sh`) method of starting. Since our Galaxy now starts the application directly with uWSGI or the "headless" `galaxy-main`, we need to install it into Galaxy's virtualenv directly.

The `galaxyprojectdotorg.galaxy` Ansible role *does* install conditional dependencies. An alternative option would be to modify `job_conf.xml` as described in the next part and then rerun the Ansible playbook from the second exercise in the Ansible section.

Assuming we will install DRMAA Python ourselves, we must:

1. Become the `galaxy` user.
2. Run `pip` from Galaxy's virtualenv in `/srv/galaxy/venv`
3. Install the `drmaa` package from PyPI.

We can do this with a single command: `sudo -H -u galaxy /srv/galaxy/venv/bin/pip install drmaa`:

```console
$ sudo -H -u galaxy /srv/galaxy/venv/bin/pip install drmaa
Collecting drmaa
  Downloading drmaa-0.7.6-py2.py3-none-any.whl
Installing collected packages: drmaa
Successfully installed drmaa-0.7.6
You are using pip version 8.1.2, however version 9.0.0 is available.
You should consider upgrading via the 'pip install --upgrade pip' command.
$
```

**Part 2 - Configure Galaxy**

At the top of the stack sits Galaxy. Galaxy must now be configured to use the cluster we've just set up. The DRMAA Python documentation (and Galaxy's own documentation) instruct that you should set the `$DRMAA_LIBRARY_PATH` environment variable so that DRMAA Python can find `libdrmaa.so` (aka slurm-drmaa). Because Galaxy is now being started under supervisor, the environment that Galaxy starts under is controlled by the `environment` option in `/etc/supervisor/conf.d/galaxy.conf`. The `[program:handler]` should thus be updated to refer to the path to slurm-drmaa, which is `/usr/lib/slurm-drmaa/lib/libdrmaa.so.1`:

```ini
environment     = VIRTUAL_ENV="/srv/galaxy/venv",PATH="/srv/galaxy/venv/bin:%(ENV_PATH)s",DRMAA_LIBRARY_PATH="/usr/lib/slurm-drmaa/lib/libdrmaa.so.1"
```

This change is not read until `supervisord` is notified with `sudo supervisorctl update`, but we'll wait to do that until after we've updated Galaxy's job configuration.

We need to modify `job_conf.xml` to instruct Galaxy's job handlers to load the Slurm job runner plugin, and set the Slurm job submission parameters. This file was installed by Ansible and can be found in `/srv/galaxy/config` (remember, it's owned by the `galaxy` user so you'll need to use `sudo` to edit it). A job runner plugin definition must have the `id`, `type`, and `load` attributes. The entire `<plugins>` tag group should look like:

```xml
    <plugins workers="4">
        <plugin id="local" type="runner" load="galaxy.jobs.runners.local:LocalJobRunner"/>
        <plugin id="slurm" type="runner" load="galaxy.jobs.runners.slurm:SlurmJobRunner"/>
    </plugins>
```

Next, we need to add a new destination for the Slurm job runner. This is a basic destination with no parameters, Galaxy will do the equivalent of submitting a job as `sbatch /path/to/job_script.sh`. Note that we also need to set a default destination now that more than one destination is defined. In a `<destination>` tag, the `id` attribute is a unique identifier for that destination and the `runner` attribute must match the `id` of defined plugin:

```xml
    <destinations default="slurm">
        <destination id="slurm" runner="slurm"/>
        <destination id="local" runner="local"/>
    </destinations>
```

To reread the job config, Galaxy must be restarted. You can do this with `sudo supervisorctl restart all`. Technically these changes only require restarting the handlers (if we were changing a tool-to-handler mapping it'd require restarting the web server as well) so `sudo supervisorctl restart gx:handler0 gx:handler1` would suffice. However, the handlers will be restarted when we update supervisor to reread the config that we changed earlier. Do this now with `sudo supervisorctl update`

Before you restart, you can follow the handler log files using `tail`: `tail -f /srv/galaxy/log/handler?.log`.

```xml
$ sudo supervisorctl update
gx: stopped
gx: updated process group
```

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

**Part 2 - Go!**

You should now be able to run a Galaxy job through Slurm. The simplest way to test is using the upload tool to upload some text. If you're not still following the log files with `tail`, do so now.

Then, upload to Galaxy to create a new job:

1. Click the upload button at the top of the tool panel (on the left side of the Galaxy UI).
2. In the resulting modal dialog, click the "Paste/Fetch data" button.
3. Type some random characters into the text field that has just appeared.
4. Click "Start" and then "Close"

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

## So, what did we learn?

Hopefully, you now understand:
- How the various DRM and DRMAA pieces fit together to allow Galaxy to interface with a cluster
- Some basic Slurm inspection and usage: commands for other DRMs are different but the process is similar

## Further Reading

- [Galaxy's cluster documentation](https://wiki.galaxyproject.org/Admin/Config/Performance/Cluster) describes in detail alternative cluster configurations
- [The job_conf.xml documentation](https://wiki.galaxyproject.org/Admin/Config/Jobs) fully describes the syntax of the job configuration file.
- The [Distributed Resource Management Application API (DRMAA)](https://www.drmaa.org/) page contains the DRMAA specification as well as documentation for various implementations. It also includes a list of DRMs supporting DRMAA.
- The [Slurm documentation](http://slurm.schedmd.com/) is extensive and covers all the features and myriad of ways in which you can configure slurm.
- [PSNC slurm-drmaa](http://apps.man.poznan.pl/trac/slurm-drmaa)'s page includes documentation and the SVN repository, which has a few minor fixes since the last released version. PSNC also wrote the initial implementations of the DRMAA libraries for PBSPro and LSF, so all three are similar.
- [My own fork of slurm-drmaa](http://github.com/natefoo/slurm-drmaa) includes support for Slurms `-M`/`--clusters` multi-cluster functionality.
- [Slurm Accounting documentation](http://slurm.schedmd.com/accounting.html) explains how to set up SlurmDBD.

## Notes

<sup>1. The package and config directory name oddities are due to an unrelated `slurm` package existing in Debian before Slurm was added to Debian.
`slurm-llnl` refers to Lawrence Livermore National Laboratory, where Slurm was originally developed,
but the package was later renamed to `slurm-wlm` (for **W**ork**L**oad **M**anager) when the Slurm authors quit LLNL.</sup>



### Section 4 - Statically map a tool to a job destination

**Part 1 - Create a tool**

We don't want to overload our training VMs trying to run real tools, so to demonstrate how to map a multicore tool to a multicore destination, we'll create a fake tool. Since most of these operations are performed as the `galaxy` user it's probably easiest to open a shell as that user before starting:

```console
ubuntu$ sudo -su galaxy
galaxy$
```

If you prefer (as I do) to do less user switching, you can open a second shell to your VM to do this.

Now, open a new file at `/srv/galaxy/server/tools/multi.xml` and add the contents:

```xml
<tool id="multi" name="Multicore Tool">
    <command>
        echo "Running with '\${GALAXY_SLOTS:-1}' threads" &gt; "$output1"
    </command>
    <inputs>
        <param name="input1" type="data" format="txt" label="Input Dataset"/>
    </inputs>
    <outputs>
        <data name="output1" format="txt" />
    </outputs>
</tool>
```

Of course, this tool doesn't actually *use* the allocated number of cores. In a real tool, you would call the tools's underlying command with whatever flag that tool provides to control the number of threads or processes it starts, such as `foobar -t \${GALAXY_SLOTS:-1} ...`.

Up until now we've been using the default tool panel config file, located at `/srv/galaxy/server/config/tool_conf.xml.sample`. Copy this to `/srv/galaxy/config/tool_conf.xml` in the same directory as the sample and open it up with an editor. We need to add the entry for our new tool. This can go anywhere, but I suggest putting it at the very top, between the opening `<toolbox>` and first `<section>`, so that it appears right at the top of the toolbox.

```console
$ cp /srv/galaxy/server/config/tool_conf.xml.sample /srv/galaxy/config/tool_conf.xml
$ vim /srv/galaxy/config/tool_conf.xml
```

The tag to add is:

```xml
  <tool file="multi.xml"/>
```

Galaxy needs to be instructed to read `tool_conf.xml` instead of `tool_conf.xml.sample`. Normally it does this automatically if `tool_conf.xml` exists, but the Ansible role we used to install Galaxy explicitly instructed Galaxy to load `tool_conf.xml.sample`.

Edit `/srv/galaxy/config/galaxy.ini` and modify the value of `tool_config_file` accordingly:

```ini
tool_config_file = /srv/galaxy/config/tool_conf.xml,/srv/galaxy/config/shed_tool_conf.xml
```

Finally, in order to read the toolbox changes, Galaxy should be restarted. You'll need to return to the `ubuntu` user to do this (since the `galaxy` user does not have `sudo` privileges). It is, as usual, `sudo supervisorctl restart all`.

Reload Galaxy in your browser and the new tool should now appear in the tool panel. If you have not already created a dataset in your history, upload a random text dataset. Once you have a dataset, click the tool's name in the tool panel, then click Execute. When your job completes, the output should be:

```
Running with '1' threads
```

**Part 2 - Create a destination and map the tool**

We want our tool to run with more than one core. To do this, we need to instruct Slurm to allocate more cores for this job. This is done in the job configuration file.

As the `galaxy` user, open up `/srv/galaxy/config/job_conf.xml` and add the following new destination:

```xml
        <destination id="slurm-2c" runner="slurm">
            <param id="nativeSpecification">--nodes=1 --ntasks=2</param>
        </destination>
```

Then, map the new tool to the new destination using the tool ID (`<tool id="multi">`) and destination id (`<destination id="slurm-2c">`) by adding a new section to the job config, `<tools>`:

```xml
    <tools>
        <tool id="multi" destination="slurm-2c"/>
    </tools>
```

And finally, restart Galaxy with `sudo supervisorctl restart all` (as the `ubuntu` user).

Now, click the rerun button on the last history item, or click **Multicore Tool** in the tool panel, and then click the tool's Execute button. If successful, your tool's output should now be:

```
Running with '2' threads
```

### Section 2 - Dynamically map a tool to a job destination

**Part 1 - Write a Dynamic Tool Destination**

Dynamic tool destinations utilize the dynamic job runner to provide dynamic job mapping functionality without having to explicitly write code to perform the mapping. The mapping functionality is mostly limited to input sizes, but often input size is the most important factor in deciding what resources to allocate for a job.

Dynamic tool destinations are configured via a YAML file at `/srv/galaxy/config/tool_destinations.yml`. As before, we'll use a fake example. Create the file with the following contents:

```yaml
---
tools:
  multi:
    rules:
      - rule_type: file_size
        lower_bound: 16
        upper_bound: Infinity
        destination: slurm-2c
    default_destination: slurm
default_destination: local
verbose: True
```

The rule says:
- If the tool has ID `multi`:
  - If the input dataset is >=16 bytes, run on the destination `slurm-2c`
  - If the input dataset is <16 bytes, run on the destination `slurm`
- Else, run on the destination `local`

We also need to inform Galaxy of the path to the file we've just created, which is done using the `tool_destinations_config_file` in `galaxy.ini`:

```ini
tool_destinations_config_file = /srv/galaxy/config/tool_destinations.yml
```

Once the dynamic tool definition has been written, we need to update Galaxy's job configuration to use this rule. Open `/srv/galaxy/config/job_conf.xml` and add a DTD destination:

```xml
        <destination id="dtd" runner="dynamic">
            <param id="type">dtd</param>
        </destination>
```

Also, comment out the previous `<tool>` definition for the `multi` tool, and replace it with a mapping to the dtd destination like so:

```xml
    <tools>
<!--
        <tool id="multi" destination="slurm-2c"/>
-->
        <tool id="multi" destination="dtd"/>
    </tools>
```

Then, restart Galaxy with `sudo supervisorctl restart all`.

**Part 2 - Verify**

Our rule specified that any invocation of the `multi` tool with an input dataset with size <16 bytes would run on the 1 core destination, whereas any with >= 16 bytes would run on the 2 core destination. To verify, create a dataset using the upload paste tool of just a few (<16) characters, and another with >16 characters and run the Multicore Tool on each. The former will run "with '1' thread" whereas the latter will run "with '2' threads".

## Section 5 - Implement a job resource selector

**Part 1 - Define the resource selector**

You may find that certain tools can benefit from having form elements added to them to allow for controlling certain job parameters, so that users can select based on their own knowledge. For example, a user might know that a particular set of parameters and inputs to a certain tool needs a larger memory allocation than the standard amount for a given tool. This of course assumes that your users are well behaved enough not to choose the maximum whenever available, although such concerns can be mitigated somewhat by the use of concurrency limits on larger memory destinations.

Such form elements can be added to tools without modifying each tool's configuration file through the use of the **job resource parameters configuration file**, `/srv/galaxy/config/job_resource_params_conf.xml`. Create this file and add the following contents:

```xml
<parameters>
    <param label="Cores" name="cores" type="select" help="Number of cores to run job on.">
        <option value="1">1 (default)</option>
        <option value="2">2</option>
    </param>
  <param label="Time" name="time" type="integer" size="3" min="1" max="24" value="1" help="Maximum job time in hours, 'walltime' value (1-24). Leave blank to use default value." />
</parameters>
```

This defines two resource fields, a select box where users can choose between 1 and 2 cores, and a text entry field where users can input an integer value from 1-24 to set the walltime for a job.

As usual, we need to instruct Galaxy of where to find this file in `galaxy.ini` using the `job_resource_params_file` option:

```ini
job_resource_params_file = /srv/galaxy/config/job_resource_params_conf.xml
```

**Part 2 - Configure Galaxy to use the resource selector**

Next, we define a new section in `job_conf.xml`: `<resources>`. This groups together parameters that should appear together on a tool form. Add the following section to `/srv/galaxy/server/job_conf.xml`:

```xml
    <resources>
        <group id="multi_resources">cores,time</group>
    </resources>
```

The group ID will be used to map a tool to job resource parameters, and the text value of the `<group>` tag is a comma-separated list of `name`s from `job_resource_params_conf.xml` to include on the form of any tool that is mapped to the defined `<group>`.

Finally, in `job_conf.xml`, move the previous `<tool>` definition for the `multi` tool into the comment and define a new `<tool>` that defines the `resources` for the tool:

```xml
    <tools>
<!--
        <tool id="multi" destination="slurm-2c"/>
        <tool id="multi" destination="dtd"/>
-->
        <tool id="multi" destination="dynamic_cores_time" resources="multi_resources"/>
    </tools>
```

We have assigned the `multi` tool to a new destination: `dynamic_cores_time`, but this destination does not exist. We need to create it. Add the following destination:

```xml
        <destination id="dynamic_cores_time" runner="dynamic">
            <param id="type">python</param>
            <param id="function">dynamic_cores_time</param>
        </destination>
```

This is a **Python function dynamic destination**. Galaxy will load a function from `/srv/galaxy/server/lib/galaxy/jobs/rules/*.py` named `dynamic_cores_time` and that function will determine the job destination for this tool.

**Part 3 - Python function dynamic rule**

Lastly, we need to write the rule that will read the value of the job resource parameter form fields and decide how to submit the job. But first, let's see the fruits of our labor thus far. Restart Galaxy with `sudo supervisorctl restart all`, then click the **Multicore Tool**. You should see that a "Job Resource Parameters" select box has been added to the bottom of the tool form. If this is switched to "**Specify job resource parameters**", the fields that were defined in `job_resource_params_conf.xml` are displayed.

We need to write a Python function that will process these rules. Such rules live, by default, in `/srv/galaxy/server/lib/galaxy/jobs/rules/*.py` (although this can be configured). Create a new file, `/srv/galaxy/server/lib/galaxy/jobs/rules/cores_time.py`. This file should contains the function that we named in `job_conf.xml`: `dynamic_cores_time`:

```python
import logging
from galaxy.jobs.mapper import JobMappingException

log = logging.getLogger(__name__)

DESTINATION_IDS = {
    1 : 'slurm',
    2 : 'slurm-2c'
}
FAILURE_MESSAGE = 'This tool could not be run because of a misconfiguration in the Galaxy job running system, please report this error'


def dynamic_cores_time(app, tool, job, user_email):
    destination = None
    destination_id = 'slurm'

    # build the param dictionary
    param_dict = dict( [ ( p.name, p.value ) for p in job.parameters ] )
    param_dict = tool.params_from_strings( param_dict, app )

    # handle job resource parameters
    if '__job_resource' in param_dict:
        if param_dict['__job_resource']['__job_resource__select'] == 'yes':
            try:
                # validate params
                cores = int(param_dict['__job_resource']['cores'])
                time = int(param_dict['__job_resource']['time'])
                assert cores in (1, 2), "Invalid value for core selector"
                assert time in range(1, 25), "Invalid value for time selector"
            except (TypeError, AssertionError) as exc:
                log.exception(exc)
                log.error('(%s) param_dict was: %s', job.id, param_dict)
                raise JobMappingException( FAILURE_MESSAGE )
            # params validated
            destination_id = DESTINATION_IDS[cores]
            destination = app.job_config.get_destination(destination_id)
            # set walltime
            if 'nativeSpecification' not in destination.params:
                destination.params['nativeSpecification'] = ''
            destination.params['nativeSpecification'] += ' --time=%s:00:00' % time
        elif param_dict['__job_resource']['__job_resource__select'] != 'no':
            # someone's up to some shennanigans
            log.error('(%s) resource selector not yes/no, param_dict was: %s', job.id, param_dict)
            raise JobMappingException( FAILURE_MESSAGE )
    else:
        # resource param selector not sent with tool form, job_conf.xml misconfigured
        log.warning('(%s) did not receive the __job_resource param, keys were: %s', job.id, param_dict.keys())
        raise JobMappingException( FAILURE_MESSAGE )

    if destination is not None and 'nativeSpecification' in destination.params:
        log.info("native specification: %s", destination.params['nativeSpecification'])
    log.info('returning destination: %s', destination_id)
    return destination or destination_id
```

It is important to note that **you are responsible for parameter validation, including the job resource selector**. This function only handles the job resource parameter fields, but it could do many other things - examine inputs, job queues, other tool paramters, etc.

Once written, restart Galaxy with `sudo supervisorctl restart all`.

**Part 4 - Verify***

Run the **Multicore Tool* with various resource parameter selections:
- Use default job resource parameters
- Specify job resource parameters:
  - 1 core
  - 2 cores
  - Some value for walltime from 1-24

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

### So, what did we learn?

Hopefully, you now understand:
- The various ways in which tools can be mapped to destinations, both statically and dynamically
- How to write a dynamic tool destination (DTD)
- How to write a dynamic python function destination
- How to use the job resource parameter selection feature

### Further Reading

- The [sample dynamic tool destination config file](https://github.com/galaxyproject/galaxy/blob/dev/config/tool_destinations.yml.sample) fully describes the configuration language
- [Dynamic destination documentation](https://wiki.galaxyproject.org/Admin/Config/Jobs#Dynamic_Destination_Mapping)
- Job resource parameters are not as well documented as they could be, but the [sample configuration file](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/test/files/galaxy/config/job_resource_params_conf.xml) shows some of the possibilities.
- [usegalaxy.org's job_conf.xml](https://github.com/galaxyproject/usegalaxy-playbook/blob/master/env/main/templates/galaxy/config/job_conf.xml.j2) is publicly available for reference.



