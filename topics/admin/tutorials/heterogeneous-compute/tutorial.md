![galaxy logo](../../docs/shared-images/galaxy_logo_25percent_transparent.png)

### Galaxy Administrators Course

# Running jobs on remote resources using Pulsar - Exercise

#### Authors: Nate Coraor. 2017, Marius van den Beek. 2018, Simon Gladman. 2019

## Learning Outcomes

By the end of this tutorial, you should:

1. Have an understanding of what Pulsar is and how it works
2. Install and configure a Pulsar server on a remote linux machine
3. Be able to get Galaxy to send jobs to a remote Pulsar server

## Introduction

A slideshow presentation on this subject can be found [here](https://galaxyproject.github.io/dagobah-training/2019-pennstate/17-heterogeneous/heterogeneous.md)

More details on Pulsar can be found at:

* Pulsar Read-the-docs
    * [https://pulsar.readthedocs.io/en/latest/index.html](https://pulsar.readthedocs.io/en/latest/index.html)
* Pulsar on galaxyproject.org
    * [https://galaxyproject.org/admin/config/pulsar/](https://galaxyproject.org/admin/config/pulsar/)
* Pulsar Github
    * [https://github.com/galaxyproject/pulsar](https://github.com/galaxyproject/pulsar)
* Pulsar Ansible
    * [https://github.com/galaxyproject/ansible-pulsar](https://github.com/galaxyproject/ansible-pulsar)

Pulsar is the Galaxy Project's remote job running system. It was written by John Chilton (@jmchilton) of the Galaxy Project. It is a python server application that can accept jobs from a Galaxy server, submit them to a local resource and then send the results back to the originating Galaxy server.

Transport of data, tool information and other metadata can be configured as a web application via a RESTful interface or using a message passing system such as RabbitMQ.

At the Galaxy end, it is configured within the `job_conf.xml` file and uses one of two special Galaxy job runners.
* `galaxy.jobs.runners.pulsar:PulsarRESTJobRunner` for the RESTful interface
* `galaxy.jobs.runners.pulsar:PulsarMQJobRunner` for the message passing interface.

In this tutorial, we will:

* Install and configure a Pulsar server on a remote linux machine using ansible
    * We will configure the Pulsar server to run via the RESTful interface
* Configure our Galaxy servers to run a job there
* Run a job remotely

You have been assigned a remote linux instance for you to install Pulsar onto. (located in Sydney, Australia - it's actually remote.. :) ) These instances have 2 vcpus, 8GB RAM, a 10GB root disk and a 60GB disk located at `/mnt`. We will install Pulsar onto the `/mnt` disk.

The ip address of your instance is shown in the column **Pulsar Host** of the original spreadsheet where the instances for the course were assigned.

The demonstrators/instructors will pass around a ssh private key to allow you to access your instance (they are not configured to use a password.)

## Section 1: Installing Pulsar

We need to create a new ansible playbook to install Pulsar. We will be using a *role* written by Nate Coraor - `galaxyproject.pulsar`

#### Step 1: Install the `galaxyproject.pulsar` ansible role

From your ansible working directory, edit the `requirements.yml` file and add the following line:

```yaml
- src: galaxyproject.pulsar
```

Now install it with:

```bash
ansible-galaxy install -p roles -r requirements.yml
```

#### Step 2: Create a new group_vars file

From the galaxyproject.pulsar ansible role documentation, we need to specify some variables.

There is one required variable:

`pulsar_server_dir` - The location in which to install pulsar

Then there are a lot of optional variables. They are listed here for information. We will set some for this tutorial but not all.

| Variable Name | Description |
|---------------|-------------|
| `pulsar_pip_install` | Set to `true` to get pulsar to be sourced from pip. (Default = `false`) |
| `pulsar_yaml_config` | a YAML dictionary whose contents will be used to create Pulsar's `app.yml` |
| `pulsar_git_repo` | (default: https://github.com/galaxyproject/pulsar): Upstream git repository from which Pulsar should be cloned. |
| `pulsar_changeset_id` | (default: master): A changeset id, tag, branch, or other valid git identifier for which changeset Pulsar should be updated to. This is also possible when installing from pip (pulsar will be installed using the git+https:// pip scheme). |
| `pulsar_venv_dir` | (default: `<pulsar_server_dir>/venv` if installing via pip, `<pulsar_server_dir>/.venv` if not): The role will create a [virtualenv] from which Pulsar will run, this controls where the virtualenv will be placed. |
| `pulsar_config_dir` | (default: `<pulsar_server_dir>/config` if installing via pip, `<pulsar_server_dir>` if not): Directory that will be used for Pulsar configuration files. |
| `pulsar_optional_dependencies` | (default: None): List of optional dependency modules to install. Whether or not you need these depends on what features you are enabling. |
| `pulsar_install_environments` | Installing dependencies may require setting certain environment variables to compile successfully. |


Additional options from Pulsar's server.ini are configurable via the following variables (these options are explained in the Pulsar documentation and server.ini.sample):

| Variable | Description |
|----------|-------------|
| `pulsar_host` | (default: `localhost`) Though in our situation, we will set it to `0.0.0.0` This is the interface pulsar will listen to. |
| `pulsar_port` | (default: `8913`) |
| `pulsar_uwsgi_socket` | (default: if unset, uWSGI will be configured to listen for HTTP requests on pulsar_host port pulsar_port): If set, uWSGI will listen for uWSGI protocol connections on this socket. |
| `pulsar_uwsgi_options` | (default: empty hash): Hash (dictionary) of additional uWSGI options to place in the [uwsgi] section of server.ini |


Some of the other options we will be using are:

* We are going to run in RESTful mode so we will need to specify a `private_token` variable so we can "secure" the connection. (For a given value of "secure".)
* We will be using the uwsgi web server to host the RESTful interface.
* We will set the tool dependencies to rely on **conda** for tool installs.

Create a new file in `group_vars` called `pulsarservers.yml` and set some of the above variables as well as some others.

**Replace <some_really_long_string_here> with a long randomish (or not) string.**

```yaml
pulsar_server_dir: /mnt/pulsar/server
pulsar_venv_dir: /mnt/pulsar/venv
pulsar_config_dir: /mnt/pulsar/config
pulsar_staging_dir: /mnt/pulsar/staging
pulsar_pip_install: true

pulsar_host: 0.0.0.0
pulsar_port: 8913

private_token: '<some_really_long_string_here>'

pulsar_optional_dependencies:
  - pyOpenSSL
  - psutil
  - pycurl
  - requests
  - poster

pulsar_yaml_config:
  dependency_resolvers_config_file : dependency_resolvers_conf.xml
  conda_auto_init : True
  conda_auto_install : True
  staging_directory: "{{ pulsar_staging_dir }}"
  private_token: "{{ private_token }}"
```

#### Step 3: Add `[pulsarservers]` to the `hosts`  file

Add the following lines to your `hosts` file:

```ini
[pulsarservers]
<ip_address of your pulsar server>
```

#### Step 4: Write a playbook for the Pulsar install

We will now write a new playbook for the pulsar installation similar to the one we did for the CVMFS installation earlier in the week.

We need to include a couple of pre-tasks to install python-virtualenv, python-pip and git etc.

Create a `pulsar_playbook.yml` file with the following contents:

```yaml
- hosts: pulsarservers
  pre_tasks:
    - name: Install some packages
      apt:
        name: "{{ item }}"
        state: installed
        update_cache: yes
      become: yes
      with_items:
        - build-essential
        - vim
        - git
        - python-dev
        - libcurl4-openssl-dev
        - libssl-dev
        - virtualenv
    - name: chown the /mnt dir to ubuntu
      file:
        path: /mnt
        owner: ubuntu
        group: ubuntu
        mode: 0755
      become: yes
  roles:
    - galaxyproject.pulsar
```

There are a couple of *pre-tasks* here. This is because we need to install some base packages on these very vanilla ubuntu instances as well as give ourselves ownership of the directory we are installing into.

#### Step 5: Create the `dependency_resolvers_conf.xml` file

We also need to create the dependency resolver file so pulsar knows how to find and install dependencies for the tools we ask it to run. The simplest method which covers 99% of the use cases is to use conda auto installs similar to how Galaxy works. We need to create the file and put it where the `galaxyproject.pulsar` role can find it.

Create a `templates` directory in your working directory.

```bash
mkdir templates
```

Create a `dependency_resolvers_conf.xml.j2` file inside the `templates` directory with the following contents:

```xml
<dependency_resolvers>
    <conda auto_install="True" auto_init="True"/>
</dependency_resolvers>
```

This tells pulsar to **only** look for dependencies in conda.

**NOTE: If the tool you want to run on Pulsar doesn't have a conda package, you will need to make alternative arrangements! This is complex and beyond our scope here. See the Pulsar documentation for details.**

#### Step 6: Run the playbook

Run the playbook we just created. We will need to supply the `ansible-playbook` command with the private key for the connection. Your demonstrator/instructor will be able to give you the key.

Therefore, we will need the `--private_key` switch.

```bash
ansible-playbook --private_key pulsar-demo-key.pem -i hosts pulsar_playbook.yml
```

After the script has run, pulsar will be installed on the remote machines!

Log in to the machines using your private key and have a look in the `/mnt/pulsar` directory. You will see the venv and config directories. All the config files created can be perused.

#### Step 7: Start Pulsar

Pulsar can now be started. We do that by activating the virtualenv, then running it and telling it which webserver to use (Paster) and where the config files are.

ssh into your Pulsar instance as the ubuntu user with your private key.

Then run the Pulsar server.

```bash
cd /mnt/pulsar
source venv/bin/activate
pulsar -m paster -c config/
```

A log will now start scrolling, showing the startup of pulsar. You'll notice that it will be initializing and installing conda. Once it is complete and listening on the assigned port, you can stop the server with Ctrl-C.

This is a pretty gross way of running Pulsar and it can run in daemon mode. To do that, add `--daemon` to the end of the above command line.

Stopping the daemon is as simple as re-running the command with `--stop-daemon` instead.

Finally: Make sure your Pulsar is running in `--daemon` mode.

## Section 2: Configuring Galaxy

Now we have a Pulsar server up and running, we need to tell our Galaxy about it.

Galaxy talks to the Pulsar server via it's `job_conf.xml` file. We need to let Galaxy know about Pulsar there and make sure Galaxy has loaded the requisite job runner, and has a destination set up.

There are three things we need to do here:

* We will need to create a job runner which uses the  `galaxy.jobs.runners.pulsar:PulsarRESTJobRunner` code.
* Create job destination which references the above job runner.
* Tell Galaxy which tools to send to the job destination.

(We will use bwa-mem)

## Step 1: Edit the job_conf

In your `job_conf.xml` file add the following job runner to the `<plugins>` section:

```xml
<plugin id="pulsar_runner" type="runner" load="galaxy.jobs.runners.pulsar:PulsarRESTJobRunner" />
```

Then add the Pulsar destination. We will need the ip address of your pulsar server and the private_token string you used when you created it.

Add the following to the `<destinations>` section of your `job_conf.xml` file:

```xml
        <destination id="pulsar" runner="pulsar_runner" >
            <param id="url">http://your_ip_address_here:8913/</param>
            <param id="private_token">your_private_token_here</param>
        </destination>
```

Finally we need to tell Galaxy which tools to send to Pulsar. We will tell it to send bwa-mem jobs to it. We use the `<tools>` section of the `job_conf.xml` file.

We need to know the full id of the tool in question, we can get this out of the `integrated_tool_panel.xml` file in the `mutable-config` directory. Then we tell Galaxy which destination to send it to (pulsar).

Add the following to the end of the `job_conf.xml` file (inside the `<tools>` section if it exists or create it if it doesn't.)

```xml
    <tools>
        <tool id="toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.1" destination="pulsar"/>
    </tools>
```

That's it. Now we need to restart Galaxy and run a job.

#### Step 2: Restart Galaxy and run a job

Restart Galaxy:

```bash
sudo supervisorctl restart galaxy
```

#### Step 3: Upload some data and run bwa-mem

Now we will upload a small set of data to run bwa-mem with.

Upload the following files from zenodo.

| File URL | filetype |
|----------|----------|
| `https://zenodo.org/record/582600/files/mutant_R1.fastq` | fastqsanger |
| `https://zenodo.org/record/582600/files/mutant_R2.fastq` | fastqsanger |
| `https://zenodo.org/record/582600/files/wildtype.fna` | fasta |

Now run the bwa-mem tool with the `wildtype.fna` file as the reference, and the two fastq files as the paired end reads. Leave everything default.

As soon as you press *execute* Galaxy will send the job to the pulsar server. You can watch the log in Galaxy using:

```
sudo supervisorctl tail -f galaxy stderr
```

You can watch the log in Pulsar by ssh'ing to it and tailing the log file with:

```
tail -f /mnt/pulsar/paster.log
```

You'll notice that the Pulsar server has received the job (all the way in Sydney!) and now should be installing bwa-mem via conda. Once this is complete (which may take a while - first time only) the job will run and the results will be returned to Galaxy!

How awesome is that? :)
