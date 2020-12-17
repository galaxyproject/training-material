---
layout: tutorial_hands_on
redirect_from:
- /topics/admin/tutorials/heterogeneous-compute/tutorial

title: "Running Jobs on Remote Resources with Pulsar"
questions:
  - How does pulsar work?
  - How can I deploy it?
objectives:
  - Have an understanding of what Pulsar is and how it works
  - Install and configure a RabbitMQ message queueing server
  - Install and configure a Pulsar server on a remote linux machine
  - Be able to get Galaxy to send jobs to a remote Pulsar server
time_estimation: "60m"
key_points:
  - Pulsar allows you to easily add geographically distributed compute resources into your Galaxy instance
  - It also works well in situations where the compute resources cannot share storage pools.
contributors:
  - natefoo
  - slugger70
  - mvdbeek
  - hexylena
subtopic: features
tags:
  - jobs
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - connect-to-compute-cluster
      - cvmfs
  - title: "A server/VM on which to deploy Pulsar"
    type: "none"
---


# Overview
{:.no_toc}

Pulsar is the Galaxy Project's remote job running system. It was written by John Chilton (@jmchilton) of the Galaxy Project. It is a python server application that can accept jobs from a Galaxy server, submit them to a local resource and then send the results back to the originating Galaxy server.

More details on Pulsar can be found at:

- [Pulsar's Documentation](https://pulsar.readthedocs.io/en/latest/index.html)
- [Pulsar's Github Repository](https://github.com/galaxyproject/pulsar)
- [Pulsar Ansible Role](https://github.com/galaxyproject/ansible-pulsar)


Transport of data, tool information and other metadata can be configured as a web application via a RESTful interface or using a message passing system such as RabbitMQ.

At the Galaxy end, it is configured within the `job_conf.xml` file and uses one of two special Galaxy job runners.
* `galaxy.jobs.runners.pulsar:PulsarRESTJobRunner` for the RESTful interface
* `galaxy.jobs.runners.pulsar:PulsarMQJobRunner` for the message passing interface.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


This tutorial assumes that you have:

- A VM or machine where you will install Pulsar, and a directory in which the installation will be done. This tutorial assumes it is `/mnt`
- That you have completed the "Galaxy Installation with Ansible" and CVMFS tutorials (Job configuration tutorial is optional) and have access to the VM/computer where it is installed.

# Overview

We will be installing the RabbitMQ server daemon onto the Galaxy server to act as an intermediary message passing system between Galaxy and the remote Pulsar. The figure below shows a schematic representation of the system.

![Schematic diagram of Galaxy communicating with a Pulsar server via the RabbitMQ server](../../images/pulsar_amqp_schema.png "Schematic diagram of Galaxy communicating with a Pulsar server via the RabbitMQ server. Red arrows represent AMQP communications and blue represent file transfers (initiated by the Pulsar server.)")


**How it will work:** 

1. Galaxy will send a message to the RabbitMQ server on the Pulsar server's particular queue saying that there is a job to be run and then will monitor the queue for job status updates.
2. The Pulsar server monitors this queue and when the job appears it will take control of it.
3. The Pulsar server will then download the required data etc. from the Galaxy server using `curl`.
4. The Pulsar server will install any required tools/tool dependencies using Conda.
5. The Pulsar server will start running the job using it's local mechanism and will send a message to the "queue" stating that the job has started.
6. Once the job has finished running, the Pulsar server will send a message to the queue stating that the job has finished.
7. Pulsar then sends the output data etc. back to the Galaxy server by `curl` again.
8. The Galaxy server acknowledges the job status and closes the job.

**Some notes:**

* RabbitMQ uses the Advanced Message Queueing Protocol (AMQP) to communicate with both the Galaxy server and the remote Pulsar VM.
* Transport of files, meta-data etc. occur via `curl` from the Pulsar end.
* RabbitMQ is written in erlang and does not add much overhead to the Galaxy VM, although in larger installations, RabbitMQ is commonly installed on a separate VM to Galaxy. e.g. Galaxy Europe, Galaxy Main and Galaxy Australia.

> ### {% icon tip %} Tip: Other file transport methods for Pulsar
>
>  Pulsar can use a variety of file transport methods including:
>  * Default: Galaxy initiates file transfer and stages files to Pulsar via http transfer.
>      * This requires that a http transfer port be open on the remote Pulsar.
>  * Remote transfer: Pulsar initiates file transfer. This can use a variety of lso available and can use a variety of methods:
>      * Curl
>      * Rsync
>      * Http
>
>  We use remote transfer using **Curl** here so we don't need an open port on the Pulsar server and tranfer robustness respectively.
>
{: .tip}


> ### {% icon details %} Why are we using Pulsar in MQ mode here and not the RESTful interface?
> We are teaching you to install Pulsar and configure it in MQ mode in this tutorial. Configuring Pulsar in RESTful mode is also possible and is quite useful in certain situations. However, in the most common situation MQ mode is preferable for a number of reasons:
> * When running Pulsar in RESTful mode, all of the job control and data transfer is controlled by the Galaxy server usually using http transfers. This can place a limit on the size of files that can be transferred without constant configuring of the webserver.
> * When running in RESTful mode, you also need to run a https server such as nginx, including securing it, configuring it, getting certificates and opening ports. This can be very difficult to do if you are attempting to submit jobs to an institutional HPC where the admins probably won't let you do any of these things.
> * In MQ mode, you only need to open a port for the RabbitMQ server on a machine you are more likely to control.
>
> See the [Pulsar documentation](https://pulsar.readthedocs.io/en/latest/) for details.
{: .details}

# Install and configure a message queueing system

In this section we will install the RabbitMQ server on your Galaxy server VM.

RabbitMQ is an AMQP server that can queue messages between systems for all sorts of reasons. Here, we will be using the queue so that Galaxy and Pulsar can communicate jobs, job status and job metadata between them easily and robustly. More information on RabbitMQ can be found [here](https://www.rabbitmq.com/).

## Installing the RabbitMQ role

Firstly we will add and configure another *role* to our Galaxy playbook - a slightly modified version of `jasonroyle.rabbitmq`

> ### {% icon hands_on %} Hands-on: Install the modified `jasonroyle.rabbitmq` ansible role
>
> 1. From your ansible working directory, edit the `requirements.yml` file and add the following line:
>
>    ```yaml
>    - src: https://github.com/slugger70/ansible-role-rabbitmq
>      version: 0.0.5
>      name: jasonroyle.rabbitmq
>    ```
>
> 2. Now install it with:
>
>    ```bash
>    ansible-galaxy install -p roles -r requirements.yml
>    ```
>
{: .hands_on}

## Configuring RabbitMQ

We need to configure RabbitMQ to be able to handle Pulsar messages. To do this we will need to create some queues, Rabbit users, some queue vhosts and set some passwords. We also need to configure rabbit to listen on various interfaces and ports.

### Defining Virtual Hosts

Each set of queues in RabbitMQ are grouped and accessed via virtual hosts. We need to create one of these for the transactions between the Galaxy server and Pulsar server. They are set as an array under the `rabbitmq_vhosts` variable. 

### Defining users

Users need to be defined, given passwords and access to the various queues. We will need to create a user that can access this vhost. We will also create an admin user. The queue will need access to the Pulsar queue vhost. They are set as an array under the `rabbitmq_users` variable with the following structure:


```yaml
rabbitmq_users:
  - user: username
    password: somelongpasswordstring
    vhost: /vhostname
```

Optional: You can add tags to each user if required. e.g. For an admin user it could be useful to add in a *administrator* tag.

### RabbitMQ server config

We also need to set some RabbitMQ server configuration variables. Such as where its security certificates are and which ports to listen on (both via localhost and network).

**NOTE: We will need to make sure that the RabbitMQ default port is open and accessible on the server we are installing RabbitMQ onto. (In our case this is the Galaxy server). Default port number is: 5671**

More information about the rabbitmq ansible role can be found [here](https://github.com/Slugger70/ansible-role-rabbitmq).

## Add RabbitMQ configuration to Galaxy VM.

> ### {% icon hands_on %} Hands-on: Add RabbitMQ settings to Galaxy VM groupvars file.
>
> 1. Create or edit the file `group_vars/all.yml` and set your private token:
>
>    ```yaml
>    rabbitmq_password_galaxy_au: areallylongpasswordhere
>    ```
>
>    This is going in a special file because all (two) of our services need it. Both Galaxy in the job configuration, and Pulsar in its configuration. The `group_vars/all.yml` is included for every playbook run, no matter which group a machine belong to.
>
>    Replace `areallylongpasswordhere` with a long randomish (or not) string.
>
> 2. From your ansible working directory, edit the `group_vars/galaxy.yml` file and add the following lines:
>
>    ```yaml
>    rabbitmq_admin_password: somereallylongpasswordhere
>
>    rabbitmq_version: 3.8.9-1
>    rabbitmq_plugins: rabbitmq_management
>
>    rabbitmq_config:
>      - ssl_listeners:
>        - "'0.0.0.0'": 5671
>        - "'127.0.0.1'": 5671
>      - ssl_options:
>         - cacertfile: /etc/ssl/certs/fullchain.pem
>         - certfile: /etc/ssl/certs/cert.pem
>         - keyfile: /etc/ssl/user/privkey-rabbitmq.pem
>         - fail_if_no_peer_cert: 'false'
>
>    rabbitmq_vhosts:
>      - /pulsar/galaxy_au
>
>    rabbitmq_users:
>      - user: admin
>        password: "{{ rabbitmq_admin_password }}"
>        tags: administrator
>        vhost: /
>      - user: galaxy_au
>        password: "{{ rabbitmq_password_galaxy_au }}"  #This password is set in group_vars/all.yml
>        vhosts: /pulsar/galaxy_au        
>        
>    ```
>
> 3. Update the Galaxy playbook to include the *jasonroyle.rabbitmq* role.
>
> 4. Run the playbook.
>
> The rabbitmq server daemon will have been installed on your Galaxy VM. You can check it's running with `systemctl status rabbitmq-server`.
>
{: .hands_on}



# Installing and configuring Pulsar on a remote machine

Now that we have a message queueing system running on our Galaxy VM, we need to install and configure Pulsar on our remote compute VM. To do this we need to create a new ansible playbook to install Pulsar. 

## Install the Pulsar ansible role

We will be using a *role* developed by the Galaxy community - `galaxyproject.pulsar`

> ### {% icon hands_on %} Hands-on: Install the `galaxyproject.pulsar` ansible role
>
> 1. From your ansible working directory, edit the `requirements.yml` file and add the following line:
>
>    ```yaml
>    - src: galaxyproject.pulsar
>      version: 1.0.6
>    ```
>
> 2. Now install it with:
>
>    ```bash
>    ansible-galaxy install -p roles -r requirements.yml
>    ```
>
{: .hands_on}


## Configuring Pulsar

From the [`galaxyproject.pulsar` ansible role documentation](https://github.com/galaxyproject/ansible-pulsar#role-variables), we need to specify some variables.

There is one required variable:

`pulsar_server_dir` - The location in which to install pulsar

Then there are a lot of optional variables. They are listed here for information. We will set some for this tutorial but not all.

 Variable Name                 | Description                                                                                        | Default
---------------                | -------------                                                                                      | ---
`pulsar_yaml_config`           | a YAML dictionary whose contents will be used to create Pulsar's `app.yml`                         |
`pulsar_venv_dir`              | The role will create a virtualenv from which Pulsar will run                                       | `<pulsar_server_dir>/venv` if installing via pip, `<pulsar_server_dir>/.venv` if not.
`pulsar_config_dir`            | Directory that will be used for Pulsar configuration files.                                        | `<pulsar_server_dir>/config` if installing via pip, `<pulsar_server_dir>` if not
`pulsar_optional_dependencies` | List of optional dependency modules to install, depending on which features you are enabling.      | `None`
`pulsar_install_environments`  | Installing dependencies may require setting certain environment variables to compile successfully. |
`pulsar_create_user`           | Should a user be created for running pulsar?                                                       |
`pulsar_user`                  | Define the user details                                                                            |

Some of the other options we will be using are:

* We will set the tool dependencies to rely on **conda** for tool installs.

* You will need to know the FQDN or IP address of the Galaxy server VM that you installed RabbitMQ on.

> ### {% icon hands_on %} Hands-on: Configure pulsar group variables
>
>
> 2. Create a new file in `group_vars` called `pulsarservers.yml` and set some of the above variables as well as some others.
>
>    {% raw %}
>    ```yaml
>    galaxy_server_url: #please put the your Galaxy server's fqdn or ip address here (or the fqdn or ip address of the RabbitMQ server).
>
>    pulsar_root: /mnt/pulsar
>    pulsar_server_dir: /mnt/pulsar/server
>    pulsar_venv_dir: /mnt/pulsar/venv
>    pulsar_config_dir: /mnt/pulsar/config
>    pulsar_staging_dir: /mnt/pulsar/files/staging
>    pulsar_persistence_dir: /mnt/pulsar/files/persistent_data
>    pulsar_dependencies_dir: /mnt/pulsar/deps
>
>    pulsar_pip_install: true
>    pulsar_pycurl_ssl_library: openssl
>    pulsar_systemd: true
>    pulsar_systemd_runner: webless
>
>    pulsar_create_user: true
>    pulsar_user: {name: ubuntu, shell: /bin/bash}
>
>    pulsar_optional_dependencies:
>      - pyOpenSSL
>      # For remote transfers initiated on the Pulsar end rather than the Galaxy end
>      - pycurl
>      # drmaa required if connecting to an external DRM using it.
>      - drmaa
>      # kombu needed if using a message queue
>      - kombu
>      # psutil and pylockfile are optional dependencies but can make Pulsar
>      # more robust in small ways.
>      - psutil
>
>    pulsar_yaml_config:
>      dependency_resolvers_config_file: dependency_resolvers_conf.xml
>      conda_auto_init: True
>      conda_auto_install: True
>      staging_directory: "{{ pulsar_staging_dir }}"
>      persistence_directory: "{{ pulsar_persistence_dir }}"
>      # The following are the settings for the pulsar server to contact the message queue with related timeouts etc.
>      message_queue_url: "pyamqp://galaxy_au:{{ rabbitmq_password_galaxy_au }}@{{ galaxy_server_url }}:5671/pulsar/galaxu_au?ssl=1"
>      min_polling_interval: 0.5
>      amqp_publish_retry: True
>      amqp_publish_retry_max_retries: 5
>      amqp_publish_retry_interval_start: 10
>      amqp_publish_retry_interval_step: 10
>      amqp_publish_retry_interval_max: 60
>
>    ```
>    {% endraw %}
>
> 3. Add the following lines to your `hosts` file:
>
>    ```ini
>    [pulsarservers]
>    <ip_address of your pulsar server>
>    ```
>
{: .hands_on}

We will now write a new playbook for the pulsar installation as we are going to install it on a separate VM. We will also install the CVMFS client and the Galaxy CVMFS repos on this machine so Pulsar has the same access to reference data that Galaxy does.

We need to include a couple of pre-tasks to install virtualenv, git, etc.

> ### {% icon hands_on %} Hands-on: Creating the playbook
>
> 1. Create a `pulsar-playbook.yml` file with the following contents:
>
>    {% raw %}
>    ```yaml
>    - hosts: pulsarservers
>      pre_tasks:
>        - name: Install some packages
>          package:
>            name:
>              - build-essential
>              - git
>              - python3-dev
>              - libcurl4-openssl-dev
>              - libssl-dev
>              - virtualenv
>            state: present
>            update_cache: yes
>          become: yes
>        - name: chown the /mnt dir to ubuntu
>          file:
>            path: /mnt
>            owner: ubuntu
>            group: ubuntu
>            mode: 0755
>          become: yes
>      roles:
>        - role: galaxyproject.cvmfs
>          become: yes
>        - role: galaxyproject.nginx
>          become: yes
>        - galaxyproject.pulsar
>    ```
>    {% endraw %}
>
>    There are a couple of *pre-tasks* here. This is because we need to install some base packages on these very vanilla ubuntu instances as well as give ourselves ownership of the directory we are installing into.
>
>
{: .hands_on}

We also need to create the dependency resolver file so pulsar knows how to find and install dependencies for the tools we ask it to run. The simplest method which covers 99% of the use cases is to use conda auto installs similar to how Galaxy works. We need to create the file and put it where the `galaxyproject.pulsar` role can find it.


> ### {% icon hands_on %} Hands-on: Creating dependency resolver configuration
>
> 1. Create a `templates` directory in your working directory.
>
>    ```bash
>    mkdir templates
>    ```
>
> 2. Create a `dependency_resolvers_conf.xml.j2` file inside the `templates` directory with the following contents:
>
>    ```xml
>    <dependency_resolvers>
>        <conda auto_install="True" auto_init="True"/>
>    </dependency_resolvers>
>    ```
>
>    This tells pulsar to **only** look for dependencies in conda.
>
>
{: .hands_on}


> ### {% icon details %} Running non-conda tools
> If the tool you want to run on Pulsar doesn't have a conda package, you will need to make alternative arrangements! This is complex and beyond our scope here. See the [Pulsar documentation](https://pulsar.readthedocs.io/en/latest/) for details.
{: .details}


> ### {% icon hands_on %} Hands-on: Run the Playbook
>
> 1. Run the playbook. If your remote pulsar machine uses a different key, you may need to supply the `ansible-playbook` command with the private key for the connection using the `--private-key key.pem` option.
>
>    ```bash
>    ansible-playbook pulsar-playbook.yml
>    ```
>
>    After the script has run, pulsar will be installed on the remote machines!
>
> 2. Log in to the machines and have a look in the `/mnt/pulsar` directory. You will see the venv and config directories. All the config files created by Ansible can be perused.
>
> 3. Run `journalctl -f -u pulsar`
>
>    A log will now start scrolling, showing the startup of pulsar. You'll notice that it will be initializing and installing conda. Once this is completed, Pulsar will be listening on the assigned port.
>
{: .hands_on}

# Configuring Galaxy to use Pulsar as a job destination

Now we have a Pulsar server up and running, we need to tell our Galaxy about it.

Galaxy talks to the Pulsar server via it's `job_conf.xml` file. We need to let Galaxy know about Pulsar there and make sure Galaxy has loaded the requisite job runner, and has a destination set up.

There are three things we need to do here:

* We will need to create a job runner which uses the  `galaxy.jobs.runners.pulsar:PulsarMQJobRunner` code.
* Create job destination which references the above job runner.
* Tell Galaxy which tools to send to the job destination: We will use `bwa-mem`

> ### {% icon hands_on %} Hands-on: Configure Galaxy
>
> 1. In your `templates/galaxy/config/job_conf.xml.j2` file add the following job runner to the `<plugins>` section:
>
>    ```xml
>    <plugin id="pulsar_runner" type="runner" load="galaxy.jobs.runners.pulsar:PulsarMQJobRunner" >
>        <param id="amqp_url">pyamqp://galaxy_au:{{ rabbitmq_password_galaxy_au }}@localhost:5671/pulsar/galaxy_au?ssl=1</param>
>        <param id="amqp_ack_republish_time">1200</param>
>        <param id="amqp_acknowledge">True</param>
>        <param id="amqp_consumer_timeout">2.0</param>
>        <param id="amqp_publish_retry">True</param>
>        <param id="amqp_publish_retry_max_retries">60</param>
>        <param id="galaxy_url">https://your_galaxy_ip_address_or_fqdn_here</param>
>        <param id="manager">_default_</param>
>    </plugin>
>    ```
>
> **Make sure you replace *your_galaxy_ip_address_or_fqdn_here* with your Galaxy servers IP adress or FQDN**
>
>    Add the following to the `<destinations>` section of your `job_conf.xml` file:
>
>    {% raw %}
>    ```xml
>    <destination id="pulsar" runner="pulsar_runner" >
>        <param id="default_file_action">remote_transfer</param>
>        <param id="dependency_resolution">remote</param>
>        <param id="jobs_directory">/mnt/pulsar/files/staging</param>
>        <param id="persistence_directory">/mnt/pulsar/files/persisted_data</param>
>        <param id="remote_metadata">False</param>
>        <param id="rewrite_parameters">True</param>
>        <param id="transport">curl</param>
>    </destination>
>    ```
>    {% endraw %}
>
> 2. Finally we need to tell Galaxy which tools to send to Pulsar. We will tell it to send bwa-mem jobs to it. We use the `<tools>` section of the `job_conf.xml` file.
>    We need to know the full id of the tool in question, we can get this out of the `integrated_tool_panel.xml` file in the `mutable-config` directory. Then we tell Galaxy which destination to send it to (pulsar).
>
>    Add the following to the end of the `job_conf.xml` file (inside the `<tools>` section if it exists or create it if it doesn't.)
>
>    ```xml
>    <tools>
>        <tool id="bwa" destination="pulsar"/>
>        <tool id="bwa_mem" destination="pulsar"/>
>    </tools>
>    ```
>
>    You can use the full tool ID here (toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa/0.7.17.4), or the short version. By using the full version, we restrict to only running that specific version in pulsar.
>
> 3. Run the Galaxy playbook in order to deploy the updated job configuration, and to restart Galaxy.
>
{: .hands_on}


# Testing Pulsar

Now we will upload a small set of data to run bwa-mem with.


> ### {% icon hands_on %} Hands-on: Testing the Pulsar destination
>
> 1. Upload the following files from zenodo.
>
>    ```
>    https://zenodo.org/record/582600/files/mutant_R1.fastq
>    https://zenodo.org/record/582600/files/mutant_R2.fastq
>    ```
>
> 2. **Map with BWA-MEM** {% icon tool %} with the following parameters
>    - *"Will you select a reference genome from your history or use a built-in index"*: `Use a built-in genome index`
>    - *"Using reference genome"*: `Escherichia coli (str. K-12 substr MG1655): eschColi_K12`
>    - *"Single or Paired-end reads"*: `Paired end`
>    - {% icon param-file %} *"Select first set of reads"*: `mutant_R1.fastq`
>    - {% icon param-file %} *"Select second set of reads"*: `mutant_R2.fastq`
>
>    As soon as you press *execute* Galaxy will send the job to the pulsar server. You can watch the log in Galaxy using:
>
>    ```
>    journalctl -fu galaxy
>    ```
>
>    You can watch the log in Pulsar by ssh'ing to it and tailing the log file with:
>
>    ```
>    journalctl -fu pulsar
>    ```
>
{: .hands_on}

You'll notice that the Pulsar server has received the job (all the way in Australia!) and now should be installing bwa-mem via conda. Once this is complete (which may take a while - first time only) the job will run. When it starts running it will realise it needs the *E. coli* genome from CVMFS and fetch that, and then results will be returned to Galaxy!

How awesome is that? Pulsar in another continent with reference data automatically from CVMFS :)

# Pulsar in Production

If you want to make use of Pulsar on a Supercomputer, you only need access to a submit node, and you will need to run Pulsar there. We recommend that if you need to run a setup with Pulsar, that you deploy an AMQP server (e.g. RabbitMQ) alongside your Galaxy. That way, you can run Pulsar on any submit nodes, and it can connect directly to the AMQP and Galaxy. Other Pulsar deployment options require exposing ports wherever Pulsar is running, and this requires significant more coordination effort.

You can also create multiple queues on your RabbitMQ server for multiple Pulsar servers. On Galaxy Australia, we run 5 different Pulsar servers spread out all around the country. They all communicate with Galaxy via the one RabbitMQ server. Galaxy Main (USA) and Galaxy Europe have similar setups with Galaxy Europe having over 7 Pulsars each configured for different roles!

For each new Pulsar server, you will need to add:
  1. In the RabbitMQ config:
      * A vhost
      * A user - configured with a password and the new vhost
  2. In the Galaxy job_conf.xml:
      * A new job runner with the new connection string
      * A new destination or multiple destinations for the new runner.

Pulsar servers can be the head node of a cluster. You can create a cluster and use your favourite job scheduler such as Slurm or PBS to schedule jobs. You can have many destinations in your Galaxy job_conf.xml file that change the number of cpus, amount of RAM etc. It can get quite complex and flexible if you like.

