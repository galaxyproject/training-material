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

In this tutorial, we will:

* Install and configure a Pulsar server on a remote linux machine using ansible
    * We will configure the Pulsar server to run via the RESTful interface
* Configure our Galaxy servers to run a job there
* Run a job remotely

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


This tutorial assumes that you have:

- A VM or machine where you will install Pulsar, and a directory in which the installation will be done. This tutorial assumes it is `/mnt`

# Installing the Pulsar Role

We need to create a new ansible playbook to install Pulsar. We will be using a *role* developed by the Galaxy community - `galaxyproject.pulsar`

> ### {% icon hands_on %} Hands-on: Install the `galaxyproject.pulsar` ansible role
>
> 1. From your ansible working directory, edit the `requirements.yml` file and add the following line:
>
>    ```yaml
>    - src: galaxyproject.pulsar
>      version: 1.0.2
>    ```
>
> 2. Now install it with:
>
>    ```bash
>    ansible-galaxy install -p roles -r requirements.yml
>    ```
>
{: .hands_on}


# Configuring Pulsar

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

Additional options from Pulsar's server.ini are configurable via the following variables (these options are explained in the Pulsar documentation and server.ini.sample):

Variable               | Description                                                                                                                                                              | Default
----------             | -------------                                                                                                                                                            | ------
`pulsar_host`          | This is the interface pulsar will listen to, we will set it to `0.0.0.0`                                                                                                 | localhost
`pulsar_port`          |                                                                                                                                                                          | 8913
`pulsar_uwsgi_socket`  | If unset, uWSGI will be configured to listen for HTTP requests on pulsar_host port pulsar_port; If set, uWSGI will listen for uWSGI protocol connections on this socket. | `unset`
`pulsar_uwsgi_options` | Hash (dictionary) of additional uWSGI options to place in the [uwsgi] section of server.ini                                                                              | `{}`


Some of the other options we will be using are:

* We are going to run in RESTful mode so we will need to specify a `private_token` variable so we can "secure" the connection. (For a given value of "secure".)
* We will be using the uwsgi web server to host the RESTful interface.
* We will set the tool dependencies to rely on **conda** for tool installs.

> ### {% icon hands_on %} Hands-on: Configure pulsar group variables
>
> 1. Create a new file in `group_vars` called `pulsarservers.yml` and set some of the above variables as well as some others.
>
>    {% raw %}
>    ```yaml
>    pulsar_server_dir: /mnt/pulsar/server
>    pulsar_venv_dir: /mnt/pulsar/venv
>    pulsar_config_dir: /mnt/pulsar/config
>    pulsar_staging_dir: /mnt/pulsar/staging
>    pulsar_systemd: true
>
>    pulsar_host: 0.0.0.0
>    pulsar_port: 8913
>
>    private_token: your_private_token_here
>
>    pulsar_create_user: true
>    pulsar_user: {name: pulsar, shell: /bin/bash}
>
>    pulsar_optional_dependencies:
>      - pyOpenSSL
>      # For remote transfers initiated on the Pulsar end rather than the Galaxy end
>      - pycurl
>      # uwsgi used for more robust deployment than paste
>      - uwsgi
>      # drmaa required if connecting to an external DRM using it.
>      - drmaa
>      # kombu needed if using a message queue
>      - kombu
>      # requests and poster using Pulsar remote staging and pycurl is unavailable
>      - requests
>      # psutil and pylockfile are optional dependencies but can make Pulsar
>      # more robust in small ways.
>      - psutil
>
>    pulsar_yaml_config:
>      dependency_resolvers_config_file: dependency_resolvers_conf.xml
>      conda_auto_init: True
>      conda_auto_install: True
>      staging_directory: "{{ pulsar_staging_dir }}"
>      private_token: "{{ private_token }}"
>
>    # NGINX
>    nginx_selinux_allow_local_connections: true
>    nginx_servers:
>      - pulsar-proxy
>    nginx_enable_default_server: false
>    nginx_conf_http:
>      client_max_body_size: 5g
>    ```
>    {% endraw %}
>
> 2. Replace `your_private_token_here` with a long randomish (or not) string.
>
> 3. Add the following lines to your `hosts` file:
>
>    ```ini
>    [pulsarservers]
>    <ip_address of your pulsar server>
>    ```
>
> 4. Create the file `templates/nginx/pulsar-proxy.j2` with the following contents:
>
>    ```nginx
>    server {
>        # Listen on 80, you should secure your server better :)
>        listen 80 default_server;
>        listen [::]:80 default_server;
>
>        location / {
>            proxy_redirect off;
>            proxy_set_header Host $host;
>            proxy_set_header X-Real-IP $remote_addr;
>            proxy_pass http://localhost:8913;
>        }
>    }
>    ```
>
{: .hands_on}

We will now write a new playbook for the pulsar installation similar to the one we did for the CVMFS installation earlier in the week.

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
>        - role: galaxyproject.nginx
>          become: yes
>        - galaxyproject.pulsar
>    ```
>    {% endraw %}
>
>    There are a couple of *pre-tasks* here. This is because we need to install some base packages on these very vanilla ubuntu instances as well as give ourselves ownership of the directory we are installing into.
>
>    > ### {% icon comment %} Why NGINX?
>    Additionally we install NGINX, you might not have expected this! We used to use Pulsar's webserving directly via uWSGI, but in Python 3 Galaxy, the requests that are sent to Pulsar are chunked, a transfer encoding that is not part of the wsgi spec and unsupported. *Our recommendation*: avoid all of this weirdness and use RabbitMQ as the transport instead. Unfortunately that is currently outside of the scope of this tutorial. [The documentation](https://pulsar.readthedocs.io/en/latest/galaxy_with_rabbitmq_conf.html) covers it in detail.
>    {: .comment}
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

# Configuring Galaxy

Now we have a Pulsar server up and running, we need to tell our Galaxy about it.

Galaxy talks to the Pulsar server via it's `job_conf.xml` file. We need to let Galaxy know about Pulsar there and make sure Galaxy has loaded the requisite job runner, and has a destination set up.

There are three things we need to do here:

* We will need to create a job runner which uses the  `galaxy.jobs.runners.pulsar:PulsarRESTJobRunner` code.
* Create job destination which references the above job runner.
* Tell Galaxy which tools to send to the job destination: We will use `bwa-mem`

> ### {% icon hands_on %} Hands-on: Configure Galaxy
>
> 1. In your `files/galaxy/config/job_conf.xml` file add the following job runner to the `<plugins>` section:
>
>    ```xml
>    <plugin id="pulsar_runner" type="runner" load="galaxy.jobs.runners.pulsar:PulsarRESTJobRunner" />
>    ```
>
>    Then add the Pulsar destination. We will need the ip address of your pulsar server and the private_token string you used when you created it.
>
>    Add the following to the `<destinations>` section of your `job_conf.xml` file:
>
>    ```xml
>    <destination id="pulsar" runner="pulsar_runner" >
>        <param id="url">http://your_ip_address_here:80/</param>
>        <param id="private_token">your_private_token_here</param>
>    </destination>
>    ```
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
>    https://zenodo.org/record/582600/files/wildtype.fna
>    ```
>
> 2. **Map with BWA-MEM** {% icon tool %} with the following parameters
>    - *"Will you select a reference genome from your history or use a built-in index"*: `Use a genome from history and build index`
>    - {% icon param-file %} *"Use the following dataset as the reference genome"*: `wildtype.fna`
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

You'll notice that the Pulsar server has received the job (all the way in Sydney!) and now should be installing bwa-mem via conda. Once this is complete (which may take a while - first time only) the job will run and the results will be returned to Galaxy!

How awesome is that? :)
