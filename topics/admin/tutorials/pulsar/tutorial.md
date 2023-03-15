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
  - gmauro
subtopic: jobs
tags:
  - ansible
  - jobs
  - git-gat
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - connect-to-compute-cluster
      - job-destinations
      - cvmfs
  - title: "A server/VM on which to deploy Pulsar"
    type: "none"
---


Pulsar is the Galaxy Project's remote job running system. It was written by John Chilton ([@jmchilton](https://github.com/jmchilton)) of the Galaxy Project. It is a python server application that can accept jobs from a Galaxy server, submit them to a local resource and then send the results back to the originating Galaxy server.

More details on Pulsar can be found at:

- [Pulsar's Documentation](https://pulsar.readthedocs.io/en/latest/index.html)
- [Pulsar's Github Repository](https://github.com/galaxyproject/pulsar)
- [Pulsar Ansible Role](https://github.com/galaxyproject/ansible-pulsar)

Transport of data, tool information and other metadata can be configured as a web application via a RESTful interface or using a message passing system such as RabbitMQ.

At the Galaxy end, it is configured within the `job_conf.yml` file and uses one of two special Galaxy job runners.
* `galaxy.jobs.runners.pulsar:PulsarRESTJobRunner` for the RESTful interface
* `galaxy.jobs.runners.pulsar:PulsarMQJobRunner` for the message passing interface.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="pulsar" %}

**This tutorial assumes that:**

- You have a VM or machine where you will install Pulsar, and a directory in which the installation will be done. This tutorial assumes it is `/mnt`
- You have completed the "Galaxy Installation with Ansible", "Connecting Galaxy to a Compute Cluster", and the "CVMFS" tutorials
- You have access to the VM/computer where it is installed.

> <tip-title>This is NOT intended as a standalone Pulsar guide</tip-title>
> This tutorial is not intended to be a standalone Pulsar setup guide. If you read carefully and understand Ansible, it is likely you can figure out which portions are required to just setup Pulsar.
{: .tip}

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

> <tip-title>Other file transport methods for Pulsar</tip-title>
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


> <details-title>Why are we using Pulsar in MQ mode here and not the RESTful interface?</details-title>
> We are teaching you to install Pulsar and configure it in MQ mode in this tutorial. Configuring Pulsar in RESTful mode is also possible and is quite useful in certain situations. However, in the most common situation MQ mode is preferable for a number of reasons:
> * When running Pulsar in RESTful mode, all of the job control and data transfer is controlled by the Galaxy server usually using http transfers. This can place a limit on the size of files that can be transferred without constant configuring of the webserver.
> * When running in RESTful mode, Pulsar also needs to have an https server such as nginx, including securing it, configuring it, getting certificates and opening ports. This can be very difficult to do if you are attempting to submit jobs to an institutional HPC where the admins probably won't let you do any of these things.
> * In MQ mode, you only need to open a port for the RabbitMQ server on a machine you are more likely to control. The HPC side running Pulsar can just connect back to you.
>
> See the [Pulsar documentation](https://pulsar.readthedocs.io/en/latest/) for details.
{: .details}

# Install and configure a message queueing system

In this section we will install the RabbitMQ server on your Galaxy server VM.

RabbitMQ is an AMQP server that can queue messages between systems for all sorts of reasons. Here, we will be using the queue so that Galaxy and Pulsar can communicate jobs, job status and job metadata between them easily and robustly. More information on RabbitMQ can be found [on their website](https://www.rabbitmq.com/).

## Installing the roles

Firstly we will add and configure another *role* to our Galaxy playbook - we maintain a slightly modified version of `jasonroyle.rabbitmq` to support python3 and other minor updates. Additionally we will use the Galaxy community role for deploying Pulsar

> <hands-on-title>Install the Ansible roles</hands-on-title>
>
> 1. From your ansible working directory, edit the `requirements.yml` file and add the following lines:
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -24,3 +24,7 @@
>       version: 0.0.2
>     - src: galaxyproject.slurm
>       version: 0.1.3
>    +- name: usegalaxy_eu.rabbitmq
>    +  version: 0.1.0
>    +- src: galaxyproject.pulsar
>    +  version: 1.0.8
>    {% endraw %}
>    ```
>    {: data-commit="Add requirements"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Now install it with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
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

Optional: You can add tags to each user if required. e.g. For an admin user it could be useful to add in a *administrator* tag. These tags allow you to grant permissions to every user with a specific tag.

### RabbitMQ server config

We also need to set some RabbitMQ server configuration variables. Such as where its security certificates are and which ports to listen on (both via localhost and network).

> <tip-title>Port accessibility is important!</tip-title>
> We will need to make sure that the RabbitMQ default port is open and accessible on the server we are installing RabbitMQ onto. (In our case this is the Galaxy server). Default port number is: `5671`
{: .tip}

More information about the rabbitmq ansible role can be found [in the repository](https://github.com/usegalaxy-eu/ansible-role-rabbitmq).

## Add RabbitMQ configuration to Galaxy VM.

> <hands-on-title>Add RabbitMQ settings to Galaxy VM groupvars file.</hands-on-title>
>
> 1. Edit your `group_vars/secret.yml` and define some random passwords:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault edit group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    ```yaml
>    vault_rabbitmq_password_vhost: "a-really-long-password-here"
>    vault_rabbitmq_admin_password: "a-different-really-long-password"
>    ```
>
>    <!-- Ignore this, just for the gat-automation. Vaults are ugly to work with :(
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/secret.yml
>    +++ b/group_vars/secret.yml
>    @@ -1,7 +1,13 @@
>     $ANSIBLE_VAULT;1.1;AES256
>    -32653961383866636531396135663630386630346237333333653633313436663439643535323964
>    -6363626330336430363332643638646262316338313937320a666566306539373462386266383166
>    -30326165393863633463353234613561393939326164376432633732316264636464313061383161
>    -3532373937656138320a616361343664353264613332616236623231326137316635323465623562
>    -66656539346130353639623736633034653932373438663330646436656336666637313933666264
>    -3636313438626533633831323239373461373538646635613637
>    +62346261323266656232393034396134316636376533376139666437363535393562663838613938
>    +6336666266633563346337623265353935646361326337610a393834333233313461346439376438
>    +63383338346530656561636631666134373238366364363164313166346461383736613162653237
>    +3461363334323431370a656132303965653262386130353332623937376261396530393761353834
>    +38336565666437666436643163363831633331333766653266356163613138393734656465323634
>    +39366362383433366437353534663134313330316337393335383962613961386665633261616237
>    +35366635373063313631323939396164336330356361393464326636353037336461323531336434
>    +35613933303333623031353936393265636130363335376533393335663266313863376135383338
>    +36613464373231623938373434306266373234633036343636633963353361356631363533353066
>    +39323064336237646432323530313065303331326636353334343862373330313133326363363063
>    +38383564636161396435666164643334656435393533643163393434623434656238633631633939
>    +33353232666432376661
>    {% endraw %}
>    ```
>    {: data-commit="Add rabbitmq passwords to the vault"}
>
>    -->
>
>    This is going in the vault as they are secrets we need to set. Both of our services, Galaxy and Pulsar, need these variables, so we'll need to make sure they're in both playbooks. Both Galaxy in the job configuration, and Pulsar in its configuration.
>
>    Replace both with long random (or not) string.
>
> 2. From your ansible working directory, edit the `group_vars/galaxyservers.yml` file and add make the following changes.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -123,8 +123,10 @@ certbot_environment: staging
>     certbot_well_known_root: /srv/nginx/_well-known_root
>     certbot_share_key_users:
>       - nginx
>    +  - rabbitmq
>     certbot_post_renewal: |
>         systemctl restart nginx || true
>    +    systemctl restart rabbitmq-server || true
>     certbot_domains:
>      - "{{ inventory_hostname }}"
>     certbot_agree_tos: --agree-tos
>    @@ -180,6 +182,34 @@ slurm_config:
>       SelectType: select/cons_res
>       SelectTypeParameters: CR_CPU_Memory  # Allocate individual cores/memory instead of entire node
>     
>    +# RabbitMQ
>    +rabbitmq_version: 3.8.35-1
>    +rabbitmq_plugins: rabbitmq_management
>    +
>    +rabbitmq_config:
>    +- rabbit:
>    +  - tcp_listeners:
>    +    - "'127.0.0.1'": 5672
>    +  - ssl_listeners:
>    +    - "'0.0.0.0'": 5671
>    +  - ssl_options:
>    +     - cacertfile: /etc/ssl/certs/fullchain.pem
>    +     - certfile: /etc/ssl/certs/cert.pem
>    +     - keyfile: /etc/ssl/user/privkey-rabbitmq.pem
>    +     - fail_if_no_peer_cert: 'false'
>    +
>    +rabbitmq_vhosts:
>    +  - /pulsar/galaxy_au
>    +
>    +rabbitmq_users:
>    +  - user: admin
>    +    password: "{{ vault_rabbitmq_admin_password }}"
>    +    tags: administrator
>    +    vhost: /
>    +  - user: galaxy_au
>    +    password: "{{ vault_rabbitmq_password_vhost }}"
>    +    vhost: /pulsar/galaxy_au
>    +
>     # TUS
>     galaxy_tusd_port: 1080
>     tusd_instances:
>    {% endraw %}
>    ```
>    {: data-commit="Configure RabbitMQ"}
>
>    > <tip-title>RabbitMQ installation errors?</tip-title>
>    > RabbitMQ depends on specific Erlang versions. If the Erlang version has been updated, you may need to change the value of `rabbitmq_version:` in the configuration above. [Information on the RabbitMQ Erlag version requirements.](https://www.rabbitmq.com/which-erlang.html)
>    {: .tip}
>
> 3. Update the Galaxy playbook to include the *usegalaxy_eu.rabbitmq* role.
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -29,6 +29,7 @@
>         - role: uchida.miniconda
>           become: true
>           become_user: "{{ galaxy_user.name }}"
>    +    - usegalaxy_eu.rabbitmq
>         - galaxyproject.nginx
>         - galaxyproject.tusd
>         - galaxyproject.cvmfs
>    {% endraw %}
>    ```
>    {: data-commit="Add role"}
>
>    > <tip-title>Why is this at the end?</tip-title>
>    > This is one of the constant problems with Ansible, how do you order everything correctly? Does an ordering exist such that a single run of the playbook will have everything up and working? We encounter one such instance of this problem now.
>    >
>    > Here are the dependencies between the roles:
>    >
>    > From             | To               | Purpose
>    > ---------------- | --               | -------
>    > nginx            | galaxy           | The nginx templates depend on variables only available after the Galaxy role is run
>    > SSL certificates | nginx            | A running nginx is required
>    > RabbitMQ         | SSL certificates | RabbitMQ will silently start with incorrect configuration if SSL certificates are not present at boot time.
>    > Galaxy           | RabbitMQ         | Galaxy needs the RabbitMQ available to submit jobs.
>    >
>    > And as you can see there is a circular dependency. Galaxy requires RabbitMQ, but RabbitMQ depends on a long chain of things that depends finally on Galaxy.
>    >
>    > There are some mitigating factors, some software will start with incomplete configuration. We can rely on Galaxy retrying access to RabbitMQ if it isn't already present. Additionally on first run, Galaxy is restarted by a handler which runs at the end. (Except that the nginx role triggers all pending handlers as part of the SSL certificate deployment.)
>    >
>    > We try to present the optimal version here but due to these interdependencies and Ansible specifics, sometimes it is not possible to determine a good ordering of roles, and multiple runs might be required.
>    >
>    {: .tip}
>
> 4. Run the playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> The rabbitmq server daemon will have been installed on your Galaxy VM. Check that it's running now:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > systemctl status rabbitmq-server
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    >
>    > ```ini
>    > ● rabbitmq-server.service - RabbitMQ broker
>    >      Loaded: loaded (/lib/systemd/system/rabbitmq-server.service; enabled; vendor preset: enabled)
>    >      Active: active (running) since Fri 2020-12-18 13:52:14 UTC; 8min ago
>    >    Main PID: 533733 (beam.smp)
>    >      Status: "Initialized"
>    >       Tasks: 163 (limit: 19175)
>    >      Memory: 105.3M
>    >      CGroup: /system.slice/rabbitmq-server.service
>    >              ├─533733 /usr/lib/erlang/erts-11.1.4/bin/beam.smp -W w -K true -A 128 -MBas ageffcbf -MHas ageffcbf -MBlmbcs 512 -MHlmbcs 512 -MMmcs 30 -P 1048576 -t 5000000 -stbt db -zdbbl 1280>
>    >              ├─533923 erl_child_setup 32768
>    >              ├─533969 /usr/lib/erlang/erts-11.1.4/bin/epmd -daemon
>    >              ├─534002 inet_gethost 4
>    >              └─534003 inet_gethost 4
>    >
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   ##########  Licensed under the MPL 2.0. Website: https://rabbitmq.com
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   Doc guides: https://rabbitmq.com/documentation.html
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   Support:    https://rabbitmq.com/contact.html
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   Tutorials:  https://rabbitmq.com/getstarted.html
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   Monitoring: https://rabbitmq.com/monitoring.html
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   Logs: /var/log/rabbitmq/rabbit@gat-0.log
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:         /var/log/rabbitmq/rabbit@gat-0_upgrade.log
>    > Dec 18 13:52:10 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   Config file(s): (none)
>    > Dec 18 13:52:14 gat-0.training.galaxyproject.eu rabbitmq-server[533733]:   Starting broker... completed with 0 plugins.
>    > Dec 18 13:52:14 gat-0.training.galaxyproject.eu systemd[1]: Started RabbitMQ broker.
>    > ```
>    {: .code-out.code-max-300}
>
>    But this doesn't tell the whole story, so run the diagnostics command to
>    check that the interfaces are setup and listening. RabbitMQ has a bad
>    habit of silently failing when processing the configuration, without any
>    logging information If RabbitMQ has any problem reading the configuration
>    file, it falls back to the default configuration (listens *without* ssl on
>    `tcp/5672`) so be sure to check that everything is OK before continuing.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > sudo rabbitmq-diagnostics status
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    >
>    > ```ini
>    > ...
>    >
>    > Listeners
>    >
>    > Interface: [::], port: 15672, protocol: http, purpose: HTTP API
>    > Interface: [::], port: 25672, protocol: clustering, purpose: inter-node and CLI tool communication
>    > Interface: [::], port: 5672, protocol: amqp, purpose: AMQP 0-9-1 and AMQP 1.0
>    > Interface: 0.0.0.0, port: 5671, protocol: amqp/ssl, purpose: AMQP 0-9-1 and AMQP 1.0 over TLS
>    > ```
>    {: .code-out.code-max-300}
>
>    But wait! There are more ways it can go wrong. To be extra sure, run a quick `curl` command.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > curl http://localhost:5672
>    > curl -k https://localhost:5671
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    >
>    > These should *both* report the same response:
>    >
>    > ```console
>    > curl: (1) Received HTTP/0.9 when not allowed
>    > ```
>    >
>    > if they don't, consider the following debugging steps:
>    >
>    > 1. Restarting RabbitMQ
>    > 2. Check that the configuration looks correct (ssl private key path looks valid)
>    > 3. Check that the private key is shared correctly with the rabbitmq user
>    {: .code-out.code-max-300}
>
{: .hands_on}


# Installing and configuring Pulsar on a remote machine

Now that we have a message queueing system running on our Galaxy VM, we need to install and configure Pulsar on our remote compute VM. To do this we need to create a new ansible playbook to install Pulsar.

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

> <hands-on-title>Configure pulsar group variables</hands-on-title>
>
>
> 2. Create a new file in `group_vars` called `pulsarservers.yml` and set some of the above variables as well as some others.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/pulsarservers.yml
>    @@ -0,0 +1,51 @@
>    +galaxy_server_hostname: "" # Important!!!
>    +# Put your Galaxy server's fully qualified domain name (FQDN) (or the FQDN of the RabbitMQ server) above.
>    +
>    +pulsar_root: /mnt/pulsar
>    +
>    +pulsar_pip_install: true
>    +pulsar_pycurl_ssl_library: openssl
>    +pulsar_systemd: true
>    +pulsar_systemd_runner: webless
>    +
>    +pulsar_create_user: true
>    +pulsar_user: {name: pulsar, shell: /bin/bash}
>    +
>    +pulsar_optional_dependencies:
>    +  - pyOpenSSL
>    +  # For remote transfers initiated on the Pulsar end rather than the Galaxy end
>    +  - pycurl
>    +  # drmaa required if connecting to an external DRM using it.
>    +  - drmaa
>    +  # kombu needed if using a message queue
>    +  - kombu
>    +  # amqp 5.0.3 changes behaviour in an unexpected way, pin for now.
>    +  - 'amqp==5.0.2'
>    +  # psutil and pylockfile are optional dependencies but can make Pulsar
>    +  # more robust in small ways.
>    +  - psutil
>    +
>    +pulsar_yaml_config:
>    +  staging_directory: "{{ pulsar_staging_dir }}"
>    +  persistence_directory: "{{ pulsar_persistence_dir }}"
>    +  tool_dependency_dir: "{{ pulsar_dependencies_dir }}"
>    +  # The following are the settings for the pulsar server to contact the message queue with related timeouts etc.
>    +  message_queue_url: "pyamqp://galaxy_au:{{ vault_rabbitmq_password_vhost }}@{{ galaxy_server_hostname }}:5671//pulsar/galaxy_au?ssl=1"
>    +  min_polling_interval: 0.5
>    +  amqp_publish_retry: True
>    +  amqp_publish_retry_max_retries: 5
>    +  amqp_publish_retry_interval_start: 10
>    +  amqp_publish_retry_interval_step: 10
>    +  amqp_publish_retry_interval_max: 60
>    +  # We also need to create the dependency resolvers configuration so pulsar knows how to find and install dependencies
>    +  # for the tools we ask it to run. The simplest method which covers 99% of the use cases is to use conda auto installs
>    +  # similar to how Galaxy works.
>    +  dependency_resolution:
>    +    resolvers:
>    +      - type: conda
>    +        auto_init: true
>    +        auto_install: true
>    +
>    +# Pulsar should use the same job metrics plugins as Galaxy. This will automatically set `job_metrics_config_file` in
>    +# `pulsar_yaml_config` and create `{{ pulsar_config_dir }}/job_metrics_conf.yml`.
>    +pulsar_job_metrics_plugins: "{{ galaxy_job_metrics_plugins }}"
>    {% endraw %}
>    ```
>    {: data-commit="Add pulsar group variables"}
>
>    > <details-title>Running non-conda tools</details-title>
>    > If the tool you want to run on Pulsar doesn't have a conda package, you will need to make alternative arrangements! This is complex and beyond our scope here. See the [Pulsar documentation](https://pulsar.readthedocs.io/en/latest/) for details.
>    {: .details}
>
> 3. Add the following lines to your `hosts` file:
>
>    {% raw %}
>    ```diff
>    --- a/hosts
>    +++ b/hosts
>    @@ -1,2 +1,4 @@
>     [galaxyservers]
>     gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +[pulsarservers]
>    +gat-0.au.training.galaxyproject.eu ansible_user=ubuntu
>    {% endraw %}
>    ```
>    {: data-commit="Add pulsar host"}
>
{: .hands_on}

We will now write a new playbook for the pulsar installation as we are going to install it on a separate VM. We will also install the CVMFS client and the Galaxy CVMFS repos on this machine so Pulsar has the same access to reference data that Galaxy does.

We need to include a couple of pre-tasks to install virtualenv, git, etc.

> <hands-on-title>Creating the playbook</hands-on-title>
>
> 1. Create a `pulsar.yml` file with the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/pulsar.yml
>    @@ -0,0 +1,20 @@
>    +- hosts: pulsarservers
>    +  vars_files:
>    +    - group_vars/secret.yml
>    +  pre_tasks:
>    +    - name: Install some packages
>    +      package:
>    +        name:
>    +          - build-essential
>    +          - git
>    +          - python3-dev
>    +          - libcurl4-openssl-dev
>    +          - libssl-dev
>    +          - virtualenv
>    +        state: present
>    +        update_cache: yes
>    +      become: yes
>    +  roles:
>    +    - role: galaxyproject.cvmfs
>    +      become: yes
>    +    - galaxyproject.pulsar
>    {% endraw %}
>    ```
>    {: data-commit="Add pulsar playbook"}
>
>    There are a couple of *pre-tasks* here. This is because we need to install some base packages on these very vanilla ubuntu instances as well as give ourselves ownership of the directory we are installing into.
>
{: .hands_on}

> <hands-on-title>Run the Playbook</hands-on-title>
>
> 1. Run the playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook pulsar.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    After the script has run, pulsar will be installed on the remote machines!
>
>    > <tip-title>Connection issues?</tip-title>
>    > If your remote pulsar machine uses a different key, you may need to supply the `ansible-playbook` command with the private key for the connection using the `--private-key key.pem` option.
>    {: .tip}
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

Galaxy talks to the Pulsar server via it's `job_conf.yml` file. We need to let Galaxy know about Pulsar there and make sure Galaxy has loaded the requisite job runner, and has a destination set up.

There are three things we need to do here:

* Create a job runner which uses the  `galaxy.jobs.runners.pulsar:PulsarMQJobRunner` code.
* Create a job destination referencing the above job runner.
* Tell Galaxy which tools to send to this job destination.

For this tutorial, we will configure Galaxy to run the BWA and BWA-MEM tools on Pulsar.

> <hands-on-title>Configure Galaxy</hands-on-title>
>
> 1. In your `templates/galaxy/config/job_conf.yml.j2` file add the following job runner to the `<plugins>` section:
>
>    {% raw %}
>    ```diff
>    --- a/templates/galaxy/config/job_conf.yml.j2
>    +++ b/templates/galaxy/config/job_conf.yml.j2
>    @@ -4,6 +4,16 @@ runners:
>         workers: 4
>       slurm:
>         load: galaxy.jobs.runners.slurm:SlurmJobRunner
>    +  pulsar_runner:
>    +    load: galaxy.jobs.runners.pulsar:PulsarMQJobRunner
>    +    amqp_url: "pyamqp://galaxy_au:{{ vault_rabbitmq_password_vhost }}@localhost:5671/{{ rabbitmq_vhosts[0] }}?ssl=1"
>    +    amqp_acknowledge: true
>    +    amqp_ack_republish_time: 1200
>    +    amqp_consumer_timeout: 2
>    +    amqp_publish_retry: true
>    +    amqp_publish_retry_max_retries: 60
>    +    galaxy_url: "https://{{ inventory_hostname }}"
>    +    manager: _default_
>     
>     execution:
>       default: slurm
>    {% endraw %}
>    ```
>    {: data-commit="Add pulsar plugin"}
>
>    Add the following to the `<destinations>` section of your `job_conf.yml` file:
>
>    {% raw %}
>    ```diff
>    --- a/templates/galaxy/config/job_conf.yml.j2
>    +++ b/templates/galaxy/config/job_conf.yml.j2
>    @@ -20,6 +20,16 @@ execution:
>       environments:
>         local_dest:
>           runner: local_runner
>    +    pulsar:
>    +      runner: pulsar_runner
>    +      default_file_action: remote_transfer
>    +      dependency_resolution: remote
>    +      jobs_directory: /mnt/pulsar/files/staging
>    +      persistence_directory: /mnt/pulsar/files/persisted_data
>    +      remote_metadata: false
>    +      rewrite_parameters: true
>    +      transport: curl
>    +      outputs_to_working_directory: false
>         slurm:
>           runner: slurm
>           singularity_enabled: true
>    {% endraw %}
>    ```
>    {: data-commit="Add pulsar destination"}
>
>    You'll notice we need to know a lot about the configuration of the remote end, this is an unfortunate requirement with pulsar. Changes to e.g. the staging directory need to be coordinated between Pulsar and Galaxy. That's fine if both are under your administration, but for a completely remote Pulsar it can be difficult.
>
>    Notably we also override `outputs_to_working_directory`, as this option is incompatible with running Pulsar, and, unnecessary. Pulsar already provides the same job isolation and safety that we request when we set that option by default in Galaxy's configuration.
>
> 2. Install the BWA and BWA-MEM tools, if needed.
>
>    {% snippet topics/admin/faqs/install_tool.md query="bwa" name="Map with BWA-MEM" section="Mapping" %}
>
> 3. We now need to tell Galaxy to send BWA and BWA-MEM jobs to the `pulsar` destination. We specify this in the `<tools>` section of the `job_conf.yml` file.
>
>    Add the following to the end of the `job_conf.yml` file (inside the `<tools>` section if it exists or create it if it doesn't.)
>
>    {% raw %}
>    ```diff
>    --- a/templates/galaxy/config/job_conf.yml.j2
>    +++ b/templates/galaxy/config/job_conf.yml.j2
>    @@ -86,3 +86,7 @@ tools:
>     - id: testing
>       environment: dynamic_cores_time
>       resources: testing
>    +- id: bwa
>    +  environment: pulsar
>    +- id: bwa_mem
>    +  environment: pulsar
>    {% endraw %}
>    ```
>    {: data-commit="Send bwa and bwa-mem to pulsar"}
>
>    Note that here we are using the short tool IDs. If you want to run only a specific version of a tool in Pulsar, you have to use the full tool ID (e.g. `toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa/0.7.17.4`) instead. The full tool ID can be found inside the `integrated_tool_panel.xml` file in the `mutable-config` directory.
>
> 4. Finally run the Galaxy playbook in order to deploy the updated job configuration, and to restart Galaxy.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on}

> ```bash
> 1-pulsar.sh
> ```
> {: data-test="true"}
{: .hidden}


# Testing Pulsar

Now we will upload a small set of data to run bwa-mem with.

> <hands-on-title>Testing the Pulsar destination</hands-on-title>
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

> ```bash
> 2-run-job.sh
> ```
> {: data-test="true"}
{: .hidden}

You'll notice that the Pulsar server has received the job (all the way in Australia!) and now should be installing bwa-mem via conda. Once this is complete (which may take a while - first time only) the job will run. When it starts running it will realise it needs the *E. coli* genome from CVMFS and fetch that, and then results will be returned to Galaxy!

How awesome is that? Pulsar in another continent with reference data automatically from CVMFS :)

{% snippet topics/admin/faqs/missed-something.md step=9 %}

# Retries of the staging actions

When the staging actions are carried out by the Pulsar server itself (like in the case when driving Pulsar by message queue), there are some parameters that can be tweaked to ensure reliable communication between the Galaxy server and the remote Pulsar server.
The aim of these parameters is to control the retrying of staging actions in the event of a failure.

For each action (preprocess/input or postprocess/output), you can specify:
```text
 - *_action_max_retries    - the maximum number of retries before giving up
 - *_action_interval_start - how long start sleeping between retries (in seconds)
 - *_action_interval_step  - by how much the interval is increased for each retry (in seconds)
 - *_action_interval_max   - the maximum number of seconds to sleep between retries
```
substitute the * with `preprocess` or `postprocess`

In the following box, as an example, we have collected the values adopted in a Pulsar site with an unreliable network connection:

```yaml
preprocess_action_max_retries: 30
preprocess_action_interval_start: 2
preprocess_action_interval_step: 10
preprocess_action_interval_max: 300
postprocess_action_max_retries: 30
postprocess_action_interval_start: 2
postprocess_action_interval_step: 10
postprocess_action_interval_max: 300

```
In this case, for both actions, Pulsar will try to carry out the staging action 30 times, sleeping 2 secs after the first retry and adding 10 secs more to each next retries, until a maximum of 300 seconds between retries.

We hope you never have to experience a situation like this one, but if needed just adapt the numbers to your case and add the parameters in the `pulsar_yaml_config` section of your `pulsarservers.yml` file.

# Pulsar in Production

If you want to make use of Pulsar on a Supercomputer, you only need access to a submit node, and you will need to run Pulsar there. We recommend that if you need to run a setup with Pulsar, that you deploy an AMQP server (e.g. RabbitMQ) alongside your Galaxy. That way, you can run Pulsar on any submit nodes, and it can connect directly to the AMQP and Galaxy. Other Pulsar deployment options require exposing ports wherever Pulsar is running, and this requires significant more coordination effort.

For each new Pulsar server, you will need to add:
  1. In the RabbitMQ config:
      * A vhost
      * A user - configured with a password and the new vhost
  2. In the Galaxy job_conf.yml:
      * A new job runner with the new connection string
      * A new destination or multiple destinations for the new runner.

Pulsar servers can be the head node of a cluster. You can create a cluster and use your favourite job scheduler such as Slurm or PBS to schedule jobs. You can have many destinations in your Galaxy job_conf.yml file that change the number of cpus, amount of RAM etc. It can get quite complex and flexible if you like.

## Australia

You can also create multiple queues on your RabbitMQ server for multiple Pulsar servers. On Galaxy Australia, we run 5 different Pulsar servers spread out all around the country. They all communicate with Galaxy via the one RabbitMQ server.

![Map of australia with 6 pulsar nodes marked around the country.](../../images/pulsar_australia.png)

## Europe

Galaxy Europe has taken Pulsar and built [The Pulsar Network](https://pulsar-network.readthedocs.io/en/latest). This provides a framework for easily deploying Pulsar clusters in the cloud, something needed to support compute centers which might not have as much experience. This way they get an easy package they can deploy and the European Galaxy team can manage.

![Map of europe with pulsar nodes marked in many countries. An inset shows australia with a node there too.](https://pulsar-network.readthedocs.io/en/latest/_images/nodes.png)

The main purpose of this network is to support the workload of the UseGalaxy.eu instance by distributing it across several European data centers and clusters. If you're interested in setting up something similar, they [provide documentation](https://pulsar-network.readthedocs.io/en/latest/introduction.html) on how to install and configure a Pulsar network endpoint on a cloud infrastructure and how to connect it to your server.

# Conclusion

You're ready to ship your Galaxy jobs around the world! Now wherever you have compute space, you know how to setup a Pulsar node and connect it to Galaxy. Let us know if you come up with creative places to run your Galaxy jobs (coworker's laptops, your IoT fridge, the sky is the limit if it's x86 and has python)
