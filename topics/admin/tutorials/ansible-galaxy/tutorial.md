---
layout: tutorial_hands_on

title: "Galaxy Installation with Ansible"
questions:
- How does the Galaxy Ansible module work internally?
- How can I install a Galaxy server with Ansible
objectives:
- Have an understanding of how Galaxy's Ansible roles are structured and interact with one another
- Be able to use an Ansible playbook to install different flavors of Galaxy for different purposes
time_estimation: "2h30m"
key_points:
- Basic deployment with Ansible is surprisingly easy
- Complexity can grow over time as your organisation does, no need to start with playbooks like UseGalaxy.org
contributions:
  authorship:
  - hexylena
  - natefoo
  - slugger70
  editing:
  - shiltemann
  - nsoranzo
  - gmauro
  - mvdbeek
  - martenson
  - jmchilton
  - davebx
  - nsoranzo
  - lecorguille
  - abretaud
  testing:
  - mira-miracoli
tags:
  - ansible
  - deploying
  - git-gat
subtopic: core
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
  - type: "none"
    title: "A VM with at least 2 vCPUs and 4 GB RAM, preferably running Ubuntu 18.04 - 20.04."
abbreviations:
    DVCS: Distributed Version Control System
---

> <warning-title>22.05 not yet released</warning-title>
> We updated the tutorials ahead of GCC, preparing for the GCC Galaxy Admin Training. However, 22.05 is not yet released and we still have some bugs in the roles to work our, so, until that time we recommend you follow the [GTN Training Archive](https://training.galaxyproject.org/archive/2022-05-01/topics/admin/tutorials/ansible-galaxy/tutorial.html) materials for this tutorial (and all others in this series.)
{: .warning}

This tutorial assumes you have some familiarity with [Ansible](https://www.ansible.com/resources/get-started) and are comfortable with writing and running playbooks. Here we'll see how to install a Galaxy server using an Ansible playbook. The Galaxy Project has decided on Ansible for all of its deployment recipes. For our project, Ansible is even more fitting due to its name:

> An ansible is a category of fictional device or technology capable of instantaneous or faster-than-light communication. It can send and receive messages to and from a corresponding device over any distance or obstacle whatsoever with no delay, even between star systems (Source: [Wikipedia](https://en.wikipedia.org/wiki/Ansible))
{: .quote}

We want to give you a comprehensive understanding of how the Galaxy installation occurs, but we want to avoid you having to write a "custom" Galaxy installation playbook which you would eventually throw away, in order to use the official playbooks. Given these goals, we will go through the playbook in depth first, and then move to a hands-on portion later. If you are not interested in the inner workings, you can [skip to that section now](#installing-galaxy).

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="ansible-galaxy" %}

{% snippet topics/admin/faqs/admin-testing.md %}

# Playbook Overview


## Configuration

We'll be using the [official Galaxy role](https://github.com/galaxyproject/ansible-galaxy) to install and manage Galaxy. This role is found in [Ansible Galaxy](https://galaxy.ansible.com/) (no relation - it is Ansible's system for sharing reusable Ansible roles) as [galaxyproject.galaxy](https://galaxy.ansible.com/galaxyproject/galaxy).

The official role is extremely configurable, everything that you want to change is exposed as a variable, and then tasks will change behaviour based on that. The [role documentation](https://github.com/galaxyproject/ansible-galaxy#role-variables) is the most up-to-date source of documentation for the variables. You should take a minute and read over the variables listed there.

The important variables for this tutorial are:

- `galaxy_root`
- `galaxy_commit_id`
- `galaxy_config`
- `galaxy_server_dir` (automatically set based on `galaxy_root`)

These are largely self explanatory: a directory for all of Galaxy's code and configuration, which commit should be installed, and the Galaxy configuration. We will not explain Galaxy configuration variables in detail as they are covered sufficiently in the `galaxy.yml` sample file or the [online documentation](https://docs.galaxyproject.org/en/master/admin/config.html#configuration-options).

The official recommendation is that you should have a variables file such as a `group_vars/galaxyservers.yml` for storing all of the Galaxy configuration.

## Tasks

As with every role, the entry point for execution is the `tasks/main.yml` file. [In the case of the galaxyproject.galaxy role](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/main.yml), this includes a few groups of important tasks:

- [Clone (or Download) Galaxy](#cloning-galaxy)
- [Managing Configuration](#managing-configuration)
- [Fetching Dependencies](#dependencies)
- [Managing Mutable Setup](#mutable-setup)
- [Managing the Database](#managing-the-database)
- [Building the Client](#building-the-client)

The actions each set of tasks performs are described below.

### Cloning Galaxy

The [clone](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/clone.yml) task is the one which is primarily interesting to us, it downloads Galaxy, using git, at a specific commit (or more generally, any [git reference](https://git-scm.com/book/en/v2/Git-Internals-Git-References)).

1. Ansible tries to update Galaxy, cloning it if it is missing, or otherwise attempting to update to the correct commit (or latest commit of the given branch).
2. Any change is reported.
3. The virtualenv is set up:
    1. An empty virtualenv is created.
    2. Pip is updated within the virtualenv.
4. Any `.pyc` files are removed, as this can occasionally result in Python loading the cached code, even if the corresponding `.py` file is no longer present at the checked-out commit. For safety, all of these are removed.

With that Galaxy is cloned to disk and is ready to be configured by the next task.

### Managing Configuration

The [static configuration setup](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/static_setup.yml) is relatively straightforward:

1. The directories for Galaxy configuration data and for the shed tools are created
2. Any config files are copied over
3. Any templates are copied over
4. The `galaxy.yml` (or `.ini`) is deployed

The setup for deploying extra Galaxy configuration templates and files is a little bit non-standard by Ansible best practices. Here you are expected to provide your own templates and static config files, and then describe them as a list of files and where they should be deployed to.

Using the [UseGalaxy.eu](https://github.com/usegalaxy-eu/infrastructure-playbook/blob/02ca578211bfee45044facf36635d28208e5dbb3/group_vars/galaxy.yml#L578) configuration as an example, we have something like:


{% raw %}
```yaml
galaxy_config_templates:
  - src: templates/galaxy/config/builds.txt
    dest: "{{ galaxy_config.galaxy.builds_file_path }}"
  - src: templates/galaxy/config/data_manager_conf.xml
    dest: "{{ galaxy_config.galaxy.data_manager_config_file }}"
  - src: templates/galaxy/config/datatypes_conf.xml
    dest: "{{ galaxy_config.galaxy.datatypes_config_file }}"
  - src: templates/galaxy/config/dependency_resolvers_conf.xml
    dest: "{{ galaxy_config.galaxy.dependency_resolvers_config_file }}"
  - src: templates/galaxy/config/disposable_email_blocklist.conf
    dest: "{{ galaxy_config.galaxy.blocklist_file }}"
```
{% endraw %}

The configuration here is a bit different, it references the `galaxy_config`, which is structured like:

{% raw %}
```yaml
galaxy_config:
  galaxy:
    builds_file_path: "{{ galaxy_config_dir  }}/builds.txt"
    datatypes_config_file: "{{ galaxy_config_dir  }}/datatypes_conf.xml"
```
{% endraw %}

So the references in `galaxy_config_templates` to `galaxy_config` are done to ensure that the setting for e.g. "location of the datatypes config file" is the same between where we have configured Galaxy to looking for it, and where the file has been deployed, without requiring us to make variables changes in numerous places.

{% snippet topics/admin/faqs/ansible_define-once.md %}

### Dependencies

Now that Galaxy is available on disk, Ansible is ready to start processing [dependencies](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/dependencies.yml) of Galaxy.

1. The virtualenv is updated with data from the `galaxy_requirements_file`, by default pointing to the requirements file in the codebase: {% raw %}`{{ galaxy_server_dir }}/lib/galaxy/dependencies/pinned-requirements.txt`.{% endraw %}
2. Any necessary conditional dependencies of Galaxy are [collected by processing the config file](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/dependencies/__init__.py)
3. and then installed to the virtualenv.

### Mutable Setup

[This task](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/mutable_setup.yml) creates a directory and initializes "mutable" (written/managed by Galaxy itself) configuration files. It also deploys any hand-managed mutable config files, but it is unlikely that you want to manage these directly, as Galaxy does a sufficient job. Any changes you make to Galaxy, for example installing some tools, would result in the tools being "forgotten about", if you re-ran the playbook and overwrote the `shed_tool_conf.xml` mutable config file with a hand-managed one.

### Managing the Database

The [database management tasks](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/database.yml) are extremely convenient; any time you run the playbook to update Galaxy, this will automatically run the database schema migration as needed.

1. Galaxy first obtains the current DB version and the maximum possible DB version based on the codebase.
2. If needed, the database is created.
3. Both numbers are reported for the runner of the playbook.
4. If the numbers are different, then Ansible runs the command to upgrade the database to the latest version.

As an administrator who often forgot to run the upgrade, and would only notice it once Galaxy crashed during startup, having this process completely automated is extremely nice.

### Building the Client

Galaxy is a modern web application that includes both a server (written in Python) and a client (written in Javascript). After the server is installed and its database prepared, the next step is to build the client application. This means fetching its dependencies, bundling components, creating minified copies of static content, etc. This process ensures that the smallest possible amount of data is transferred to the user when accessing Galaxy over the web, which is important for the performance of the website.

This process can be lengthy and resource intensive. Future plans for Galaxy include pre-building the client so that the build process is not necessary as long as you run "release" versions of Galaxy.

The client lives in the Galaxy code under the [client/](https://github.com/galaxyproject/galaxy/blob/dev/client/) directory, and the build process deploys it to the [static/](https://github.com/galaxyproject/galaxy/tree/dev/static/) directory, which we will configure a webserver to serve in this tutorial.

## Handlers

A number of the tasks that are executed will trigger a restart of Galaxy. Currently the only auto-magic implementation of this restart [handler](https://docs.ansible.com/ansible/2.9/user_guide/playbooks_intro.html#handlers-running-operations-on-change) is when Galaxy is started and run using systemd. Thankfully, this tutorial uses systemd, but if systemd is not an option wherever you're deploying Galaxy, you will have to do something custom that fits for your setup. If this is the case for you, the role provides a way to reference your own handler by setting `galaxy_restart_handler_name`.

## Defaults

As with other roles, numerous [default values](https://github.com/galaxyproject/ansible-galaxy/blob/master/defaults/main.yml) are provided, but these are useful mostly as reference, and not to go through individually.

## Summary

Installation of Galaxy with the playbook follows generally the steps you would expect:

- Galaxy is cloned (or updated)
- A virtualenv is created if it doesn't exist
- Configuration files are installed
- Any missing dependencies are installed
- Any database updates are applied
- The client application is built and deployed

It would not be difficult to write a role that does this yourself, but by using
the `galaxyproject.galaxy` role, you know that you're getting all of the Galaxy
best practices and knowledge from previous admins codified for you.

# Installing Galaxy

With the necessary background in place, you are ready to install Galaxy with Ansible. The playbooks will start simple, and grow over time. We will start with the minimal Galaxy playbook which only requires setting the `galaxy_root` and expand from there. First, however, we need a database for Galaxy to connect to, so we will do that now.

To proceed from here it is expected that:

> <comment-title>Requirements for Running This Tutorial</comment-title>
>
> 1. You have [Ansible installed](https://docs.ansible.com/ansible/2.9/installation_guide/intro_installation.html) on the machine where you will install Galaxy
>
>    > <comment-title>Running Ansible on remote machine</comment-title>
>    > It is possible to have Ansible installed on your laptop/local machine and run it against some remote hosts as well. We will **not** do that in this training.
>    {: .comment}
>
> 2. Your `ansible` version is `>=2.10`, you can check this by running `ansible --version`
> 3. You have an [inventory file](../ansible/tutorial.html#inventory-file) with the VM or host specified where you will deploy Galaxy. We will refer to this group of hosts as "galaxyservers."
> 4. Your VM has a public DNS name: this tutorial sets up SSL certificates from the start and as an integral part of the tutorial.
> 5. Your VM has `python3` installed.
>
>    If you follow the official Ansible instructions to install Ansible on a blank machine, you will probably end up with py2 installed. You will need to install python3 in addition.
>
> 6. In your inventory file, you have written the full DNS hostname that has been provided, and **not** `localhost`, as we will be requesting SSL certificates.
>
> 7. You have the following ports exposed:
>
>    - 22 for SSH, this can be a different port or via VPN or similar.
>    - 80 for HTTP, this needs to be available to the world if you want to follow the LetsEncrypt portion of the tutorial.
>    - 443 for HTTPs, this needs to be available to the world if you want to follow the LetsEncrypt portion of the tutorial.
>    - 5671 for AMQP for Pulsar, needed if you plan to setup Pulsar for remote job running.
>
{: .comment}


> <comment-title>Ubuntu or Debian, CentOS or RHEL?</comment-title>
> The training tutorial is only tested on Ubuntu, as that is the training platform most often used. The training *should* work on RHEL/CentOS, and the roles we use are definitely compatible with both. If any of the variable values differ between Ubuntu and CentOS, we try to note it in the tutorial. Any places we don't note it are bugs.
>
> There are known issues with CentOS7 and python3, which is used in this tutorial. If you use this setup, you should have python2 and python3 coexisting, and use python2 for Ansible (i.e. do not set `interpreter_python` in your `ansible.cfg`) and python3 for Galaxy. This setup requires numerous changes that you will need to discover, which are not covered in this tutorial.
>
> Both `python-psycopg2` and `python3-psycopg2` need to be installed.
>
> CentOS8 (and RHEL in general) have significantly different package names for some modules, beware!
{: .comment}


> <warning-title>22.05 not yet released</warning-title>
> We updated the tutorials ahead of GCC, preparing for the GCC Galaxy Admin Training. However, 22.05 is not yet released and we still have some bugs in the roles to work our, so, until that time we recommend you follow the [GTN Training Archive](https://training.galaxyproject.org/archive/2022-05-01/topics/admin/tutorials/ansible-galaxy/tutorial.html) materials for this tutorial (and all others in this series.)
{: .warning}

## Requirements

We have codified all of the dependencies you will need into a YAML file that `ansible-galaxy` can install.

> <hands-on-title>Installing roles</hands-on-title>
>
> 1. Create a new directory `galaxy` in your home folder, and `cd` into that directory
>
> 2. Create a new file in your working directory called `requirements.yml` and include the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/requirements.yml
>    @@ -0,0 +1,14 @@
>    +- src: galaxyproject.galaxy
>    +  version: 0.10.4
>    +- src: galaxyproject.nginx
>    +  version: 0.7.0
>    +- src: galaxyproject.postgresql
>    +  version: 1.0.3
>    +- src: galaxyproject.postgresql_objects
>    +  version: 1.1.0
>    +- src: geerlingguy.pip
>    +  version: 2.0.0
>    +- src: uchida.miniconda
>    +  version: 0.3.0
>    +- src: usegalaxy_eu.certbot
>    +  version: 0.1.5
>    {% endraw %}
>    ```
>    {: data-commit="Add requirements"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
>    > <details-title>What do each of these roles do?</details-title>
>    > We'll cover it in more detail as we use each of the roles but briefly:
>    >
>    >  | Role | Purpose |
>    >  |------|---------|
>    >  |`galaxyproject.galaxy` | Installs and configures the Galaxy application|
>    >  |`galaxyproject.nginx` | Sets up a webserver|
>    >  |`galaxyproject.postgresql` | Installs our database, PostgreSQL|
>    >  |`galaxyproject.postgresql_objects` | Creates users and databases within PostgreSQL|
>    >  |`geerlingguy.pip` | Ensures that pip is available|
>    >  |`uchida.miniconda` | Installs miniconda, which is used by Galaxy|
>    >  |`usegalaxy_eu.certbot` | Installs certbot and requests SSL certificates|
>    {: .details}
>
>    > <details-title>Role version vs Galaxy version?</details-title>
>    > Q: Is there a correspondence between galaxy role versions and galaxy versions?
>    > A: They  are correlated (because generally new Galaxy versions could introduce e.g. new configuration options that the galaxy role would then need to adopt), but the Galaxy role can install older versions of Galaxy.
>    {: .details}
>
> 3. In the same directory, run:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    This will install all of the required modules for this training into the `roles/` folder. We choose to install to a folder to give you easy access to look through the different roles when you have questions on their behaviour.
>
>    > <tip-title>Using `git`?</tip-title>
>    > If you're using git, it can make sense to tell git to ignore the `roles`
>    > directory, as all of that data can be perfectly recreated from our
>    > variables in the requirements.yml file.
>    >
>    > Simply create a file `.gitignore` with the following content:
>    >
>    > ```console
>    > roles/
>    > .vault-password.txt
>    > ```
>    >
>    > When you run `git status` you'll notice that the `roles/` folder is not listed among the other "Untracked files". We're adding the `.vault-password.txt` file as well, ahead of time, because this should *never* ever be committed to the repository.
>    >
>    > Now you can do `git add .` to add all of the files in the current directory, and not worry about committing generated artifacts!
>    {: .tip}
>
> 4. Inspect the contents of the newly created `roles` directory in your working directory.
{: .hands_on}

> <hands-on-title>Configuration files</hands-on-title>
>
> 1. Create a `ansible.cfg` file (next to your playbook) to [configure settings](https://docs.ansible.com/ansible/2.9/reference_appendices/config.html) like the inventory file (and save ourselves some typing!), or the Python interpreter to use:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/ansible.cfg
>    @@ -0,0 +1,6 @@
>    +[defaults]
>    +interpreter_python = /usr/bin/python3
>    +inventory = hosts
>    +retry_files_enabled = false
>    +stdout_callback = yaml # Use the YAML callback plugin.
>    +bin_ansible_callbacks = True # Use the stdout_callback when running ad-hoc commands.
>    {% endraw %}
>    ```
>    {: data-commit="Add ansible.cfg"}
>
>    > <tip-title>CentOS7</tip-title>
>    > As mentioned in the "Ubuntu or Debian, CentOS or RHEL?" comment above, if you are using CentOS7 do not set `interpreter_python` in `ansible.cfg` .
>    {: .tip}
>
>    There is an additional useful option that you might want to add to your `ansible.cfg` file if you are connecting over SSH:
>
>    ```ini
>    [ssh_connection]
>    pipelining = true
>    ```
>
>    Pipelining will make [Ansible run faster](https://docs.ansible.com/ansible/2.9/reference_appendices/config.html#ansible-pipelining) by significantly reducing the number of new SSH connections that must be opened.
>
> 2. Create the `hosts` inventory file if you have not done so yet, defining a `[galaxyservers]` group with the address of the host where you want to install Galaxy. If you are running Ansible on the same machine where Galaxy will be installed to, you should set the `ansible_connection=local` variable. Lastly, you should explicitly set the `ansible_user` variable to the username to use when connecting to the server. Ansible has changed its behaviour over time regarding whether or not `ansible_user` is defined, and it is most effective to define it explicitly even when it can sometimes be inferred.
>
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > cat hosts
>    > > ```
>    > {: .code-in}
>    >
>    > > <code-out-title>Bash</code-out-title>
>    > >
>    > > Your hostname is probably different:
>    > >
>    > > {% raw %}
>    > > ```diff
>    > > --- /dev/null
>    > > +++ b/hosts
>    > > @@ -0,0 +1,2 @@
>    > > +[galaxyservers]
>    > > +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    > > {% endraw %}
>    > > ```
>    > > {: data-commit="Add hosts"}
>    > {: .code-out}
>    {: .code-2col}
>
{: .hands_on}

## PostgreSQL

Galaxy is capable of talking to multiple databases through SQLAlchemy drivers. SQLite is the development database, but PostgreSQL (/ˈpoʊstɡrɛs ˌkjuː ˈɛl/, POHST-gres kyoo el) is recommended in production. MySQL is a possibility, but does not receive the same testing or bugfixes from the main development team as PostgreSQL, so we will only show installation with PostgreSQL.

PostgreSQL maintains its own user database apart from the system user database. By default, PostgreSQL uses the "peer" authentication method which allows access for system users with matching PostgreSQL usernames (other authentication mechanisms are available, see the [PostgreSQL Client Authentication documentation](https://www.postgresql.org/docs/current/static/client-authentication.html).

For this tutorial, we will use the default "peer" authentication, so we need to create a PostgreSQL user matching the system user under which Galaxy will be running, i.e. `galaxy`. This is normally done with the PostgreSQL `createuser` command, and it must be run as the `postgres` user. In our case, we will use the `galaxyproject.postgresql_objects` role to handle this step. Additionally we're setting a couple of variables to control the automatic backups, they'll be placed in the `/data/backups` folder next to our user uploaded Galaxy data.

> <hands-on-title>Installing PostgreSQL</hands-on-title>
>
> 1. Create and edit `group_vars/galaxyservers.yml` and add some variables to configure PostgreSQL:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/galaxyservers.yml
>    @@ -0,0 +1,15 @@
>    +---
>    +# Python 3 support
>    +pip_virtualenv_command: /usr/bin/python3 -m virtualenv # usegalaxy_eu.certbot, usegalaxy_eu.tiaas2, galaxyproject.galaxy
>    +certbot_virtualenv_package_name: python3-virtualenv    # usegalaxy_eu.certbot
>    +pip_package: python3-pip                               # geerlingguy.pip
>    +
>    +# PostgreSQL
>    +postgresql_objects_users:
>    +  - name: galaxy
>    +postgresql_objects_databases:
>    +  - name: galaxy
>    +    owner: galaxy
>    +# PostgreSQL Backups
>    +postgresql_backup_dir: /data/backups
>    +postgresql_backup_local_dir: "{{ '~postgres' | expanduser }}/backups"
>    {% endraw %}
>    ```
>    {: data-commit="Add initial group variables file"}
>
>    > <tip-title>Using postgres via the network</tip-title>
>    > To use postgres via another machine, or via the network, you can add lines like the following:
>    >
>    > ```yaml
>    > postgresql_pg_hba_conf:
>    > - host all all 127.0.0.1/32 md5
>    > ```
>    >
>    > Here you should either add multiple lines per IP address that needs access to the Postgres server (with netmask `/32`), or a less specific IP range defined also [via netmask.](https://en.wikipedia.org/wiki/Classless_Inter-Domain_Routing#IPv4_CIDR_blocks)
>    >
>    > Additionally you'll need to change the `postgresql_objects_users` statement to include a password (maybe stored in a vault, discussed later!)
>    >
>    > ```yaml
>    > postgresql_objects_users:
>    >   - name: galaxy
>    >     password: super-secret-password-that-should-be-in-vault
>    > ```
>    {: .tip}
>
> 2. Create and open `galaxy.yml` which will be our playbook. Add the following:
>
>    - Add a pre-task to install the necessary dependencies: `python3-psycopg2`, `acl`
>    - A role for `galaxyproject.postgresql`. This will handle the installation of PostgreSQL.
>    - A role for `galaxyproject.postgresql_objects`, run as the postgres user. (You will need `become`/`become_user`.) This role allows for managing users and databases within postgres.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/galaxy.yml
>    @@ -0,0 +1,13 @@
>    +---
>    +- hosts: galaxyservers
>    +  become: true
>    +  become_user: root
>    +  pre_tasks:
>    +    - name: Install Dependencies
>    +      package:
>    +        name: ['acl', 'python3-psycopg2']
>    +  roles:
>    +    - galaxyproject.postgresql
>    +    - role: galaxyproject.postgresql_objects
>    +      become: true
>    +      become_user: postgres
>    {% endraw %}
>    ```
>    {: data-commit="Add initial galaxy playbook"}
>
>    > <question-title></question-title>
>    >
>    > What should we expect to see in our folder?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > - a playbook
>    > > - a hosts file
>    > > - a requirements file
>    > > - Some group variables
>    > > - all of our roles
>    > >
>    > > > <code-in-title>Bash</code-in-title>
>    > > > Let's investigate our directory.
>    > > > ```bash
>    > > > tree -L 2
>    > > > ```
>    > > {: .code-in}
>    > >
>    > > > <code-out-title>Bash</code-out-title>
>    > > >
>    > > > Your directory should look approximately like this:
>    > > >
>    > > > ```ini
>    > > > .
>    > > > ├── ansible.cfg
>    > > > ├── galaxy.yml
>    > > > ├── group_vars
>    > > > │   └── galaxyservers.yml
>    > > > ├── hosts
>    > > > ├── requirements.yml
>    > > > └── roles
>    > > >     ├── galaxyproject.galaxy
>    > > >     ├── galaxyproject.nginx
>    > > >     ├── galaxyproject.postgresql
>    > > >     ├── geerlingguy.pip
>    > > >     ├── galaxyproject.postgresql_objects
>    > > >     ├── uchida.miniconda
>    > > >     └── usegalaxy_eu.certbot
>    > > >
>    > > > 9 directories, 5 files
>    > > > ```
>    > > {: .code-out.code-max-300}
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>    {% snippet topics/admin/faqs/ansible_role-prefix.md %}
>
>    {% snippet topics/admin/faqs/ansible_yaml.md %}
>
> 3. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    >
>    > ```ini
>    > PLAY [galaxyservers] ************************************************************************************
>    >
>    > TASK [Gathering Facts] **********************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [Install Dependencies] *****************************************************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : include_tasks] *********************************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.postgresql/tasks/debian.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.postgresql : Install pgdg package signing key (Debian/pgdg)] ************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Install pgdg repository (Debian/pgdg)] *********************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Install PostgreSQL (Debian)] *******************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Get installed version] *************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set version fact] ******************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set version fact] ******************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set OS-specific variables] *********************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set pgdata fact] *******************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set conf dir fact] *****************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : include_tasks] *********************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Create conf.d] *********************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Check for conf.d include in postgresql.conf] ***************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set conf.d include in postgresql.conf] *********************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Include 25ansible_postgresql.conf in postgresql.conf] ******************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set config options] ****************************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Install pg_hba.conf] ***************************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : include_tasks] *********************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Ensure PostgreSQL is running] ******************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql_objects : Revoke extra privileges] *********************************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Drop databases] ******************************************************
>    > skipping: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy', 'owner': 'galaxy'})
>    >
>    > TASK [galaxyproject.postgresql_objects : Create and drop users] ***********************************************
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy'})
>    > [WARNING]: Module did not set no_log for no_password_changes
>    >
>    > TASK [galaxyproject.postgresql_objects : Create groups] *******************************************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Add or remove users from groups] *************************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Drop groups] *********************************************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Create databases] ****************************************************
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy', 'owner': 'galaxy'})
>    >
>    > TASK [galaxyproject.postgresql_objects : Grant user database privileges] **************************************
>    > skipping: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy'})
>    >
>    > TASK [galaxyproject.postgresql_objects : Grant extra privileges] **********************************************
>    >
>    > RUNNING HANDLER [galaxyproject.postgresql : Reload PostgreSQL] ******************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > PLAY RECAP **********************************************************************************************
>    > gat-0.training.galaxyproject.eu : ok=17   changed=7    unreachable=0    failed=0    skipped=14   rescued=0    ignored=0
>    >
>    > ```
>    {: .code-out.code-max-300}
>
>    > <comment-title>Error: `postgresql_version is version_compare('9.3', '>=')' failed`</comment-title>
>    > You might see an error like this, if you're running with `--check` mode. Unfortunately here `--check` mode meets real life: not all modules support it because some rely on running command line tools to obtain version numbers, something that cannot be mocked with check mode. You can solve this by running it without `--check`.
>    >
>    > ```
>    > fatal: [gat-34.us.galaxy.training]: FAILED! => {"msg": "The conditional check 'postgresql_version is version_compare('9.3', '>=')' failed. The error was: Input version value cannot be empty\n\nThe error appears to be in '/home/ubuntu/galaxy/roles/galaxyproject.postgresql/tasks/main.yml': line 42, column 3, but may\nbe elsewhere in the file depending on the exact syntax problem.\n\nThe offending line appears to be:\n\n# etc.). So check for a match first and then add if there's no match.\n- name: Check for conf.d include in postgresql.conf\n ^ here\n"}
>    > ```
>    {: .comment}
>
>    > <comment-title>When running Ansible</comment-title>
>    > Always pay close attention to tasks reported as **changed** and ensure that the changes were expected!
>    {: .comment}
>
>    > <tip-title>Why didn't we use `-i` in our ansible command?</tip-title>
>    > In our `ansible.cfg` file we specified the inventory was stored in a file called `hosts`:
>    >
>    > ```ini
>    > [defaults]
>    > interpreter_python = /usr/bin/python3
>    > inventory = hosts
>    > retry_files_enabled = false
>    > ```
>    {: .tip}
>
>    {% snippet topics/admin/faqs/ansible_error-no-hosts-matched.md %}
>
> 4. Inspect the changes that have been made on your Galaxy server. Places to look include:
>
>    - `/etc/postgresql`
>    - Databases and users in PostgreSQL.
>
>    > > <code-in-title>Bash</code-in-title>
>    > > List available databases.  You can now login and access the database, but only as the `postgres` user. You will need to `sudo -iu postgres` first, and then you can run `psql galaxy`.
>    > > ```bash
>    > > sudo -iu postgres psql -l
>    > > ```
>    > {: .code-in}
>    >
>    > > <code-out-title>Bash</code-out-title>
>    > >
>    > > ```ini
>    > >                               List of databases
>    > >    Name    |  Owner   | Encoding | Collate |  Ctype  |   Access privileges
>    > > -----------+----------+----------+---------+---------+-----------------------
>    > >  galaxy    | galaxy   | UTF8     | C.UTF-8 | C.UTF-8 |
>    > >  postgres  | postgres | UTF8     | C.UTF-8 | C.UTF-8 |
>    > >  template0 | postgres | UTF8     | C.UTF-8 | C.UTF-8 | =c/postgres          +
>    > >            |          |          |         |         | postgres=CTc/postgres
>    > >  template1 | postgres | UTF8     | C.UTF-8 | C.UTF-8 | =c/postgres          +
>    > >            |          |          |         |         | postgres=CTc/postgres
>    > > (4 rows)
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    > <code-in-title>Bash</code-in-title>
>    > Access postgres and list some users
>    > ```bash
>    > sudo -iu postgres psql
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    > The database will currently be empty (check relations with \d and list of roles with \dg) as Galaxy has not yet connected to it. Once you install Galaxy in the next step, the database will be populated. (Quit database with \q)
>    > ```ini
>    > psql (10.12 (Ubuntu 10.12-0ubuntu0.18.04.1))
>    > Type "help" for help.
>    >
>    > postgres=# \d
>    > Did not find any relations.
>    > postgres=# \dg
>    >                                    List of roles
>    >  Role name |                         Attributes                         | Member of
>    > -----------+------------------------------------------------------------+-----------
>    >  galaxy    |                                                            | {}
>    >  postgres  | Superuser, Create role, Create DB, Replication, Bypass RLS | {}
>    >
>    > postgres=# \q
>    > ```
>    {: .code-out}
>
>    > > <code-in-title>Bash</code-in-title>
>    > > Inspect the postgres directory
>    > > ```bash
>    > > ls /etc/postgresql/*/main/
>    > > ```
>    > {: .code-in}
>    >
>    > > <code-out-title>Bash</code-out-title>
>    > >
>    > > ```ini
>    > > conf.d
>    > > environment
>    > > pg_ctl.conf
>    > > pg_hba.conf
>    > > pg_hba.conf.3350.2020-07-08@11:51:47~
>    > > pg_ident.conf
>    > > postgresql.conf
>    > > start.conf
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
{: .hands_on}


## Galaxy

Next we will dive right in to deploying a copy of Galaxy onto our server, but it will just be a static copy of the code without anything running.

For a normal Galaxy instance there are a few configuration changes you make very early during deployment:

- Changing the database connection
- Configuring the admin user list
- Changing the "brand"

Additionally we'll go ahead and set up the production-ready Gunicorn which will handle processing Galaxy jobs. Gunicorn launches as many workers as you request, and then they take turns placing a lock, accepting a job, releasing that lock, and then going on to process that job.

Finally, best admin practices are to not run Galaxy as a user with `sudo` access, like your login user probably has. Additionally, it is best to install the Galaxy code and configs as a separate user, for security purposes. So we will instruct the `galaxyproject.galaxy` role to create a new user account specifically to run Galaxy under.

> <details-title>Mules are not the only option</details-title>
>
> Galaxy can be run in a [couple of other configurations](https://docs.galaxyproject.org/en/master/admin/scaling.html#deployment-options) depending on your needs. Mules are generally a good solution for most production needs.
>
{: .details}

The configuration is quite simple thanks to the many sensible defaults that are provided in the Ansible roles.

> <hands-on-title>Minimal Galaxy Playbook</hands-on-title>
>
> 1. Open `galaxy.yml` with your text editor and set the following:
>
>    - Amend the [package installation](https://docs.ansible.com/ansible/2.9/modules/package_module.html#package-module) pre-task to install some additional necessary dependencies: `bzip2`, `git`, `make`, `tar`, and `virtualenv`.
>    - Add the roles `geerlingguy.pip`, `galaxyproject.galaxy` and `uchida.miniconda` (in this order) at the end, with `uchida.miniconda` run as the `galaxy` user.
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -5,9 +5,14 @@
>       pre_tasks:
>         - name: Install Dependencies
>           package:
>    -        name: ['acl', 'python3-psycopg2']
>    +        name: ['acl', 'bzip2', 'git', 'make', 'python3-psycopg2', 'tar', 'virtualenv']
>       roles:
>         - galaxyproject.postgresql
>         - role: galaxyproject.postgresql_objects
>           become: true
>           become_user: postgres
>    +    - geerlingguy.pip
>    +    - galaxyproject.galaxy
>    +    - role: uchida.miniconda
>    +      become: true
>    +      become_user: "{{ galaxy_user.name }}"
>    {% endraw %}
>    ```
>    {: data-commit="Add pip, miniconda, galaxy to playbook"}
>
>    > <tip-title>Miniconda fails to work</tip-title>
>    > The Galaxy user is created to separate privileges. Then we add `uchida.miniconda`, which is run as the Galaxy user.
>    >
>    > The miniconda role attempts to install tar and bzip2, even when the user doesn't have permissions to do this. This issue has been addressed in the codebase, but, the role has not seen a release which would address this. You do *not* need to run this role, Galaxy will attempt to install conda when it is missing. We added it to prevent a possible race condition between the two mules if both attempt to install conda at the same time. So, remove the role from your playbook, and carry on, if you have issues.
>     {: .tip}
>
> 2. Edit your group variables file for your group (`group_vars/galaxyservers.yml`).
>
>    We need to set the following variables at the top level:
>
>    Variable                        | Value                                                            | Purpose
>    ---                             | -----                                                            | ---
>    `galaxy_create_user`            | `true`                                                           | Instruct the role to create a Galaxy user
>    `galaxy_separate_privileges`    | `true`                                                           | Enable separation mode to install the Galaxy code as `root` but run the Galaxy server as `galaxy`
>    `galaxy_manage_paths`           | `true`                                                           | Instruct the role to create the needed directories.
>    `galaxy_layout`                 | `root-dir`                                                       | This enables the `galaxy_root` Galaxy deployment layout: all of the code, configuration, tools, and mutable-data (like caches, location files, etc.) folders will live by default beneath `galaxy_root`. User data is stored under `file_path`, a variable we will set later.
>    `galaxy_root`                   | `/srv/galaxy`                                                    | This is the root of the Galaxy deployment.
>    `galaxy_user`                   | `{name: galaxy, shell: /bin/bash}`                               | The user that Galaxy will run as.
>    `galaxy_commit_id`              | `release_22.05`                                                  | The git reference to check out, which in this case is the branch for Galaxy Release 22.05
>    `galaxy_force_checkout`         | `true`                                                           | If we make any modifications to the Galaxy codebase, they will be removed. This way we know we're getting an unmodified Galaxy and no one has made any unexpected changes to the codebase.
>    `miniconda_prefix`              | {% raw %}`"{{ galaxy_tool_dependency_dir }}/_conda"`{% endraw %} | We will manually install conda as well. Normally Galaxy will attempt to auto-install this, but since we will set up a production-ready instance with multiple handlers, there is the chance that they can get stuck.
>    `miniconda_version`             | `4.7.12`                                                         | Install a specific miniconda version, the latest one at the time of writing that was tested and working.
>    `miniconda_manage_dependencies` | `false`                                                          | Specify whether to install the miniconda installer dependencies.
>
>    > <tip-title>Different Galaxy Releases!</tip-title>
>    > In the time between this tutorial was last updated ({{ page.last_modified_at | date: "%Y-%m-%d" }}), and when you are now reading it, one or more new releases of Galaxy may have occured.
>    >
>    > If you wish to use a different version of Galaxy than the one listed, you can do so, but there is no guarantee this tutorial will work without modifications.
>    >
>    > Please ensure you consult the [releases documentation](https://docs.galaxyproject.org/en/master/releases/index.html)
>    {: .tip}
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -13,3 +13,16 @@ postgresql_objects_databases:
>     # PostgreSQL Backups
>     postgresql_backup_dir: /data/backups
>     postgresql_backup_local_dir: "{{ '~postgres' | expanduser }}/backups"
>    +
>    +# Galaxy
>    +galaxy_create_user: true
>    +galaxy_separate_privileges: true
>    +galaxy_manage_paths: true
>    +galaxy_layout: root-dir
>    +galaxy_root: /srv/galaxy
>    +galaxy_user: {name: galaxy, shell: /bin/bash}
>    +galaxy_commit_id: release_22.05
>    +galaxy_force_checkout: true
>    +miniconda_prefix: "{{ galaxy_tool_dependency_dir }}/_conda"
>    +miniconda_version: 4.7.12
>    +miniconda_manage_dependencies: false
>    {% endraw %}
>    ```
>    {: data-commit="Configure miniconda and galaxy"}
>
>    > <tip-title>Following this training outside of a GAT course?</tip-title>
>    > Consider updating the Galaxy Commit ID to the latest version of Galaxy available, this will probably give better results (given that dependencies are always updating) than using a potentially outdated version.
>    {: .tip}
>
> 3. Again edit the group variables file and add a variable for `galaxy_config`. It will be a hash with one key, `galaxy` which will also be a hash. Inside here you can place all of your Galaxy configuration.
>
>    So the structure looks like:
>    ```yaml
>    galaxy_config:
>       galaxy:
>           key: value
>    ```
>
>    Now you should set:
>    1. `admin_users` to the email address you will use with this Galaxy.
>    2. `brand` to something fun! (E.g. 🧬🔬🚀)
>    3. `database_connection` to point to the database you setup earlier (`postgresql:///galaxy?host=/var/run/postgresql`).
>    4. `file_path` to a place to store data, `/data` is fine for this lesson which sets up a single-node Galaxy. If you have separate compute machines, this will normally need to be storage shared between the Galaxy node and compute nodes.
>    5. `check_migrate_tools` must be set to `false` due to a new installation of Galaxy.
>    6. `tool_data_path` to {% raw %}`{{ galaxy_mutable_data_dir }}/tool-data`{% endraw %}, so that when tools are installed, due to privilege separation, this will happen in a directory Galaxy can actually write into.
>    7. `object_store_store_by` to `uuid`, this is a better way of storing files that will ensure better filesystem balancing than the older system.
>    8. `id_secret` to {% raw %}`{{ vault_id_secret }}`{% endraw %}, we'll define this variable next but it will be used to encode the IDs used in Galaxy URLs and for securing session cookies.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -26,3 +26,14 @@ galaxy_force_checkout: true
>     miniconda_prefix: "{{ galaxy_tool_dependency_dir }}/_conda"
>     miniconda_version: 4.7.12
>     miniconda_manage_dependencies: false
>    +
>    +galaxy_config:
>    +  galaxy:
>    +    brand: "🧬🔬🚀"
>    +    admin_users: admin@example.org
>    +    database_connection: "postgresql:///galaxy?host=/var/run/postgresql"
>    +    file_path: /data
>    +    check_migrate_tools: false
>    +    tool_data_path: "{{ galaxy_mutable_data_dir }}/tool-data"
>    +    object_store_store_by: uuid
>    +    id_secret: "{{ vault_id_secret }}"
>    {% endraw %}
>    ```
>    {: data-commit="Configure galaxy config"}
>
>    > <tip-title>Data storage</tip-title>
>    > Galaxy datasets cannot be separated by user or other attribute currently, but you can spread data unintelligently across 1 or more storage pools.
>    {: .tip}
>
>    > <tip-title>PostgreSQL connection string</tip-title>
>    > If you want to run your database on a different machine, you will need to change the connection string. In your hosts file, place the hostname of the machine you're installing on. `ansible_connection` can be left off entirely and it will connect over ssh. You may need to set `ansible_user` to the username of the admin user (who can run sudo).
>    >
>    > Here are some examples of connection strings:
>    >
>    > ```
>    > sqlite:///./database/universe.sqlite?isolation_level=IMMEDIATE
>    > postgres://<name>:<password>@localhost:5432/galaxy
>    > postgresql:///galaxy?host=/var/run/postgresql
>    > ```
>    >
>    {: .tip}
>
>    > <tip-title>Using postgres via the network</tip-title>
>    > Here you'll need to re-use the connection details you specified during `postgresql_objects_users`. You can reference it like the following, for example. `localhost` may need to change if you're hosting the database on another host.
>    > ```
>    > +    database_connection: "postgres://{{ postgresql_objects_users[0].name }}:{{ postgresql_objects_users[0].password }}@localhost:5432/{{ postgresql_objects_databases[0].name }}"
>    > ```
>    {: .tip}
>
>    > <comment-title>Ansible Variable Templating</comment-title>
>    > In this step we use some templated variables. These are seen in our group variables, among other places, and look like {% raw %}`miniconda_prefix: "{{ galaxy_tool_dependency_dir  }}/_conda"`{% endraw %}.
>    >
>    > When Ansible runs:
>    >
>    > 1. It collects variables defined in group variables and other places
>    > 2. The first task for each machine is the [`setup` module](https://docs.ansible.com/ansible/2.9/modules/setup_module.html) which gathers facts about the host, which are added to the available variables
>    > 3. As roles are executed:
>    >    1. Their defaults are added to the set of variables (the group variables having precedence over these variables)
>    >    2. They can also dynamically define more variables which may not be set until that role is run
>    > 4. Before use (in templates, commands, etc.), variables are resolved to their final value
>    > So it is not always easy to tell what variables will be set, or what their finaly value will be, without running the playbook. It is possible, but non trivial.
>    >
>    {: .comment}
>
> 4. In order to enable the Gunicorn + Webless strategy, we need to edit the Gravity configuration of Galaxy. This has a default value, but we will have to override it. Add the following configuration as a child of the `galaxy_config` variable:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -37,3 +37,23 @@ galaxy_config:
>         tool_data_path: "{{ galaxy_mutable_data_dir }}/tool-data"
>         object_store_store_by: uuid
>         id_secret: "{{ vault_id_secret }}"
>    +  gravity:
>    +    galaxy_root: "{{ galaxy_root }}/server"
>    +    app_server: gunicorn
>    +    gunicorn:
>    +      # listening options
>    +      bind: "unix:{{ galaxy_mutable_config_dir }}/gunicorn.sock"
>    +      # performance options
>    +      workers: 2
>    +      # Other options that will be passed to gunicorn
>    +      extra_args: '--forwarded-allow-ips="*"'
>    +      preload: true
>    +    celery:
>    +      concurrency: 2
>    +      loglevel: DEBUG
>    +    handlers:
>    +      handler:
>    +        processes: 3
>    +        pools:
>    +          - job-handler
>    +          - workflow-scheduler
>    {% endraw %}
>    ```
>    {: data-commit="Configure gravity"}
>
>
>    > <tip-title>Options explanations.</tip-title>
>    > 1`workers`: Controls the number of Galaxy application processes Gunicorn will spawn. Increased web performance can be attained by increasing this value. If Gunicorn is the only application on the server, a good starting value is the number of CPUs * 2 + 1. 4-12 workers should be able to handle hundreds if not thousands of requests per second.
>    > 2 `extra_args`: You can specify additional arguments to pass to gunicorn here.
>    {: .tip}
>
> 5. Let's set up our vault to store the secrets for these tutorials.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > openssl rand -base64 24 > .vault-password.txt
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 6. Tell Ansible where to find the decryption key. Edit your file `ansible.cfg` and add the `vault_password_file` variable.
>
>    {% raw %}
>    ```diff
>    --- a/ansible.cfg
>    +++ b/ansible.cfg
>    @@ -4,3 +4,4 @@ inventory = hosts
>     retry_files_enabled = false
>     stdout_callback = yaml # Use the YAML callback plugin.
>     bin_ansible_callbacks = True # Use the stdout_callback when running ad-hoc commands.
>    +vault_password_file = .vault-password.txt
>    {% endraw %}
>    ```
>    {: data-commit="Setup the vault file"}
>
> 7. Create the vault:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault create group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    This will open the editor. Within that file, define your `vault_id_secret` to be a long random value.
>
>    ```
>    vault_id_secret: BxI6zlQVhoHLPVf3gqQ
>    ```
>    > <tip-title>How to get a good random value?</tip-title>
>    > You can use the same command we ran before: `openssl rand -base64 24` to get a good, secure `id_secret` for your Galaxy
>    {: .tip}
>
>    Your vault should probably look *something* like this, it's safe to commit to GitHub (or another public {DVCS}) now that you've got the secrets safely encrypted.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/secret.yml
>    @@ -0,0 +1,7 @@
>    +$ANSIBLE_VAULT;1.1;AES256
>    +32653961383866636531396135663630386630346237333333653633313436663439643535323964
>    +6363626330336430363332643638646262316338313937320a666566306539373462386266383166
>    +30326165393863633463353234613561393939326164376432633732316264636464313061383161
>    +3532373937656138320a616361343664353264613332616236623231326137316635323465623562
>    +66656539346130353639623736633034653932373438663330646436656336666637313933666264
>    +3636313438626533633831323239373461373538646635613637
>    {% endraw %}
>    ```
>    {: data-commit="Add id_secret to the vault" data-comment="The password is password"}
>
>    > <tip-title>Using `git`?</tip-title>
>    > If you're using git, and storing all of these steps in your git history,
>    > you'll find that working with Vault secrets can be rather annoying as
>    > the opaque blobs are impossible to read. There's a way to make this
>    > *much* easier though:
>    >
>    > Create a file `.gitattributes` with the following content:
>    >
>    > ```
>    > group_vars/secret.yml diff=ansible-vault merge=binary
>    > ```
>    >
>    > You can add this file to your repository with `git add .gitattributes`
>    > to ensure colleagues get a copy of the file too. Just **be sure**
>    > `.vault-password.txt` is listed in your `.gitignore` file!
>    >
>    > If you have more vault secrets, you can adjust this line (or add more,
>    > wildcards are supported) to list all of your secret files. This tells
>    > `git` to use `ansible-vault` to diff the two files, as you can see in
>    > the following real (redacted) snippet from a vault stored in a public
>    > github project.
>    >
>    > ```console
>    > $ git log -p group_vars/all/secret.yml
>    > commit a137d7d6aa4ed374c29545ac4728837815c460aa
>    > Author: Helena Rasche <hxr@hx42.org>
>    > Date:   Thu Dec 2 12:52:55 2021 +0100
>    >
>    >     fix automation password
>    >
>    > diff --git a/group_vars/all/secret.yml b/group_vars/all/secret.yml
>    > index 3fe4444..11831ae 100644
>    > --- a/group_vars/all/secret.yml
>    > +++ b/group_vars/all/secret.yml
>    > @@ -1,6 +1,6 @@
>    >  ---
>    >  id_secret: "..."
>    > -jenkins_password: "old-secret-value"
>    > +jenkins_password: "new-secret-value"
>    >  grafana_admin_pass:    "..."
>    >  secret_tiaas_admin_pw: "..."
>    >
>    > ```
>    {: .tip}
>
> 8. Load the secrets in the playbook
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -2,6 +2,8 @@
>     - hosts: galaxyservers
>       become: true
>       become_user: root
>    +  vars_files:
>    +    - group_vars/secret.yml
>       pre_tasks:
>         - name: Install Dependencies
>           package:
>    {% endraw %}
>    ```
>    {: data-commit="Load the vault in the playbook"}
>
> 9. Run the playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    >
>    > ```ini
>    > PLAY [galaxyservers] ***************************************************************************
>    >
>    > TASK [Gathering Facts] *************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [Install Dependencies] ********************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : include_tasks] ************************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.postgresql/tasks/debian.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.postgresql : Install pgdg package signing key (Debian/pgdg)] ***************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Install pgdg repository (Debian/pgdg)] ************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Install PostgreSQL (Debian)] **********************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Get installed version] ****************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set version fact] *********************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set version fact] *********************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set OS-specific variables] ************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set pgdata fact] **********************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set conf dir fact] ********************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : include_tasks] ************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Create conf.d] ************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Check for conf.d include in postgresql.conf] ******************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set conf.d include in postgresql.conf] ************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    >
>    > TASK [galaxyproject.postgresql : Include 25ansible_postgresql.conf in postgresql.conf] *********
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Set config options] *******************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Install pg_hba.conf] ******************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : include_tasks] ************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql : Ensure PostgreSQL is running] *********************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.postgresql_objects : Revoke extra privileges] ************************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Drop databases] *********************************************
>    > skipping: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy', 'owner': 'galaxy'})
>    >
>    > TASK [galaxyproject.postgresql_objects : Create and drop users] **************************************
>    > ok: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy'})
>    > [WARNING]: Module did not set no_log for no_password_changes
>    >
>    > TASK [galaxyproject.postgresql_objects : Create groups] **********************************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Add or remove users from groups] ****************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Drop groups] ************************************************
>    >
>    > TASK [galaxyproject.postgresql_objects : Create databases] *******************************************
>    > ok: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy', 'owner': 'galaxy'})
>    >
>    > TASK [galaxyproject.postgresql_objects : Grant user database privileges] *****************************
>    > skipping: [gat-0.training.galaxyproject.eu] => (item={'name': 'galaxy'})
>    >
>    > TASK [galaxyproject.postgresql_objects : Grant extra privileges] *************************************
>    >
>    > TASK [geerlingguy.pip : Ensure Pip is installed.] **********************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [geerlingguy.pip : Ensure pip_install_packages are installed.] ****************************
>    >
>    > TASK [galaxyproject.galaxy : Ensure that mutually exclusive options are not set] ***************
>    > ok: [gat-0.training.galaxyproject.eu] => {
>    >     "changed": false,
>    >     "msg": "All assertions passed"
>    > }
>    >
>    > TASK [galaxyproject.galaxy : Set privilege separation default variables] ***********************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include layout vars] **********************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Set any unset variables from layout defaults] *********************
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_venv_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_server_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_config_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_mutable_data_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_mutable_config_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_shed_tools_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_cache_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy_local_tools_dir)
>    >
>    > TASK [galaxyproject.galaxy : Check that any explicitly set Galaxy config options match the values of explicitly set variables] ***
>    > skipping: [gat-0.training.galaxyproject.eu] => (item=tool_dependency_dir)
>    > skipping: [gat-0.training.galaxyproject.eu] => (item=file_path)
>    > skipping: [gat-0.training.galaxyproject.eu] => (item=job_working_directory)
>    > skipping: [gat-0.training.galaxyproject.eu] => (item=shed_tool_config_file)
>    >
>    > TASK [galaxyproject.galaxy : Set any unset variables corresponding to Galaxy config options from galaxy_config or layout defaults] ***
>    > ok: [gat-0.training.galaxyproject.eu] => (item=tool_dependency_dir)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=file_path)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=job_working_directory)
>    >
>    > TASK [galaxyproject.galaxy : Include user creation tasks] **************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/user.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Create Galaxy group] **********************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Create Galaxy user] ***********************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Create Galaxy privilege separation user] **************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include path management tasks] ************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/paths.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Get group IDs for Galaxy users] ***********************************
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=root)
>    >
>    > TASK [galaxyproject.galaxy : Get group names for Galaxy users] *********************************
>    > ok: [gat-0.training.galaxyproject.eu] => (item=galaxy)
>    > ok: [gat-0.training.galaxyproject.eu] => (item=root)
>    >
>    > TASK [galaxyproject.galaxy : Set Galaxy user facts] ********************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Determine whether to restrict to group permissions] ***************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Create galaxy_root] ***********************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Create additional privilege separated directories] ****************
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/venv)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/server)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/config)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/local_tools)
>    >
>    > TASK [galaxyproject.galaxy : Create additional directories] ************************************
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/var)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/var/config)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/var/cache)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/var/shed_tools)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/var/dependencies)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/data)
>    > changed: [gat-0.training.galaxyproject.eu] => (item=/srv/galaxy/jobs)
>    >
>    > TASK [galaxyproject.galaxy : Include clone tasks] **********************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/clone.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Update Galaxy to specified ref] ***********************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Report Galaxy version change] *************************************
>    > changed: [gat-0.training.galaxyproject.eu] => {
>    >     "msg": "Galaxy version changed from '' to 'e75a792fc631ff216de6664d399f6bbe371086cd'"
>    > }
>    >
>    > TASK [galaxyproject.galaxy : Create Galaxy virtualenv] *****************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Ensure pip is the latest release] *********************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Remove orphaned .pyc files and compile bytecode] ******************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include download tasks] *******************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include static config setup tasks] ********************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/static_setup.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Ensure Galaxy version is set] *************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/_inc_galaxy_version.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Collect Galaxy version file] **************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Determine Galaxy version] *****************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Install additional Galaxy config files (static)] ******************
>    >
>    > TASK [galaxyproject.galaxy : Install additional Galaxy config files (template)] ****************
>    >
>    > TASK [galaxyproject.galaxy : Install local tools] **********************************************
>    >
>    > TASK [galaxyproject.galaxy : Install local_tool_conf.xml] **************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Append local_tool_conf.xml to tool_config_file Galaxy config option] ***
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Append shed_tool_conf.xml to tool_config_file Galaxy config option] ***
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Ensure dynamic job rules paths exists] ****************************
>    >
>    > TASK [galaxyproject.galaxy : Install dynamic job rules] ****************************************
>    >
>    > TASK [galaxyproject.galaxy : Ensure dynamic rule __init__.py's exist] **************************
>    >
>    > TASK [galaxyproject.galaxy : Create Galaxy configuration file] *********************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include dependency setup tasks] ***********************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/dependencies.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Create Galaxy virtualenv] *****************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Ensure pip is the latest release] *********************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Install Galaxy base dependencies] *********************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Collect Galaxy conditional dependency requirement strings] ********
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Install Galaxy conditional dependencies] **************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include mutable config setup tasks] *******************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/mutable_setup.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Ensure Galaxy version is set] *************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Instantiate mutable configuration files] **************************
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_data_manager_conf.xml', 'dest': '/srv/galaxy/var/config/shed_data_manager_conf.xml'})
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_tool_data_table_conf.xml', 'dest': '/srv/galaxy/var/config/shed_tool_data_table_conf.xml'})
>    >
>    > TASK [galaxyproject.galaxy : Instantiate mutable configuration templates] **********************
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_tool_conf.xml.j2', 'dest': '/srv/galaxy/var/config/migrated_tools_conf.xml'})
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_tool_conf.xml.j2', 'dest': '/srv/galaxy/var/config/shed_tool_conf.xml'})
>    >
>    > TASK [galaxyproject.galaxy : Include database management tasks] ********************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/database.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Get current Galaxy DB version] ************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Get maximum Galaxy DB version] ************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Report current and max Galaxy database] ***************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Upgrade Galaxy DB] ************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include client build tasks] ***************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/client.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Ensure client_build_hash.txt exists] ******************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Get current client commit id] *************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Install Galaxy base dependencies] *********************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Collect Galaxy conditional dependency requirement strings] ********
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Install Galaxy conditional dependencies] **************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include mutable config setup tasks] *******************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/mutable_setup.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Ensure Galaxy version is set] *************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Instantiate mutable configuration files] **************************
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_data_manager_conf.xml', 'dest': '/srv/galaxy/var/config/shed_data_manager_conf.xml'})
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_tool_data_table_conf.xml', 'dest': '/srv/galaxy/var/config/shed_tool_data_table_conf.xml'})
>    >
>    > TASK [galaxyproject.galaxy : Instantiate mutable configuration templates] **********************
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_tool_conf.xml.j2', 'dest': '/srv/galaxy/var/config/migrated_tools_conf.xml'})
>    > changed: [gat-0.training.galaxyproject.eu] => (item={'src': 'shed_tool_conf.xml.j2', 'dest': '/srv/galaxy/var/config/shed_tool_conf.xml'})
>    >
>    > TASK [galaxyproject.galaxy : Include database management tasks] ********************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/database.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Get current Galaxy DB version] ************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Get maximum Galaxy DB version] ************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Report current and max Galaxy database] ***************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Upgrade Galaxy DB] ************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include client build tasks] ***************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/client.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Ensure client_build_hash.txt exists] ******************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Get current client commit id] *************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Set client build version fact] ************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Set client build version fact] ************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Report client version mismatch] ***********************************
>    > changed: [gat-0.training.galaxyproject.eu] => {
>    >     "msg": "Galaxy client is out of date: new-unbuilt != e75a792fc631ff216de6664d399f6bbe371086cd"
>    >
>    > }
>    >
>    > TASK [galaxyproject.galaxy : Ensure galaxy_node_version is set] ********************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/_inc_node_version.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Collect Galaxy Node.js version file] ******************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Set Galaxy Node.js version fact] **********************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Report preferred Node.js version] *********************************
>    > ok: [gat-0.training.galaxyproject.eu] => {
>    >     "galaxy_node_version": "10.15.3"
>    >
>    > }
>    >
>    > TASK [galaxyproject.galaxy : Install node] *****************************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Install yarn] *****************************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include client build process] *************************************
>    > included: /home/ubuntu/galaxy/roles/galaxyproject.galaxy/tasks/_inc_client_build_make.yml for gat-0.training.galaxyproject.eu
>    >
>    > TASK [galaxyproject.galaxy : Build client] ****************************************************************************************************************************************************************************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Fetch client version] ********************************************************************************************************************************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Set client build version fact] ***********************************************************************************************************************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Ensure that client update succeeded] *****************************************************************************************************************************************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [galaxyproject.galaxy : Include error document setup tasks] ******************************************************************************************************************************************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : miniconda installer is downloaded] ***********************************************************************************************************************************************************************************
>    > [WARNING]: Module remote_tmp /home/galaxy/.ansible/tmp did not exist and was created with a mode of 0700, this may cause issues when running as another user. To avoid this, create the remote_tmp dir with the correct permissions manually
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : directory /srv/galaxy/var/dependencies exists] ***********************************************************************************************************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : tar is installed] ****************************************************************************************************************************************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : bzip2 is installed] **************************************************************************************************************************************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : miniconda is installed] **********************************************************************************************************************************************************************************************
>    > changed: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : miniconda is up-to-date] *********************************************************************************************************************************************************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : conda environment file /tmp/{{ miniconda_env.name }}-environment.yml is created] *************************************************************************************************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : conda environment {{ miniconda_env.name }} is created] ***************************************************************************************************************************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > TASK [uchida.miniconda : conda environment {{ miniconda_env.name }} is up-to-date] ************************************************************************************************************************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > RUNNING HANDLER [galaxyproject.galaxy : default restart galaxy handler] ***********************************************************************************************************************************************************************
>    > ok: [gat-0.training.galaxyproject.eu] => {
>    >     "msg": "RESTARTER NOT IMPLEMENTED - Please restart Galaxy manually. You can define your own handler and enable it with `galaxy_restart_handler_name`"
>    > }
>    >
>    > RUNNING HANDLER [galaxyproject.galaxy : email administrator with commit id] *******************************************************************************************************************************************************************
>    > skipping: [gat-0.training.galaxyproject.eu]
>    >
>    > PLAY RECAP ************************************************************************************************************************************************************************************************************************************
>    > gat-0.training.galaxyproject.eu : ok=79   changed=20   unreachable=0    failed=0    skipped=38   rescued=0    ignored=0
>    > ```
>    {: .code-out.code-max-300}
>
>    > <comment-title>Free knowledge</comment-title>
>    > You will notice the significant volume of output from this playbook! There are
>    > many, many steps for deploying a production Galaxy, and there is a vast
>    > amount of knowledge encoded into these playbooks by the developers and
>    > admins who use them. While we recommend this, choosing to use Ansible
>    > and the official playbooks means you are opting in to receive all of this
>    > knowledge applied automatically to your systems. You are welcome to use
>    > a different Configuration Management system which works for your infrastructure, but it might
>    > mean re-creating a large amount of existing work that is maintained by a
>    > large community of admins.
>    {: .comment}
>
>    > <tip-title>Slow Deployment</tip-title>
>    > The deployment can be slowed down by the client build.
>    > The client is only re-built when there are changes in the files needed for the Galaxy user interface (JavaScript, CSS). Because we are tracking a release branch, we’ll receive updates that are published to that branch during the training since the last time the playbook was run.
>    {: .tip}
>
>    > <tip-title>ERROR: Failed to set permissions on the temporary files</tip-title>
>    > Did you get an error message like this?
>    > ```
>    > fatal: [localhost]: FAILED! => {"msg": "Failed to set permissions on the temporary files Ansible needs to create when becoming an unprivileged user (rc: 1, err: chown: changing ownership of '/var/tmp/ansible-tmp-1607430009.739602-32983298209838/': Operation not permitted\nchown: changing ownership of '/var/tmp/ansible-tmp-1607430009.739602-32983298209838/source': Operation not permitted\n}). For information on working around this, see https://docs.ansible.com/ansible/become.html#becoming-an-unprivileged-user"}
>    > ```
>    >
>    > You're missing the `setfacl` command provided by the `acl` package. [reference](https://github.com/georchestra/ansible/issues/55#issuecomment-588313638)
>    >
>    {: .tip}
>
> 6. Explore what has been set up for you.
>
>    - Galaxy has been deployed to `/srv/galaxy/server`
>    - The configuration lives in `/srv/galaxy/config/galaxy.yml` - be sure to look through it to see what options have been set for you
>    - Note the ownership and permissions of the contents of `/srv/galaxy`
>    - Some config files that Galaxy maintains itself, such as `shed_tool_conf.xml`, which controls what tools that you have installed from the Tool Shed will be loaded, have been instantiated in `/srv/galaxy/var/config`
>    - A Python virtualenv - an isolated Python environment - with all of the Galaxy framework's dependencies has been installed in `/srv/galaxy/venv`
>
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > tree -L 1 /srv/galaxy/
>    > > ```
>    > > - `config` holds all static Galaxy config files
>    > > - `jobs` is where all job files will go (e.g. temporary working dirs, job scripts, etc.)
>    > > - `local_tools` is a directory for custom, non-ToolShed tools managed by the Ansible playbook
>    > > - `server` contains all of the Galaxy server code
>    > > - `var` is a directory for all files created by Galaxy, e.g. whenever tools are installed from the ToolShed, the Galaxy-managed cache, and the integrated tool panel file.
>    > > - `venv` contains the Galaxy virtual environment and all dependencies, like Gunicorn.
>    > {: .code-in}
>    >
>    > > <code-out-title>Bash</code-out-title>
>    > >
>    > > ```ini
>    > > /srv/galaxy/
>    > > ├── config
>    > > ├── jobs
>    > > ├── local_tools
>    > > ├── server
>    > > ├── var
>    > > └── venv
>    > >
>    > > 6 directories, 0 files
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > cat /srv/galaxy/config/galaxy.yml
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    > You'll notice that the file is significantly different from the configuration you have set up in your group variables. The Ansible role adds a significant number of additional configuration options which all require a path, and templates the appropriate paths into all of them.
>    > ```yaml
>    > ---
>    > ##
>    > ## This file is managed by Ansible.  ALL CHANGES WILL BE OVERWRITTEN.
>    > ##
>    > gravity:
>    >     app_server: gunicorn
>    >     celery:
>    >         concurrency: 2
>    >         loglevel: DEBUG
>    >     galaxy_root: /srv/galaxy/server
>    >     gunicorn:
>    >         bind: unix:/srv/galaxy/var/config/gunicorn.sock
>    >         extra_args: --forwarded-allow-ips="*"
>    >         preload: true
>    >         workers: 2
>    >     handlers:
>    >         handler:
>    >             pools:
>    >             - job-handler
>    >             - workflow-scheduler
>    >             processes: 3
>    >
>    > uwsgi:
>    >     buffer-size: 16384
>    >     die-on-term: true
>    >     enable-threads: true
>    >     farm: job-handlers:1,2
>    >     hook-master-start: unix_signal:2 gracefully_kill_them_all
>    >     hook-master-start: unix_signal:15 gracefully_kill_them_all
>    >     socket: 127.0.0.1:5000
>    >     master: true
>    >     module: galaxy.webapps.galaxy.buildapp:uwsgi_app()
>    >     mule: lib/galaxy/main.py
>    >     mule: lib/galaxy/main.py
>    >     offload-threads: 2
>    >     processes: 1
>    >     py-call-osafterfork: true
>    >     pythonpath: /srv/galaxy/server/lib
>    >     static-map: /static=/srv/galaxy/server/static
>    >     static-map: /favicon.ico=/srv/galaxy/server/static/favicon.ico
>    >     static-safe: client/galaxy/images
>    >     threads: 4
>    >     thunder-lock: true
>    >     virtualenv: /srv/galaxy/venv
>    >
>    > galaxy:
>    >     admin_users: admin@example.org
>    >     brand: 🧬🔬🚀
>    >     builds_file_path: /srv/galaxy/server/tool-data/shared/ucsc/builds.txt.sample
>    >     citation_cache_data_dir: /srv/galaxy/var/cache/citations/data
>    >     citation_cache_lock_dir: /srv/galaxy/var/cache/citations/locks
>    >     data_manager_config_file: /srv/galaxy/server/config/data_manager_conf.xml.sample
>    >     database_connection: postgresql:///galaxy?host=/var/run/postgresql
>    >     datatypes_config_file: /srv/galaxy/server/config/datatypes_conf.xml.sample
>    >     external_service_type_config_file: /srv/galaxy/server/config/external_service_types_conf.xml.sample
>    >     file_path: /data
>    >     integrated_tool_panel_config: /srv/galaxy/var/config/integrated_tool_panel.xml
>    >     job_working_directory: /srv/galaxy/jobs
>    >     mulled_resolution_cache_data_dir: /srv/galaxy/var/mulled/data
>    >     mulled_resolution_cache_lock_dir: /srv/galaxy/var/mulled/lock
>    >     new_file_path: /srv/galaxy/var/tmp
>    >     object_store_cache_path: /srv/galaxy/var/cache/object_store_cache
>    >     openid_config_file: /srv/galaxy/server/config/openid_conf.xml.sample
>    >     openid_consumer_cache_path: /srv/galaxy/var/cache/openid_consumer_cache
>    >     shed_data_manager_config_file: /srv/galaxy/var/config/shed_data_manager_conf.xml
>    >     shed_tool_config_file: /srv/galaxy/var/config/shed_tool_conf.xml
>    >     shed_tool_data_table_config: /srv/galaxy/var/config/shed_tool_data_table_conf.xml
>    >     template_cache_path: /srv/galaxy/var/cache/template_cache
>    >     tool_cache_data_dir: /srv/galaxy/var/cache/tool_cache
>    >     tool_config_file: /srv/galaxy/server/config/tool_conf.xml.sample
>    >     tool_data_path: /srv/galaxy/var/tool-data
>    >     tool_data_table_config_path: /srv/galaxy/server/config/tool_data_table_conf.xml.sample
>    >     tool_dependency_dir: /srv/galaxy/var/dependencies
>    >     tool_search_index_dir: /srv/galaxy/var/cache/tool_search_index
>    >     tool_sheds_config_file: /srv/galaxy/server/config/tool_sheds_conf.xml.sample
>    >     ucsc_build_sites: /srv/galaxy/server/tool-data/shared/ucsc/ucsc_build_sites.txt.sample
>    >     visualization_plugins_directory: config/plugins/visualizations
>    >     whoosh_index_dir: /srv/galaxy/var/cache/whoosh_cache
>    > ```
>    {: .code-out.code-max-300}
>
{: .hands_on}

Galaxy is now configured with an admin user, a database, and a place to store data. Additionally we've immediately configured the mules for production Galaxy serving. So we're ready to set up systemd which will manage the Galaxy processes! Launching Galaxy by hand is not a good use of your time, so we will immediately switch to a process manager for that, [systemd](https://freedesktop.org/wiki/Software/systemd/).

> <hands-on-title>systemd</hands-on-title>
>
> 1. Configure systemd in `group_vars/galaxyservers.yml` file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -57,3 +57,6 @@ galaxy_config:
>             pools:
>               - job-handler
>               - workflow-scheduler
>    +
>    +# systemd
>    +galaxy_manage_systemd: true
>    {% endraw %}
>    ```
>    {: data-commit="Setup systemd variables"}
>
>    > <tip-title>Using Galaxy 21.01?</tip-title>
>    > For versions of Galaxy older than 22.05, using `gravity` was not the default. There you'll need to set these additional variables:
>    > {% raw %}
>    > ```diff
>    > +galaxy_manage_gravity: "{{ false if __galaxy_major_version is version('22.01', '<') else true }}"
>    > +galaxy_systemd_mode: "{{ 'mule' if __galaxy_major_version is version('22.01', '<') else 'gravity' }}"
>    > ```
>    > {% endraw %}
>    {: .tip}
>
> 2. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Log in and check the status with `sudo systemctl status galaxy`
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > sudo systemctl status galaxy
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    >
>    > ```ini
>    > ● galaxy.service - Galaxy
>    >     Loaded: loaded (/etc/systemd/system/galaxy.service; enabled; vendor preset: enabled)
>    >     Active: active (running) since Mon 2022-06-27 03:46:49 UTC; 11h ago
>    >   Main PID: 11778 (galaxyctl)
>    >      Tasks: 73 (limit: 9495)
>    >     Memory: 1.7G (limit: 16.0G)
>    >        CPU: 1h 11min 33.684s
>    >     CGroup: /system.slice/galaxy.service
>    >             ├─11778 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/galaxyctl start --foreground --quiet
>    >             ├─11791 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/supervisord -c /srv/galaxy/var/gravity/supervisor/supervisord.conf --nodaemon
>    >             ├─11801 /bin/tail -f /srv/galaxy/var/gravity/supervisor/supervisord.log
>    >             ├─11802 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/celery --app galaxy.celery worker --concurrency 2 --loglevel DEBUG --pool threads --queues celery,galaxy.internal,galaxy.external
>    >             ├─11827 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/celery --app galaxy.celery beat --loglevel DEBUG --schedule /srv/galaxy/var/gravity/celery-beat-schedule
>    >             ├─11840 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/gunicorn "galaxy.webapps.galaxy.fast_factory:factory()" --timeout 300 --pythonpath lib -k galaxy.webapps.galaxy.workers.Worker -b unix:/srv/>
>    >             ├─11891 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/gunicorn "galaxy.webapps.galaxy.fast_factory:factory()" --timeout 300 --pythonpath lib -k galaxy.webapps.galaxy.workers.Worker -b unix:/srv/>
>    >             ├─11898 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/gunicorn "galaxy.webapps.galaxy.fast_factory:factory()" --timeout 300 --pythonpath lib -k galaxy.webapps.galaxy.workers.Worker -b unix:/srv/>
>    >             ├─11929 python ./lib/galaxy/main.py -c /srv/galaxy/config/galaxy.yml --server-name=handler_0 --attach-to-pool=job-handler --attach-to-pool=workflow-scheduler --pid-file=/srv/galaxy/var/gravity/supe>
>    >             ├─11949 python ./lib/galaxy/main.py -c /srv/galaxy/config/galaxy.yml --server-name=handler_1 --attach-to-pool=job-handler --attach-to-pool=workflow-scheduler --pid-file=/srv/galaxy/var/gravity/supe>
>    >             └─11968 python ./lib/galaxy/main.py -c /srv/galaxy/config/galaxy.yml --server-name=handler_2 --attach-to-pool=job-handler --attach-to-pool=workflow-scheduler --pid-file=/srv/galaxy/var/gravity/supe>
>    >
>    > Jun 27 03:47:58 training-test galaxyctl[11791]: 2022-06-27 03:47:58,960 INFO success: handler_1 entered RUNNING state, process has stayed up for > than 20 seconds (startsecs)
>    > Jun 27  03:47:58 training-test galaxyctl[11801]: 2022-06-27 03:47:58,960 INFO success: handler_1 entered RUNNING state, process has stayed up for > than 20 seconds (startsecs)
>    > Jun 27 03:47:58 training-test galaxyctl[11791]: 2022-06-27 03:47:58,975 INFO waiting for handler_2 to stop
>    > Jun 27 03:47:58 training-test galaxyctl[11801]: 2022-06-27 03:47:58,975 INFO waiting for handler_2 to stop
>    > Jun 27 03:48:00 training-test galaxyctl[11791]: 2022-06-27 03:48:00,264 INFO stopped: handler_2 (exit status 0)
>    > ```
>    {: .code-out.code-max-300}
>
>    > <tip-title>Did this fail?</tip-title>
>    > ```
>    > ● galaxy.service - Galaxy
>    >      Loaded: loaded (/etc/systemd/system/galaxy.service; enabled; vendor preset: enabl>
>    >      Active: failed (Result: exit-code) since Mon 2021-06-28 17:07:58 CEST; 18min ago
>    >     Process: 521705 ExecStart=/srv/galaxy/venv/bin/uwsgi --yaml /srv/galaxy/config/gal>
>    >    Main PID: 521705 (code=exited, status=1/FAILURE)
>    >         CPU: 21ms
>    > Jun 28 17:07:58 gat-14.be.training.galaxyproject.eu systemd[1]: galaxy.service: Schedu>
>    > Jun 28 17:07:58 gat-14.be.training.galaxyproject.eu systemd[1]: Stopped Galaxy.
>    > Jun 28 17:07:58 gat-14.be.training.galaxyproject.eu systemd[1]: galaxy.service: Start >
>    > Jun 28 17:07:58 gat-14.be.training.galaxyproject.eu systemd[1]: galaxy.service: Failed>
>    > Jun 28 17:07:58 gat-14.be.training.galaxyproject.eu systemd[1]: Failed to start Galaxy.
>    > ```
>    >
>    > Check your /srv/galaxy/config/galaxy.yml and ensure that it lines up exactly with what you expect.
>    {: .tip}
>
> 6. Some things to note:
>
>    1. Although the playbook will restart Galaxy upon config changes, you will sometimes need to restart it by hand, which can be done with `systemctl restart galaxy`
>    2. You can use `journalctl -fu galaxy` to see the logs of Galaxy
>    3. You may have noticed after setting up Galaxy in the previous section, that Ansible printed the message "RESTARTER NOT IMPLEMENTED - Please restart Galaxy manually. You can define your own handler and enable it with `galaxy_restart_handler_name`." You should no longer see this message on future playbook runs, as enabling systemd has automatically set `galaxy_restart_handler_name` for us - the role now knows how to restart Galaxy!
>
{: .hands_on}

> <details-title>Ansible, failures, and notifications</details-title>
>
> Sometimes Ansible tasks will fail. Usually due to misconfiguration, but occasionally due to other issues like your coworker restarted the server while you were doing maintenance, or network failures, or any other possible error. It happens. An unfortunate side effect can be observed in specific situations:
>
> Let's say you're running a playbook that updates the `galaxy.yml`, which will in turn notify the handler `Restart Galaxy`. If this change is made, and notification triggered, but a failure occurs before Ansible can reach the step where it runs the handlers. The handlers will not run during this Ansible execution.
>
> The next time you run the playbook, Ansible will not observe any configuration files changing (because they were changed in the last run.) And so the `Restart Galaxy` handler will not run.
>
> If you encounter this situation you just have to be mindful of the fact, and remember to manually run the handler. There is no general solution to this problem unfortunately. This applies mostly to development setups. In production you're probably running that playbook somewhat regularly and do not expect failures as everything is quite stable.
>
{: .details}

> ```bash
> 0-galaxy-up.sh
> ```
> {: data-test="true"}
{: .hidden}

## NGINX

With this we have:

- PostgreSQL running
- Galaxy running (managed by systemd)

Although Gunicorn can server HTTP for us directly, by moving to NGINX (pronounced "engine X" /ˌɛndʒɪnˈɛks/ EN-jin-EKS) (or another reverse proxy), it can automatically compress selected content, and we can easily apply caching headers to specific types of content like CSS or images. It is also necessary if we want to serve multiple sites at once, e.g. with a group website at `/` and Galaxy at `/galaxy`. Lastly, it can provide authentication as well, as noted in the [External Authentication]({{ site.baseurl }}/topics/admin/tutorials/external-auth/tutorial.html) tutorial.

For this, we will use NGINX. It is possible to configure Galaxy with Apache and potentially other webservers but this is not the configuration that receives the most testing. We recommend NGINX unless you have a specific need for Apache.

> <hands-on-title>NGINX</hands-on-title>
>
> 1. Add the role `galaxyproject.nginx` to the end of your playbook and have it run as root.
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -18,3 +18,4 @@
>         - role: uchida.miniconda
>           become: true
>           become_user: "{{ galaxy_user.name }}"
>    +    - galaxyproject.nginx
>    {% endraw %}
>    ```
>    {: data-commit="Add nginx to playbook"}
>
> 3. We need to configure the virtualhost. This is a slightly more complex process as we have to write the proxying configuration ourselves. This may seem annoying, but it is often the case that sites have individual needs to cater to, and it is difficult to provide a truly generic webserver configuration. Additionally, we will enable secure communication via HTTPS using SSL/TLS certificates provided by [certbot](https://certbot.eff.org/).
>
>    Add the following to your group variables file:
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -60,3 +60,55 @@ galaxy_config:
>     
>     # systemd
>     galaxy_manage_systemd: true
>    +
>    +# Certbot
>    +certbot_auto_renew_hour: "{{ 23 |random(seed=inventory_hostname)  }}"
>    +certbot_auto_renew_minute: "{{ 59 |random(seed=inventory_hostname)  }}"
>    +certbot_auth_method: --webroot
>    +certbot_install_method: virtualenv
>    +certbot_auto_renew: yes
>    +certbot_auto_renew_user: root
>    +certbot_environment: staging
>    +certbot_well_known_root: /srv/nginx/_well-known_root
>    +certbot_share_key_users:
>    +  - nginx
>    +certbot_post_renewal: |
>    +    systemctl restart nginx || true
>    +certbot_domains:
>    + - "{{ inventory_hostname }}"
>    +certbot_agree_tos: --agree-tos
>    +
>    +# NGINX
>    +nginx_selinux_allow_local_connections: true
>    +nginx_servers:
>    +  - redirect-ssl
>    +nginx_ssl_servers:
>    +  - galaxy
>    +nginx_enable_default_server: false
>    +nginx_conf_http:
>    +  client_max_body_size: 1g
>    +  # gzip: "on" # This is enabled by default in Ubuntu, and the duplicate directive will cause a crash.
>    +  gzip_proxied: "any"
>    +  gzip_static: "on"   # The ngx_http_gzip_static_module module allows sending precompressed files with the ".gz" filename extension instead of regular files.
>    +  gzip_vary: "on"
>    +  gzip_min_length: 128
>    +  gzip_comp_level: 6  # Tradeoff of better compression for slightly more CPU time.
>    +  gzip_types: |
>    +      text/plain
>    +      text/css
>    +      text/xml
>    +      text/javascript
>    +      application/javascript
>    +      application/x-javascript
>    +      application/json
>    +      application/xml
>    +      application/xml+rss
>    +      application/xhtml+xml
>    +      application/x-font-ttf
>    +      application/x-font-opentype
>    +      image/png
>    +      image/svg+xml
>    +      image/x-icon
>    +nginx_ssl_role: usegalaxy_eu.certbot
>    +nginx_conf_ssl_certificate: /etc/ssl/certs/fullchain.pem
>    +nginx_conf_ssl_certificate_key: /etc/ssl/user/privkey-nginx.pem
>    {% endraw %}
>    ```
>    {: data-commit="Configure nginx, certbot"}
>
>    > <details-title>Certbot details</details-title>
>    >
>    > This is a lot of configuration but it is not very complex to understand. We'll go through it step by step:
>    >
>    > - `certbot_auto_renew_hour/minute`: Certbot certificates are short lived, they only last 90 days. As a consequence, automated renewal is a significant part of the setup and well integrated. The certbot role installs a cron job which checks if the certificate needs to be renewed (when it has <30 days of lifetime left) and attempts to renew the certificate as needed. In order to reduce load on the certbot servers, we randomly set the time when the request will be made, so not all of the requests occur simultaneously. For training VMs this will likely never be reached. For real-life machines, this is more important.
>    > - `certbot_auth_method`: [Multiple authentication methods](https://certbot.eff.org/docs/using.html) are supported, we will use the webroot method since that integrates nicely with `galaxyproject.nginx`. This writes out a file onto the webserver's root (that we specify in `certbot_well_known_root`) which certbot's servers will check.
>    > - `certbot_auto_renew`: Automatically attempt renewing the certificate as the `certbot_auto_renew_user`
>    > - `certbot_environment`: The options here are `production` and `staging`, we will set this to staging and obtain a verified but invalid certificate as browsers are intentionally not configured to trust the certbot staging certificates. The staging environment has higher [rate limits](https://letsencrypt.org/docs/rate-limits/) and allows requesting more certificates during trainings. If you are deploying on a production machine you should set this to `production`.
>    > - `certbot_share_key_users`: This variable automatically shares the certificates with any system users that might need to access them. Here just nginx needs access.
>    > - `certbot_post_renewal`: Often services need to be notified or restarted once the certificates have been updated.
>    > - `certbot_domains`: These are the domains that are requested for verification. Any entries you place here **must** all be publicly resolvable.
>    > - `certbot_agree_tos`: We automatically agree to the certbot TOS. You can read the current one on [their website](https://letsencrypt.org/repository/)
>    {: .details}
>
>    > <details-title>Nginx details</details-title>
>    >
>    > Likewise the nginx configuration has a couple of important points:
>    > - `nginx_selinux_allow_local_connections`: Specific to CentOS hosts where Nginx will need to access Galaxy
>    > - `nginx_enable_default_server/vhost`: Most Nginx packages come with a default configuration for the webserver. We do not want this.
>    > - `nginx_conf_http`: Here we can write any extra configuration we have, `client_max_body_size: 1g` increases the POST limit to 1Gb which makes uploads easier.
>    >
>    > These control the SSL configuration
>    > - `nginx_conf_ssl_certificate/key`: Location of the certificate / private key.
>    >
>    > The configuration variables we added in our group variables file has the following variables
>    > ```yaml
>    > nginx_servers:
>    >   - redirect-ssl
>    > nginx_ssl_servers:
>    >   - galaxy
>    > ```
>    >
>    > The `galaxyproject.galaxy` role expects to find two files with these names in `templates/nginx/redirect-ssl.j2` and `templates/nginx/galaxy.j2`
>    >
>    {: .details}
>
>    > <details-title>Running this tutorial with your own SSL certificates</details-title>
>    >
>    > If you want, you can run this with your own SSL certificates. You can read about all of the variables you need to set [in the ansible NGINX role](https://github.com/galaxyproject/ansible-nginx#ssl-configuration).
>    >
>    {: .details}
>
>    > <details-title>Running this tutorial <i>without</i> SSL</details-title>
>    >
>    > If you want, you can run this tutorial without SSL. We will provide a sketch of the configuration changes needed, but this is of course not recommended for production, so we will not go into detail here:
>    >
>    > Instead of the above step you should do:
>    >
>    > {% raw %}
>    > ```diff
>    > --- a/group_vars/galaxyservers.yml
>    > +++ b/group_vars/galaxyservers.yml
>    > @@ -83,12 +83,7 @@ certbot_agree_tos: --agree-tos
>    >  # NGINX
>    >  nginx_selinux_allow_local_connections: true
>    >  nginx_servers:
>    > -  - redirect-ssl
>    > -nginx_ssl_servers:
>    >    - galaxy
>    >  nginx_enable_default_server: false
>    >  nginx_conf_http:
>    >    client_max_body_size: 1g
>    > -nginx_ssl_role: usegalaxy_eu.certbot
>    > -nginx_conf_ssl_certificate: /etc/ssl/certs/fullchain.pem
>    > -nginx_conf_ssl_certificate_key: /etc/ssl/user/privkey-nginx.pem
>    > ```
>    > {% endraw %}
>    >
>    {: .details}
>
> 4. Create the directory `templates/nginx` (staying in galaxy directory, after which groups_vars, roles will be siblings of templates), where we will place our configuration files which should be templated out to the server.
>
>    Create the `templates/nginx/redirect-ssl.j2` with the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/nginx/redirect-ssl.j2
>    @@ -0,0 +1,14 @@
>    +server {
>    +    listen 80 default_server;
>    +    listen [::]:80 default_server;
>    +
>    +    server_name "{{ inventory_hostname }}";
>    +
>    +    location /.well-known/ {
>    +        root {{ certbot_well_known_root }};
>    +    }
>    +
>    +    location / {
>    +        return 302 https://$host$request_uri;
>    +    }
>    +}
>    {% endraw %}
>    ```
>    {: data-commit="Setup http to https redirect"}
>
>    This will redirect all requests to use HTTPS.
>
> 5. Create `templates/nginx/galaxy.j2` with the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/nginx/galaxy.j2
>    @@ -0,0 +1,61 @@
>    +upstream galaxy {
>    +    server unix:{{ galaxy_mutable_config_dir }}/gunicorn.sock;
>    +
>    +    # Or if you serve galaxy at a path like http(s)://fqdn/galaxy
>    +    # Remember to set galaxy_url_prefix in the galaxy.yml file.
>    +    # server unix:{{ galaxy_mutable_config_dir }}/gunicorn.sock:/galaxy;
>    +}
>    +
>    +server {
>    +    # Listen on port 443
>    +    listen        *:443 ssl default_server;
>    +    # The virtualhost is our domain name
>    +    server_name   "{{ inventory_hostname }}";
>    +
>    +    # Our log files will go here.
>    +    access_log  syslog:server=unix:/dev/log;
>    +    error_log   syslog:server=unix:/dev/log;
>    +
>    +    # The most important location block, by default all requests are sent to gunicorn
>    +    # If you serve galaxy at a path like /galaxy, change that below (and all other locations!)
>    +    location / {
>    +        # This is the backend to send the requests to.
>    +        proxy_pass http://galaxy;
>    +
>    +        proxy_set_header Host $http_host;
>    +        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
>    +        proxy_set_header X-Forwarded-Proto $scheme;
>    +        proxy_set_header Upgrade $http_upgrade;
>    +    }
>    +
>    +    # Static files can be more efficiently served by Nginx. Why send the
>    +    # request to Gunicorn which should be spending its time doing more useful
>    +    # things like serving Galaxy!
>    +    location /static {
>    +        alias {{ galaxy_server_dir }}/static;
>    +        expires 24h;
>    +    }
>    +
>    +    # In Galaxy instances started with run.sh, many config files are
>    +    # automatically copied around. The welcome page is one of them. In
>    +    # production, this step is skipped, so we will manually alias that.
>    +    location /static/welcome.html {
>    +        alias {{ galaxy_server_dir }}/static/welcome.html.sample;
>    +        expires 24h;
>    +    }
>    +
>    +    # serve visualization and interactive environment plugin static content
>    +    location ~ ^/plugins/(?<plug_type>[^/]+?)/((?<vis_d>[^/_]*)_?)?(?<vis_name>[^/]*?)/static/(?<static_file>.*?)$ {
>    +        alias {{ galaxy_server_dir }}/config/plugins/$plug_type/;
>    +        try_files $vis_d/${vis_d}_${vis_name}/static/$static_file
>    +                  $vis_d/static/$static_file =404;
>    +    }
>    +
>    +    location /robots.txt {
>    +        alias {{ galaxy_server_dir }}/static/robots.txt;
>    +    }
>    +
>    +    location /favicon.ico {
>    +        alias {{ galaxy_server_dir }}/static/favicon.ico;
>    +    }
>    +}
>    {% endraw %}
>    ```
>    {: data-commit="Configure nginx vhost for galaxy"}
>
>    You'll notice that we have set a 24 hour cache timeout on static content served by NGINX. [Google's PageSpeed Tools](https://developers.google.com/speed/pagespeed/insights/) can identify any additional compression or caching improvements you can make.
>
>    > <details-title>Running this tutorial <i>without</i> SSL</details-title>
>    >
>    > In your `galaxy.j2` in the above step, you should change the `listen` parameter:
>    >
>    > {% raw %}
>    > ```diff
>    > --- a/templates/nginx/galaxy.j2
>    > +++ b/templates/nginx/galaxy.j2
>    > -listen        *:443 ssl default_server;
>    > +listen        *:80 default_server;
>    > ```
>    > {% endraw %}
>    >
>    {: .details}
>
>    {% snippet topics/admin/tutorials/ansible-galaxy/faqs/custom-welcome.md %}
>
> 6. Run the playbook. At the very end, you should see output like the following indicating that Galaxy has been restarted:
>
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > ansible-playbook galaxy.yml
>    > > ```
>    > > {: data-cmd="true"}
>    > {: .code-in}
>    > > <code-out-title>Output</code-out-title>
>    > > ```
>    > > ...
>    > > RUNNING HANDLER [restart galaxy] ****************************************
>    > > changed: [gat-0.training.galaxyproject.eu]
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
>    If you didn't, you might have missed the first step in this hands-on.
>
>    > <tip-title>Install nginx - Failed [emerg] socket() [::]:80 failed</tip-title>
>    > If you see this message:
>    >
>    > ```
>    > nginx: [emerg] socket() [::]:80 failed (97: Address family not supported by protocol)
>    > ```
>    >
>    > This may be due to missing IPv6 support on your machine. It is extremely unlikely your operating system doesn't support it, instead it was probably disabled. Check `/etc/sysctl.conf` and `/etc/sysctl.d/*` for `net.ipv6.conf.*.disable_ipv6=1`, or `/etc/default/grub` for `ipv6.disable=1`
>    {: .tip}
>
> 7. Check out the changes made to your server in `/etc/nginx/sites-enabled/`, particularly the directory containing the Galaxy virtualhost.
>
> 8. Check out the nginx logs with `journalctl`
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > journalctl -fu nginx
>    > ```
>    {: .code-in}
{: .hands_on}

> <details-title>"Potential Security Risk" / LetsEncrypt Staging Environment</details-title>
>
> LetsEncrypt has rate limits on requesting trusted certificates to prevent abuse of their service.
> In a training setting there is no need to request certificates that will be trusted by all browsers. So we will request a testing certificate to show how it works, and by changing `staging` to `production`, you can request browser trusted certificates.
>
> You will probably see an error like this, when trying to access your Galaxy:
>
> ![Browser warning for invalid certificate](../../images/ssl-warning0.png "A browser warning for an invalid certificate. But because we requested a staging certificate, we expected this.")
>
> If you view the details of the certificate, you can see that it is trusted, but by the Fake LE Intermediate, which browsers do not trust.
>
> ![Certificate information](../../images/ssl-warning1.png "Investigating the certificate a little, we can see that it was signed, just untrusted.")
>
> Clicking through the warnings (with full knowledge of why) we will see our secured Galaxy:
>
> ![The finally working Galaxy](../../images/working-galaxy.png "Galaxy is alive!")
>
{: .details}

> <comment-title>Role Dependencies</comment-title>
>
> Throughout the playbook we added roles in a specific order. Partially this was to mimic the original training and build up a working Galaxy server from nothing, but partially this is also because of some hidden role dependencies on each other. Some must run before others, in order to set certain variables. Looking at the dependencies in detail:
>
>  Role                         | Role-Role Dependencies
>  ----                         | ------------
>  `galaxyproject.postgresql`   | None
>  `galaxyproject.postgresql_objects` | None
>  `geerlingguy.pip`            | None
>  `galaxyproject.galaxy`       | None
>  `uchida.miniconda`           | In our group variables, we define the path to {% raw %}`{{ galaxy_tool_dependency_dir }}/_conda`{% endraw %}, so Galaxy needs to have set those variables
>  `galaxyproject.nginx`        | This requires Galaxy variables to find the static assets.
{: .comment}

## Log in to Galaxy

Now that your production-ready Galaxy is running, try registering a new user and logging in!

In order to be the administrator user, you will need to register an account with the same email address you used in the group variables under the `admin_users` setting.

## Job Configuration

One of the most important configuration files for a large Galaxy server is the `job_conf.yml` file. This file tells Galaxy where to run all of the jobs that users execute. If Galaxy can't find a job conf file or none has been specified in the `galaxy.yml` file, it will use a default configuration, [`job_conf.sample_advanced.yml`](https://github.com/galaxyproject/galaxy/blob/fd8826904be3f3201f7e7a3b6ffbadd27f268896/test/unit/app/jobs/job_conf.sample_advanced.yml).

The job configuration file allows Galaxy to run jobs in multiple locations using a variety of different mechanisms. Some of these mechanisms include:

* Local - Galaxy runs jobs on the same computer that Galaxy itself is running on.
* DRMAA - Galaxy can connect to a cluster and submit jobs via a distributed resource manager such as Slurm, Condor, PBS Torque or Sun Grid Engine.
* Pulsar - Galaxy can also send jobs to remote compute hosts over the internet using Pulsar.

### The job conf file - basics

The `job_conf.yml` file has three basic sections:

* **runners** - This section lists the types of job management systems that this Galaxy server is configured to use, and tells Galaxy to load the drivers for each type.
* **execution** - This section lists the different locations, queues, etc. that Galaxy can send jobs to. Each one has a name and uses a *plugin* to communicate with that location. They can specify things like the number of CPUs, amount of RAM to be allocated, etc. for DRMAA locations. Usually, one of the destinations is set to be the default.
* **tools** - This section lists the various tools that you would like to send to a non-default *destination*. Each line in this section pairs up a tool in Galaxy with a particular job *destination*. Every time Galaxy gets a job for that particular tool, it is always sent to that *destination*.

The basic `job_conf.yml` file looks like this:

{% raw %}
```yml
runners:
  local_runner:
    load: galaxy.jobs.runners.local:LocalJobRunner
    workers: 4
execution:
  default: local_dest
  environments:
    local_dest:
      runner: local_runner
tools:
- id: bwa
  destination: local_dest
```
{% endraw %}

The above `job_conf.yml` file defines a *runner* and an *execution* to allow Galaxy to run user jobs on the local computer (i.e. The computer that Galaxy is running on.)

Firstly, the plugins section contains a plugin called "local" which is of type "runner" and then loads the python code module for supporting local jobs. Next the destinations section contains a destination called "local" using the runner "local". As this is the only destination specified, it is also the default. So now everytime a user clicks "Execute" on a tool form, Galaxy will run the corresponding job locally using the python code specified.

> <tip-title>Want to use something else?</tip-title>
> There are a lot of other plugins available for Galaxy for using other resources such as docker containers, kubernetes hosts, Pulsar destinations and HPC clusters to name a few. See the Galaxy documentation on [job configuration](https://docs.galaxyproject.org/en/master/admin/jobs.html) for more details on these plugins and their configuration.
{: .tip}

> <hands-on-title>Job Conf</hands-on-title>
>
> 1. If the folder does not exist, create `templates/galaxy/config` next to your `galaxy.yml` playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > mkdir -p templates/galaxy/config/
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 2. Create `templates/galaxy/config/job_conf.yml.j2` with the following contents (note that we have changed the names of the plugin and destination from the basic sample file to provide a bit more clarity):
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/galaxy/config/job_conf.yml.j2
>    @@ -0,0 +1,14 @@
>    +runners:
>    +  local_runner:
>    +    load: galaxy.jobs.runners.local:LocalJobRunner
>    +    workers: 4
>    +
>    +execution:
>    +  default: local_dest
>    +  environments:
>    +    local_dest:
>    +      runner: local_runner
>    +
>    +tools:
>    +- class: local # these special tools that aren't parameterized for remote execution - expression tools, upload, etc
>    +  environment: local_dest
>    {% endraw %}
>    ```
>    {: data-commit="Add job conf"}
>
>    > <tip-title>workers=4</tip-title>
>    > In the local runner, `workers: 4` means "number of jobs that can be running at one time". For every other job runner, it means the number of threads that are created to start/manage/finish jobs. E.g. if you are in a class and 50 people submit jobs, then there are four threads that can handle these jobs at once. But additional job handlers can be more useful as well.
>    {: .tip}
>
> 3. Inform the `galaxyproject.galaxy` role of where you would like the `job_conf.yml` to reside, by setting it in your `group_vars/galaxyservers.yml`:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -37,6 +37,7 @@ galaxy_config:
>         tool_data_path: "{{ galaxy_mutable_data_dir }}/tool-data"
>         object_store_store_by: uuid
>         id_secret: "{{ vault_id_secret }}"
>    +    job_config_file: "{{ galaxy_config_dir }}/job_conf.yml"
>       gravity:
>         galaxy_root: "{{ galaxy_root }}/server"
>         app_server: gunicorn
>    {% endraw %}
>    ```
>    {: data-commit="Load job conf in galaxy config"}
>
>    And then deploy the new config file using the `galaxy_config_templates` var (also from the `galaxyproject.galaxy` role) in your group vars:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -59,6 +59,10 @@ galaxy_config:
>               - job-handler
>               - workflow-scheduler
>     
>    +galaxy_config_templates:
>    +  - src: templates/galaxy/config/job_conf.yml.j2
>    +    dest: "{{ galaxy_config.galaxy.job_config_file }}"
>    +
>     # systemd
>     galaxy_manage_systemd: true
>     
>    {% endraw %}
>    ```
>    {: data-commit="Deploy job conf to config dir"}
>
> 4. Run the playbook. At the very end, you should see output like the following indicating that Galaxy has been restarted:
>
>    > > <code-in-title>Bash</code-in-title>
>    > > ```bash
>    > > ansible-playbook galaxy.yml
>    > > ```
>    > > {: data-cmd="true"}
>    > {: .code-in}
>    > > <code-out-title>Output</code-out-title>
>    > > ```
>    > > ...
>    > > RUNNING HANDLER [restart galaxy] ****************************************
>    > > changed: [gat-0.training.galaxyproject.eu]
>    > > ```
>    > {: .code-out}
>    {: .code-2col}
>
> 5. Checkout the new job_conf.yml file.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > cat /srv/galaxy/config/job_conf.yml
>    > ```
>    {: .code-in}
>
>    > <code-out-title>Bash</code-out-title>
>    > ```yml
>    > runners:
>    >   local_runner:
>    >     load: galaxy.jobs.runners.local:LocalJobRunner
>    >     workers: 4
>    > execution:
>    >   default: local_dest
>    >   environments:
>    >     local_dest:
>    >       runner: local_runner
>    > ```
>    {: .code-out.code-max-300}
>
{: .hands_on}

## Productionising Galaxy

This is a fantastic base Galaxy installation but there are numerous additional options we should be setting for a real production Galaxy.

> <hands-on-title>More Features!</hands-on-title>
>
> 1. Make the following changes to your `group_vars/galaxyservers.yml`:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -38,6 +38,28 @@ galaxy_config:
>         object_store_store_by: uuid
>         id_secret: "{{ vault_id_secret }}"
>         job_config_file: "{{ galaxy_config_dir }}/job_conf.yml"
>    +    # SQL Performance
>    +    database_engine_option_server_side_cursors: true
>    +    slow_query_log_threshold: 5
>    +    enable_per_request_sql_debugging: true
>    +    # File serving Performance
>    +    nginx_x_accel_redirect_base: /_x_accel_redirect
>    +    # Automation / Ease of Use / User-facing features
>    +    watch_job_rules: 'auto'
>    +    allow_path_paste: true
>    +    enable_quotas: true
>    +    allow_user_deletion: true
>    +    show_welcome_with_login: true
>    +    expose_user_name: true
>    +    expose_dataset_path: true
>    +    expose_potentially_sensitive_job_metrics: true
>    +    # NFS workarounds
>    +    retry_job_output_collection: 3
>    +    # Debugging
>    +    cleanup_job: onsuccess
>    +    allow_user_impersonation: true
>    +    # Tool security
>    +    outputs_to_working_directory: true
>       gravity:
>         galaxy_root: "{{ galaxy_root }}/server"
>         app_server: gunicorn
>    {% endraw %}
>    ```
>    {: data-commit="Add production facing vars"}
>
>    > <tip-title>What do these do?</tip-title>
>    > Check out the full details in the [Galaxy documentation](https://docs.galaxyproject.org/en/master/admin/config.html#configuration-options), but we'll discuss a couple briefly:
>    > - `nginx_x_accel_redirect_base`: This is required to have NGINX serve user files. You don't want Galaxy to waste time reading a 100GB fastq file a user has asked for, so you offload that to NGINX. The request is passed through to Galaxy, so permissions checks still occur, but Galaxy instead replies to NGINX just the path to the file that it should send to the requesting user.
>    > - `enable_quotas`: You definitely want to set a default quota for your users!
>    > - `expose_user_name`: This exposes usernames in the history and dataset sharing forms which makes life easier for your users.
>    > - `expose_dataset_path`: This shares the entire dataset path and enable the user to have access to the command line that the wrapper launched.
>    > - `cleanup_job`: This prevents a job from being cleaned up if it failed, allowing you to debug and see what went wrong.
>    > - `allow_user_impersonation`: Users submit horrible bug reports (often screenshots of unrelated things), impersonation means you can just see their history and look at what went wrong.
>    {: .tip}
>
> 2. Make the following changes to your `templates/nginx/galaxy.j2`:
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -58,4 +58,14 @@ server {
>         location /favicon.ico {
>             alias {{ galaxy_server_dir }}/static/favicon.ico;
>         }
>    +
>    +    location /_x_accel_redirect {
>    +        internal;
>    +        alias /;
>    +    }
>    +
>    +    # Support click-to-run in the GTN-in-Galaxy Webhook
>    +    location /training-material/ {
>    +        proxy_pass https://training.galaxyproject.org/training-material/;
>    +    }
>     }
>    {% endraw %}
>    ```
>    {: data-commit="Add nginx x-accel-redir and g-i-g webhook config to nginx"}
>
>    > <tip-title>What do these do?</tip-title>
>    > The `_x_accel_redirect` is required for the NGINX file serving discussed above. For information on the GTN-in-Galaxy Webhook, see the [Galaxy Documentation](https://docs.galaxyproject.org/en/master/admin/special_topics/gtn.html?highlight=gtn%20galaxy). It's a very cool feature which helps your users access training materials directly in Galaxy.
>    {: .tip}
>
> 3. Run the playbook
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
> 1-galaxy-up.sh
> ```
> {: data-test="true"}
{: .hidden}

## Disaster Strikes! (Optional)

Because you're an admin, you need to be prepared for any situation, including the worst case scenarios. So we're going to simulate a disaster and show you how you can recover from it. It'll be fun!

For this "disaster", we will pretend that:

1. Your database is on another machine
2. Your datasets are on an NFS server or some other remote machine.

> <hands-on-title>Summon the Apocalypse</hands-on-title>
>
> So let's have a "worst case scenario", where your Galaxy server gets destroyed
>
> 1. Log on to your machine.
> 2. **Carefully**, as root, `rm -rf /srv/galaxy`, completely wipe out your Galaxy home directory.
>
{: .hands_on}

Your entire Galaxy server is gone! 😱

You were a responsible admin and had your user data and database stored on a separate system (and backed up), so at least those survived. Nevertheless, this is when most of us start feeling really bad; bosses start yelling, we start crying or reaching for bad coping habits.

But not you! You spent the day writing this Ansible playbook that describes your environment completely; all of the software that was installed, all of the configuration changes you have made. It leverages many community maintained roles and can be used to completely rebuild the server! With minimal effort on your part.

> <hands-on-title>Revert the Apocalypse</hands-on-title>
>
> 1. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 2. Get a hot beverage of your choice while it rebuilds.
>
> 3. Congratulate yourself for saving the day.
>
{: .hands_on}

And with that, Galaxy should be up and running again. If you log in, you should see the results of any jobs you ran earlier, you should still be able to log in with your old account, and everything should just work.

Ansible can save you from some really bad scenarios, *if and only if*:

- You can replace the hardware or find somewhere new to re-deploy
- You've made absolutely certain that every change made to a system is recorded within your playbooks and roles (i.e. no manual package installation)

Then you can potentially use it to recover.

> <comment-title>We have experience</comment-title>
>
> We can tell you this, we can repeat it over and over, but unless you really have a disaster happen to you, it is hard to appreciate how important it is that machines are completely controlled in terms of configuration and software deployment.
>
> We've experienced these incidents and we know how horribly stressful it can be if an important service like Galaxy goes down and you cannot immediately replace it with another instance. We hope you will immediately apply the lessons from this training material, it can potentially save you a lot of stress and worry.
>
{: .comment}

# Production & Maintenance

The time required to maintain a production Galaxy instance depends on the number of users and their specific needs, but a smallish server (<= 25 users) will typically require a day or two per month of maintenance. Large public servers like usegalaxy.org and usegalaxy.eu are largely full time jobs (although even their admins do find time to do other things).

## Keeping Galaxy Updated

If you have set your `galaxy_commit_id` group variable to a branch name like `release_22.05`, then all you need to do to keep Galaxy up to date (e.g. for security and bug fixes) is to run the playbook regularly. The `git` module in Ansible checks if you are on the latest commit of a given branch, and will update the clone of the repository if it is not.

## Upgrading Galaxy (Optional)

With Ansible, upgrading Galaxy to a new release is incredibly easy. Here is a commit from UseGalaxy.eu's upgrade:

```diff
--- a/group_vars/galaxyservers.yml
+++ b/group_vars/galaxyservers.yml
@@ -345,7 +345,7 @@ galaxy_instance_hostname: usegalaxy.eu
 
 galaxy_repo: 'https://github.com/usegalaxy-eu/galaxy.git'
-galaxy_commit_id: 'release_19.05'
+galaxy_commit_id: 'release_19.09'
 galaxy_force_checkout: true # discard any modified files
```

This is all that is required, changing the `galaxy_commit_id` and re-running the playbook. Everything else (building the client, database migrations, etc.) are taken care of for you.

It is recommend that you also do the following during an upgrade:

- [**Read the release notes**](https://docs.galaxyproject.org/en/latest/releases/index.html) for the new version of Galaxy
- Check out the latest [`galaxy.yml.sample`](https://github.com/galaxyproject/galaxy/blob/master/config/galaxy.yml.sample) for the new release, and see if there are any variables you want to set for your server
- Compare the [other configuration files](https://github.com/galaxyproject/galaxy/tree/master/config) to see if there are new features you want to take advantage of (e.g. new job runner options or metrics you wish to capture.)

When you've read the documentation and checked out the new features, you can then run the playbook and you're done!

## User Support

There are many user support resources available to you and your users online. [help.galaxyproject.org](https://help.galaxyproject.org) is the primary landing point for Galaxy users. For helping individual users, you might find it useful to impersonate them if they fail to send adequate bug reports.

### Impersonating

You can activate user impersonation with:

```
allow_user_impersonation: true
```

It is recommended to ask permission or consent before impersonating a user: "Hey, mind if we look at your history?"
Also, since Galaxy is not normally a real-time collaborative activity, you should avoid running jobs in their history, which can confuse users when datasets show up unexpectedly.
Additionally you can automatically send failing job error reports, even if users do not submit one, and maybe proactively address those issues (depending on the number of your users).

## Running on a cluster

If you need to run on a cluster with a shared file system, you will need to expose several directories to your cluster:

- `galaxy_shed_tools_dir`
- `galaxy_tool_dependency_dir`
- `galaxy_file_path`
- `galaxy_job_working_directory`
- `galaxy_server_dir`
- `galaxy_venv_dir`

But this isn't strictly true, not every one of these directories needs to be exported over the network. Instead, there are a couple of those directories that can be recreated in different locations (e.g. locally to nodes) using playbook tasks or simply by copying them to the relevant locations. Sometimes admins choose to do this for performance reasons when NFS can be slow or simply unnecessary:

- `galaxy_server_dir`
- `galaxy_venv_dir`

Both of those directories can be re-created by running parts of the `galaxyproject.galaxy` role on different machines. As neither of these directories is written to during system operations, it is fine for them to be copies living on different machines, rather coming from a single source of truth like an NFS server. For the other directories (e.g. job working directory or tools directory) these need to be written to, and read from, simultaneously from different machines, and should be consistent across the network.

Most of us use NFS, those who are using something more exotic (ceph, gluster, etc) have some reason for that like "my uni provided it" or "we really wanted to try something shiny". But NFS in most cases is decent and well tested and can be used. For larger deployments, a single NFS node may be insufficient, but at these scales it is common that your university or organisation provides some managed NFS service, e.g. Isilon which has load balancing built into its NFS service.

## Other software

But what about your other software, things that are deployed along with Galaxy? Things without an ansible role or are quite weird and require "manual tricks" to deploy?

You can write roles for that! Sometimes they are really ugly roles, but it at least keeps it documented + in place. E.g. UseGalaxy.eu has a custom role for rewriting users and it’s ugly and untested and should not be used by anyone else in case it breaks their site. But it's one of these manual tricks or bits of glue code, but we can encapsulate it as ansible. You can include tarballs in your role to be deployed and so on.

It may seem daunting to use ansible, but you don't have to do everything in ansible! You can just do a little bit, for managing just Galaxy, and manage the rest of your stack separately. Whatever fits best for your deployment.

# Loving Ansible? Convert your own servers!

A common question we get is:

> Is it possible to transform a Galaxy server which was installed by hand, into an ansible-managed one?
{: .quote}

Because of the great variance between sites and deployments it is not really possible to produce a guide for doing so. However, if you back your current deployment up and set all of the variables that control paths appropriately, it can be done!

It's very possible to do this piecemeal, taking one component at a time to ansibilise. nginx is a very easy place to start, the database is as well. Galaxy can be rebuilt in a new location, many of us ran a hand-managed galaxy setup somewhere like `/home/janedoe/work/projects/galaxy` that suddenly became a permanent project, this is a good chance to keep that but rebuild an identical one in `/srv/galaxy` or similar, and then switch over the traffic to that new, production Galaxy. Relocating data is a bit tougher and can require rewriting entries in the database.

# Final Notes

If you've been following along you should have a production-ready Galaxy, secured, everything ready to go.

{% snippet topics/admin/faqs/missed-something.md step=1 %}

![diagram of galaxy setup with postgres, galaxy, nginx, gunicorn, storage, and compute shown attached to each other.](../../images/ansible-galaxy-intro-2.png "This approximately represents what you have setup today, NGINX proxies requests to Gunicorn which runs the Galaxy App. That talks to Postgres, and Storage, and the Compute handlers which interact with the storage. As we setup mules, the compute handlers are the same processes which respond to web requests for nginx.")
