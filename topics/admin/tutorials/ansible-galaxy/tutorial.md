---
layout: tutorial_hands_on

title: "Galaxy Installation with Ansible"
questions:
- How does the Galaxy Ansible module work internally?
- How can I install a Galaxy server with Ansible
objectives:
- Have an understanding of how Galaxy's Ansible roles are structured and interact with one another
- Be able to use an Ansible playbook to install different flavors of Galaxy for different purposes
time_estimation: "2h"
key_points:
- Basic deployment with Ansible is surprisingly easy
- Complexity can grow over time as your organisation does, no need to start with playbooks like UseGalaxy.org
contributors:
  - erasche
  - natefoo
  - slugger70
  - shiltemann
tags:
  - ansible
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
  - title: "A server/VM to deploy Galaxy on"
    type: "none"
  - type: "external"
    title: Ansible setup on your local machine
    link: "https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html"
  - title: Comfort with Ansible roles and tasks
    type: "none"
---

# Overview
{:.no_toc}

This tutorial assumes you have some familiarity with [Ansible](https://www.ansible.com/resources/get-started) and are comfortable with writing and running playbooks. Here we'll see how to install a Galaxy server using an Ansible playbook. Galaxy Project has decided on Ansible for all of its deployment recipes. For our project, Ansible is even more fitting due to its name:

> An ansible is a category of fictional device or technology capable of instantaneous or faster-than-light communication. It can send and receive messages to and from a corresponding device over any distance or obstacle whatsoever with no delay, even between star systems ([wikipedia](https://en.wikipedia.org/wiki/Ansible))

We want to give you a comprehensive understanding of how the Galaxy installation occurs, but we want to avoid you having to write a "custom" Galaxy installation playbook which you would eventually throw away, in order to use the official playbooks. Given these goals, we will go through the playbook in depth first, and then move to a hands-on portion later. If you are not interested in the inner workings, you can [skip to that section now](#installing-galaxy).

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Playbook Overview


## Configuration

We'll be using the [official Galaxy role](https://github.com/galaxyproject/ansible-galaxy) to install and manage Galaxy. This role is found in [Ansible Galaxy](https://galaxy.ansible.com/) (no relation - it is Ansible's ) as [galaxyproject.galaxy](https://galaxy.ansible.com/galaxyproject/galaxy).

The official role is extremely configurable, everything that you want to change is exposed as a variable, and then tasks will change behaviour based on that. The [role documentation](https://github.com/galaxyproject/ansible-galaxy#role-variables) is the most up-to-date source of documentation for the variables. You should take a minute and read over the variables listed there.

The important variables for this tutorial are:

- `galaxy_root`
- `galaxy_server_dir`
- `galaxy_commit_id`
- `galaxy_config`

These are largely self explanatory: a directory for all of Galaxy's code and configuration, which commit should be installed, and the Galaxy configuration. We will not explain Galaxy configuration variables in detail as they are covered sufficiently in the `galaxy.yml` sample file or the [online documentation](https://docs.galaxyproject.org/en/master/admin/config.html#configuration-options).

The official recommendation is that you should have a variables file such as a `group_vars/galaxy.yml` for storing all of the Galaxy configuration.

## Tasks

As with every role, the entry point for execution is the `tasks/main.yml` file. For the [ansible-galaxy](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/main.yml) file, this includes a few groups of important tasks:

- [Clone (or Download) Galaxy](#cloning-galaxy)
- [Managing Configuration](#managing-configuration)
- [Fetching Dependencies](#dependencies)
- [Managing Mutable Setup](#mutable-setup)
- [Managing the Database](#managing-the-database)

### Cloning Galaxy

The [clone](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/clone.yml) task is the one which is primarily interesting to us, it downloads Galaxy, using git, to a specific commit.

1. Ansible tries to update Galaxy, cloning it if it is missing, or otherwise attempting to update to the correct commit (or latest commit of that branch.)
2. Any change is reported.
3. The virtualenv is set up:
    1. An empty virtualenv is created.
    2. Pip is updated within the virtualenv.
4. Any `.pyc` files are removed, as this can occasionally result in Python loading the cached code, even if the corresponding `.py` file is no more present at the checked-out commit. For safety, all of these are removed.

With that Galaxy is cloned to disk and is ready to be configured by the next task.

### Managing Configuration

The [static configuration setup](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/static_setup.yml) is relatively straightforward:

1. The directories for Galaxy configuration data and for the shed tools are created
2. Any config files are copied over
3. Any templates are copied over
4. The `galaxy.yml` (or `.ini`) is deployed

The setup for deploying templates and configuration files is a little bit non-standard by Ansible standards. Here you are expected to provide your own templates and static config files, and then describe them as a list of files and where they should be deployed to.

Using the [UseGalaxy.eu](https://github.com/usegalaxy-eu/infrastructure-playbook/blob/02ca578211bfee45044facf36635d28208e5dbb3/group_vars/galaxy.yml#L578) configuration as an example, we have something like:


{% raw %}
```yaml
galaxy_config_files:
  - src: files/galaxy/config/builds.txt
    dest: "{{ galaxy_config['galaxy']['builds_file_path'] }}"
  - src: files/galaxy/config/data_manager_conf.xml
    dest: "{{ galaxy_config['galaxy']['data_manager_config_file'] }}"
  - src: files/galaxy/config/datatypes_conf.xml
    dest: "{{ galaxy_config['galaxy']['datatypes_config_file'] }}"
  - src: files/galaxy/config/dependency_resolvers_conf.xml
    dest: "{{ galaxy_config['galaxy']['dependency_resolvers_config_file'] }}"
  - src: files/galaxy/config/disposable_email_blacklist.conf
    dest: "{{ galaxy_config['galaxy']['blacklist_file'] }}"
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

So the references in `galaxy_config_file` to `galaxy_config` are done to ensure that the setting for e.g. "location of the blacklist file" is the same between where we have configured Galaxy to looking for it, and where the file has been deployed, without requiring us to make variables changes in numerous places.

### Dependencies

Now that Galaxy is available on disk, Ansible is ready to start processing [dependencies](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/dependencies.yml) of Galaxy.

1. The virtualenv is updated with data from the `galaxy_requirements_file`, by default pointing to the requirements file in the codebase: {% raw %}`{{ galaxy_server_dir  }}/lib/galaxy/dependencies/pinned-requirements.txt`.{% endraw %}
2. Any necessary conditional dependencies of Galaxy are [collected by processing the config file](https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy/dependencies/__init__.py)
2. and then installed to the virtualenv.

### Mutable Setup

[This task](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/mutable_setup.yml) creates a directory and deploys any hand-managed mutable configuration files. It is unlikely that you want to manage these, as Galaxy does a sufficient job. Any changes you make to Galaxy like installing tools would result in the tools being "forgotten about", if you re-ran the playbook and overwrote that file.


### Managing the Database

The [database management tasks](https://github.com/galaxyproject/ansible-galaxy/blob/master/tasks/database.yml) are extremely convenient; any time you run the playbook to update Galaxy, this will automatically run the database schema migration as needed.

1. Galaxy first obtains the current DB version and the maximum possible DB version based on the codebase.
2. If needed, the database is created.
3. Both numbers are reported for the runner of the playbook.
4. If the numbers are different, then Ansible runs the command to upgrade the database to the latest version.

As an administrator who often forgot to run the upgrade, and would only notice it once Galaxy crashed during startup, having this process completely automated is extremely nice.

## Handlers

A number of the tasks that are executed will trigger a restart of Galaxy. Currently there is no auto-magic implementation of this, and you will have to do something that fits for your setup. The role provides a way to reference your own [handler](https://docs.ansible.com/ansible/latest/user_guide/playbooks_intro.html#handlers-running-operations-on-change), which we will do in this exercise. As Galaxy continues to standardise on setup, something will be implemented directly in the role to automatically restart the correct processes.

## Defaults

As with other roles, numerous [default values](https://github.com/galaxyproject/ansible-galaxy/blob/master/defaults/main.yml) are provided, but these are useful mostly as reference, and not to go through individually.

## Summary

Installation of Galaxy with the playbook follows generally the steps you would expect:

- Galaxy is cloned (or updated)
- A virtualenv is created if it doesn't exist
- Configuration files are installed
- Any missing dependencies are installed
- Any database updates are applied

It would not be difficult to write a role that does this yourself, but by using
the galaxyproject.galaxy role, you know that you're getting all of the Galaxy
best practices and knowledge from previous admins codified for you.

# Installing Galaxy

With the necessary background in place, you're ready to install Galaxy with Ansible. The playbooks will start simple, and grow over time. We will start with the minimal Galaxy playbook which only requires setting the `galaxy_server_dir` and expand from there. First, however, we need a database for Galaxy to connect to, so we will do that now.

To proceed from here it is expected that:

- You have [Ansible installed](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html) on your local machine

  > ### {% icon tip %} Tip: Running Ansible on your remote machine
  > It is possible to have Ansible installed on the remote machine and run it there as well. You will need to update your hosts file to point to localhost, and pass the `-c local` parameter. Be **certain** that the playbook that you're building is stored somewhere safe like your user home directory. We will remove data at one point during this tutorial and would not want you to lose your progress.
  {: .tip}

- You have a [hosts file](../ansible/tutorial.html#hosts-file) with the VM or host specified where you will deploy galaxy. We will refer to this group of hosts as "galaxyservers." You can use a different name if you prefer or are working on an existing playbook, just be sure to update all references later on.

## Requirements

We have codified all of the dependencies you will need into a yaml file that `ansible-galaxy` can install


> ### {% icon hands_on %} Hands-on: Minimal Galaxy Playbook
>
> 1. If you haven't done so already, create a new directory for your playbook, `playbook`, and `cd` to that directory
>
> 2. Create a new file in your working directory called `requirements.yml` and include the following contents:
>
>    ```yaml
>    - galaxyproject.repos
>    - galaxyproject.galaxy
>    - geerlingguy.nginx
>    - natefoo.postgresql_objects
>    - galaxyproject.postgresql
>    - galaxyproject.proftpd
>    - geerlingguy.pip
>    - src: https://github.com/usegalaxy-eu/ansible-role-supervisor
>      name: usegalaxy-eu.supervisor
>    - src: https://github.com/usegalaxy-eu/ansible-certbot
>      name: usegalaxy-eu.certbot
>    ```
>
> 3. In the same directory, run `ansible-galaxy install -p roles -r requirements.yml`. This will install all of the required modules for this training into the `roles/` folder. We choose to install to a folder to give you easy access to look through the different roles when you have questions on their behaviour.
>
> 4. Create the hosts file if you have not done so, include a group for `[galaxyservers]` with the address of the host where you will install Galaxy.
>
> 5. Inspect the contents of the newly created `roles` directory in your working directory.
>
{: .hands_on}

## PostgreSQL

Galaxy is capable of talking to multiple databases through SQLAlchemy drivers. SQLite is the development database, but PostgreSQL is recommended in production. MySQL is a possibility, but does not receive the same testing or bugfixes from the main development team as PostgreSQL, so we will only show installation with PostgreSQL.

PostgreSQL maintains its own user database apart from the system user database. By default, PostgreSQL uses the "peer" authentication method which allows access for system users with matching PostgreSQL usernames (other authentication mechanisms are available, see the [PostgreSQL Client Authentication documentation](https://www.postgresql.org/docs/current/static/client-authentication.html).

For this tutorial, we will use the default "peer" authentication, so we need to create a PostgreSQL user matching the system user under which Galaxy will be running, i.e. `galaxy`. This is normally done with the PostgreSQL `createuser` command, and it must be run as the `postgres` user. In our case, we will use the `natefoo.postgresql_objects` role to handle this step.

> ### {% icon hands_on %} Hands-on: Installing PostgreSQL
>
> 1. Create and edit `group_vars/galaxyservers.yml` and add some variables to configure PostgreSQL:
>
>    ```yaml
>    postgresql_objects_users:
>      - name: galaxy
>    postgresql_objects_databases:
>      - name: galaxy
>        owner: galaxy
>    ```
>
> 2. Create and open `playbook.yml` with your text editor and add the following:
>
>    - Add a pre-task to install the necessary dependency, `python-psycopg2`
>    - A role for `galaxyproject.repos`. This will add the additional repositories that are needed by Galaxy in various places.
>    - A role for `galaxyproject.postgresql`. This will handle the installation of PostgreSQL.
>    - A role for `natefoo.postgresql_objects`, run as the postgres user. (You will need `become`/`become_user`.) This role allows for managing users and databases within postgres.
>
>    > ### {% icon question %} Question
>    >
>    > How does your current playbook look?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > ```yaml
>    > > - hosts: galaxyservers
>    > >   become: true
>    > >   pre_tasks:
>    > >     - name: Install Dependencies
>    > >       package:
>    > >         name: 'python-psycopg2'
>    > >   roles:
>    > >     - galaxyproject.repos
>    > >     - galaxyproject.postgresql
>    > >     - role: natefoo.postgresql_objects
>    > >       become: true
>    > >       become_user: postgres
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 3. Run the playbook (`ansible-playbook -i hosts playbook.yml`)
>
>    > ### {% icon tip %} Tip: When running Ansible
>    > Always pay close attention to tasks reported as **changed** and ensure that the changes were expected!
>    {: .tip}
>
> 4. Inspect the changes that have been made on your Galaxy server. Places to look include:
>     - `/etc/postgresql`
>     - Databases and users in PostgreSQL (connection details are below). Hint: try the commands `\du` and `\l`.
>
{: .hands_on}

You can now login and access the database, but only as the `postgres` user. You will need to `sudo -iu postgres` first, and then you can run `psql galaxy`. The database will currently be empty as Galaxy has never connected to it yet. Once you install Galaxy in the next step, the database will be populated.


## Galaxy

Next we will dive right in to deploying a copy of Galaxy onto our server, but it will just be a static copy of the code without anything running.

For a normal Galaxy instance there are a few configuration changes you make very early during deployment:

- Changing the database connection
- Configuring the admin user list
- Changing the "brand"

Additionally we'll go ahead and set up the production-ready [uWSGI Mules](https://uwsgi-docs.readthedocs.io/en/latest/Mules.html) which will handle processing Galaxy jobs. With Mules, uWSGI launches as many as you request, and then they take turns placing a lock, accepting a job, releasing that lock, and then going on to process that job.

Finally, best admin practices are to not run Galaxy as a user with `sudo` access, like your login user probably has. Additionally, it is best to install the Galaxy code and configs as a separate user, for security purposes. So we will instruct the `galaxyproject.galaxy` role to create a new user account specifically to run Galaxy under.

> ### {% icon tip %} Tip: Mules are not the only option
>
> Galaxy can be run in a [couple of other configurations](https://docs.galaxyproject.org/en/master/admin/scaling.html#deployment-options) depending on your needs. Mules are generally a good solution for most production needs.
>
{: .tip}

The configuration is quite simple thanks to the many sensible defaults that are provided in the Ansible roles.

> ### {% icon hands_on %} Hands-on: Minimal Galaxy Playbook
>
> 1. Open `playbook.yml` with your text editor and set the following:
>
>    - Amend the dependency installation pre-task to install additional necessary dependencies: git, python-virtualenv, make
>    - Add the role `galaxyproject.galaxy` to the roles to be executed, at the end
>
>    > ### {% icon question %} Question
>    >
>    > How does your final configuration look?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > ```yaml
>    > > - hosts: galaxyservers
>    > >   pre_tasks:
>    > >     - name: Install Dependencies
>    > >       package:
>    > >         name: ['python-psycopg2', 'git', 'python-virtualenv', 'make']
>    > >   roles:
>    > >     - galaxyproject.repos
>    > >     - galaxyproject.postgresql
>    > >     - role: natefoo.postgresql_objects
>    > >       become: true
>    > >       become_user: postgres
>    > >     - galaxyproject.galaxy
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 2. Edit your group variables file for your group (`group_vars/galaxyservers.yml`).
>
>    We need to set the following variables:
>
>    {% raw %}
>    Variable                     | Value                                      | Purpose
>    ---                          | -----                                      | ---
>    `galaxy_create_user`         | `true`                                     | Instruct the role to create a Galaxy user
>    `galaxy_separate_privileges` | `true`                                     | Enable separation mode to install the Galaxy code as `root` but run the Galaxy server as `galaxy`
>    `galaxy_manage_paths`        | `true`                                     | Instruct thre role to create the needed directories.
>    `galaxy_layout`              | `root-dir`                                 | This enables the `galaxy_root` Galaxy deployment layout:all of the code, configuration, and data folders will live beneath `galaxy_root`.
>    `galaxy_root`                | `/srv/galaxy`                              | This is the root of the Galaxy deployment.
>    `galaxy_file_path`           | `/data`                                    | The directory where Galaxy datasets (user data) will be stored. On a real deployment, this would likely be a mounted network filesystem.
>    `galaxy_user`                | `{name: galaxy, shell: /bin/bash}`         | The user that Galaxy will run as.
>    `galaxy_commit_id`           | `release_18.09`                            | The git reference to check out, which in this case is the <br>branch for Galaxy Release 18.09.
>    `galaxy_config_style`        | `yaml`                                     | We want to opt-in to the new style YAML configuration.
>    `galaxy_force_checkout`      | `true`                                     | If we make any modifications to the Galaxy codebase, they will be removed. This way we know we're getting an unmodified Galaxy and no one has made any unexpected changes to the codebase.
>    `check_migrate_tools`        | `false`                                    | Must be set to false in this case due to a new install of Galaxy
>    {% endraw %}
>
>
> 3. Again edit the group variables file and add a variable for `galaxy_config`. It will be a hash with one key, `galaxy` which will also be a hash. Inside here you can place all of your Galaxy configuration.
>
>    So tthe structure looks like:
>    ```yaml
>    galaxy_config:
>       galaxy:
>           key: value
>    ```
>
>    Now you should set:
>    1. `admin_users` to the email address you will use with this Galaxy.
>    2. `brand` to something fun!
>    3. `database_connection` to point to the database you setup earlier (`postgresql:///galaxy?host=/var/run/postgresql`)
>    4. `file_path` to a place to store data, `/data` for this lesson.
>
>    > ### {% icon question %} Question
>    >
>    > How does your current group variables file look?
>    >
>    > > ### {% icon solution %} Solution
>    > > {% raw %}
>    > > ```yaml
>    > > ---
>    > > # PostgreSQL
>    > > postgresql_objects_users:
>    > >   - name: galaxy
>    > >     password: null
>    > > postgresql_objects_databases:
>    > >   - name: galaxy
>    > >     owner: galaxy
>    > >
>    > > # Galaxy
>    > > galaxy_create_user: true
>    > > galaxy_separate_privileges: true
>    > > galaxy_manage_paths: true
>    > > galaxy_layout: root-dir
>    > > galaxy_root: /srv/galaxy
>    > > galaxy_file_path: /data
>    > > galaxy_user:
>    > >   name: galaxy
>    > >   shell: /bin/bash
>    > > galaxy_commit_id: release_18.09
>    > > galaxy_config_style: yaml
>    > > galaxy_force_checkout: true
>    > >
>    > > galaxy_config:
>    > >   galaxy:
>    > >     brand: "My Galaxy"
>    > >     admin_users: admin@example.org
>    > >     database_connection: "postgresql:///galaxy?host=/var/run/postgresql"
>    > >     file_path: /data
>    > >     check_migrate_tools: false
>    > > ```
>    > > {% endraw %}
>    > {: .solution }
>    >
>    {: .question}
>
> 4. In order to use mule messaging, we need to edit the uWSGI configuration of Galaxy. This has a default value, but we will have to override it. Add the following configuration as a child of the `galaxy_config` variable:
>
>    {% raw %}
>    ```yaml
>    uwsgi:
>      # Default values
>      http: 0.0.0.0:8080
>      buffer-size: 16384
>      processes: 1
>      threads: 4
>      offload-threads: 2
>      static-map:
>        - /static/style={{ galaxy_server_dir }}/static/style/blue
>        - /static={{ galaxy_server_dir }}/static
>      master: true
>      virtualenv: "{{ galaxy_venv_dir }}"
>      pythonpath: "{{ galaxy_server_dir }}/lib"
>      module: galaxy.webapps.galaxy.buildapp:uwsgi_app()
>      thunder-lock: true
>      die-on-term: true
>      hook-master-start:
>        - unix_signal:2 gracefully_kill_them_all
>        - unix_signal:15 gracefully_kill_them_all
>      py-call-osafterfork: true
>      enable-threads: true
>      # Our additions
>      mule:
>        - lib/galaxy/main.py
>        - lib/galaxy/main.py
>      farm: job-handlers:1,2
>    ```
>    {% endraw %}
>
>    > ### {% icon question %} Question
>    >
>    > How does your current group variables file look?
>    >
>    > > ### {% icon solution %} Solution
>    > > {% raw %}
>    > > ```yaml
>    > > ...
>    > > galaxy_config:
>    > >   galaxy:
>    > >     brand: "My Galaxy"
>    > >     admin_users: admin@example.org
>    > >     database_connection: "postgresql:///galaxy?host=/var/run/postgresql"
>    > >     file_path: /data
>    > >     check_migrate_tools: false
>    > >   uwsgi:
>    > >     ...
>    > >     mule:
>    > >       - lib/galaxy/main.py
>    > >       - lib/galaxy/main.py
>    > >     farm: job-handlers:1,2
>    > > ```
>    > > {% endraw %}
>    > {: .solution }
>    >
>    {: .question}
>
>
> 5. Run the playbook.
>
>    `ansible-playbook -i hosts playbook.yml`
>
> 6. SSH into your server and explore what has been set up for you.
>     - Galaxy has been deployed to `/srv/galaxy/server`
>     - The configuration lives in `/srv/galaxy/config/galaxy.yml` - be sure to look through it to see what default options have been set for you
>     - Note the permissions of the contents of `/srv/galaxy`
>     - Some config files that Galaxy maintains itself, such as `shed_tool_conf.xml`, which controls what tools that you have installed from the Tool Shed will be loaded, have been instantiated in `/srv/galaxy/var/config`
>     - A Python virtualenv - an isolated Python environment - with all of the Galaxy framework's dependencies has been installed in `/srv/galaxy/venv`
>
{: .hands_on}

> ### {% icon tip %} Tip: ansible.cfg
> Typing `-i hosts` every time can be a bit repetitive, you can save having to type this flag by creating an `ansible.cfg` file (next to your playbook) with the following contents:
>
> ```ini
> [defaults]
> inventory = hosts
> ```
>
> There are some additional useful options that you might want to add to your `ansible.cfg` file:
>
> ```ini
> [ssh_connection]
> pipelining = true
> [defaults]
> retry_files_enabled = false
> ```
>
> Pipelining will make [ansible run faster](https://docs.ansible.com/ansible/latest/reference_appendices/config.html#ansible-pipelining) by significantly reducing the number of new SSH connections that must be opened. Setting `retry_files_enabled = false` will prevent Ansible from creating `playbook.retry` files whenever a playbook crashes before finishing. These are rarely useful for the cases in which we run Ansible.
>
> For users running with the local connection, you can specify this in your `hosts` file:
>
> ```ini
> [galaxyservers]
> localhost ansible_connection=local
> ```
{: .tip}

Galaxy is now configured with an admin user, a database, and a place to store data. Additionally we've immediately configured the mules for production Galaxy serving. So we're ready to set up supervisord which will manage the Galaxy processes!

> ### {% icon hands_on %} Hands-on: (Optional) Launching uWSGI by hand
>
> 1. SSH into your server
> 2. Switch user to Galaxy account (`sudo -iu galaxy`)
> 3. Change directory into `/srv/galaxy/server`
> 4. Activate virtualenv (`. ../venv/bin/activate`)
> 5. `uwsgi --yaml ../config/galaxy.yml`
> 6. Accces at port `<ip address>:8080` once the server has started
{: .hands_on}

### Supervisord

Launching Galaxy by hand is not a good use of your time, so we will immediately switch to a process manager for that, [supervisord](http://supervisord.org/). If you're familiar with systemd, supervisord does many of the same things. We use supervisord instead of the native init system as it supports some of Galaxy's use cases better and was fully featured long before SystemD became common.

> ### {% icon hands_on %} Hands-on: Supervisord
>
> 1. Add the roles `geerlingguy.pip` and `usegalaxy-eu.supervisor` to your playbook, these need to install things and should run as the root user.
>
> 2. Supervisor defines `programs` which should be executed with additional metadata like whether or not they should be restarted, what user they should run as, etc. Just like any other init system like SystemD. We will define a program for Galaxy which will directly invoke uWSGI, rather than run.sh, as run.sh does some additional tasks we do not need to do on every Galaxy start (e.g. rebuilding the client). For some setups like zerglings it is required that you use supervisord as you need to start multiple processes.
>
>    {% raw %}
>    ```yaml
>    supervisor_programs:
>      - name: galaxy
>        state: present
>        command: uwsgi --yaml {{ galaxy_config_dir }}/galaxy.yml
>        configuration: |
>          autostart=true
>          autorestart=true
>          startretries=1
>          startsecs=10
>          user=galaxy
>          umask=022
>          directory={{ galaxy_server_dir }}
>          environment=HOME={{ galaxy_mutable_data_dir }},VIRTUAL_ENV={{ galaxy_venv_dir }},PATH={{ galaxy_venv_dir }}/bin:%(ENV_PATH)s
>    ```
>    {% endraw %}
>
>    Here we've defined a `galaxy` command that should be `present`. It will run the command `uwsgi ...` and is set to automatically start when supervisord starts and restart if it crashes, with 1 second between the retries. It will wait 10 seconds to see if the program has not crashed, and if it reaches this threshold it will be marked as `running`. It starts as the galaxy user with a umask of `022` (files created will be world readable by default). Its working directory on startup is the root of the Galaxy (cloned) code, and will run with the defined environment variables set.
>
> 3. Now that we have defined a process manager for Galaxy, we can also instruct `galaxyproject.galaxy` to use Supervisor to restart it when Galaxy is upgraded or config changes are made. To do so, open `playbook.yml` and add a `handlers:` section at the same level as `pre_tasks:` and `roles:`, and add a task to restart Galaxy using the [supervisorctl Ansible module](https://docs.ansible.com/ansible/latest/modules/supervisorctl_module.html). Handlers are structured just like tasks:
>
>    ```yaml
>    - hosts: galaxyservers
>      pre_tasks:
>        - name: Install Dependencies
>          package:
>            name: ['python-psycopg2', 'git', 'python-virtualenv', 'make']
>      handlers:
>        - name: Restart Galaxy
>          supervisorctl:
>            name: galaxy
>            state: restarted
>      roles:
>        ...
>    ```
>
> 4. Open your group variables file and we'll add some variables for supervisor. Supervisor communicates over a unix or tcp socket; we will use the unix socket without password authentication, instead using user/group authentication. We will thus need to set a couple of variables to allow our Galaxy user to access this. Additionally, we need to inform `galaxyproject.galaxy` what the name of our new handler is. Add the following to your group variables file:
>
>    ```yaml
>    supervisor_socket_user: galaxy
>    supervisor_socket_chown: galaxy
>    galaxy_restart_handler_name: Restart Galaxy
>    ```
>
> 5. Log in and check the status with `supervisorctl status` (remember to change to the Galaxy user)
>
>    > ### {% icon question %} Question
>    >
>    > How does the output look?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > If everything went correctly you should see something like
>    > >
>    > > ```yaml
>    > > galaxy                  RUNNING   pid 2246972, uptime 0:02:00
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>    Take a look at the supervisor configs that have been written in `/etc/supervisor` and `/etc/supervisor/conf.d`.
>
> 6. Some things to note:
>
>    1. Refreshing the page before Galaxy has restarted will hang, which is a nice feature of uWSGI
>    2. Although the playbook will restart Galaxy upon config changes, you will sometimes need to restart it by hand, which can be done with `supervisorctl restart galaxy`
>    3. You can use `supervisorctl tail -f galaxy` and `supervisorctl tail -f galaxy stderr` to see the logs of Galaxy
>
{: .hands_on}

Galaxy should now be accessible over port :8080, again try connecting to your VM now and checking that Galaxy is working. Note that the welcome page is broken, this is a known issue, and a good reminder to write your own :)

### NGINX

With this we have:

- PostgreSQL running
- Galaxy running (managed by supervisord)

When we first configured Galaxy, we used the setting `http: 0.0.0.0:8080`, which instructed uWSGI to handle the serving of Galaxy, and to process the HTTP requests itself. This has some overhead and is not as efficient as is desired in production. So we will set up a reverse proxy to handle the HTTP processing, and translate this into the more efficient uWSGI protocol. Additionally it can handle serving static files for us without the requests going through uWSGI, allowing it to spend more time on useful tasks like processing jobs.

For this, we will use NGINX. It is possible to configure Galaxy with Apache and potentially other webservers but this is not the configuration that receives the most testing.

> ### {% icon hands_on %} Hands-on: NGINX
>
> 1. Edit your `group_vars/galaxyservers.yml`, we will update the line that `http: 0.0.0.0:8080` to be `socket: 127.0.0.1:8080`. This will cause uWSGI to only respond to uWSGI protocol, and only to requests originating on localhost.
>
> 2. Add the role `geerlingguy.nginx` to your playbook and have it run as root.
>
> 3. In our group variables file, we will set a convenience variable, `hostname`, to refer to our VM. If you have a DNS entry for your server, add `hostname: <your server name>` to your group variables. Otherwise, add {% raw %}`hostname: "{{ ansible_hostname }}"`{% endraw %}.
>
> 4. We need to configure the virtualhost. This is a slightly more complex process as we have to write the proxying configuration ourselves. This may seem annoying, but it is often the case that sites have individual needs to cater to, and it is difficult to provide a truly generic webserver configuration.
>
>    Add the following to your group variables file:
>    {% raw %}
>    ```yaml
>    nginx_package_name: nginx-full # nginx-galaxy on RHEL/CentOS
>    nginx_remove_default_vhost: true
>    nginx_server_names_hash_bucket_size: "128"
>    nginx_vhosts:
>      - listen: "80"
>        server_name: "{{ hostname }}"
>        root: "/var/www/{{ hostname }}"
>        index: "index.html"
>        access_log: "/var/log/nginx/access.log"
>        error_log: "/var/log/nginx/error.log"
>        state: "present"
>        filename: "{{ hostname }}.conf"
>        extra_parameters: |
>            client_max_body_size 10G; # aka max upload size, defaults to 1M
>            uwsgi_read_timeout 2400;
>
>            location / {
>                uwsgi_pass      127.0.0.1:8080;
>                uwsgi_param UWSGI_SCHEME $scheme;
>                include         uwsgi_params;
>            }
>
>
>            location /static {
>                    alias {{ galaxy_server_dir }}/static;
>                    expires 24h;
>            }
>
>            location /static/style {
>                    alias {{ galaxy_server_dir }}/static/style/blue;
>                    expires 24h;
>            }
>
>            location /static/scripts {
>                    alias {{ galaxy_server_dir }}/static/scripts;
>                    expires 24h;
>            }
>
>            location /robots.txt {
>                    alias {{ galaxy_server_dir }}/static/robots.txt;
>            }
>
>            location /favicon.ico {
>                    alias {{ galaxy_server_dir }}/static/favicon.ico;
>            }
>
>            location /static/welcome.html {
>                    alias {{ galaxy_server_dir }}/static/welcome.html.sample;
>            }
>    ```
>    {% endraw %}
>
>    This is a lot of configuration but it is not very complex to understand. We'll go through it step by step:
>
>    This snippet says that the role should install the nginx-full package rather than the normal `nginx` package. Earlier we added the role `galaxyproject.repos`, this role gave us access to the Galaxy Project's APT and YUM repositories.
>    {% raw %}
>    ```yaml
>    nginx_package_name: nginx-full # nginx-galaxy on RHEL/CentOS
>    nginx_remove_default_vhost: true
>    nginx_server_names_hash_bucket_size: "128"
>    ```
>    {% endraw %}
>
>
>    Here we define a virtual host, usually used to respond on multiple hostnames with different content. In our case we only have a single virtualhost.
>
>    The value of `server_name` does not need to be a valid domain name in our case, because we only have a single virtualhost. Normally NGINX will try and match the hostname to the server_name, but in the case there is no match it will route it to the default server for that port, which is the first one defined.
>
>    Additionally we specify log files and a directory for static files which will never be used.
>
>    {% raw %}
>    ```yaml
>    nginx_vhosts:
>      - listen: "80"
>        server_name: "{{ hostname }}"
>        root: "/var/www/{{ hostname }}"
>        index: "index.html"
>        access_log: "/var/log/nginx/access.log"
>        error_log: "/var/log/nginx/error.log"
>        state: "present"
>        filename: "{{ hostname }}.conf"
>    ```
>    {% endraw %}
>
>    Lastly we specify the "extra parameters" for this virtualhost which are included within the configuration file.
>
>    {% raw %}
>    ```yaml
>        extra_parameters: |
>            client_max_body_size 10G; # aka max upload size, defaults to 1M
>            uwsgi_read_timeout 2400;  # Wait up to 40 minutes for requests to timeout.
>
>            location / {
>                uwsgi_pass      127.0.0.1:8080;    # Proxy our requests to port 8080, just like we defined in our Galaxy group variables for the uWSGI configuration
>                uwsgi_param UWSGI_SCHEME $scheme;  # Passes the scheme (http or https) to the uWSGI application. Used in generating links.
>                include         uwsgi_params;      # Include some default uwsgi settings
>            }
>
>            # Proxy any requests for static content directly to their folders on disk
>            location /static {
>                    alias {{ galaxy_server_dir }}/static;
>                    expires 24h;
>            }
>
>            location /static/style {
>                    alias {{ galaxy_server_dir }}/static/style/blue;
>                    expires 24h;
>            }
>
>            location /static/scripts {
>                    alias {{ galaxy_server_dir }}/static/scripts;
>                    expires 24h;
>            }
>
>            location /robots.txt {
>                    alias {{ galaxy_server_dir }}/static/robots.txt;
>            }
>
>            location /favicon.ico {
>                    alias {{ galaxy_server_dir }}/static/favicon.ico;
>            }
>
>            location /static/welcome.html {
>                    alias {{ galaxy_server_dir }}/static/welcome.html.sample;
>            }
>    ```
>    {% endraw %}
>
> 5. Run the playbook. At the very end, you should see output like the following indicating that Galaxy has been restarted:
>
>    ```
>    RUNNING HANDLER [restart galaxy] ****************************************
>    changed: [galaxy.example.org]
>    ```
>
>
> 6. Check out the changes made to your server in `/etc/nginx`, particularly the directory containing the Galaxy virtualhost.
>
> 7. Your Galaxy should now be accessible and served efficiently! Try registering (using the admin email from earlier) and maybe executing a couple of jobs. The author's favourite tool (speaking as an admin) is the `secure hash digest` tool, it's perfect for testing.
>
{: .hands_on}

## ProFTPD

With a large Galaxy instance, users will often request FTP access in order to upload large data sets. The [ProFTPD](http://proftpd.org/) server works well for Galaxy, it can leverage an SQL authentication module, to allow users to login to the FTP server with their normal Galaxy username and password.

> ### {% icon hands_on %} Hands-on: ProFTPD
>
> 1. Add the role `galaxyproject.proftpd` to your playbook and have it run as root.
>
> 2. Add a `pre_task` to [create the FTP directory](https://docs.ansible.com/ansible/2.5/modules/file_module.html#file-module), owned by galaxy:galaxy, with mode `0750`.
>
> 3. Edit the group variables file, we will define some variables for the proftpd role:
>
>    {% raw %}
>    ```yaml
>    galaxy_ftp_upload_dir: "{{ galaxy_root }}/ftp"
>    proftpd_display_connect: |
>      Unauthorized access is prohibited
>    proftpd_galaxy_auth: yes
>    __galaxy_user_name: galaxy
>    proftpd_options:
>      - User: galaxy
>      - Group: galaxy
>    proftpd_sql_user: galaxy
>    proftpd_sql_db: galaxy@/var/run/postgresql
>
>    proftpd_virtualhosts:
>      - id: galaxy
>        address: "{{ hostname }}"
>        options:
>          - ServerAdmin: admin@example.org
>          - ServerName: Galaxy FTP
>          - TransferLog: /var/log/proftpd/xfer.log
>          - MaxLoginAttempts: 3
>          - RequireValidShell: no
>          - AllowOverwrite: yes
>    ```
>    {% endraw %}
>
> 4. And additionally two changes to the Galaxy configuration portion:
>
>    {% raw %}
>    ```yaml
>    ...
>    galaxy_config:
>      galaxy:
>        ...
>        ftp_upload_dir: "{{ galaxy_ftp_upload_dir }}"
>        ftp_upload_site: "ftp://{{ hostname}}"
>        ...
>      uwsgi: ...
>    ```
>    {% endraw %}
>
> 5. Apply the playbook.
{: .hands_on}

With this, users can see an "Upload via FTP" button in their upload interface, and they can upload their data using Filezilla or other preferred FTP client.


> ### {% icon tip %} Tip: Active vs Passive
> If connecting does not work, try forcing your client to use active FTP, university firewalls often have issues with passive FTP.
{: .tip}

## Disaster Strikes!

Because you're an admin, you need to be prepared for any situation, including the worst case scenarios. So we're going to simulate a disaster and show you how you can recover from it. It'll be fun!

For this "disaster", we will pretend that:

1. Your database is on another machine
2. Your datasets are on an NFS server or some other remote machine.

> ### {% icon hands_on %} Hands-on: Summon the Apocalypse
>
> So let's have a "worst case scenario", where your Galaxy server gets destroyed
>
> 1. Log on to your machine.
> 2. **Carefully**, as root, `rm -rf /srv/galaxy`, completely wipe out your Galaxy home directory.
>
{: .hands_on}

Your entire Galaxy server is gone. You were a responsible admin and had your user data and database stored on a separate system (and backed up), so at least those survived. Nevertheless, this is when most of us start feeling really bad; bosses start yelling, we start crying or reaching for bad coping habits.

But not you! You spent the day writing this Ansible playbook that describes your environment completely; all of the software that was installed, all of the configuration changes you've made. It leverages many community maintained roles and can be used to completely rebuild the server! With minimal effort on your part.

> ### {% icon hands_on %} Hands-on: Revert the Apocalypse
>
> 1. `ansible-playbook -i hosts playbook.yml`
>
{: .hands_on}

And with that, Galaxy should be up and running again. If you log in, you should see the results of any jobs you ran earlier, you should still be able to login with your old account, everything should just work.

Ansible can save you from some really bad scenarios, *if and only if*:

- You can replace the hardware or find somewhere new to re-deploy
- You've made absolutely certain that every change made to a system is recorded within your playbooks and roles (i.e. no manual package installation)

Then you can potentially use it to recover.

> ### {% icon tip %} Tip: We have experience
>
> We can tell you this, we can repeat it over and over, but unless you really have a disaster happen to you, it is hard to appreciate how important it is that machines are completely controlled in terms of configuration and software deployment.
>
> We've experienced these incidents and we know how horribly stressful it can be if an important service like Galaxy goes down and you cannot immediately replace it with another instance. We hope you will immediately apply the lessons from this training material, it can potentially save you a lot of stress and worry.
>
{: .tip}

## (Optional) Securing your Instance

Now that you've gotten through the worst case scenario, we'll attack the next worst case scenario, auditors! They've shown up and demanded that everything have valid SSL certificates.

This step uses [Let's Encrypt](https://letsencrypt.org/) for generating certificates, so it assumes that:

1. Your machine is publicly accessible
2. It has a publicly resolvable DNS entry

If you do not meet these requirements, you should read through them to see the changes that are required.

> ### {% icon hands_on %} Hands-on: SSL Certificates with Let's Encrypt (LE)!
>
> 1. We need to configure NGINX to proxy requests to LE, because LE needs to authenticate on 443 or 80, and nginx is already bound to 80. Add this configuration to the top of your `nginx_vhosts` > `extra_parameters` section in your group variables file:
>
>    ```yaml
>    location /.well-known/ {
>        proxy_set_header           Host $host:$server_port;
>        proxy_set_header           X-Real-IP $remote_addr;
>        proxy_set_header           X-Forwarded-For $proxy_add_x_forwarded_for;
>        proxy_set_header           X-Forwarded-Proto $scheme;
>        proxy_pass                 http://127.0.0.1:8118;
>        proxy_pass_request_headers on;
>    }
>    ```
>
> 2. Run the playbook
>
> 3. Add the role `usegalaxy-eu.certbot` to your playbook and have it run as root.
>
> 4. We need to add some variables for this:
>
>    {% raw %}
>    ```yaml
>    certbot_auto_renew: yes
>    certbot_auto_renew_user: root
>    certbot_auto_renew_hour:    # Pick a number between 0 and 23, inclusive
>    certbot_auto_renew_minute:  # Pick a number between 0 and 59, inclusive
>    certbot_auto_renew_extra: "--preferred-challenges http-01 --http-01-port 8118"
>    certbot_environment: staging
>    certbot_domains:
>      - "{{ hostname }}"
>    certbot_agree_tos: --agree-tos
>    certbot_admin_email: security@usegalaxy.eu
>    certbot_share_key_users:
>      - nginx
>    certbot_post_renewal: |
>        systemctl restart nginx || true
>    ```
>    {% endraw %}
>
>    Be sure to pick a actual value for `certbot_auto_renew_minute` and `certbot_auto_renew_hour`
>    as noted in the commands above.
>
> 5. Run the playbook
>
{: .hands_on}

And now we should have valid SSL certificates! We just need to go back and update our various services to support this.

> ### {% icon tip %} Tip: Issues on CentOS
> Check that setenforce is permissive, or allow nginx to connect back to localhost on a non-standard port.
{: .tip}


> ### {% icon hands_on %} Hands-on: Securing NGINX
>
> 1. Open your group variables, we need to re-structure the NGINX configuration to look like the following. As before we will show you the necessary configuration.
>
>
>    {% raw %}
>    ```yaml
>    nginx_vhosts:
>      - listen: "80"
>        server_name: "{{ hostname }}"
>        return: "301 https://{{ hostname }}$request_uri"
>        filename: "{{ hostname }}.80.conf"
>
>      - listen: "443 ssl"
>        server_name: "{{ hostname }}"
>        root: "/var/www/{{ hostname }}"
>        index: "index.html"
>        access_log: "/var/log/nginx/access.log"
>        error_log: "/var/log/nginx/error.log"
>        state: "present"
>        filename: "{{ hostname }}.conf"
>        extra_parameters: |
>            client_max_body_size 10G; # aka max upload size, defaults to 1M
>            uwsgi_read_timeout 2400;
>
>            location / {
>                uwsgi_pass      127.0.0.1:8080;
>                uwsgi_param UWSGI_SCHEME $scheme;
>                include         uwsgi_params;
>            }
>
>            location /.well-known/ {
>                proxy_set_header           Host $host:$server_port;
>                proxy_set_header           X-Real-IP $remote_addr;
>                proxy_set_header           X-Forwarded-For $proxy_add_x_forwarded_for;
>                proxy_set_header           X-Forwarded-Proto $scheme;
>                proxy_pass                 http://127.0.0.1:8118;
>                proxy_pass_request_headers on;
>            }
>
>            location /static {
>                    alias {{ galaxy_server_dir }}/static;
>                    expires 24h;
>            }
>
>            location /static/style {
>                    alias {{ galaxy_server_dir }}/static/style/blue;
>                    expires 24h;
>            }
>
>            location /static/scripts {
>                    alias {{ galaxy_server_dir }}/static/scripts;
>                    expires 24h;
>            }
>
>            location /robots.txt {
>                    alias {{ galaxy_server_dir }}/static/robots.txt;
>            }
>
>            location /favicon.ico {
>                    alias {{ galaxy_server_dir }}/static/favicon.ico;
>            }
>
>            location /static/welcome.html {
>                    alias {{ galaxy_server_dir }}/static/welcome.html.sample;
>            }
>            ssl_certificate /etc/ssl/certs/fullchain.pem;
>            ssl_certificate_key /etc/ssl/private/privkey-nginx.pem;
>
>            ssl_protocols TLSv1.2;# Requires nginx >= 1.13.0 else use TLSv1.2
>            ssl_prefer_server_ciphers on;
>            ssl_ciphers ECDHE-RSA-AES256-GCM-SHA512:DHE-RSA-AES256-GCM-SHA512:ECDHE-RSA-AES256-GCM-SHA384:DHE-RSA-AES256-GCM-SHA384:ECDHE-RSA-AES256-SHA384;
>            ssl_ecdh_curve secp384r1;
>            ssl_session_timeout  10m;
>            ssl_session_cache shared:SSL:10m;
>            ssl_session_tickets off;
>            #ssl_stapling on;
>            #ssl_stapling_verify on;
>
>            resolver 8.8.8.8 8.8.4.4 valid=300s;
>            resolver_timeout 5s;
>            add_header X-Content-Type-Options nosniff;
>            add_header X-XSS-Protection "1; mode=block";
>            add_header X-Robots-Tag none;
>    ```
>    {% endraw %}
>
> 2. Add the following proftpd configuration:
>
>    ```yaml
>    proftpd_tls_cipher_suite: AES128+EECDH:AES128+EDH
>    proftpd_tls_protocol: TLSv1.2
>    proftpd_conf_ssl_certificate: /etc/ssl/certs/cert.pem
>    proftpd_conf_ssl_certificate_key: /etc/ssl/private/privkey.pem
>    proftpd_conf_ssl_ca_certificate: /etc/ssl/certs/fullchain.pem
>    ```
>
> 3. Run the playbook
>
> 4. You are now running with staging certificates. You will need production certificates
>    for most applications, so now change the `certbot_environment` group variable to
>    `production`, wipe out the existing certifacts with `sudo rm -rf /etc/letsencrypt`,
>    and re-run the playbook.
>
{: .hands_on}

# Final Notes

If you've been following along you should have a production-ready Galaxy, secured, everything ready to go.

If you missed any steps, you can compare against a reference [playbook.yml](./playbook.txt), and [group_vars/galaxyservers.yml](./galaxyservers.txt)
