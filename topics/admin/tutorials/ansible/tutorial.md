---
layout: tutorial_hands_on

title: "Ansible"
zenodo_link: ""
questions:
  - Why Ansible?
  - How and when to use Ansible?
  - How to write a role?
  - How to leverage community build roles?
objectives:
  - Learn Ansible basics
  - Write a simple role
  - Install a role from Ansible Galaxy
time_estimation: "60m"
key_points:
  - Ansible lets you do system administration at scale
  - Many system administration, software installation, and software management tasks are already available as ansible tasks or roles
contributors:
  - erasche
  - shiltemann
subtopic: core

---

# Overview
{:.no_toc}

In this tutorial we have briefly cover what Ansible is and how to understand what it does. This guide is not meant to make you an expert on Ansible, but perhaps give you enough that you can debug broken roles and modify them to suit your needs. Or maybe to contribute to the [Galaxyproject Ansible roles](https://github.com/galaxyproject?q=ansible).

This will be a very practical training with emphasis on looking at examples from modules and becoming self sufficient.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}


# What is Ansible?

Ansible runs commands on local or remote computers. It can move files around, create files from templates, and run command line tools. Primarily used for system administration tasks at scale. It has a push model rather than a pull model like puppet. If you've used Puppet, Ansible doesn't evaluate what changes need to be made and make those, it just runs through all of commands every time.

Some terms that you should know first:

Inventory file
:    An Ansible-specific file that defines the systems ("hosts") and groups of hosts on which Ansible should operate.

Ansible module
:    A piece of Python code that converts some parameters into an invocation. An example would be the `command` module which converts parameters like `command: ls` into a command line that is executed. There are pre-built modules for just about everything.

task
:    A call to an Ansible module that should be executed and the configuration for this module.

role
:    A folder containing some tasks, templates, files, and default values for variables. People share roles on ["Ansible Galaxy"](https://galaxy.ansible.com/).

playbook
:    a YAML file listing a set of tasks and/or roles that should be applied to a group of hosts.

vault
:    An encrypted YAML file. You put your secrets here and then you can use them in tasks/roles/playbooks.

Looking at each of these briefly:

## Inventory file

```ini
[webservers]
web_a
web_b

[databases]
db_1.example.org ansible_user=root
```

Here we've defined two groups of computers, `webservers` and `databases`. `ansible_user` is used to specify which user to connect with.

> ### {% icon details %} Ansible Inventory Documentation
> For more advanced features of the inventory file, check out [the official documentation on this topic](https://docs.ansible.com/ansible/latest/user_guide/intro_inventory.html).
{: .details}

## Roles

We will look at [ansible-cvmfs](https://github.com/galaxyproject/ansible-cvmfs) as our example for the layout of a role:

```
├── defaults
│   └── main.yml
├── files
│   ├── cvmfs_wipecache.c
│   ├── cvmfs_wipecache.centos_6
│   └── cvmfs_wipecache.centos_7
├── handlers
│   └── main.yml
├── meta
│   └── main.yml
├── tasks
│   └── main.yml
└── templates
    └── stratum1_squid.conf.j2
```

These are the folders that are included in many complex roles. Simpler roles will often not need all of the folders.

Folder    | Usage
--------- | ------
defaults  | Default variables values go here (e.g. "version of software to install").
files     | These are files which should be copied as-is over to the remote location.
handlers  | This is for restarting processes usually.
meta      | Only needed if you publish your role to Ansible Galaxy.
tasks     | **Always start reading here**. This is the most important folder and the best place to start when trying to understand what an unfamiliar role does. Anything that is loaded will be referenced here, e.g. variables to load, handlers, files, templates.
templates | Files that are templated out with variables before being copied.

> ### {% icon details %} Ansible Role Documentation
>
> For more information check out [the official documentation on this topic](https://docs.ansible.com/ansible/latest/user_guide/playbooks_reuse_roles.html).
{: .details}

## Modules and Tasks

A `tasks/main.yml` file calls multiple Ansible modules to accomplish its goal. A typical tasks file looks like:

```yaml
---
- name: Install CernVM-FS client (yum)
  yum:
    name: cvmfs
    state: {% raw %}"{{ 'latest' if cvmfs_upgrade_client else 'present' }}"{% endraw %}
  when: ansible_os_family == "RedHat"

- name: Ensure AutoFS is enabled + running
  service:
    name: autofs
    enabled: yes
    state: started
```

Here we have two tasks. Each has a `name` that will be shown to the person running the playbook.

The first invokes the `yum` module with the arguments `name: cvmfs, state: ...`. This will use yum to install the package named `cvmfs`. The state parameter uses a [Jinja2](http://jinja.pocoo.org/) template to evaluate the value of the variable `cvmfs_upgrade_client`. We can [grep through](https://github.com/galaxyproject/ansible-cvmfs/search?q=cvmfs_upgrade_client&unscoped_q=cvmfs_upgrade_client) the repository and see that `defaults/main.yml` sets that to `false` by default. We can override this if we need, we'll come back to that later. The first task also has a `when` condition to ensure it only runs on RHEL family machines, so RedHat or CentOS. It is better to use the OS-generic [`package`](http://docs.ansible.com/ansible/latest/package_module.html) module, if possible.

The second invokes the [`service`](http://docs.ansible.com/ansible/latest/service_module.html) module. The arguments to this one are quite legible and the functionality can be inferred from the names for the most part: The service `name: autofs` will be `enabled` and its `state` should be `started`.

[Many modules](https://docs.ansible.com/ansible/latest/modules/modules_by_category.html) are available for you to use.

### Stylistic Choices

Ansible accepts two ways to format tasks:

YAML style:

```yaml
- package:
    name: "{% raw %}{{ package_name }}{% endraw %}"
    state: "{% raw %}{{ package_state }}{% endraw %}"
```

And an inline style.

```yaml
{% raw %}- package name={{ package_name }} state={{ package_state }}{% endraw %}
```

Some groups prefer one style or another. You can mix both of these but you probably shouldn't. The inline style does not require quoting of templated values.

## Playbooks

```yaml
- name: CVMFS
  hosts: webservers
  vars:
    cvmfs_numfiles: 4096
  roles:
    - galaxyproject.cvmfs
```

This is a quite minimal playbook. It selects a `hosts` group named `webservers`, overrides the variable `cvmfs_numfiles`, and then says the following set of roles will be executed for this group of hosts. Ansible makes it easy to collect tasks that should apply to a group of hosts and run a playbook for all of those hosts. Some good uses of this are things like ensuring a certain set of users are installed on all of your managed machines, or using one of the package autoupdating roles to make sure your machines are up-to-date.

> ### {% icon details %} Ansible Playbook Documentation
>
> For more information check out [the official documentation on this topic](https://docs.ansible.com/ansible/latest/user_guide/playbooks_intro.html).
{: .details}

### Philosophies

Different groups use playbooks differently. Here is a non-exhaustive list of ways that the author has seen playbooks used:

Strategy                                | Use Cases
--------------------------------------- | ---
A single playbook that does everything  | Works well in cloud use cases where you would rather terminate a VM and just restart from nothing. If 100% of the tasks are idempotent (can be re-run without causing problems) then this can also work well for ensuring all machines that are managed are meeting your compliance standards.
Playbooks that execute one-off commands | Works well when you have commands that cannot safely be re-run (e.g. you are managing a repository on a remote machine and executing a 'commit' action.) Also works if you have multiple people who are responsible for managing a machine but you want to ensure exactly the same commands are run each time a task needs to be done, with no variation.

## Variables

There are a bunch of places variables can be set. [The list is ridiculous.](https://docs.ansible.com/ansible/latest/user_guide/playbooks_variables.html#variable-precedence-where-should-i-put-a-variable) Try and find a convention within your group and stick to that to help manage the chaos.

There are some special places you can put variables that will be loaded automatically:

- `group_vars/<group_name>.yml`: if you run a playbook which has `hosts: webservers`, and `group_vars/webservers.yml` exists, it will be loaded.
- `host_vars/<host_name>.yml`: if you run a playbook and it affects host (e.g.) `api01`, and `host_vars/api01.yml` exists, it will be loaded automatically.

For a real-life example, UseGalaxy.eu [generally attempts](https://github.com/usegalaxy-eu/infrastructure-playbook) to put variables in an appropriately named `group_vars/` file.

## Vaults

[`ansible-vault`](https://docs.ansible.com/ansible/latest/user_guide/vault.html) is a super useful tool that lets you encrypt secrets for your group using a pre-shared key. This allows you to commit ALL of your Ansible playbooks and variables to source control, without the concern of leaking secrets. E.g. [UseGalaxy.eu](https://github.com/usegalaxy-eu/infrastructure-playbook/blob/master/secret_group_vars/all.yml)'s vault file. The encrypted version is not very interesting to look at, but it is mostly to show that we confidently place an encrypted copy of our secrets online, under configuration management. This has made our life a lot easier.


# Your First Playbook and First Role

The above introduction was certainly not enough for you to feel confident in Ansible, so we will now build a small role and a playbook to run this role.

> ### {% icon warning %} Safety First
>
> Many of the things you can do with Ansible can be quite dangerous. As dangerous as normally being at the Linux command line, but scaled across N machines. Be very careful with the changes you plan to make.
> Ansible provides some flags which can help you identify changes before they're made to production systems:
>
> **`--diff`**
> :    This will show the difference whenever a file is changed.
>
> **`--check`**
> :    Will do a dry-run of the playbook, attempting to show you which tasks will execute, which files will be updated. Not all actions can be predicted because commands are not run, if a downstream task depends on the result of a command execution task, Ansible has no way of knowing whether or not it will execute.
>
> In this tutorial we will write to files in `/tmp` as that is a *relatively* safe thing to do. The training material community does not have the resources to test this tutorial across all of the platforms you might want to run it on. Additionally we do not want to be responsible if you accidentally cause permanent damage by following this tutorial.
>
{: .warning}


## A Basic Role

> ### {% icon hands_on %} Hands-on: Setting up our workspace
>
> 1. In this training we will run Ansible on the machine it will modify. This is not best practice, but it is convenient for trainings. You should probably run this in a VM either on the Cloud or in VirtualBox or similar
>
>    All of the steps are the same, no matter which machine Ansible will manage and where you run it. The only difference is the connection setup
>
> 2. [Install Ansible.](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html) where you will run it
>
> 3. Create an directory named `intro` and `cd` into it
>
> 4. Create your inventory file (named `hosts`) in this folder
>
>    1. We will call our group `my_hosts`
>
>    2. Create [an inventory file](https://docs.ansible.com/ansible/latest/user_guide/intro_inventory.html) with the group `my_hosts` and `localhost ansible_connection=local`, which tells ansible to not use SSH, and just use the local connection.
>
>       > ### {% icon solution %} Solution
>       > The file should look like:
>       >
>       > ```ini
>       > [my_hosts]
>       > localhost ansible_connection=local
>       > ```
>       {: .solution }
>
>    > ### {% icon details %} Running Ansible on a different machine
>    >
>    > The best practice is to have playbooks in git, then it doesn't matter whether you run it locally or remotely. The advantage of running remotely is that you can manage dozens of machines simultaneously, rather than just the local machine. This scaling out to N machine is one of the strengths of Ansible.
>    >
>    > In order to run remotely:
>    >
>    > 1. Make sure you can SSH into it. (Test it now)
>    >
>    > 2. We will call our group `my_hosts`
>    >
>    > 3. Create a hosts file with the group `my_hosts` and your host.
>    >
>    >    > ### {% icon solution %} Solution
>    >    > The file should look like:
>    >    >
>    >    > ```ini
>    >    > [my_hosts]
>    >    > your.host.fqdn.or.ip
>    >    > ```
>    >    > Remember that if you SSH in with a username different than your current local user account's name, you will need to specify `ansible_ssh_user=remote-user-name`
>    >    {: .solution }
>    >
>    {: .details}
>
> 5. Create the roles directory, your role, and the tasks folder: `mkdir -p roles/my-role/tasks/`
>
> 6. Create a YAML file in that directory, `roles/my-role/tasks/main.yml` and open it for editing
>
> 7. Define a `copy` task like below:
>
>    ```yaml
>    ---
>    - name: Copy a file to the remote host
>      copy:
>        src: test.txt
>        dest: /tmp/test.txt
>    ```
>
>    You can read about all of the parameters available to the [`copy`](http://docs.ansible.com/ansible/latest/copy_module.html) module on Ansible's documentation.
>
>    > ### {% icon details %} Ansible Module Documentation
>    > You can usually find a module that will represent most commands you will run at the linux cli. Usually by searching the internet for "ansible $do-some-action" e.g. "ansible copy file to server" or "ansible restart service." If you cannot find a module that does it, there is the [`command`](http://docs.ansible.com/ansible/latest/command_module.html) module, but this should be avoided if possible. Expect to have a browser session with 10-30 different Ansible module documentation tabs if you work with Ansible regularly, no one remembers what arguments are available to every module.
>    >
>    {: .details }
>
> 8. Create a `roles/my-role/files` folder, and within it a file named `test.txt`, containing the content "Hello, World"
>
> 9. This is a complete role by itself and will copy the file `test.txt` from the `roles/my-role/files/` folder over to the remote server and place it in `/tmp`.
>
> 10. Open `playbook.yml` for editing in the root folder. Place the following content in there:
>
>     ```yaml
>     ---
>     - hosts: my_hosts
>       roles:
>         - my-role
>     ```
>
>     > ### {% icon question %} Question
>     >
>     > How does your file tree look now? Use `find` or `tree`.
>     >
>     > > ### {% icon solution %} Solution
>     > >
>     > > ```
>     > > .
>     > > ├── hosts
>     > > ├── playbook.yml
>     > > └── roles
>     > >     └── my-role
>     > >         ├── files
>     > >         │   └── test.txt
>     > >         └── tasks
>     > >             └── main.yml
>     > > ```
>     > >
>     > {: .solution }
>     {: .question}
>
> 11. Run one of the following command, whichever is appropriate:
>
>     - Real remote host: `ansible-playbook -i hosts playbook.yml`
>     - Localhost: `ansible-playbook -i hosts -c local playbook.yml`
>
>     Even local users can run the 'real remote host' command, Ansible will just issue a warning. Running with `-c local` silences this warning.
>
>     > ### {% icon question %} Question
>     >
>     > How does the output look?
>     >
>     > > ### {% icon solution %} Solution
>     > >
>     > > The important thing is `failed=0`
>     > >
>     > > ```
>     > > $ ansible-playbook -i hosts playbook.yml -c local
>     > > PLAY [my_hosts] *********************************
>     > > TASK [Gathering Facts] *************************
>     > > ok: [localhost]
>     > > TASK [my-role : Copy] **************************
>     > > changed: [localhost]
>     > > PLAY RECAP *************************************
>     > > localhost                  : ok=2    changed=1    unreachable=0    failed=0    skipped=0    rescued=0    ignored=0
>     > > ```
>     > >
>     > > You can re-run this and it should say `changed=0`
>     > {: .solution }
>     {: .question}
>
>
> 12. Login to the appropriate host and `cat /tmp/test.txt` to see that the change was made.
>
{: .hands_on}

Now that you've done this, here are some starting points for exploration:

- Add more hosts, watch as Ansible executes over all of them in parallel.
- Identify a task you do regularly, e.g. restarting a service. Find the Ansible service module and add that to your playbook.

> ### {% icon comment %} Too Many Cows?
> If you've installed the `cowsay` tool, Ansible (for some reason) will take advantage of that to output a lot of the output with cowsay. To disable this you can `export ANSIBLE_NOCOWS=1` (Remember that exporting will only last as long as the current invocation of your terminal does, so consider adding this to your user profile if you wish to keep cowsay installed and still have legible output.)
>
> ```
>  ________________
> < PLAY [my_hosts] >
>  ----------------
>         \   ^__^
>          \  (oo)\_______
>             (__)\       )\/\
>                 ||----w |
>                 ||     ||
>  ________________________
> < TASK [Gathering Facts] >
>  ------------------------
>         \   ^__^
>          \  (oo)\_______
>             (__)\       )\/\
>                 ||----w |
>                 ||     ||
> ```
{: .comment}


## Facts

In the last step of the last hands on, you ran your playbook. The first `TASK` that was executed was not one that you had written, it read:

```
TASK [Gathering Facts] *************************
ok: [localhost]
```

The [`setup`](https://docs.ansible.com/ansible/latest/modules/setup_module.html) module runs by default for every host and gathers facts about the host.

> ### {% icon hands_on %} Hands-on: The Setup Module
>
> 1. Run the command `ansible -i hosts -c local -m setup my_hosts`.
>
>    The `my_hosts` at the end of the command refers to the group we defined in our `hosts` inventory file.
>
> 2. Investigate the output. See what sort of information is made available to you.
>
>    > ### {% icon question %} Question
>    >
>    > 1. What variable stores the OS name?
>    > 2. Where are all of the places that you can find your machine's IP (v4) addresses?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > 1. The OS name is stored in `ansible_distribution`. We saw `ansible_os_family` used above in the `ansible-cvmfs` role. You can use these variables if you are writing a generic role but packages or commands are named different on different operating systems.
>    > >
>    > > 2. Ansible stores network information in quite a few places, sometimes one place is more convenient or more correct to use:
>    > >
>    > >    - `ansible_all_ipv4_addresses`
>    > >    - `ansible_default_ipv4`
>    > >    - `ansible_<interface_name>.ipv4`
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}


## Templates

Templates give you greater control over the files you are deploying to the remote system. If you need to deploy a file to multiple hosts, but configure it differently on each host, you should use templates. For instance deploying a service that should only listen on the correct IP address for that host would be a good use case for templates. All of the facts you discovered in the previous hands on are available to you to use in templates, `when` statements (like the [ansible-cvmfs example we saw earlier](#modules-and-tasks)). Additionally all of the variables you've defined are available as well.

> ### {% icon details %} Template Syntax
> Templates use Jinja2 syntax. If you are not familiar with it, you should [read about it](http://jinja.pocoo.org/docs/2.10/templates/) first, before moving on with the tutorial.
{: .details}

> ### {% icon hands_on %} Hands-on: Variables and Templates
>
> 1. Create the directories: `roles/my-role/templates`, `roles/my-role/defaults`
>
> 2. Create and edit a file for your role's variables, `roles/my-role/defaults/main.yml`
>
> 3. Insert the following content:
>
>    ```yaml
>    ---
>    server_name: Cats!
>    ```
>
>    This will define a variable `server_name` that can then be used in templates.
>
> 4. Create and edit `roles/my-role/templates/test.ini.j2`
>
> 5. Insert the following content:
>
>    {% raw %}
>    ```ini
>    [example]
>    server_name = {{ server_name }}
>    listen = {{ ansible_default_ipv4.address }}
>    ```
>    {% endraw %}
>
> 6. Edit `roles/my-role/tasks/main.yml` and append a new task to the end to template this file:
>
>    ```yaml
>    - name: Template the configuration file
>      template:
>        src: test.ini.j2
>        dest: /tmp/test.ini
>    ```
>
> 7. Run the playbook again.
>
> 8. Check the contents of `/tmp/test.ini`
>
>    > ### {% icon question %} Question
>    >
>    > How does it look?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The file should look like:
>    > >
>    > > ```ini
>    > > [example]
>    > > server_name = Cats!
>    > > listen = 192.168.0.2
>    > > ```
>    > >
>    > > Where the last line has the machine's IP address.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>    Now that this has worked successfully, we will setup a `group_vars` folder
>    to show how a person using `my-role` would override the `server_name` variable.
>
> 9. Create the folder `group_vars/` (in the root of your directory)
>
> 10. Create and edit `group_vars/my_hosts.yml`
>
> 11. Insert the following:
>
>     ```yaml
>     ---
>     server_name: Dogs!
>     ```
>
> 12. Run the playbook again, but imagine you are worried about this change, and supply the `--check --diff` flag to see what changes are made before committing to make them.
>
>     > ### {% icon question %} Question
>     >
>     > How does the output look?
>     >
>     > > ### {% icon solution %} Solution
>     > >
>     > > ```
>     > > $ ansible-playbook -i hosts playbook.yml -c local --check --diff
>     > > PLAY [my_hosts] ******************************************
>     > > TASK [Gathering Facts] **********************************
>     > > ok: [localhost]
>     > > TASK [my-role : Copy] ***********************************
>     > > ok: [localhost]
>     > > TASK [my-role : Template the configuration file] ********
>     > > --- before: /tmp/test.ini
>     > > +++ after: /home/hxr/.ansible/tmp/ansible-local-1906887dr2u6j8n/tmptx9pdelg/test.ini.j2
>     > > @@ -1,3 +1,3 @@
>     > >  [example]
>     > > -server_name = Cats!
>     > > +server_name = Dogs!
>     > >  listen = 192.168.0.25
>     > > changed: [localhost]
>     > > PLAY RECAP **********************************************
>     > > localhost                  : ok=3    changed=1    unreachable=0    failed=0
>     > > ```
>     > >
>     > > Here you can see that the server_name value will be changed. Despite Ansible reporting `changed=1`, no changes have actually been applied to the system.
>     > {: .solution }
>     {: .question}
>
> 13. Run the playbook again, without the `--check` flag to apply your changes.
{: .hands_on}

> ### {% icon comment %} Ansible Variable Templating
> In this step we use some templated variables. We defined them in a template, but they are also commonly used in group variables file. Our templated variable looked like: {% raw %}`listen = {{ ansible_default_ipv4.address }}`{% endraw %}.
>
> It is common to see things like this in Ansible roles:
>
> ```yaml
> root_dir = /opt/my-app
> config_dir = "{{ root_dir }}/config"
> ```
>
> When Ansible runs:
>
> 1. It collects variables defined in group variables and other places
> 2. The first task for each machine is the [`setup` module](https://docs.ansible.com/ansible/latest/modules/setup_module.html) which gathers facts about the host, which are added to the available variables
> 3. When multiple roles execute in a playbook:
>    1. Their defaults are added to the set of variables (the group variables having precedence over these variables)
>    2. They can also dynamically define more variables which may not be set until that role is run
> 4. Before use (in templates, commands, etc.), variables are resolved to their final value
>
>
{: .comment}

# Ansible Galaxy

Now that you've built a small role, you can imagine that building real roles that manage the full installation of a piece of software are not simple things. Ansible Galaxy is the answer here. Many roles for common administration tasks, and software installation and setup are readily available on Ansible Galaxy.

**Warning**: This will install git on the remote machine.

> ### {% icon hands_on %} Hands-on: Installing a module using ansible-galaxy
>
> 1. Run the command `ansible-galaxy install -p roles/ geerlingguy.git`
>
>    This will install the new role into your `roles` folder, alongside your own role.
>
> 2. Edit your playbook.yml and add the role `geerlingguy.git` at the bottom, after `my-role`
>
> 3. Run the playbook
>
>    > ### {% icon question %} Question
>    >
>    > Did something go wrong?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > Since you have been running the playbook as a non-root user (or at least you should have been!), the step to install a package fails.
>    > > The solution to this is to set `become: true`. Edit your playbook.yml and add `become: true` just below `hosts: my_hosts`.
>    > >
>    > > `become` causes Ansible to attempt to become a different user (using sudo/su/whatever is appropriate), by default this is `root`. If you want to become a different user, just set `become_user`. Beware, the user should be able to privilege escalate without a password prompt. Otherwise when you execute the playbook you should set `--ask-become-pass`, using the privilege escalation password for that host.
>    > >
>    > > > ### {% icon details %} Ansible Become
>    > > > See the [documentation](https://docs.ansible.com/ansible/latest/user_guide/become.html) if you need to control this behaviour differently. `become` can be set either at the task level or the playbook level.
>    > > >
>    > > {: .details}
>    > >
>    > {: .solution }
>    {: .question}
>
{: .hands_on}

## Choosing a Role

Picking the best role for a task from Ansible Galaxy is not always a trivial task. Sometimes there will only be a single role doing what you need. Other times you'll have to choose between 20 different roles that all look more or less the same. Here are some tips to guide you in identifying appropriate and well-written roles:

- The name should match the software you are using (I.e. ignore a role named `stackstorm` when you are trying to set up `rabbitmq`. Ansible Galaxy does not have perfect search.)
- `geerlingguy` wrote a huge number of roles for many pieces of standard software. If there is a role from this person, this is usually a safe choice.
- Check the GitHub readme of each role you consider using. Look for:
  - Extensive documentation of all of the variables, their default values, and how they behave.
  - An example playbook using the role
- A large number of downloads is a good sign that other people found it useful

These are usually good proxies for quality, but do not treat them as strict rules. For an example of a role meeting many of these qualities, [`ansible-cvmfs`](https://github.com/galaxyproject/ansible-cvmfs) is good; the variables are well documented and there are example playbooks that you can (more or less) copy-and-paste and run.

Sometimes a role will accomplish 95% of what you need to do, but not everything. Once you have installed the role with `ansible-galaxy install`, you can edit it locally to make any changes. In an ideal world you would contribute this back, but this is not always a high priority. Many projects copy roles directly into their repositories, e.g. [galaxyproject](https://github.com/galaxyproject/infrastructure-playbook/tree/master/roles) and [usegalaxy.eu](https://github.com/usegalaxy-eu/infrastructure-playbook/tree/master/roles)

# (Optional) Ansible Vault

Now that you have a small role built up, you might start thinking about deploying larger and more complex services and infrastructure. One last common task we want to cover here is the inclusion of secrets.

> ### {% icon hands_on %} Hands-on: Setting up secrets
>
> 0. Run `mkdir -p secret_group_vars`
>
> 1. Now we'll create a strong password: `openssl rand -base64 24 > vault-password.txt`
>
> 2. Run `ansible-vault create secret_group_vars/all.yml --vault-password-file=vault-password.txt`
>
>    This will open `secret_group_vars/all.yml` in your text editor.
>
> 3. In this file, enter the following contents and save it:
>
>    ```yaml
>    apikey: super-secret-api-key-wow!
>    ```
>
>    > ### {% icon question %} Question
>    >
>    > How does your file look? Is it readable? Run `cat secret_group_vars/all.yml`
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The file will look like this, it is encrypted by Ansible Vault with AES256 encryption.
>    > > ```
>    > > $ cat secret_group_vars/all.yml
>    > > $ANSIBLE_VAULT;1.1;AES256
>    > > 64373665366130333437393639343534653134346538636239393363373062393830653333323966
>    > > 3134333366363130326139323162323131643763336236320a393262303938316262643764323862
>    > > 36393161666663353231366336613838633866323230303031313465646333613862363264323139
>    > > 3263383530626262370a666139666462663938343531656432353239346532316630366165376566
>    > > 34313765353766666330366632303836353863396430343264303032363739666139383830323565
>    > > 6133663637356331613062353834646561653366386665623930
>    > > ```
>    > >
>    > {: .solution}
>    {: .question}
>
> 4. Use the new variable in our `.ini` file from earlier. Edit `roles/my-role/templates/test.ini.j2` and add the line `{% raw %}apikey = {{ apikey }}{% endraw %}`
>
> 5. Add encrypted variable file to `playbook.yml`
>
>    ```yaml
>     - hosts: my_hosts
>       ...
>       vars_files:
>        - secret_group_vars/all.yml
>     ```
>
> 6. Tell ansible where to find the decryption file. Create file `ansible.cfg` with content
>
>     ```yaml
>     [default]
>     vault_password_file=vault-password.txt
>     ```
> 7. Run the playbook
>
> 8. Check the contents of `/tmp/test.ini`
>
>    > ### {% icon question %} Question
>    >
>    > How does it look?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > The file should look like:
>    > >
>    > > ```ini
>    > > [example]
>    > > server_name = Cats!
>    > > listen = 192.168.0.2
>    > > apikey = super-secret-api-key-wow!
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}

Ansible Vault is really useful to include encrypted secrets in your playbook repository. In real life scenarios where you are sharing your playbooks publicly, be sure to encrypt all secrets from the start (or fix/remove the git history if you ever did.) If you are storing your vault password in a file, remember to add it to your `.gitignore` (or VCS appropriate file.)

# Other Stuff

Ansible has a huge array of features and we can't cover them all. Some commonly used features are documented below:

## With Items

Duplicating tasks ten times to install ten packages is not efficient, so Ansible provides `with_items` construct

```yaml
- name: Install stuff
  package:
    name: {% raw %}"{{ item }}"{% endraw %}
    state: installed
  with_items:
    - htop
    - git
    - vim
```

This works for any task, not just package installation if you have things you'd like to repeat.

## When Changed

Doing something only when a task is in the "changed" state is a common pattern. This is often used for reloading a service when some configuration files have changed.

```yaml
- name: Copy configuration file
  copy:
    src: main.conf
    dest: "/etc/service/main.conf
    owner: "root"
    group: "root"
    mode: "0644"
  register: service_conf

- name: Restart the service
  service:
    name: service
    enabled: yes
    state: restarted
  when: service_conf.changed
```

## Notifying Handlers

Often you want to restart a service whenever something has changed, like above. But if you need to restart the service whenever one of many different tasks has changed, there is a simpler way than writing `when: a.changed or b.changed or c.changed or ...`

First, move the restarting or reloading of the service into the `handlers/main.yml`

```
- name: restart service
  service:
    name: service
    enabled: yes
    state: restarted

- name: reload httpd
  service:
    name: httpd
    enabled: yes
    state: reloaded
```

Now you can change your task definitions:

```
- name: Copy configuration file
  copy:
    src: main.conf
    dest: "/etc/service/main.conf
    owner: "root"
    group: "root"
    mode: "0644"
  notify: 'restart service'
```

We no longer `register` the command, instead the `notify` attribute automatically marks that the appropriate handler should be called *at the end of the playbook run*. Additionally tasks from outside of this role can use the same `notify`. So if you use a community built role for managing apache, and then have a custom role that builds the apache configuration, you can use changes there to reload the service using their definition of the reload/restart.
