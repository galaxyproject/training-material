---
layout: tutorial_hands_on

title: "Galaxy Backups & Cleanup"
questions:
- How can I back up my Galaxy?
- What data should be included?
- How can I ensure jobs get cleaned up appropriately?
objectives:
- Setup postgres backups
- Setup cleanups
time_estimation: "30m"
key_points:
- Most of this is handled for you
- But it is important to be aware of proper backup procedures
contributions:
  authorship:
  - hexylena
  - natefoo
tags:
  - ansible
  - deploying
  - git-gat
subtopic: maintenance
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible-galaxy
  - type: "none"
    title: "A VM with at least 2 vCPUs and 4 GB RAM, preferably running Ubuntu 18.04 - 20.04."
abbreviations:
    WORM: Write Once Read Many
---

Keeping your Galaxy cleaned up is an important way to retain space, especially since for many groups that is the
limiting factor in their deployment.

Additionally, backups are necessary to ensure that if you ever experience system level failures, you can safely recover
from these.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="backup-cleanup" %}

# Cleanups

There are two kinds of data that are produced when running a Galaxy: files users create and then delete or purge, and
then files Galaxy creates itself. Both of these can be cleaned to save space.

## User Created Files

You can use `gxadmin` to cleanup user created files. `gxadmin` is covered in more detail in [it's own dedicated
tutorial]({% link topics/admin/tutorials/gxadmin/tutorial.md %}).

> <hands-on-title>Installing gxadmin with Ansible</hands-on-title>
>
> 1. Edit your `requirements.yml` and add the following:
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -28,3 +28,5 @@
>       version: 0.1.0
>     - src: galaxyproject.pulsar
>       version: 1.0.8
>    +- src: galaxyproject.gxadmin
>    +  version: 0.0.8
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement"}
>
> 2. Install the role with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Add the role to your playbook:
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -39,3 +39,4 @@
>         - galaxyproject.nginx
>         - galaxyproject.tusd
>         - galaxyproject.cvmfs
>    +    - galaxyproject.gxadmin
>    {% endraw %}
>    ```
>    {: data-commit="Add the gxadmin role"}
>
> 3. Setup a cleanup task to run regularly:
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -39,3 +39,4 @@
>    post_tasks:
>      - name: Setup gxadmin cleanup task
>        ansible.builtin.cron:
>          name: "Cleanup Old User Data"
>          user: galaxy # Run as the Galaxy user
>          minute: "0"
>          hour: "0"
>          job: "GALAXY_LOG_DIR=/tmp/gxadmin/ GALAXY_ROOT={{ galaxy_root }}/server /usr/bin/gxadmin galaxy cleanup 60"
>    {% endraw %}
>    ```
>    {: data-commit="Configure gxadmin to cleanup data"}
>
>    This will cause datasets deleted for more than 60 days to be purged.
>
> 4. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on}

Whenever `gxadmin` runs, it will create logs you can read in `/tmp/gxadmin` which you can check later.

## Galaxy Created Files

Before we begin backing up our Galaxy data, let's set up automated cleanups to ensure we backup the minimal required set of data.

> <hands-on-title>Configuring PostgreSQL Backups</hands-on-title>
>
> 1. Edit `galaxy.yml` to install `tmpwatch` (if using RHEL/CentOS/Rocky) and `tmpreaper` if using Debian/Ubuntu
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -8,6 +8,12 @@
>         - name: Install Dependencies
>           package:
>             name: ['acl', 'bzip2', 'git', 'make', 'python3-psycopg2', 'tar', 'virtualenv']
>    +    - name: Install RHEL/CentOS/Rocky specific dependencies
>    +      package:
>    +        name: ['tmpwatch']
>    +    - name: Install Debian/Ubuntu specific dependencies
>    +      package:
>    +        name: ['tmpreaper']
>       roles:
>         - galaxyproject.postgresql
>         - role: galaxyproject.postgresql_objects
>    {% endraw %}
>    ```
>    {: data-commit="Install tmpwatch/tmpreaper"}
>
> 1. Edit `group_vars/galaxyservers.yml` and add some variables to configure PostgreSQL:
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -15,6 +15,7 @@ postgresql_objects_databases:
>     galaxy_create_user: true
>     galaxy_separate_privileges: true
>     galaxy_manage_paths: true
>    +galaxy_manage_cleanup: true
>     galaxy_layout: root-dir
>     galaxy_root: /srv/galaxy
>     galaxy_user: {name: galaxy, shell: /bin/bash}
>    {% endraw %}
>    ```
>    {: data-commit="Configure automated cleanup"}
>
> 2. > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Check out the cleanup task which has been generated in: `/etc/cron.d/ansible_galaxy_tmpclean`
>
{: .hands_on}

This will setup `tmpwatch` to cleanup a few folders:

- the job working directory, important if you set `cleanup: onsuccess`, to cleanup old failed jobs once you're done debugging their failures.
- the new file upload path, to catch uploaded temporary files that are no longer necessary.

# Backups

There are two important things to back up with your Ansible Galaxy:

- Galaxy
- The Database
- The Data

## Galaxy

By using Ansible, as long as you are storing your playbooks on another system, you are generally safe from failues of
the Galaxy node, and you'll be able to re-run your playbook at a later date.

However, playbooks often do not include:

- Which tools you've installed (have you ever installed a tool outside of ephemeris? This might be lost!)
- Conda environments, which will not always resolve identically over time. If strong guarantees of reproducibility are
  important, then consider backing these up as well.

## Database Backups

We're setting a couple of variables to control the automatic backups, they'll be placed in the `/data/backups` folder next to our user uploaded Galaxy data.

> <hands-on-title>Configuring PostgreSQL Backups</hands-on-title>
>
> 1. Edit `group_vars/galaxyservers.yml` and add some variables to configure PostgreSQL:
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -11,6 +11,10 @@ postgresql_objects_databases:
>       - name: galaxy
>         owner: galaxy
>     
>    +# PostgreSQL Backups
>    +postgresql_backup_dir: /data/backups
>    +postgresql_backup_local_dir: "{{ '~postgres' | expanduser }}/backups"
>    +
>     # Galaxy
>     galaxy_create_user: true
>     galaxy_separate_privileges: true
>    {% endraw %}
>    ```
>    {: data-commit="Add backups"}
>
>    <!-- TODO: expanduser tip here -->
{: .hands_on}

This will setup our backups to run as a cron job. 

<!-- TODO: explore cron jobs -->

## Data Backup

With Galaxy it is *technically* only necessary to backup your inputs, as the downstream files *should*, *in theory* be re-createable due to the reproducibility of Galaxy.

In practice, some groups either choose to not backup, or to backup everything, often to extremely cheap and slow storage like Glacier or a tape library.

Most groups choose to implement this as a custom cron job, e.g.

```yaml
post_tasks:
  - name: Setup backup cron job
    ansible.builtin.cron:
      name: "Backup User Data"
      minute: "0"
      hour: "5,2"
      job: "rsync -avr /data/galaxy/ backup@backup.example.org:/backups/$(date -I)/"
```

> <tip-title>This isn't a backup strategy!</tip-title>
> People who, let's say, care strongly about backups will often insist that you need to version files. This is of course unnecessary in the Galaxy case as files are essentially {WORM}s, which is a really good file storage practice. Files can get removed so it isn't a true {WORM} strategy that you'd use for e.g. audit logs, but it is close.
> That said, since files never get changed, keeping multiple versions is unnecesary.
>
> Please consider communicating **very well** with your users what the data backup policy is.
{: .tip}
