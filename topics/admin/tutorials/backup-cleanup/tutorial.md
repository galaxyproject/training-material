---
layout: tutorial_hands_on

title: "Server Maintenance: Cleanup, Backup, and Restoration"
questions:
- How can I back up my Galaxy?
- What data should be included?
- How can I ensure jobs get cleaned up appropriately?
- How do I maintain a Galaxy server?
- What happens if I lose everything?
objectives:
- Learn about different maintenance steps
- Setup postgres backups
- Setup cleanups
- Learn what to back up and how to recover
time_estimation: "30m"
key_points:
  - Use configuration management (e.g. Ansible)
  - Store configuration management in git
  - Back up the parts of Galaxy that can't be recreated
contributions:
  authorship:
  - hexylena
  - lldelisle
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

edam_ontology:
- topic_3489 # Database Management
- topic_0605 # Informatics
- topic_3071 # Data Management
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

You can use `gxadmin` to cleanup user created files. `gxadmin` is covered in more detail in [its own dedicated
tutorial]({% link topics/admin/tutorials/gxadmin/tutorial.md %}).

> <hands-on-title>Installing gxadmin with Ansible</hands-on-title>
>
> 1. Edit your `requirements.yml` and add the following:
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -11,3 +11,6 @@
>       version: 0.3.1
>     - src: usegalaxy_eu.certbot
>       version: 0.1.11
>    +# gxadmin (used in cleanup, and later monitoring.)
>    +- src: galaxyproject.gxadmin
>    +  version: 0.0.12
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
>    @@ -27,3 +27,4 @@
>           become: true
>           become_user: "{{ galaxy_user_name }}"
>         - galaxyproject.nginx
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
>    @@ -28,3 +28,11 @@
>           become_user: "{{ galaxy_user_name }}"
>         - galaxyproject.nginx
>         - galaxyproject.gxadmin
>    +  post_tasks:
>    +    - name: Setup gxadmin cleanup task
>    +      ansible.builtin.cron:
>    +        name: "Cleanup Old User Data"
>    +        user: galaxy # Run as the Galaxy user
>    +        minute: "0"
>    +        hour: "0"
>    +        job: "SHELL=/bin/bash source {{ galaxy_venv_dir }}/bin/activate &&  GALAXY_LOG_DIR=/tmp/gxadmin/ GALAXY_ROOT={{ galaxy_root }}/server GALAXY_CONFIG_FILE={{ galaxy_config_file }} /usr/local/bin/gxadmin galaxy cleanup 60"
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
>    @@ -21,6 +21,14 @@
>         - name: Install Dependencies
>           package:
>             name: ['acl', 'bzip2', 'git', 'make', 'tar', 'python3-venv', 'python3-setuptools']
>    +    - name: Install RHEL/CentOS/Rocky specific dependencies
>    +      package:
>    +        name: ['tmpwatch']
>    +      when: ansible_os_family == 'RedHat'
>    +    - name: Install Debian/Ubuntu specific dependencies
>    +      package:
>    +        name: ['tmpreaper']
>    +      when: ansible_os_family == 'Debian'
>       roles:
>         - galaxyproject.galaxy
>         - role: galaxyproject.miniconda
>    {% endraw %}
>    ```
>    {: data-commit="Install tmpwatch/tmpreaper"}
>
> 1. Edit `group_vars/galaxyservers.yml` and add some variables to configure PostgreSQL:
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -2,6 +2,7 @@
>     galaxy_create_user: true # False by default, as e.g. you might have a 'galaxy' user provided by LDAP or AD.
>     galaxy_separate_privileges: true # Best practices for security, configuration is owned by 'root' (or a different user) than the processes
>     galaxy_manage_paths: true # False by default as your administrator might e.g. have root_squash enabled on NFS. Here we can create the directories so it's fine.
>    +galaxy_manage_cleanup: true
>     galaxy_layout: root-dir
>     galaxy_root: /srv/galaxy
>     galaxy_user: {name: "{{ galaxy_user_name }}", shell: /bin/bash}
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

There are a few important things to back up with your Ansible Galaxy:

- Galaxy
	- The Galaxy-managed config files
	- The playbooks
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
>    --- a/group_vars/dbservers.yml
>    +++ b/group_vars/dbservers.yml
>    @@ -5,3 +5,7 @@ postgresql_objects_users:
>     postgresql_objects_databases:
>       - name: "{{ galaxy_db_name }}"
>         owner: "{{ galaxy_user_name }}"
>    +
>    +# PostgreSQL Backups
>    +postgresql_backup_dir: /data/backups
>    +postgresql_backup_local_dir: "{{ '~postgres' | expanduser }}/backups"
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

{% snippet topics/admin/faqs/missed-something.md step=2 %}

# Restoration

Sometimes failures happen! We're sorry you have to read this section.

## Restoring the Database

This procedure is more complicated, you can [read about the restoration procedure](https://github.com/galaxyproject/ansible-postgresql/pull/30#issuecomment-963600656) in the associated PR.

This step assumes you have pre-existing backups in place, you must check this first:

```console
ls /data/backups/
```

If you have backups, you're ready to restore:

```console
# Stop Galaxy, you do NOT want galaxy to connect mid-restoration in case it
# tries to modify the database.
sudo systemctl stop galaxy

# Stop the database
sudo systemctl stop postgresql
# Ensure that it is stopped
sudo systemctl status postgresql

# Begin the backup procedure by becoming postgres:
sudo su - postgres

# Move the current, live database to a backup location just in case:
mkdir /tmp/test/

# ====
# NOTE THAT THIS NUMBER MAY BE DIFFERENT FOR YOU!
# You will need to change 12 to whatever version of postgres you're running
# in every subsequent command
# ====
mv /var/lib/postgresql/12/main/* /tmp/test/

# Add backup
rsync -av /data/backups/YOUR_LATEST_BACKUP/ /var/lib/postgresql/12/main
# Add the restore_command, to your backup file:
# restore_command = 'cp "/tmp/backup/current/wal/%f" "%p"'
$EDITOR ./12/main/postgresql.auto.conf

# Touch a recovery file
touch /var/lib/postgresql/12/main/recovery.signal

# As $username (with sudo right)
sudo systemctl restart postgresql
sudo systemctl status postgresql
# Restart Galaxy
sudo systemctl start galaxy
```

If you encounter issues, we suggest reading [Lucille's log of her experiences restoring](https://github.com/lldelisle/galaxyduboule-infrastructure/blob/778e5b056f03970a1af5fc2f9dd133b4e4791cac/testUpgradeOnVM.sh#L53-L130) as you might encounter similar issues.

## Restoring Galaxy

Restoring Galaxy is easy via Ansible (maybe ensuring users cannot login by disabling the routes in nginx)

```console
ansible-playbook galaxy.yml
```

And if you are following best practices, you probably have your tools stored in a YAML file to use with [Ephemeris]({% link topics/admin/tutorials/tool-management/tutorial.md %}):

```console
shed-tools install -g https://galaxy.example.org -a <api-key> -t our_tools.yml
```

## Restoring User Data

This should simply be `rsync`ing your data from the backup location back into `/data/galaxy`.

{% snippet topics/admin/faqs/git-gat-path.md tutorial="backup-cleanup" %}
