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



> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="backup-cleanup" %}

# Cleanups

Before we begin backing up our Galaxy data, let's set up automated cleanups to ensure we backup the minimal required set of data.

> <hands-on-title>Configuring PostgreSQL Backups</hands-on-title>
>
> 1. Edit `group_vars/galaxyservers.yml` and add some variables to configure PostgreSQL:
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/galaxyservers.yml
>    @@ -0,0 +1,15 @@
>    +pretasks: install tmpwatch (rhel) tmpreaper (deb)
>    +galaxy_manage_cleanup: yes
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

- The Database
- The Data

## Database Backups

We're setting a couple of variables to control the automatic backups, they'll be placed in the `/data/backups` folder next to our user uploaded Galaxy data.

> <hands-on-title>Configuring PostgreSQL Backups</hands-on-title>
>
> 1. Edit `group_vars/galaxyservers.yml` and add some variables to configure PostgreSQL:
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/galaxyservers.yml
>    @@ -0,0 +1,15 @@
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


