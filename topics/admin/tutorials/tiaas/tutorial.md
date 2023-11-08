---
layout: tutorial_hands_on

title: "Training Infrastructure as a Service (TIaaS)"
zenodo_link: ""
questions:
  - How to deploy EU's TIaaS
objectives:
  - Setup TIaaS
  - Request and manage trainings
  - Join a training
time_estimation: "30m"
key_points:
  - TIaaS is an additional service you can deploy which can help you provide a better service to your users
contributors:
  - hexylena
  - shiltemann
subtopic: features
tags:
  - ansible
  - training
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
      - pulsar
abbreviations:
  TIaaS: Training Infrastructure as a Service
---

Galaxy is widely used for teaching. In order to facilitate instructors, the Galaxy Project has developed {TIaaS}.
Workshop instructors can apply for {TIaaS}, and on the day of their workshop, their participants will be placed in a special group and use dedicated
resources, thus reducing queue times on the day of the training.

![TIaaS concept](../../images/tiaas/tiaas_intro.png "With TIaaS, all of your users visit the same server. In the background, the scheduler recognises which users are training users, and directs their jobs to special resources. In the EU deployment of TIaaS jobs preferentially use private resources, but can spill over to the main queue if there is not enough space available."){: width="70%"}

This tutorial will go cover how to set up such a service on your own Galaxy server.


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="tiaas" %}


# Setting up TIaaS

> <hands-on-title>Setup TIaaS</hands-on-title>
>
> 1. In your `requirements.yml` add the {TIaaS} ansible role:
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -51,3 +51,6 @@
>       version: 6f6fdf7f5ead491560783d52528b79e9e088bd5b
>     - src: cloudalchemy.grafana
>       version: 0.14.2
>    +# Training Infrastructure as a Service
>    +- src: galaxyproject.tiaas2
>    +  version: 2.1.5
>    {% endraw %}
>    ```
>    {: data-commit="Add tiaas2 requirement"}
>
>    And run the install step:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. In your `galaxyservers` group variables file, add the following:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -346,3 +346,8 @@ telegraf_plugins_extra:
>           - timeout = "10s"
>           - data_format = "influx"
>           - interval = "15s"
>    +
>    +# TIaaS setup
>    +tiaas_dir: /srv/tiaas
>    +tiaas_admin_user: admin
>    +tiaas_admin_pass: changeme
>    {% endraw %}
>    ```
>    {: data-commit="Configure tiaas"}
>
> 2. In the `galaxyservers` group variables file, we also need to set the database permissions correctly for TIaaS. It needs to be able to access some Galaxy tables, and we will carefully define only the ones we really need:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/dbservers.yml
>    +++ b/group_vars/dbservers.yml
>    @@ -3,6 +3,7 @@
>     postgresql_objects_users:
>       - name: "{{ galaxy_user_name }}"
>       - name: telegraf
>    +  - name: tiaas
>     postgresql_objects_databases:
>       - name: "{{ galaxy_db_name }}"
>         owner: "{{ galaxy_user_name }}"
>    @@ -11,7 +12,26 @@ postgresql_objects_privileges:
>         roles: telegraf
>         privs: SELECT
>         objs: ALL_IN_SCHEMA
>    -
>    +  - database: galaxy
>    +    roles: tiaas
>    +    objs: galaxy_user,galaxy_session,job,history,workflow,workflow_invocation
>    +    type: table
>    +    privs: SELECT
>    +  - database: galaxy
>    +    roles: tiaas
>    +    objs: user_group_association,galaxy_group,role,group_role_association
>    +    type: table
>    +    privs: SELECT,INSERT
>    +  - database: galaxy
>    +    roles: tiaas
>    +    objs: group_role_association
>    +    type: table
>    +    privs: DELETE
>    +  - database: galaxy
>    +    roles: tiaas
>    +    objs: role_id_seq,galaxy_group_id_seq,group_role_association_id_seq,user_group_association_id_seq
>    +    type: sequence
>    +    privs: USAGE,SELECT
>     
>     # PostgreSQL Backups
>     postgresql_backup_dir: /data/backups
>    {% endraw %}
>    ```
>    {: data-commit="Add database privileges for TIaaS"}
>
>    > <tip-title>Why does TIaaS get `DELETE` privileges on Galaxy's Database?</tip-title>
>    > The `DELETE` privilege is limited in scope to one table: `group_role_association`. This allows TIaaS to
>    > disassociate training groups from roles in the Galaxy database after the training event date has passed, so that
>    > users who participated in a training return to using normal (non-training) resources after the training ends.
>    >
>    > The `galaxyproject.tiaas2` role will create a [cron](https://manpages.debian.org/stable/cron/cron.8.en.html) job
>    > to perform this process every night at midnight. You can control when this runs (or disable it) using
>    > [the tiaas_disassociate_training_roles variable](https://github.com/galaxyproject/ansible-tiaas2/blob/d5be2a064c49e010f67bfcea18e36812da23d7d8/defaults/main.yml#L20).
>    >
>    {: .tip}
>
>    > <tip-title>Running the playbook from scratch</tip-title>
>    > This is one of the few statements we've provided that presents difficulties when running the playbook completely from scratch on a blank machine. Setting postgresql roles is one of the first steps in our playbook, but the rules we've provided above depend on the Galaxy tables existing in that database. If those tables aren't there, it will fail. If you do someday run this from scratch, you'll find that you need to comment out those roles.
>    {: .tip}
>
> 3. We need to add the `galaxyproject.tiaas2` role before the `nginx` role, as TIaaS defines variables that Nginx needs.
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -47,6 +47,7 @@
>         - galaxyproject.nginx
>         - geerlingguy.docker
>         - usegalaxy_eu.rabbitmqserver
>    +    - galaxyproject.tiaas2
>         - galaxyproject.gxadmin
>         - galaxyproject.cvmfs
>         - dj-wasabi.telegraf
>    {% endraw %}
>    ```
>    {: data-commit="Add TIaaS role to the Galaxy playbook"}
>
> 4. Lastly we should add the routes for TIaaS to the NGINX template for Galaxy. TIaaS provides a set of default nginx routes that can be used.
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -113,4 +113,6 @@ server {
>     		proxy_set_header Host $http_host;
>     	}
>     
>    +	{{ tiaas_nginx_routes }}
>    +
>     }
>    {% endraw %}
>    ```
>    {: data-commit="Add nginx routes for TIaaS"}
>
> 5. Run the playbook
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
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}


{TIaaS} should be available now! The following routes on your server are now configured (we will run through these in the next section)


| URL                                                                                       | Use                                                                         | Audience                 |
| :----                                                                                     | ----                                                                        | -----                    |
| [/tiaas/new/](https://my.gat.galaxy.training/?path=/tiaas/new/)                           | Request a new TIaaS training                                                | Instructors              |
| [/tiaas/admin/](https://my.gat.galaxy.training/?path=/tiaas/admin/)                       | Approve and Manage requests                                                 | Admin                    |
| [/tiaas/stats/](https://my.gat.galaxy.training/?path=/tiaas/stats/)                       | Overall TIaaS statistics ([EU Stats](https://usegalaxy.eu/tiaas/stats/))    | Admins, Funding Agencies |
| [/tiaas/calendar/](https://my.gat.galaxy.training/?path=/tiaas/calendar/)                 | Calendar of trainings ([EU Calendar](https://usegalaxy.eu/tiaas/calendar/)) | Admins, Funding Agencies |
| [/join-training/ID](https://my.gat.galaxy.training/?path=/join-training/ID)               | Join an TIaaS training                                                      | Participants             |
| [/join-training/ID/status](https://my.gat.galaxy.training/?path=/join-training/ID/status) | Dashboard with job states of trainees.                                      | Instructors              |


Let's see it in action!

> <hands-on-title>Using TIaaS</hands-on-title>
>
> 1. **Create a new TIaaS request**
>    - Go to [/tiaas/new](https://my.gat.galaxy.training/?path=/tiaas/new)
>    - Here you will find the request form users will fill in to request TIaaS:
>      ![TIaaS request form](../../images/tiaas/tiaas_request_form.png)
>    - For *"Training Identifier"*, fill in `gat`
>      - This is the `<training-id>` used in the URLs listed above used for:
>        1. Workshop participants to join the tiaas group
>        2. Workshop instructors to monitor the progress of their participants.
>    - Fill in the rest of the form as you like
>    - Submit the form and you should see a confirmation dialog:
>      ![TIaaS requested successfully](../../images/tiaas/tiaas_request_form_submitted.png)
>
> 2. **Approve TIaaS request**
>    - Next, the request will have to be approved by an admin
>    - Go to [/tiaas/admin](https://my.gat.galaxy.training/?path=/tiaas/admin)
>    - **Log in** using the values you configured `tiaas_admin_user` and `tiaas_admin_pass` in your group variables file
>      - Default values were `admin:changeme`
>    - You should now see the admin panel:
>      ![TIaaS admin console](../../images/tiaas/tiaas_admin_console.png)
>    - Click on **Trainings**, you should see the TIaaS request listed here:
>      ![TIaaS request list](../../images/tiaas/tiaas_request_list.png)
>    - **Approve the request**
>      - Click on the training
>      - Scroll down to the bottom
>      - Change *"Processed"* to `Approved` and **Save**
>        ![Approve TIaaS](../../images/tiaas/tiaas_request_approve.png)
>    - At this point, you would likely email the person who made the request to inform them of approval
>
> 3. **Join TIaaS Training**
>    - Make sure you are logged in to Galaxy
>    - On the day of the workshop, participants will visit a following URL to join the TIaaS group
>      - [/join-training/gat](https://my.gat.galaxy.training/?path=/join-training/gat)
>      - A confirmation dialog should appear if all went well:
>        ![Join TIaaS](../../images/tiaas/tiaas_join_training.png)
>
> 4. **Monitor TIaaS status**
>    - This is very useful for instructors to monitor the job state of their participants
>    - Go to [/join-training/gat](https://my.gat.galaxy.training/?path=/join-training/gat)
>    - In the Dasboard you should see that one user (you) has joined the training \
>    - Run some jobs to see the dashboard in action
>      ![TIaaS dashboard](../../images/tiaas/tiaas_dashboard.png)
>    - Scroll down to get some more information on a per-user level (anonymized)
>      - Every user designated by their own identifier and colour, but no personal information
>      ![TIaaS dashboard](../../images/tiaas/tiaas_dashboard2.png)
>
{: .hands_on}


> <comment-title>Note: GDPR assistance</comment-title>
>
> Since this setup tracks additional personal information (submitter name & email, users in the queue view), TIaaS includes some always-on features to assist with your GDPR compliance.
>
>  - Users in public status dashboard are only visible by an anonymized identifier and colour
>  - Email addressses in the TIaaS admin panel will be automatically expunged 60 days after a training event
>
>  Of course you need to review any GDPR compliance concerns with your group's legal representative(s), this only attempts to ensure some protections exist for the users of the system.
>
{: .comment}


# Job Configuration

While observability for teachers or trainers is already a huge benefit, one of the primary benefits of {TIaaS} is that your jobs get sent to dedicated compute resources, which won't be used by anyone else, during the period of the training. We will send all of the training jobs to pulsar if you have completed that tutorial, or one of the slurm destinations from the job configuration training.

In order to achieve this, we first need some way to sort the jobs of the training users into these private queues, while letting the other jobs continue on. So let's create a *traffic controller* to figure out where jobs belong.


> <hands-on-title>Writing a dynamic job destination</hands-on-title>
>
> 1. This destination will check that the `user_email` is in a training group (role starting with `training-`).
>
>    {% raw %}
>    ```diff
>    --- a/files/galaxy/config/tpv_rules_local.yml
>    +++ b/files/galaxy/config/tpv_rules_local.yml
>    @@ -35,6 +35,15 @@ tools:
>           require:
>             - pulsar
>     
>    +roles:
>    +  training.*:
>    +    max_cores: 2
>    +    max_mem: max_cores * 3.8  # TODO check multiplier
>    +    scheduling:
>    +      require:
>    +        - slurm
>    +        - training
>    +
>     destinations:
>       local_env:
>         runner: local_runner
>    @@ -62,6 +71,19 @@ destinations:
>         max_mem: 8
>         params:
>           native_specification: --nodes=1 --ntasks=1 --cpus-per-task={cores} --time={params['walltime']}:00:00
>    +  slurm-training:
>    +    inherits: singularity
>    +    runner: slurm
>    +    max_accepted_cores: 12
>    +    max_accepted_mem: 120
>    +    max_cores: 2 # Limit the cores
>    +    max_mem: 8 # Limit the memory
>    +    params:
>    +      native_specification: --nodes=1 --ntasks=1 --mem={round(mem*1024)} --cpus-per-task={cores} --time=00:30:00
>    +    scheduling:
>    +      require:
>    +        - slurm
>    +        - training
>     
>       pulsar:
>         runner: pulsar_runner
>    {% endraw %}
>    ```
>    {: data-commit="Add to list of deployed rules"}
>
> 6. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 7. Ensure your user is joined to a training
>
> 8. Run a job and observe the logs to see where it goes (`journalctl -u galaxy -f`)
>
{: .hands_on}

Congratulations! you have now set up {TIaaS} on your Galaxy server.

> ```bash
> 2.sh
> ```
> {: data-test="true"}
{: .hidden}

{% snippet topics/admin/faqs/git-commit.md page=page %}

{% snippet topics/admin/faqs/missed-something.md step=14 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="tiaas" %}
