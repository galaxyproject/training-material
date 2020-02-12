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
  - training
  - jobs
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - connect-to-compute-cluster
      - heterogeneous-compute
---


# Overview
{:.no_toc}


> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introduction

Galaxy is widely used for teaching. In order to facilitate instructors, [Galaxy Europe](https://usegalaxy.eu) has developed Training Infrastructure as a Service (TIaaS).
Workshop instructors can apply for TIaaS, and on the day of their workshop, their participants will be placed in a special group and use dedicated
resources, thus reducing queue times on the day of the training.

![TIaaS idea](../../../instructors/images/workshop-intro/tiaas.png){: width="60%"}

This tutorial will go cover how to set up such a service on your own Galaxy server.


# Setting up TIaaS

> ### {% icon hands_on %} Hands-on: Setup TIaaS
>
> 1. In your `requirements.yml` add the TIaaS ansible role:
>
>    ```yml
>    - src: usegalaxy_eu.tiaas2
>      version: 0.0.1
>    ```
>
>    And run the install step:
>
>    ```
>    ansible-galaxy install -p roles -r requirements.yml
>    ```
>
> 2. In your `galaxyservers` group variables file, add the following:
>
>    ```yml
>    # TIaaS setup
>    tiaas_dir: /opt/tiaas
>    tiaas_user: tiaas
>    tiaas_group: tiaas
>    tiaas_version: master
>    tiaas_admin_user: admin
>    tiaas_admin_pass: changeme
>    ```
>
> 2. In the `galaxyservers` group variables file, we also need to set the database permissions correctly for TIaaS. It needs to be able to access some Galaxy tables, and we will carefully define only the ones we really need:
>
>    ```diff
>    +++ group_vars/galaxyservers.yml
>     postgresql_objects_users:
>       - name: galaxy
>         password: null
>    +  - name: tiaas
>    +    password: null
>     postgresql_objects_databases:
>       - name: galaxy
>         owner: galaxy
>    +postgresql_objects_privileges:
>    +- database: galaxy
>    +  roles: tiaas
>    +  objs: galaxy_user,galaxy_session,job
>    +  type: table
>    +  privs: SELECT
>    +- database: galaxy
>    +  roles: tiaas
>    +  objs: user_group_association,galaxy_group,role,group_role_association
>    +  type: table
>    +  privs: SELECT,INSERT
>    +- database: galaxy
>    +  roles: tiaas
>    +  objs: role_id_seq,galaxy_group_id_seq,group_role_association_id_seq,user_group_association_id_seq
>    +  type: sequence
>    +  privs: USAGE,SELECT
>    ```
>
>
> 3. We need to add the `usegalaxy_eu.tiaas2` role to the end of the playbook (`galaxy.yml`)
>
> 4. Lastly we should add the routes for TIaaS to the NGINX template for Galaxy:
>
>    ```nginx
>
>    location /tiaas {
>        uwsgi_pass 127.0.0.1:5000;
>        uwsgi_param UWSGI_SCHEME $scheme;
>        include uwsgi_params;
>    }
>
>    location /tiaas/static {
>        alias /opt/tiaas/static;
>    }
>
>    location /join-training {
>        uwsgi_pass 127.0.0.1:5000;
>        uwsgi_param UWSGI_SCHEME $scheme;
>        include uwsgi_params;
>    }
>
>    ```
>
> 5. Run the playbook (`ansible-playbook -i hosts galaxy.yml`)
>
> 6. TIaaS should be available now! Check out the following routes on your server:
>
>
>    |URL | Use
>    |----|----
>    |https://server/tiaas/new/ | Request a new training, send this URLs to trainers who will use your instance |
>    |https://server/tiaas/admin/ | You can login here with the values you set for `tiaas_admin_user` and `tiaas_admin_pass` in your group variables file, and then approve the training you requested|
>    |https://server/join-training/\<training-id\>/ | Join an existing training, it must be status=approved in the admin interface, and you must be logged in to Galaxy|
>    |https://server/join-training/\<training-id\>/status | See a training overview, the state of jobs in the queue of trainees.|
>
{: .hands_on}


# Job Configuration

While observability for teachers or trainers is already a huge benefit, one of the primary benefits of TIaaS from UseGalaxy.eu is that your jobs get sent to dedicated compute resources, which won't be used by anyone else, during the period of the training. We will send all of the training jobs to pulsar.


> ### {% icon hands_on %} Hands-on: Writing a dynamic job destination
>
> 1. Create and open `files/galaxy/dynamic_job_rules/hogwarts.py`
>
>    ```python
>    from galaxy.jobs import JobDestination
>    from galaxy.jobs.mapper import JobMappingException
>    import os
>
>    def sorting_hat(app, user):
>        # Check if the user is in a training group
>        user_roles = [role.name for role in user.all_roles() if not role.deleted]
>        if any([role.startswith('training-') for role in user_roles]):
>            # This is a training user, we will send their jobs to pulsar
>            return JobDestination(runner="pulsar")
>        else:
>            return JobDestination(runner="slurm")
>    ```
>
>    This destination will check that the `user_email` is in a training group (role starting with `training-`).
>
> 2. As usual, we need to instruct Galaxy of where to find this file:
>
>    - Edit your group variables file and add the following:
>
>      ```diff
>       galaxy_dynamic_job_rules:
>         - my_rules.py
>      +  - hogwarts.py
>      ```
>
> 3. We next need to configure this plugin in our job configuration (`files/galaxy/config/job_conf.xml`):
>
>    ```xml
>    <destination id="sorting_hat" runner="dynamic">
>        <param id="type">python</param>
>        <param id="function">sorting_hat</param>
>    </destination>
>    ```
>
>    This is a **Python function dynamic destination**. Galaxy will load all python files in the {% raw %}`{{ galaxy_dynamic_rule_dir }}`{% endraw %}, and all functions defined in those will be available to be used in the `job_conf.xml`
>
> 4. Finally, in `job_conf.xml`, update the top level `<destinations>` definition and point it to the sorting hat:
>
>    ```xml
>    <destinations default="sorting_hat">
>    ...
>    </destinations>
>    ```
>
> 5. Run the playbook
>
> 6. Ensure your user is joined to a training
>
> 7. Run a job and observe the logs to see where it goes (`journalctf -u galaxy -f`)
>
{: .hands_on}
