---
layout: tutorial_hands_on

title: "Training Infrastructure as a Service"
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
subtopic: features
tags:
  - training
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---


# Overview
{:.no_toc}


> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

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
> 3. We need to add the `usegalaxy_eu.tiaas2` role to the end of the playbook
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
> 5. Run the playbook
>
> 6. TIaaS should be available now! Check out the following routes on your server:
>
>    URL | Use
>    --- | ---
>    https://server/tiaas/add | Request a new training, send this URLs to trainers who will use your instance
>    https://server/tiaas/admin/ | You can login here with the credentials you set for `tiaas_admin_user` and `tiaas_admin_pass`, and then approve trainings
>    https://server/join-training/<training-id>/ | Join an existing training, it must be status=approved in the admin interface, and you must be logged in to Galaxy
>    https://server/join-training/<training-id>/status | See a training overview, the state of jobs in the queue of trainees.
>
{: .hands_on}
