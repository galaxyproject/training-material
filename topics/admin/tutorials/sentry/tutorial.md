---
layout: tutorial_hands_on

title: "Monitoring Galaxy and Pulsar with Sentry"
zenodo_link: ""
questions:
objectives:
  - Have an understanding of Sentry
  - Install Sentry
  - Configure Galaxy and Pulsar to send errors to Sentry
  - Monitor performance with Sentry
time_estimation: "1h"
key_points:
contributions:
  authorship:
  - mvdbeek
  editing:
  - hexylena
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - pulsar
subtopic: data
tags:
  - ansible
  - git-gat
---

# Overview

Sentry is an error tracking software that helps admins and developers monitor and diagnose issues in their applications. It provides real-time alerts for errors and allows users to capture context information about each error, such as stack traces and user feedback. It is often possible to find and fix errors before users report them. Galaxy and Pulsar can log issues and failing tool runs to Sentry.


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="sentry" %}

We're going to set up a local Sentry instance using docker-compose and connect Galaxy and Pulsar to that Sentry instance. Alternatively, you can use the hosted Sentry at https://sentry.io/.

# Installing and Configuring

To proceed from here it is expected that:

> <comment-title>Requirements for Running This Tutorial</comment-title>
>
> 1. You have set up a working Galaxy instance as described in the [ansible-galaxy](../ansible-galaxy/tutorial.html) tutorial.
>

{: .comment}

# Installing and Configuring

First we need to add our new Ansible role to `requirements.yml`:

> <hands-on-title>Set up Sentry with Ansible</hands-on-title>
>
> 1. In your working directory, add the roles to your `requirements.yml`
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -48,3 +48,5 @@
>       src: https://github.com/Paprikant/ansible-role-beacon
>     - name: paprikant.beacon-importer
>       src: https://github.com/Paprikant/ansible-role-beacon_importer
>    +- name: mvdbeek.sentry_selfhosted
>    +  src: https://github.com/mvdbeek/ansible-role-sentry/archive/main.tar.gz
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement" data-ref="add-req"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Install the roles with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true" data-ref="req-install"}
>    {: .code-in}
>
> 3. Create a new playbook, `sentry.yml` with the following:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/sentry.yml
>    @@ -0,0 +1,7 @@
>    +- hosts: sentryservers
>    +  become: true
>    +  pre_tasks:
>    +    - pip:
>    +        name: docker-compose
>    +  roles:
>    +    - mvdbeek.sentry_selfhosted
>    {% endraw %}
>    ```
>    {: data-commit="Setup the sentry playbook"}
>
>    During this tutorial we will install everything on the same host, but often one keeps the monitoring infrastructure (Grafana, InfluxDB, Sentry) on a separate host.
>
> 4. Edit the inventory file (`hosts`) an add a group for sentry like:
>
>    {% raw %}
>    ```diff
>    --- a/hosts
>    +++ b/hosts
>    @@ -15,3 +15,6 @@ beacon_server
>     gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>     [beacon_import]
>     gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +
>    +[sentryservers]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    {% endraw %}
>    ```
>    {: data-commit="Add the monitoring host"}
>
>    **Ensure that the hostname is the full hostname of your machine.**
>
>    Sentry requires its own (sub)domain. For the admin training we have set up the sentry.gat-N.eu.galaxy.training subdomain. If you run this tutorial outside of the training and you cannot obtain a domain or subdomain for sentry you can use the free [Duck DNS](https://www.duckdns.org/) service to map an IP address to a domain name.
>
> 5. Edit the file `group_vars/sentryservers.yml` and set the following variables:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/sentryservers.yml
>    @@ -0,0 +1,7 @@
>    +sentry_version: 23.3.1
>    +sentry_url: "https://{{ sentry_domain }}"
>    +sentry_docker_compose_project_folder: /srv/sentry
>    +sentry_superusers:
>    +  - email:  admin@example.com
>    +    password: "{{ vault_sentry_password }}"
>    {% endraw %}
>    ```
>    {: data-commit="Configure Sentry"}
>
> 6. Add the nginx routes
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -125,3 +125,24 @@ server {
>                    proxy_set_header Host $host;
>            }
>     }
>    +
>    +server {
>    + # Listen on port 443
>    + listen        *:443 ssl;
>    + # The virtualhost is our domain name
>    + server_name   "{{ sentry_domain }}";
>    +
>    + # Our log files will go here.
>    + access_log  syslog:server=unix:/dev/log;
>    + error_log   syslog:server=unix:/dev/log;
>    +
>    + location / {
>    +     # This is the backend to send the requests to.
>    +     proxy_pass "http://localhost:9000";
>    +
>    +     proxy_set_header Host $http_host;
>    +     proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
>    +     proxy_set_header X-Forwarded-Proto $scheme;
>    +     proxy_set_header Upgrade $http_upgrade;
>    +  }
>    +}
>    {% endraw %}
>    ```
>    {: data-commit="Add nginx server "}
>
> 7. Run the sentry playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook sentry.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 8. Generate a project for Galaxy in Sentry
>  Go to the domain you configured for your Sentry instance. You need to log in with the username and admin you've set up in `group_vars/sentryservers.yml`. Click "continue" on the next page. Click "Projects", "Create Project", "Python", select "I'll create my own alerts later", and set "galaxy" as the Project Name. You'll see your project dsn that will look like `https://b0022427ee5345a8ad4cb072c73e62f4@sentry.gat-N.eu.galaxy.training/2`. We will need this string to let Galaxy know where to send data to. To avoid requesting an additional certificate for communication between Galaxy and Sentry we've set up communication via localhost:9000, so you can manually change the @ portion to localhost:9000.
>
> 9. We will add the project dsn to the vault. Edit your `group_vars/secret.yml` and add the sentry dsn.
>
>    ><code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault edit group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    ```yaml
>    vault_sentry_dsn: 'https://b0022427ee5345a8ad4cb072c73e62f4@localhost:9000/2'
>    ```
>
> 9. Edit `group_vars/galaxyservers.yml` to reference the new vault secret:
>
>    This will let Galaxy know that captured logs should be sent to our Sentry instance.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -106,6 +106,7 @@ galaxy_config:
>         # Monitoring
>         statsd_host: localhost
>         statsd_influxdb: true
>    +    sentry_dsn: "{{ vault_sentry_dsn }}"
>         # FTP
>         ftp_upload_dir: /data/uploads
>         ftp_upload_site: "{{ inventory_hostname }}"
>     ```
>    {% endraw %}
>    {: data-commit="Configure Galaxy to report to Sentry"}
>
> 10. Run the galaxy playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>

# Monitor Sentry for Errors

Let's see what kind of errors Sentry is logging right now:

> <hands-on-title>Open Galaxy Project in Sentry</hands-on-title>
> 1. Go to your Sentry instance and click on issues. ...
{: .hands_on}

{% snippet topics/admin/faqs/missed-something.md step=12 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="sentry" %}
