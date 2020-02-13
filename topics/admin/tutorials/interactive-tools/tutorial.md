---
layout: tutorial_hands_on

title: "Galaxy Interactive Tools"
zenodo_link: ""
questions:
objectives:
- Have an understanding of what Galaxy Interactive Tools are and how they work
- Configure your Galaxy to serve Interactive Tools using an Ansible Playbook
time_estimation: "1h"
key_points:
contributors:
  - natefoo
  - slugger70
  - hexylena
tags:
  - ansible
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

> ### {% icon warning %} Evolving Topic
> Galaxy Interactive Tools are a **new feature** and there are some rough edges and ongoing work to improve the experience of deploying and using them.
{: .warning}

# Overview
{:.no_toc}

Galaxy Interactive Tools (GxITs) are a method to run containerized tools that are interactive in nature. Interactive Tools typically run a persistent service accessed on a specific port and run until terminated by the user. One common example of such a tool is [Jupyter Notebook][jupyter]. Galaxy Interactive Tools are similar in purpose to [Galaxy Interactive Environments][gie-docs] (GIEs), but are implemented in a significantly different manner. Most notably, instead of directly invoking containers on the Galaxy server, dedicated Docker node, or as a Docker Swarm service (as is done for GIEs), Interactive Tools are submitted through Galaxy's job management system and thus are scheduled the same as any other Galaxy tool - on a Slurm cluster, for instance. Galaxy Interactive Tools were introduced in Galaxy Release 19.09.

There are two sections to this exercise. The first shows you how to use Ansible to setup and configure Galaxy Interactive Tools. The second shows you how to do everything manually. It is recommended that you use the Ansible method. The manual method is included here mainly for a more in depth understanding of what is happening.

[jupyter]: https://jupyter.org/
[gie-docs]: https://docs.galaxyproject.org/en/release_19.09/admin/special_topics/interactive_environments.html

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Configuring Galaxy Interactive Tools using Ansible

If the terms "Ansible," "role," and "playbook" mean nothing to you, please checkout [the Ansible introduction slides]({{ site.baseurl }}{% link topics/admin/tutorials/ansible/slides.html %}) and [the Ansible introduction tutorial]({{ site.baseurl }}{% link topics/admin/tutorials/ansible/tutorial.md %}).

**This section of the tutorial builds upon the work in the Ansible introduction tutorial, please ensure that you have completed that tutorial first.**

> ### {% icon comment %} Ansible Best Practices
> If you've set up your Galaxy server using the [Galaxy Installation with Ansible]({{ site.baseurl }}{% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, you will have created a `galaxyservers` group in your inventory file, `hosts`, and placed your variables in `group_vars/galaxyservers.yml`. Although for the purposes of this tutorial, the Galaxy server and cluster node are one and the same, in a real world deployment they are very likely to be different hosts. We will continue to use the `galaxyservers` group for simplicity, but in your own deployment you should consider creating an additional group for cluster nodes.
{: .comment}

## Installing Ansible Roles

We will use several Ansible roles for this tutorial. In order to avoid repetetively adding them to `requirements.yml` and installing them, we can simply install them all before getting started. Each role will be discussed in further detail later in the tutorial.

> ### {% icon hands_on %} Hands-on: Installing New Ansible Roles
>
> 1. In your working directory, add the docker role to your `requirements.yml`:
>
>    ```yaml
>    - src: geerlingguy.docker
>      version: 2.6.0
>    - src: usegalaxy_eu.gie_proxy
>      version: 0.0.1
>    ```
>
> 2. Install the requirements with `ansible-galaxy`:
>
>    ```console
>    ansible-galaxy role install -p roles -r requirements.yml
>    ```

## Installing Docker

Currently, Galaxy Interactive Tools must be run in Docker containers. It may be possible to run them in Singularity or other types of containers in the future. Thus, the first step is ensuring that the *nodes* where Galaxy will run have Docker installed. Both the Galaxy Project and Galaxy Project EU organizations have their own docker roles, but these are not published to Ansible Galaxy because they were mostly developed for internal purposes. For now, we will use the [docker role][geerlingguy-docker] by the prolific Ansible Galaxy publisher, [Jeff Geerling (geerlingguy)][geerlingguy]. Have a look at the geerlingguy.docker [README][geerlingguy-docker-readme] and [defaults/main.yml][geerlingguy-docker-defaults] to get an understanding of what variables are used to control the role.

[geerlingguy-docker]: https://galaxy.ansible.com/geerlingguy/docker
[geerlingguy]: https://galaxy.ansible.com/geerlingguy
[geerlingguy-docker-readme]: https://github.com/geerlingguy/ansible-role-docker/blob/master/README.md
[geerlingguy-docker-defaults]: https://github.com/geerlingguy/ansible-role-docker/blob/master/defaults/main.yml

> ### {% icon question %} Question
>
> What variables might be relevant to using this role?
>
> > ### {% icon solution %} Solution
> >
> > The `docker_users` variable (a *list*) controls which users are able to interact with the Docker daemon, which our Galaxy user will need to do. Additionally, Docker Compose is configured by default, which we do not need, so it can be disabled with `docker_install_compose: false`.
> >
> {: .solution }
>
{: .question}

> ### {% icon hands_on %} Hands-on: Installing Docker with Ansible
>
> 1. Edit the group variables file, `group_vars/galaxyservers.yml`:
>
>    The relevant variables to set for this role are:
>
>    | Variable                 | Type            | Description                                     |
>    | ----------               | -------         | -------------                                   |
>    | `docker_users`           | list of strings | List of users to be added to the `docker` group |
>    | `docker_install_compose` | boolean         | Whether to install and configure Docker Compose |
>
>    Add the following lines to your `group_vars/galaxyservers.yml` file:
>
>    {% raw %}
>    ```yaml
>    # Interactive Tools
>    docker_install_compose: false
>    docker_users:
>      - "{{ galaxy_user.name }}"
>    ```
>    {% endraw %}
>
>    > ### {% icon question %} Question
>    >
>    > {% raw %}
>    > Why is `"{{ galaxy_user.name }}"` specified instead of just the user `galaxy`?
>    > {% endraw %}
>    >
>    > > ### {% icon solution %} Solution
>    > > Duplicating values is never a good idea. If we needed to change the Galaxy user down the line or wanted to reuse this playbook on another host where the Galaxy username was different, we would have to change the value in multiple locations.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 2. Add the new role to the list of roles under the `roles` key in your playbook, `playbook.yml`:
>
>    ```yaml
>    ---
>    - hosts: galaxyservers
>      become: true
>      roles:
>        # ... existing roles ...
>        - geerlingguy.docker
>    ```
>
> 3. Run the playbook:
>
>    ```
>    ansible-playbook -i hosts playbook.yml
>    ```
{: .hands_on}

Congratulations, you've set up Docker. Verify the installation using the `docker info` command (but keep in mind: what users did we authorize to interact with Docker?).

## Installing the Interactive Tools Proxy

When an Interactive Tool's Docker container starts, it will be assigned a random port. In order to connect clients to the Interactive Tool, Galaxy needs to determine this port (and the node on which the tool is running) and configure a *proxy* from Galaxy to the GxIT's host and port. Consider the following example of running the Jupyter Notebook Interactive Tool, shown in Figure 1 below:

- nginx listens for requests from the client on **port 443** (https)
- Requests for Galaxy are delivered from nginx to Galaxy over a UNIX domain socket (uWSGI protocol)
- Requests for Interactive Tools are delivered from nginx to the Interactive Tools Proxy over (by default) **port 8000** (http)
  - GxIT http requests are forwarded by the proxy to Docker on the node on the container's (randomly assigned) **port 32768**
  - GxIT http requests are again forwarded by Docker to Jupyter on its in-container "published" **port 8888*

![Galaxy Interactive Tools Proxy Diagram](../../images/interactive-tools/gxit-proxy-diagram.png "Galaxy Interactive Tools Proxy Diagram")

As you can see, the client only ever speaks to nginx on the Galaxy server running on the standard https port (443), never directly to the interactive tool (which may be running on a node that does not even have a public IP address). The mapping of GxIT invocation and its corresponding host/port is kept in a SQLite database known as the *Interactive Tools Session Map*, and the path to this database is important, since both Galaxy and the proxy need access to it.

The GIE Proxy is written in [Node.js][nodejs] and requires some configuration. Thankfully there is an Ansible role, [usegalaxy_eu.gie_proxy][usegalaxy_eu-gie_proxy], that can install the proxy and its dependencies, and configure it for you. As usual, have a look through the [README][usegalaxy_eu-gie_proxy-readme] and [defaults][usegalaxy_eu-gie_proxy-defaults] to investigate which variables you might need to set before continuing.

[nodejs]: https://nodejs.org/
[usegalaxy_eu-gie_proxy]: https://galaxy.ansible.com/usegalaxy_eu/gie_proxy
[usegalaxy_eu-gie_proxy-readme]: https://github.com/usegalaxy-eu/ansible-gie-proxy/blob/master/README.md
[usegalaxy_eu-gie_proxy-defaults]: https://github.com/usegalaxy-eu/ansible-gie-proxy/blob/master/defaults/main.yml

> ### {% icon hands_on %} Hands-on: Installing the Proxy  with Ansible
>
> 1. Edit the group variables file, `group_vars/galaxyservers.yml`:
>
>    The relevant variables to set for this role are:
>
>    | Variable                   | Type          | Description                                                           |
>    | ----------                 | -------       | -------------                                                         |
>    | `gie_proxy_dir`            | path (string) | Path of directory into which the proxy application will be installed  |
>    | `gie_proxy_git_version`    | string        | Git reference to clone                                                |
>    | `gie_proxy_setup_nodejs`   | string        | Whether to install Node.js, options are `package` and `nodeenv`       |
>    | `gie_proxy_nodejs_version` | string        | Version of Node.js to install if using `nodeenv` method               |
>    | `gie_proxy_virtualenv`     | path (string) | Path of virtualenv into which nodeenv/Node.js/npm will be installed   |
>    | `gie_proxy_setup_service`  | string        | Whether to configure the proxy as a service, only option is `systemd` |
>    | `gie_proxy_sessions_path`  | path (string) | Path of Interactive Tools sessions map                                |
>
>    Add the following lines to your `group_vars/galaxyservers.yml` file:
>
>    {% raw %}
>    ```yaml
>    gie_proxy_dir: /srv/galaxy/gie-proxy/proxy
>    gie_proxy_git_version: ie2
>    gie_proxy_setup_nodejs: nodeenv
>    gie_proxy_nodejs_version: "10.13.0"
>    gie_proxy_virtualenv: /srv/galaxy/gie-proxy/venv
>    gie_proxy_setup_service: systemd
>    gie_proxy_sessions_path: "{{ galaxy_mutable_data_dir }}/interactivetools_map.sqlite"
>    ```
>    {% endraw %}
>
>    Note the value of `gie_proxy_git_version` is `ie2`: this is because the default branch only works with Interactive Environments, whereas the `ie2` branch has been updated for Interactive Tools. As Interactive Tools mature, this will likely be merged back to the default branch.
>
>    We have chosen to install Node.js using [nodeenv][] because the version in the training image's package manager is fairly old.
>
>    > ### {% icon question %} Question
>    >
>    > {% raw %}
>    > Why is `"{{ galaxy_user.name }}"` specified instead of just the user `galaxy`?
>    > {% endraw %}
>    >
>    > > ### {% icon solution %} Solution
>    > > Duplicating values is never a good idea. If we needed to change the Galaxy user down the line or wanted to reuse this playbook on another host where the Galaxy username was different, we would have to change the value in multiple locations.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 2. Add the new role to `playbook.yml`:
>
>    ```yaml
>    - hosts: galaxyservers
>      become: true
>      roles:
>        # ... existing roles ...
>        - geerlingguy.docker
>        - usegalaxy_eu.gie_proxy
>    ```
>
> 3. Run the playbook:
>
>    ```
>    ansible-playbook -i hosts playbook.yml
>    ```
{: .hands_on}

[nodeenv]: https://github.com/ekalinin/nodeenv
