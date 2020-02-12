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

## Installing Docker

Currently, Galaxy Interactive Tools must be run in Docker containers. It may be possible to run them in Singularity or other types of containers in the future. Thus, the first step is ensuring that the *nodes* where Galaxy will run have Docker installed. Both the Galaxy Project and Galaxy Project EU organizations have their own docker roles, but these are not published to Ansible Galaxy because they were mostly developed for internal purposes. For now, we will use the [docker role][geerlingguy-docker] by the prolific Ansible Galaxy publisher, [Jeff Geerling (geerlingguy)][geerlingguy]. Have a look at the [geerlingguy.docker][geerlingguy-docker] README and `defaults/main.yml` to get an understanding of what variables are used to control the role.

[geerlingguy-docker]: https://galaxy.ansible.com/geerlingguy/docker
[geerlingguy]: https://galaxy.ansible.com/geerlingguy

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
> 1. In your working directory, add the docker role to your `requirements.yml`
>
>    ```yaml
>    ---
>    - src: geerlingguy.docker
>      version: 2.6.0
>    ```
>
> 2. Install the requirements with `ansible-galaxy`:
>
>    ```console
>    ansible-galaxy role install -p roles -r requirements.yml
>    ```
>
> 3. Create and edit the group variables file, `group_vars/galaxyservers.yml`.
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
>    ```yaml
>    ---
>    # Interactive Tools
>    docker_install_compose: false
>    docker_users:
>      - "{{ galaxy_user.name }}"
>    ```
>
>    > ### {% icon question %} Question
>    >
>    > Why is `"{{ galaxy_user.name }}"` specified instead of just the user `galaxy`?
>    >
>    > > ### {% icon solution %} Solution
>    > > Duplicating values is never a good idea. If we needed to change the Galaxy user down the line or wanted to reuse this playbook on another host where the Galaxy username was different, we would have to change the value in multiple locations.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 4. Create and edit a new playbook which uses the Docker role, name it `gxit_playbook.yml`
>    ```yaml
>    - hosts: galaxyservers
>      become: true
>      roles:
>        - geerlingguy.docker
>    ```
>
> 5. Run the playbook
>
>    ```
>    ansible-playbook -i hosts gxit_playbook.yml
>    ```
{: .hands_on}

Congratulations, you've set up Docker.
