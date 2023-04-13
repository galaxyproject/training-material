---
layout: tutorial_hands_on

title: "Setting up Celery Workers for Galaxy"
zenodo_link: ""
questions:
objectives:
  - Have an understanding of what Celery is and how it works
  - Install Redis
  - Configure and start Celery workers
  - Install Flower to the Galaxy venv and configure it
  - Use an Ansible playbook for all of the above.
time_estimation: "1h"
key_points:
contributors:
  - mira-miracoli
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - pulsar
voice:
  id: Olivia
  lang: en-AU
  neural: true
subtopic: data
tags:
  - ansible
  - git-gat
---

# Overview


Celery is a distributed task queue written in Python that can spawn multiple workers and enables asynchronous task processing on multiple nodes. It supports scheduling, but focuses more on real-time operations.

From the Celery website:

> "Task queues are used as a mechanism to distribute work across threads or machines.
>
>A task queue’s input is a unit of work called a task. Dedicated worker processes constantly monitor task queues for new work to perform.
>
>Celery communicates via messages, usually using a broker to mediate between clients and workers. To initiate a task the client adds a message to the queue, the broker then delivers that message to a worker.
>
>A Celery system can consist of multiple workers and brokers, giving way to high availability and horizontal scaling.
>
>Celery is written in Python, but the protocol can be implemented in any language. In addition to Python there’s node-celery and node-celery-ts for Node.js, and a PHP client.
>
>Language interoperability can also be achieved exposing an HTTP endpoint and having a task that requests it (webhooks)."
>
> -- [https://docs.celeryq.dev/en/stable/getting-started/introduction.html#what-s-a-task-queue](https://docs.celeryq.dev/en/stable/getting-started/introduction.html#what-s-a-task-queue)
{: .quote id="celery-quote"}

[A slideshow presentation on this subject is available](slides.html). 

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="celery" %}

> The agenda we're going to follow today is: We're going to enable and configure celery, install a Redis server, the Flower dashboard and start Celery workers.
{: .spoken data-visual="gtn" data-target="#agenda"}

# Installing and Configuring

To proceed from here it is expected that:

> <comment-title>Requirements for Running This Tutorial</comment-title>
>
> 1. You have set up a working Galaxy instance as described in the [ansible-galaxy](../ansible-galaxy/tutorial.html) tutorial.
>
> 2. You have a working RabbitMQ server installed and added the connection string to the galaxy configuration. (RabbitMQ is installed when doing the [Pulsar](../pulsar/tutorial.html) tutorial.)
>
> 3. Your VM has a public DNS name: this tutorial sets up SSL certificates from the start and as an integral part of the tutorial.
>
> 4. You have the following ports exposed:
>
>    - 22 for SSH, this can be a different port or via VPN or similar.
>    - 80 for HTTP, this needs to be available to the world if you want to follow the LetsEncrypt portion of the tutorial.
>    - 443 for HTTPs, this needs to be available to the world if you want to follow the LetsEncrypt portion of the tutorial.
>    - 5671 for AMQP for Pulsar, needed if you plan to setup Pulsar for remote job running.
>
{: .comment}

> Redis is a very popular key-value-store database. It is very fast and a good backend for Celery.
> If you want to learn more about Redis, visit their website: (https://redis.io/)[https://redis.io/]
{: .spoken data-visual="gtn" data-target="#preparations:redis" }

Good news: Celery is already installed in your Galaxy's virtual environment if you followed the last tutorials and completed the Galaxy installation successfully.  
Also RabbitMQ should be up and running after you completed the Pulsar tutorial.
Still we need to add a few things to out Playbooks.
 - Redis is a very popular key-value-store database. It is very fast and a good backend for Celery.
If you want to learn more about Redis, visit their website: [https://redis.io/](https://redis.io/)
Installing Redis with Galaxy-EU's Ansible role is fast and simple, too!
 - Flower is a powerful dashboard for Celery and can be installed in Galaxy's venv using our role.


# Installing and Configuring


First we need to add our new Ansible Roles to the `requirements.yml`:

If the terms "Ansible", "role" and "playbook" mean nothing to you, please checkout [the Ansible introduction slides]({% link topics/admin/tutorials/ansible/slides.html %}) and [the Ansible introduction tutorial]({% link topics/admin/tutorials/ansible/tutorial.md %})

{% snippet topics/admin/faqs/ansible_local.md %}

> Okay, so let's get started. If we go back to our
> tutorial here, it says that we need to install the roles mentioned above into our
> requirements.yml and then add it to our Ansible.
{: .spoken data-visual="gtn" data-target="#hands-on-set-up-redis-flower-systemd-and-celery-with-ansible"}

> <hands-on-title>Set up Redis, Flower, Systemd and Celery with Ansible</hands-on-title>
>
> 1. In your working directory, add the roles to your `requirements.yml`
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -36,3 +36,9 @@
>       version: 2.1.3
>     - src: galaxyproject.proftpd
>       version: 0.3.1
>    +- name: geerlingguy.redis
>    +  version: 1.8.0
>    +- name: usegalaxy_eu.flower
>    +  version: 1.0.1
>    +- name: usegalaxy_eu.galaxy_systemd
>    +  version: 2.1.0
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement" data-ref="add-req"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
>    > Okay, so the first thing I'm going to do is I'm going to add the Redis, Flower and Systemd
>    > roles to the requirements.yml.
>    > Edit requirements.yml and we need to add this to the bottom of that file. Copy. Paste. And save it.
>    {: .spoken data-visual="terminal" data-ref="add-req"}
>
> 2. Install the role with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true" data-ref="req-install"}
>    {: .code-in}
>
>    > And now install the role into our local Ansible scripts using the
>    > ansible-galaxy command. And as you can see, it's downloading the
>    > roles.
>    {: .spoken data-visual="terminal" data-ref="req-install"}
>
>    > And if we look into roles now you can see that we have them.
>    {: .spoken data-visual="terminal" data-cmd="ls roles/"}
>
>    > Right, clear the screen.
>    {: .spoken data-visual="terminal" data-cmd="clear"}
>
> 3. Let's go now through all the Roles step-by-step:
>
>     1. Redis  
>        Since we can stick to the basic default settings, we will look only at a few variables:
>
>        | Variable             | Type          | Description                                                                                                                                                                    |
>        | ----------           | -------       | -------------                                                                                                                                                                  |
>        | `redis_port`         | integer        | The port Redis should listen on. 6379 by default.                               |
>        | `redis_bind_interface`         | string | The interface Redis should listen to. 127.0.0.1 is default.  |
>        | `redis_conf_path`  | string | The path where your redis configuration will be stored. Default: /etc/redis |
>        {: id="redis-variables-table"}
>
>        Luckily we can leave them all on default and don't need to change anything for Redis in the vars.  
>        Let's add the role to our playbook then:
>        {% raw %}
>        ```diff
>        --- a/galaxy.yml
>        +++ b/galaxy.yml
>        @@ -51,6 +51,7 @@
>             - galaxyproject.tusd
>             - galaxyproject.cvmfs
>             - dj-wasabi.telegraf
>        +    - geerlingguy.redis
>           post_tasks:
>             - name: Setup gxadmin cleanup task
>               ansible.builtin.cron:
>        {% endraw %}
>        ```
>        {: data-commit="Add requirement" data-ref="add-req"}
>
>     2. RabbitMQ Users  
>        Since Flower needs it's own RabbitMQ user, we should add that to the respective part of our vars
>        Edit your `group_vars/secret.yml` and define some random passwords:
>
>        ><code-in-title>Bash</code-in-title>
>        > ```
>        > ansible-vault edit group_vars/secret.yml
>        > ```
>        {: .code-in}
>    
>        ```yaml
>        vault_rabbitmq_password_flower: "a-really-long-password-here"
>        vault_rabbitmq_password_galaxy: "a-different-really-long-password"
>        vault_flower_user_password: "another-different-really-long-password"
>        ```
>
>        <!-- Ignore this, just for the gat-automation. Vaults are ugly to work with :(
>
>        {% raw %}
>        ```diff
>        --- a/group_vars/secret.yml
>        +++ b/group_vars/secret.yml
>        @@ -1,7 +1,13 @@
>         $ANSIBLE_VAULT;1.1;AES256
>        -32653961383866636531396135663630386630346237333333653633313436663439643535323964
>        -6363626330336430363332643638646262316338313937320a666566306539373462386266383166
>        -30326165393863633463353234613561393939326164376432633732316264636464313061383161
>        -3532373937656138320a616361343664353264613332616236623231326137316635323465623562
>        -66656539346130353639623736633034653932373438663330646436656336666637313933666264
>        -3636313438626533633831323239373461373538646635613637
>        +62346261323266656232393034396134316636376533376139666437363535393562663838613938
>        +6336666266633563346337623265353935646361326337610a393834333233313461346439376438
>        +63383338346530656561636631666134373238366364363164313166346461383736613162653237
>        +3461363334323431370a656132303965653262386130353332623937376261396530393761353834
>        +38336565666437666436643163363831633331333766653266356163613138393734656465323634
>        +39366362383433366437353534663134313330316337393335383962613961386665633261616237
>        +35366635373063313631323939396164336330356361393464326636353037336461323531336434
>        +35613933303333623031353936393265636130363335376533393335663266313863376135383338
>        +36613464373231623938373434306266373234633036343636633963353361356631363533353066
>        +39323064336237646432323530313065303331326636353334343862373330313133326363363063
>        +38383564636161396435666164643334656435393533643163393434623434656238633631633939
>        +33353232666432376661
>        {% endraw %}
>        ```
>        {: data-commit="Add rabbitmq passwords to the vault"}
>
>        -->
>
>        This is going in the vault as they are secrets we need to set. Both of our services, Galaxy and Pulsar, need these variables, so we'll need to make sure they're in both playbooks. Both Galaxy in the job configuration, and Pulsar in its configuration.
>
>        Replace both with long random (or not) string.  
>        Now add new users to the RabbitMQ configuration:
>        {% raw %}
>        ```diff
>        --- a/group_vars/galaxyservers.yml
>        +++ b/group_vars/galaxyservers.yml
>        @@ -305,9 +305,16 @@ rabbitmq_users:
>        rabbitmq_vhosts:
>          - /pulsar/galaxy_au
>        + - galaxy
>        @@ .... @@
>             password: "{{ vault_rabbitmq_admin_password }}"
>             tags: administrator
>             vhost: /
>        +  - user: galaxy
>        +    password: "{{ vault_rabbitmq_password_galaxy }}"
>        +    vhost: galaxy
>           - user: galaxy_au
>             password: "{{ vault_rabbitmq_password_vhost }}"
>             vhost: /pulsar/galaxy_au
>        +  - user: flower
>        +    password: "{{ vault_rabbitmq_password_flower }}"
>        +    tags: administrator
>        +    vhost: galaxy
>         
>         # Proftpd:
>         proftpd_galaxy_auth: yes
>        {% endraw %}
>        ```
>        {: data-commit="Add requirement" data-ref="add-req"}
>     2. Flower
>        Flower has a few variables, too, for example, we need to point it to our virtual environment:
>
>        | Variable             | Type          | Description                                                                                                                                                                    |
>        | ----------           | -------       | -------------                                                                                                                                                                  |
>        | `flower_python_version`         | string        | Python version to use when installing flower to a venv. Default: python39                               |
>        | `flower_port`         | integer        | The port Flower should listen on. 5555 by default.                               |
>        | `flower_bind_interface`         | string | The interface Flower should listen to. 0.0.0.0 is default.  |
>        | `flower_conf_dir`  | string | The path where your Flower configuration will be stored. Default: /etc/flower |
>        | `flower_venv_dir`  | string | The path to the venv where Flower should be installed. Default: `/home/{{ flower_user }}/.local` |
>        | `flower_user`  | string | User that owns the flower process. Default: galaxy |
>        | `flower_group`  | string | Group that owns the flower process. Default: galaxy |
>        | `flower_ui_users`  | list of dicts | Name and password of the UI users for basic auth. |
>        | `flower_app_dir`  | string | Root directory of your Python app to run with Celery. In our case `galaxy_root` |
>        | `flower_app_name`  | string | Python module to import. In our case 'galaxy.celery' |
>        | `flower_python_path`  | string | Should point to galaxy's `server/lib` directory (default) |
>        | `flower_broker_api`  | string | URL to broker's API with login credentials. |
>        | `flower_broker_url`  | string | Flower's RabbitMQ connection string. |
>        {: id="flower-variables-table"}
>        Let's add variables to our `group_vars/galaxyservers.yml`:
>
>        {% raw %}
>        ```diff
>        --- a/group_vars/galaxyservers.yml
>        +++ b/group_vars/galaxyservers.yml
>        @@ -370,3 +370,20 @@ tusd_instances:
>               - "-upload-dir={{ galaxy_config.galaxy.tus_upload_store }}"
>               - "-hooks-http=https://{{ inventory_hostname }}/api/upload/hooks"
>               - "-hooks-http-forward-headers=X-Api-Key,Cookie"
>        +
>        +# Flower
>        +flower_python_version: python3
>        +flower_app_dir: "{{ galaxy_root }}"
>        +flower_log_file: /var/log/flower
>        +flower_python_path: server/lib
>        +flower_venv_dir: "{{ galaxy_venv_dir }}"
>        +flower_app_name: galaxy.celery
>        +flower_db_file: "{{ galaxy_root }}/var/flower.db"
>        +flower_persistent: true
>        +
>        +flower_broker_api: "https://flower:{{ vault_rabbitmq_password_flower }}@localhost:5671/api/"
>        +flower_broker_url: "amqp://flower:{{ vault_rabbitmq_password_flower }}@localhost:5671/galaxy?ssl=true"
>        +
>        +flower_ui_users:
>        +  - name: admin
>        +    password: "{{ vault_flower_user_password}}"
>        {% endraw %}
>        ```
>        {: data-commit="Add requirement" data-ref="add-req"}
>        Now we can add the Flower Role to our Playbook:
>        {% raw %}
>        ```diff
>        --- a/galaxy.yml
>        +++ b/galaxy.yml
>        @@ -52,6 +52,7 @@
>             - galaxyproject.cvmfs
>             - dj-wasabi.telegraf
>             - geerlingguy.redis
>        +    - usegalaxy_eu.flower
>           post_tasks:
>             - name: Setup gxadmin cleanup task
>               ansible.builtin.cron:
>        {% endraw %}
>        ```
>        {: data-commit="Add flower role" data-ref="add-req"}
>
> 4. Now it is time to change the `group_vars/galaxyservers.yml` and enable celery in galaxy.gravity config.
>    Add the following lines to your file:
>     {% raw %}
>     ```diff
>     --- a/group_vars/galaxyservers.yml
>     +++ b/group_vars/galaxyservers.yml
>     @@ -174,6 +174,11 @@ galaxy_config:
>            preload: true
>          celery:
>            concurrency: 2
>     +      enable_celery_beat: true
>     +      enable: true
>     +      queues: celery,galaxy.internal,galaxy.external
>     +      pool: threads
>     +      memory_limit: 2G
>            loglevel: DEBUG
>          handlers:
>            handler
>     {% endraw %}
>     ```
>     {: data-commit="Add celery" data-ref="add-req"}
>     Now add the second part, Galaxy's Celery configuration:
>     {% raw %}
>     ```diff
>     --- a/group_vars/galaxyservers.yml
>     +++ b/group_vars/galaxyservers.yml
>     @@ -191,6 +191,9 @@ galaxy_config:
>            url_prefix: /reports
>            bind: "unix:{{ galaxy_mutable_config_dir }}/reports.sock"
>            config_file: "{{ galaxy_config_dir }}/reports.yml"
>     +    celery_conf:
>     +      result_backend: "redis://localhost:6379/0"
>     +      enable_celery_tasks: true
>      
>      galaxy_config_templates:
>        - src: templates/galaxy/config/container_resolvers_conf.yml.j2
>     {% endraw %}
>     ```
>     {: data-commit="Add celery-redis" data-ref="add-req"}
{: .hands_on}

# Test Celery
Now that everything is running, we want to test celery and watch it processing tasks.
We can simply do that by starting an upload to our Galaxy.

> <hands-on-title>Test Celery and monitor tasks with Flower</hands-on-title>
> 1. First, open a new tab and type `localhost:5555` then log in with `username: admin` and you password.
>    You should see an overview with active workers.  
>    Keep that tab open
> 2. In split view, open a second browser window and open you Galaxy page.
>    Click on {% icon galaxy-upload %}Upload Data, select a file from your computer and click `upload`.
> 3. The Workers should now receive a new tasks. Click on `Succeeded` and then on the UUID of the last upload task.  
>    You should see all its details here and the info that is was successful.
>
{: .hands_on}
