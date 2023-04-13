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
  - Monitor a Celery task using the Flower dashboard
time_estimation: "1h"
key_points:
contributions:
  authorship:
  - mira-miracoli
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
{: .quote}

[A slideshow presentation on this subject is available](slides.html). 

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="celery" %}

The agenda we're going to follow today is: We're going to enable and configure celery, install a Redis server, the Flower dashboard and start Celery workers.

# Installing and Configuring

To proceed from here it is expected that:

<!--  TODO: port assumptions are not correct for GAT. -->

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

Redis is a very popular key-value-store database. It is very fast and a good backend for Celery.
If you want to learn more about Redis, visit their website: (https://redis.io/)[https://redis.io/]

Good news: Celery is already installed in your Galaxy's virtual environment if you followed the last tutorials and completed the Galaxy installation successfully.  
Also RabbitMQ should be up and running after you completed the Pulsar tutorial.
Still we need to add a few things to out Playbooks:
 - Redis is a very popular key-value-store database. It is very fast and a good backend for Celery. If you want to learn more about Redis, visit their website: [https://redis.io/](https://redis.io/)
 - Flower is a powerful dashboard for Celery and can be installed in Galaxy's venv using our role.

# Installing and Configuring

First we need to add our new Ansible Roles to the `requirements.yml`:

> <hands-on-title>Set up Redis, Flower, Systemd and Celery with Ansible</hands-on-title>
>
> 1. In your working directory, add the roles to your `requirements.yml`
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -30,3 +30,7 @@
>       version: 1.4.2
>     - src: galaxyproject.pulsar
>       version: 1.0.10
>    +- name: geerlingguy.redis
>    +  version: 1.8.0
>    +- name: usegalaxy_eu.flower
>    +  version: 1.0.2
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement" data-ref="add-req"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
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
> 3. Let's go now through all the Roles step-by-step:
>
>     1. Since we can stick to the basic default settings of Redis, we will look only at a few variables:
>
>        | Variable               | Type    | Description                                                                 |
>        | ----------             | ------- | -------------                                                               |
>        | `redis_port`           | integer | The port Redis should listen on. 6379 by default.                           |
>        | `redis_bind_interface` | string  | The interface Redis should listen to. 127.0.0.1 is default.                 |
>        | `redis_conf_path`      | string  | The path where your redis configuration will be stored. Default: /etc/redis |
>
>        Luckily we can leave them all on default and don't need to change anything for Redis in the vars.  
>        Let's add the role to our playbook then:
>        {% raw %}
>        ```diff
>        --- a/galaxy.yml
>        +++ b/galaxy.yml
>        @@ -45,6 +45,7 @@
>             - role: galaxyproject.miniconda
>               become: true
>               become_user: "{{ galaxy_user_name }}"
>        +    - geerlingguy.redis
>             - galaxyproject.nginx
>             - geerlingguy.docker
>             - usegalaxy_eu.rabbitmqserver
>        {% endraw %}
>        ```
>        {: data-commit="Add requirement" data-ref="add-req"}
>
>     2. Since Flower needs it's own RabbitMQ user, we should add that to the respective part of our vars
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
>        @@ -1,13 +1,18 @@
>         $ANSIBLE_VAULT;1.1;AES256
>        -62346261323266656232393034396134316636376533376139666437363535393562663838613938
>        -6336666266633563346337623265353935646361326337610a393834333233313461346439376438
>        -63383338346530656561636631666134373238366364363164313166346461383736613162653237
>        -3461363334323431370a656132303965653262386130353332623937376261396530393761353834
>        -38336565666437666436643163363831633331333766653266356163613138393734656465323634
>        -39366362383433366437353534663134313330316337393335383962613961386665633261616237
>        -35366635373063313631323939396164336330356361393464326636353037336461323531336434
>        -35613933303333623031353936393265636130363335376533393335663266313863376135383338
>        -36613464373231623938373434306266373234633036343636633963353361356631363533353066
>        -39323064336237646432323530313065303331326636353334343862373330313133326363363063
>        -38383564636161396435666164643334656435393533643163393434623434656238633631633939
>        -33353232666432376661
>        +39393239303639616131633130366134376431396464366131373430383435656261303633346638
>        +3965633961383235386230346561366637653561363961300a316537623964343132366163313038
>        +61663061656564393331663661643039386433326134636161636661613836396536636236336161
>        +3639636434333931380a646338383664646332343364393761313462346535623930663838313833
>        +31626164643766636564356333636364343164663562663733333261393039616535313439643264
>        +31333764393130326139306362333136656336663834316462333566313532373162353331393864
>        +37633064636364363335633331623063623639353536353263333661363166633636333833303862
>        +38343132313632336631356461323733616339646364626564323932373539343964373961313035
>        +31366532323539356532663537303763343035643066653062346238616233646439643061653537
>        +37636432343161323832336236326335626564393730663563353162303230306635393463643636
>        +33323165303061383732623066623837323037396563396263363230333634643739333262373932
>        +37633864383763616364316633366262666362643132623463346465373865373732346363353238
>        +64336230346636633637613265346630343231613161373861313932633638326431653534363933
>        +63383930313835653461633164386136343833646230656535356137626337373535373530376166
>        +38643930623339653739373132623731653536613135353236333235313137323135616636333239
>        +36663831653735396330366561346361623862626163386433626431353730323562626239373333
>        +37303031373237633934333236653361366563373332366166383662323733376335
>        {% endraw %}
>        ```
>        {: data-commit="Add rabbitmq passwords to the vault"}
>
>        -->
>
>        This is going in the vault as they are secrets we need to set. Flower needs it's own RabbitMQ user with admin access and we want a different vhost for galaxy and celery.
>
>        Replace both with long random (or not) string.  
>        Now add new users to the RabbitMQ configuration:
>        {% raw %}
>        ```diff
>        --- a/group_vars/galaxyservers.yml
>        +++ b/group_vars/galaxyservers.yml
>        @@ -230,6 +230,7 @@ rabbitmq_config:
>         
>         rabbitmq_vhosts:
>           - /pulsar/pulsar_au
>        +  - galaxy_internal
>         
>         rabbitmq_users:
>           - user: admin
>        @@ -239,6 +240,13 @@ rabbitmq_users:
>           - user: pulsar_au
>             password: "{{ vault_rabbitmq_password_vhost }}"
>             vhost: /pulsar/pulsar_au
>        +  - user: galaxy
>        +    password: "{{ vault_rabbitmq_password_galaxy }}"
>        +    vhost: galaxy_internal
>        +  - user: flower
>        +    password: "{{ vault_rabbitmq_password_flower }}"
>        +    tags: administrator
>        +    vhost: galaxy_internal
>         
>         # TUS
>         galaxy_tusd_port: 1080
>        {% endraw %}
>        ```
>        {: data-commit="Configure rabbitmq users"}
>>     2. Flower
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
>        | `flower_db_file`  | string | When Flower is in persistent mode, use this path for the database. |
>
>        Let's add variables to our `group_vars/galaxyservers.yml`:
>
>        {% raw %}
>        ```diff
>        --- a/group_vars/galaxyservers.yml
>        +++ b/group_vars/galaxyservers.yml
>        @@ -260,3 +260,22 @@ tusd_instances:
>               - "-upload-dir={{ galaxy_config.galaxy.tus_upload_store }}"
>               - "-hooks-http=https://{{ inventory_hostname }}/api/upload/hooks"
>               - "-hooks-http-forward-headers=X-Api-Key,Cookie"
>        +
>        +# Flower
>        +flower_python_version: python3
>        +flower_app_dir: "{{ galaxy_root }}"
>        +flower_python_path: "{{ galaxy_root }}/server/lib"
>        +flower_venv_dir: "{{ galaxy_venv_dir }}"
>        +flower_app_name: galaxy.celery
>        +flower_db_file: "{{ galaxy_root }}/var/flower.db"
>        +flower_persistent: true
>        +flower_broker_api: "https://flower:{{ vault_rabbitmq_password_flower }}@localhost:5671/api/"
>        +flower_broker_url: "amqp://flower:{{ vault_rabbitmq_password_flower }}@localhost:5671/galaxy?ssl=true"
>        +flower_proxy_prefix: /flower
>        +
>        +flower_ui_users:
>        +  - name: admin
>        +    password: "{{ vault_flower_user_password}}"
>        +
>        +flower_environment_variables:
>        +  GALAXY_CONFIG_FILE: "{{ galaxy_config_file }}"
>        {% endraw %}
>        ```
>        {: data-commit="Configure flower"}
>
> 1. Now we can add the Flower Role to our Playbook:
>
>        {% raw %}
>        ```diff
>        --- a/galaxy.yml
>        +++ b/galaxy.yml
>        @@ -46,6 +46,7 @@
>               become: true
>               become_user: "{{ galaxy_user_name }}"
>             - geerlingguy.redis
>        +    - usegalaxy_eu.flower
>             - galaxyproject.nginx
>             - geerlingguy.docker
>             - usegalaxy_eu.rabbitmqserver
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
>
>     Now add the second part, Galaxy's Celery configuration:
>     {% raw %}
>     ```diff
>     --- a/group_vars/galaxyservers.yml
>     +++ b/group_vars/galaxyservers.yml
>     @@ -191,6 +191,9 @@ galaxy_config:
>            url_prefix: /reports
>            bind: "unix:{{ galaxy_mutable_config_dir }}/reports.sock"
>            config_file: "{{ galaxy_config_dir }}/reports.yml"
>     +    amqp_internal_connection: "pyamqp://galaxy:{{ vault_rabbitmq_password_galaxy }}@localhost:5671/galaxy?ssl=1"
>     +    celery_conf:
>     +      result_backend: "redis://localhost:6379/0"
>     +      enable_celery_tasks: true
>      
>      galaxy_config_templates:
>        - src: templates/galaxy/config/container_resolvers_conf.yml.j2
>     {% endraw %}
>     ```
>     {: data-commit="Add celery-redis" data-ref="add-req"}
>
> 1. And finally expose it via nginx:
>
>        {% raw %}
>        ```diff
>        --- a/templates/nginx/galaxy.j2
>        +++ b/templates/nginx/galaxy.j2
>        @@ -84,4 +84,13 @@ server {
>             location /training-material/ {
>                 proxy_pass https://training.galaxyproject.org/training-material/;
>             }
>        +
>        +    location /flower {
>        +        proxy_pass http://localhost:5555;
>        +        proxy_set_header Host $host;
>        +        proxy_redirect off;
>        +        proxy_http_version 1.1;
>        +        proxy_set_header Upgrade $http_upgrade;
>        +        proxy_set_header Connection "upgrade";
>        +    }
>         }
>        {% endraw %}
>        ```
>        {: data-commit="Add nginx routes"}
>
> 5. We are done with the changes and you can enter the command to run your playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>    This should also restart Galaxy and spawn the amount of Celery workers, that we defined in the Gravity configuration.
>
{: .hands_on}

# Test Celery

Now that everything is running, we want to test celery and watch it processing tasks.
We can simply do that by starting an upload to our Galaxy.

> <hands-on-title>Test Celery and monitor tasks with Flower</hands-on-title>
> 1. First, open a new tab and enter your machines hostname followed by `/flower/dashboard` then log in with `username: admin` and you password.
>    You should see an overview with active workers.  
>    Keep that tab open
> 2. In split view, open a second browser window and open you Galaxy page.
>    Click on {% icon galaxy-upload %} Upload Data, select a file from your computer and click `upload`.
> 3. The Workers should now receive a new tasks. Click on `Succeeded` and then on the UUID of the last upload task.  
>    You should see all its details here and the info that is was successful.
{: .hands_on}

{% snippet topics/admin/faqs/missed-something.md step=11 %}
