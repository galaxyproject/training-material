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
{: .quote cite="https://docs.celeryq.dev/en/stable/getting-started/introduction.html#what-s-a-task-queue"}

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
In order to run a production ready Celery setup, we need to discuss and install some other software that works together with Celery.  
We already learned about RabbitMQ in the Pulsar tutorial. The RabbitMQ server you already installed there will be our broker for Celery.
As a backend we are going to use Redis.  
Redis is a very popular key-value-store database. It is very fast and easy to set up backend for Celery.
If you want to learn more about Redis, visit their website: [https://redis.io/](https://redis.io/)

For monitoring and debugging Celery, we use the [Flower](https://github.com/mher/flower) dashboard.  
Flower is lightweight and has a clear but powerful UI and can be installed in Galaxy's venv using our role.

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
>    @@ -38,3 +38,8 @@
>       version: 1.4.4
>     - src: galaxyproject.pulsar
>       version: 1.0.10
>    +# Celery, Redis, and Flower (dashboard)
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
> 2. Install the roles with:
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
>     2. We only need to add Redis' Python package in the `group_vars/galaxyservers.yml`:
>        {% raw %}
>        ```diff
>        --- a/group_vars/galaxyservers.yml
>        +++ b/group_vars/galaxyservers.yml
>        @@ -279,3 +279,7 @@ rabbitmq_users:
>         # TUS
>         galaxy_tusd_port: 1080
>         galaxy_tus_upload_store: /data/tus
>        +
>        +#Redis
>        +galaxy_additional_venv_packages:
>        +  - redis
>        {% endraw %}
>        ```
>        {: data-commit="Configure rabbitmq users"}
>        
>     3. Let's add the role to our playbook then:
>        {% raw %}
>        ```diff
>        --- a/galaxy.yml
>        +++ b/galaxy.yml
>        @@ -42,6 +42,7 @@
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
>     4. Because Flower needs it's own RabbitMQ user, we should add that to the respective part of our vars
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
>        @@ -1,13 +1,22 @@
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
>        +66653366633665383231303635643739396466653465626163633662623230643534346666303434
>        +3531356339646666623364383435306363373534373833370a333964306130653236373438373264
>        +34366433396135303064643932313135643064373034323865363939333565623238386161333465
>        +3932623162353034370a653663333363383938393936623663343639376430653366646662353830
>        +34383964326161383966306361633162653366353162633564623738353137333936313232363764
>        +38383962396637636136316366633936316365643038333030333932336365326335373737663834
>        +39363639303764343430336435313663313762353335366562656232646334356438343535303032
>        +35636131323163323036383163363964643237313137333131303737346662373233663562343031
>        +39396531306539623661376166313534303462623532323334393736316634616262313633626637
>        +37393066663239303564356232366466333334316565363631363230626665333039386133383933
>        +31656234363064646334623938383637393534313638643635396662643163346633346237353664
>        +34666161373764663166653737323736373261363565666633653831363833393264333666356633
>        +30373062333638366666316132623736386336383933326539663065653538386363323766313664
>        +63323737396264306439303366343834326164393763376663316366353766663131303966393039
>        +61343732313865613761623733356465336263633963376161663931653864383162393839643834
>        +34373239386461666462643162613536316166656135323332316163336635393035653164363138
>        +34336265383033663436306436663536383238366665653239313661666431666339323364646633
>        +64366361653034393962316562373835623631356339313062393933633634653735656636373031
>        +65376162633866623733623961303664663931636130373936326461623465363932373165356563
>        +36323838323865326464346337393164343765613333663830636234666663303535616532323936
>        +3230
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
>        @@ -266,6 +266,7 @@ rabbitmq_config:
>         
>         rabbitmq_vhosts:
>           - /pulsar/pulsar_au
>        +  - galaxy_internal
>         
>         rabbitmq_users:
>           - user: admin
>        @@ -275,6 +276,13 @@ rabbitmq_users:
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
>     5. Flower
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
>        @@ -291,3 +291,22 @@ galaxy_tus_upload_store: /data/tus
>         #Redis
>         galaxy_additional_venv_packages:
>           - redis
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
>        +flower_broker_url: "amqp://flower:{{ vault_rabbitmq_password_flower }}@localhost:5671/galaxy_internal?ssl=true"
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
>     6. It has a dashboard, so we need to expose that via nginx:
>
>        {% raw %}
>        ```diff
>        --- a/templates/nginx/galaxy.j2
>        +++ b/templates/nginx/galaxy.j2
>        @@ -94,4 +94,13 @@ server {
>         		proxy_set_header X-Forwarded-For   $proxy_add_x_forwarded_for;
>         		proxy_set_header X-Forwarded-Proto $scheme;
>         	}
>        +
>        +	location /flower {
>        +		proxy_pass http://localhost:5555;
>        +		proxy_set_header Host $host;
>        +		proxy_redirect off;
>        +		proxy_http_version 1.1;
>        +		proxy_set_header Upgrade $http_upgrade;
>        +		proxy_set_header Connection "upgrade";
>        +	}
>         }
>        {% endraw %}
>        ```
>        {: data-commit="Add nginx routes"}
>
>     7. Now we can add the Flower Role to our Playbook:
>    
>        {% raw %}
>        ```diff
>        --- a/galaxy.yml
>        +++ b/galaxy.yml
>        @@ -43,6 +43,7 @@
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
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -127,6 +127,11 @@ galaxy_config:
>           preload: true
>         celery:
>           concurrency: 2
>    +      enable_beat: true
>    +      enable: true
>    +      queues: celery,galaxy.internal,galaxy.external
>    +      pool: threads
>    +      memory_limit: 2
>           loglevel: DEBUG
>         tusd:
>           enable: true
>    {% endraw %}
>    ```
>    {: data-commit="Add celery" data-ref="add-req"}
>
>    Now add the second part, Galaxy's Celery configuration:
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -108,6 +108,11 @@ galaxy_config:
>         # Data Library Directories
>         library_import_dir: /libraries/admin
>         user_library_import_dir: /libraries/user
>    +    # Celery
>    +    amqp_internal_connection: "pyamqp://galaxy:{{ vault_rabbitmq_password_galaxy }}@localhost:5671/galaxy_internal?ssl=1"
>    +    celery_conf:
>    +      result_backend: "redis://localhost:6379/0"
>    +    enable_celery_tasks: true
>       gravity:
>         process_manager: systemd
>         galaxy_root: "{{ galaxy_root }}/server"
>    {% endraw %}
>    ```
>    {: data-commit="Add celery-redis" data-ref="add-req"}
>
> 6. We are done with the changes and you can enter the command to run your playbook:
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
> 2. In split view, open a second browser window and open your Galaxy page.
>    Click on {% icon galaxy-upload %} Upload Data, select a file from your computer and click `upload`.
> 3. The Workers should now receive a new tasks. Click on `Succeeded` and then on the UUID of the last upload task.  
>    You should see all its details here and the info that it was successful.
{: .hands_on}

{% snippet topics/admin/faqs/missed-something.md step=12 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="celery" %}
