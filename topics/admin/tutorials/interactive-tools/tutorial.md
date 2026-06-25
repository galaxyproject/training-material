---
layout: tutorial_hands_on

title: "Galaxy Interactive Tools"
zenodo_link: ""
questions:
  - What are Galaxy Interactive Tools?
  - How can I enable Interactive Tools on my Galaxy instance?
objectives:
  - Understand what Galaxy Interactive Tools are and how they work
  - Be aware of the security implications of Interactive Tools
  - Have a basic understanding of the Interactive Tools (GxIT/GIE) Proxy, its purpose, and configuration
  - Be familiar with wildcard SSL certificates and how to get them from Let's Encrypt
  - Configure your Galaxy to serve Interactive Tools using an Ansible Playbook
  - Start, run, and use an Interactive Tool
time_estimation: "2h"
key_points:
  - Galaxy Interactive Tools run as jobs in largely the same manner as any other Galaxy job
  - nginx routes GxIT requests to the GxIT(/GIE) Proxy, which routes them to the node/port on which the GxIT is running
  - GxITs require wildcard SSL certificates
  - GxITs expose your Galaxy server's user datasets unless configured to use Pulsar
contributions:
  authorship:
    - natefoo
    - slugger70
    - hexylena
    - abretaud
    - kysrpex
  funding:
    - eurosciencegateway
    - unimelb
    - melbournebioinformatics
    - AustralianBioCommons
tags:
  - ansible
  - interactive-tools
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - connect-to-compute-cluster
      - job-destinations
subtopic: features

recordings:
- captioners:
  - abretaud
  - hexylena
  date: '2021-02-15'
  galaxy_version: '21.01'
  length: 44M
  youtube_id: lACsIhnbTbE
  speakers:
  - abretaud

---


Galaxy Interactive Tools (GxITs) are a method to run containerized tools that are interactive in nature. Interactive Tools typically run a persistent service accessed on a specific port, until terminated by the user. One common example of such a tool is [Jupyter Notebook][jupyter]. Galaxy Interactive Tools are submitted through Galaxy's job management system and thus are scheduled the same as any other Galaxy tool - on a Slurm cluster, for instance. Galaxy Interactive Tools were introduced in Galaxy Release 19.09.

> <comment-title>Evolving Topic</comment-title>
> Galaxy Interactive Tools are a **relatively new and rapidly evolving feature** and there are some rough edges. Work to improve the experience of deploying and using them is ongoing. Please watch the [Galaxy Release Notes][galaxy-release-notes] for updates, changes, new documentation, and bug fixes.
>
> You may find extra information about Interactive Tools on the [Galaxy Documentation][galaxy-docs-interactivetools].
{: .comment}

[galaxy-release-notes]: https://docs.galaxyproject.org/en/master/releases/index.html
[galaxy-docs-interactivetools]: https://docs.galaxyproject.org/en/master/admin/special_topics/interactivetools.html

> <warning-title>Before You Continue - Wildcard DNS Certificates</warning-title>
> If you are *not* completing this tutorial as part of a [Galaxy Admin Training][gat] course, **you will need a wildcard DNS record for your Galaxy server and a method for obtaining a wildcard SSL certificate for your Galaxy server**.
>
> Galaxy Interactive Tools require a [wildcard SSL certificate][wildcard-cert]. Because the **Galaxy Installation with Ansible** tutorial fetches [Let's Encrypt][lets-encrypt] certificates, this tutorial fetches Let's Encrypt wildcard certificates. However, this process is only valid for Galaxy Admin Training courses, because Let's Encrypt wildcard certificates [can only be fetched using the DNS-01 challenge method][lets-encrypt-faq], which requires control of a [dynamic DNS][ddns] server (which we have preconfigured for use at training courses). Configuring your DNS service for dynamic updates is outside the scope of this tutorial, but it will show you how to request certificates using DNS-01, which can be adapted for your site.
>
> If you are using Let's Encrypt, [a list of available DNS plugins for Certbot][certbot-dns-plugins] can be found in the Certbot documentation. If you are not using Let's Encrypt, please consult your certificate vendor's documentation for information on how to obtain a wildcard certificate. You will need a certificate with (at least) the [subject alternative name][san]s `galaxy.example.org` and `*.ep.interactivetool.galaxy.example.org` (where `galaxy.example.org` is the hostname of your Galaxy server).
>
> You will also need a wildcard DNS `CNAME` record for `*.ep.interactivetool.galaxy.example.org`. You can verify that your Galaxy server has such a record using the `host` or `dig` command line tools like so:
>
>    ```console
>    $ host -t cname foo.ep.interactivetool.live.usegalaxy.eu
>    foo.ep.interactivetool.live.usegalaxy.eu is an alias for usegalaxy.eu.
>    $ host -t cname bar.ep.interactivetool.live.usegalaxy.eu
>    bar.ep.interactivetool.live.usegalaxy.eu is an alias for usegalaxy.eu.
>    ```
>
> Please consult your DNS server software or cloud provider's documentation for information on how to set up a wildcard record.
{: .warning}

[wildcard-cert]: https://en.wikipedia.org/wiki/Wildcard_certificate
[lets-encrypt]: https://letsencrypt.org/
[gat]: https://github.com/galaxyproject/dagobah-training
[lets-encrypt-faq]: https://letsencrypt.org/docs/faq/
[ddns]: https://en.wikipedia.org/wiki/Dynamic_DNS
[certbot-dns-plugins]: https://certbot.eff.org/docs/using.html#dns-plugins
[san]: https://en.wikipedia.org/wiki/Subject_Alternative_Name
[jupyter]: https://jupyter.org/
[gie-docs]: https://docs.galaxyproject.org/en/release_19.09/admin/special_topics/interactive_environments.html

> <warning-title>Uses Ansible!</warning-title>
> If the terms "Ansible," "role," and "playbook" mean nothing to you, please check out
> [the Ansible introduction slides]({% link topics/admin/tutorials/ansible/slides.html %}) and [the Ansible introduction tutorial]({% link topics/admin/tutorials/ansible/tutorial.md %}).
{: .warning}

> <warning-title>Have you installed Galaxy?</warning-title>
> This tutorial builds upon the work in the [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, please ensure that you have completed that tutorial first.
{: .warning}

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Installing Docker

Galaxy Interactive Tools run in Docker containers. Thus, the first step is ensuring that the *nodes* that will run interactive tools have Docker installed. We will use the [Docker role][geerlingguy-docker] by the prolific Ansible Galaxy publisher, [Jeff Geerling (geerlingguy)][geerlingguy]. Have a look at the geerlingguy.docker [README][geerlingguy-docker-readme] and [defaults/main.yml][geerlingguy-docker-defaults] to get an understanding of what variables are used to control the role.

[geerlingguy-docker]: https://galaxy.ansible.com/geerlingguy/docker
[geerlingguy]: https://galaxy.ansible.com/geerlingguy
[geerlingguy-docker-readme]: https://github.com/geerlingguy/ansible-role-docker/blob/master/README.md
[geerlingguy-docker-defaults]: https://github.com/geerlingguy/ansible-role-docker/blob/master/defaults/main.yml

> <question-title></question-title>
>
> What variables might be relevant to using this role?
>
> > <solution-title></solution-title>
> >
> > The `docker_users` variable (a *list*) controls which users are able to interact with the Docker daemon, which our Galaxy user will need to do. Additionally, Docker Compose is configured by default, which we do not need, so it can be disabled with `docker_install_compose: false`.
> >
> {: .solution }
>
{: .question}

> <comment-title>Ansible Best Practices</comment-title>
> If you've set up your Galaxy server using the [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, you will have created a `galaxyservers` group in your inventory file, `hosts`, and placed your variables in `group_vars/galaxyservers.yml`. Although for the purposes of this tutorial, the Galaxy server and cluster node are one and the same, in a real world deployment they are very likely to be different hosts. We will continue to use the `galaxyservers` group for simplicity, but in your own deployment you should consider creating an additional group for cluster nodes.
{: .comment}

> <hands-on-title>Installing Docker with Ansible</hands-on-title>
>
> 1. In your working directory, add the Docker role to your `requirements.yml`:
>
>    ```yaml
>    - src: geerlingguy.docker
>      version: 7.8.0
>    ```
>
> 2. Install the requirements with `ansible-galaxy`:
>
>    ```
>    ansible-galaxy role install -p roles -r requirements.yml
>    ```
>
> 3. Edit the group variables file, `group_vars/galaxyservers.yml`:
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
>    > <question-title></question-title>
>    >
>    > {% raw %}
>    > Why is `"{{ galaxy_user.name }}"` specified instead of just the user `galaxy`?
>    > {% endraw %}
>    >
>    > > <solution-title></solution-title>
>    > > Duplicating values is never a good idea. If we needed to change the Galaxy user down the line or wanted to reuse this playbook on another host where the Galaxy username was different, we would have to change the value in multiple locations.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 4. Add the new role to the list of roles under the `roles` key in your playbook, `galaxy.yml`:
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
> 5. Run the playbook:
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
{: .hands_on}

Congratulations, you've set up Docker. Verify the installation using the `docker info` command (but keep in mind: what users did we authorize to interact with Docker?).

# Enabling the Interactive Tools Proxy

When an Interactive Tool's Docker container starts, it will be assigned a random port. In order to connect clients to the Interactive Tool, Galaxy needs to determine this port (and the node on which the tool is running) and configure a *proxy* from Galaxy to the GxIT's host and port. Consider the following example of running the Jupyter Notebook Interactive Tool, shown in Figure 1 below:

- nginx listens for requests from the client on **port 443** (https)
- Requests for Galaxy are delivered from nginx to Galaxy over a UNIX domain socket
- Requests for Interactive Tools are delivered from nginx to the Interactive Tools Proxy over (by default) **port 8000** (http)
  - GxIT http requests are forwarded by the proxy to Docker on the node on the container's (randomly assigned) **port 32768**
  - GxIT http requests are again forwarded by Docker to Jupyter on its in-container "published" **port 8888**

[//]: # The source for this figure can be found at: https://docs.google.com/presentation/d/1_4PtfM6A4mOxOlgGh6OGWvzFcxD1bdw4CydEWtm5n8k/

![Galaxy Interactive Tools Proxy Diagram](../../images/interactive-tools/gxit-proxy-diagram.png "Galaxy Interactive Tools Proxy Diagram")

As you can see, the client only ever speaks to nginx on the Galaxy server running on the standard https port (443), never directly to the interactive tool (which may be running on a node that does not even have a public IP address). By default, the mapping of GxIT invocation and its corresponding host/port is kept in a SQLite database known as the *Interactive Tools Session Map*, and the path to this database is important, since both Galaxy and the proxy need access to it.

The GxIT Proxy is written in [Node.js][nodejs]. A straightforward method to set it up is to use Galaxy's process manager, [Gravity](https://gravity.readthedocs.io/), to install and configure it.

[nodejs]: https://nodejs.org/

> <hands-on-title>Installing and configuring the Proxy using Gravity</hands-on-title>
>
> 1. Edit the group variables file, `group_vars/galaxyservers.yml`; the GxIT Proxy configuration lives under the `gravity['galaxy_config']['gx_it_proxy']` section. The available settings are:
>
>
>    | Setting         | Type    | Description                                                                                                                                                                                                                                         |
>    |-----------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
>    | `enable`        | bool    | Set to `true` to enable the GxIT Proxy. Defaults to `false`.                                                                                                                                                                                           |
>    | `version`       | str     | Version specifier for the GxIT Proxy Node.js package. Defaults to `>=0.0.6`.                                                                                                                                                                         |
>    | `ip`            | str     | Public-facing IP of the proxy. Defaults to `localhost`.                                                                                                                                                                                             |
>    | `port`          | int     | Public-facing port of the proxy. Defaults to `4002`.                                                                                                                                                                                                |
>    | `sessions`      | str     | Path to the Interactive Tools Session Map database. Defaults to `database/interactivetools_map.sqlite`.                                                                                                                                             |
>    | `verbose`       | bool    | Include verbose messages in the GxIT Proxy. Defaults to `true`.                                                                                                                                                                                      |
>    | `forward_ip`    | str     | Forward all requests to this IP. This is an advanced option that is only needed when proxying to remote interactive tool container that cannot be reached through the local network.                                                                |
>    | `forward_port`  | int     | Forward all requests to this port. This is an advanced option that is only needed when proxying to remote interactive tool container that cannot be reached through the local network.                                                              |
>    | `reverse_proxy` | str     | Rewrite location blocks with proxy port. This is an advanced option that is only needed when proxying to remote interactive tool container that cannot be reached through the local network. Defaults to `false`.                                   |
>    | `umask`         | str     | Permissions mask under which service should be executed.                                                                                                                                                                                            |
>    | `start_timeout` | int     | Value of supervisor startsecs, systemd TimeoutStartSec. Defaults to `10`.                                                                                                                                                                           |
>    | `stop_timeout`  | int     | Value of supervisor stopwaitsecs, systemd TimeoutStopSec. Defaults to `10`.                                                                                                                                                                         |
>    | `memory_limit`  | float   | Memory limit (in GB). If the service exceeds the limit, it will be killed. Default is no limit or the value of the ``memory_limit`` setting at the top level of the Gravity configuration, if set. Ignored if ``process_manager`` is ``supervisor``. |
>    | `environment`   | mapping | Extra environment variables and their values to set when running the service. A dictionary where keys are the variable names.                                                                                                                       |
>
>    Configure the settings in the `group_vars/galaxyservers.yml` file:
>
>    {% raw %}
>    ```yaml
>    # ... other settings ... #
>    galaxy_config:
>      galaxy:
>        # ... other settings ... #
>        # Interactive Tools
>        interactivetools_enable: true
>        interactivetools_map: "{{ gxit_proxy_sessions_path }}"
>        # ... other settings ... #
>      gravity:
>        # ... other settings ... #
>        gx_it_proxy:
>          enable: true
>          port: 4002
>          sessions: "{{ gxit_proxy_sessions_path }}"
>        # ... other settings ... #
>      # ... other settings ... # 
>    # ... other settings ... #
>    gxit_proxy_sessions_path: "{{ galaxy_mutable_data_dir }}/interactivetools_map.sqlite"
>    ```
>    {% endraw %}
>
> 2. Run the playbook:
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
>
{: .hands_on}

[nodeenv]: https://github.com/ekalinin/nodeenv
[path-based_interactive_tools]: https://docs.galaxyproject.org/en/master/admin/special_topics/interactivetools.html#path-based-interactivetools

Because the proxy runs as a systemd service, you can inspect the log of the service using `journalctl`. The service name is `galaxy-gx-it-proxy`:

```console
$ sudo journalctl -eu galaxy-gx-it-proxy
Feb 14 17:38:49 gcc-4 systemd[1]: Started galaxy-gx-it-proxy.service - Galaxygx-it-proxy
...
Feb 14 17:38:49 gcc-4 node[3679]: Watching path /srv/galaxy/var/interactivetools_map.sqlite
Feb 14 17:38:49 gcc-4 node[3679]: Listening on localhost:4002
```

> <comment-title>Note</comment-title>
>
> You can ignore errors about failing to read the sessions map file for now - Galaxy will create it when it's needed.
>
{: .comment}

# Proxying the Proxy

As explained in the previous section, we will proxy the Interactive Tools Proxy with nginx so that it can serve requests on the standard HTTPS port, 443. Because we've configured nginx with Ansible, this is relatively simple.

> <hands-on-title>Configuring NGINX</hands-on-title>
>
> 1. Edit the group variables file, `group_vars/galaxyservers.yml` and add a new item to the **existing** `nginx_ssl_servers` so it matches:
>
>    {% raw %}
>    ```yaml
>    nginx_ssl_servers:
>      - galaxy
>      - galaxy-gxit-proxy
>    ```
>    {% endraw %}
>
>    The nginx configuration `galaxy-gxit-proxy` doesn't exist yet, but we'll create it in a moment.
>
> 2. Create `templates/nginx/galaxy-gxit-proxy.j2` with the following contents:
>
>    {% raw %}
>    ```nginx
>    server {
>        # Listen on port 443
>        listen       *:443 ssl;
>        # Match all requests for the interactive tools subdomain
>        server_name  *.interactivetool.{{ inventory_hostname }};
>
>        # Our log files will go here.
>        access_log  syslog:server=unix:/dev/log;
>        error_log   syslog:server=unix:/dev/log;
>
>        # Proxy all requests to the GxIT Proxy application
>        location / {
>            proxy_redirect off;
>            proxy_http_version 1.1;
>            proxy_set_header Host $host;
>            proxy_set_header X-Real-IP $remote_addr;
>            proxy_set_header Upgrade $http_upgrade;
>            proxy_set_header Connection "upgrade";
>            proxy_pass http://localhost:{{ galaxy_config.gravity.gx_it_proxy.port }};
>        }
>    }
>    ```
>    {% endraw %}
>
> 3. To enable [Path-based Interactive Tools][path-based_interactive_tools], open `templates/nginx/galaxy.j2` and add the following contents:
>
>    {% raw %}
>    ```nginx
>    server {
>        # ... existing settings ...
>
>        # Route all path-based interactive tool requests to the InteractiveTool proxy application
>        location ~* ^/(interactivetool/.+)$ {
>            proxy_redirect off;
>            proxy_http_version 1.1;
>            proxy_set_header Host $host;
>            proxy_set_header X-Real-IP $remote_addr;
>            proxy_set_header Upgrade $http_upgrade;
>            proxy_set_header Connection "upgrade";
>            proxy_pass http://localhost:{{ galaxy_config.gravity.gx_it_proxy.port }};
>        }
>
>        # ... other existing settings ...
>    }
>    ```
>    {% endraw %}
>
> 4. Run the playbook:
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
>
{: .hands_on}

# Getting a Wildcard SSL Certificate

During the [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, we acquired an SSL certificate for our Galaxy server from [Let's Encrypt][lets-encrypt]. This certificate was issued for the hostname of your Galaxy server (e.g. `galaxy.example.org`). SSL certificates are valid *only for the name to which they were issued*. This presents a problem for us due to the way that Galaxy Interactive Tools work.

In order to ensure each Interactive Tool's cookies are unique, and to provide each tool with a unique entry point, they are served from a subdomain of your Galaxy server (e.g. `<unique-id>.ep.interactivetool.galaxy.example.org`). Your SSL cert is not valid for this subdomain. Further, in order to support the random `<unique-id>` in the hostname, we need a *wildcard certificate* for `*.ep.interactivetool.galaxy.example.org`.

This process is highly dependent on your site; specifically, your SSL certificate vendor, and your DNS server software or cloud provider.

Let's Encrypt, the SSL certificate vendor we use in our tutorials, [can only generate wildcard certificates using the DNS-01 challenge method][lets-encrypt-faq], which works by issuing a [dynamic DNS][ddns] update to set the requested domain's `TXT` record.

If you are completing this tutorial as part of a [Galaxy Admin Training][gat] course, we might have precreated a dynamic DNS server that you will use for this step. The *TSIG key* that allows you to perform dynamic DNS updates will be provided to you. Your instructor will also tell you which option to follow (1 or 2), depending on the DNS provider that was chosen for this course.

As we use Let's Encrypt in staging mode, the wildcard certificates generated with either option 1 or 2 will still be invalid, and you will still see a warning in your web browser when accessing an Interactive Tool. If this warning is not a problem for you, you can just skip this section of the tutorial, and move on to "Enabling Interactive Tools in Galaxy".

> <hands-on-title>Requesting a Wildcard Certificate with Certbot using Ansible - Option 1 (rfc2136)</hands-on-title>
>
> This method uses a DNS provider hosted by the Galaxy Project.
>
> 1. Edit the group variables file, `group_vars/galaxyservers.yml`:
>
>    The relevant variables to set for this role are:
>
>    | Variable                  | Type       | Description                                                                                                      |
>    | ----------                | -------    | -------------                                                                                                    |
>    | `certbot_domains`         | list       | List of domains to include as subject alternative names (the first will also be the certificate's *common name*) |
>    | `certbot_dns_provider`    | string     | Name of [Certbot DNS plugin][certbot-dns-plugins] to use                                                         |
>    | `certbot_dns_credentials` | dictionary | Plugin-specific credentials for performing dynamic DNS updates                                                   |
>    | `certbot_expand`          | boolean    | Whether to "expand" an existing certificate (add new domain names to it)                                         |
>
>    - Add a new item to the **existing** `certbot_domains` list so it matches:
>
>      {% raw %}
>      ```yaml
>      certbot_domains:
>        - "{{ inventory_hostname }}"
>        - "*.ep.interactivetool.{{ inventory_hostname }}"
>      ```
>      {% endraw %}
>
>    - Comment out the existing `certbot_auth_method` like so:
>
>      ```yaml
>      #certbot_auth_method: --webroot
>      ```
>
>        Although this is not explicitly required (setting `certbot_dns_provider` as we do overrides this setting), doing so is less confusing in the future, since it makes it clear that the "webroot" method for Let's Encrypt WEB-01 challenges is no longer in use for this server.
>
>    - Add the following lines to your `group_vars/galaxyservers.yml` file:
>
>      ```yaml
>      certbot_dns_provider: rfc2136
>      certbot_dns_credentials:
>        server: ns-training.galaxyproject.org
>        port: 53
>        name: certbot-training.
>        secret: <SECRET PROVIDED BY INSTRUCTOR>
>        algorithm: HMAC-SHA512
>      ```
>
> 2. Run the playbook **with `certbot_expand`**:
>
>    ```
>    ansible-playbook galaxy.yml -e certbot_expand=true
>    ```
>
>    > <question-title></question-title>
>    >
>    > What is the `-e` flag to `ansible-playbook` and why did we use it?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > As per `ansible-playbook --help`:
>    > >
>    > > ```
>    > >   -e EXTRA_VARS, --extra-vars EXTRA_VARS
>    > >                         set additional variables as key=value or YAML/JSON, if
>    > >                         filename prepend with @
>    > > ```
>    > >
>    > > We used this flag because `certbot_expand` only needs to be set *once*, when we are adding a new domain to the certificate. It should not be enabled on subsequent runs of the playbook, or else we would request a new certificate on each run! Thus, it does not make sense to add it to a vars file.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
>    Be patient! The certificate request step can take time due to the time allowed for DNS propagation to occur.
>
{: .hands_on}


> <hands-on-title>Requesting a Wildcard Certificate with Certbot using Ansible - Option 2 (route53)</hands-on-title>
>
> This method uses route53, the Amazon Web Services DNS provider. To manage connection to AWS, we will first install a specific role.
>
> 1. In your working directory, add the aws_cli role to your `requirements.yml`:
>
>    ```yaml
>    - src: usegalaxy_eu.aws_cli
>      version: 0.0.1
>    ```
>
> 2. Install the requirements with `ansible-galaxy`:
>
>    ```
>    ansible-galaxy role install -p roles -r requirements.yml
>    ```
> 3. Open `galaxy.yml` with your text editor to add the role `usegalaxy_eu.aws_cli` just before the nginx role:
>
>    ```diff
>    diff --git a/galaxy.yml b/galaxy.yml
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -21,6 +21,7 @@
>           become: true
>           become_user: galaxy
>         - usegalaxy_eu.galaxy_systemd
>    +    - usegalaxy_eu.aws_cli
>         - galaxyproject.nginx
>         - geerlingguy.docker
>    ```
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 4. Edit the group variables file, `group_vars/galaxyservers.yml`:
>
>    The relevant variables to set for this role are:
>
>    | Variable                  | Type       | Description                                                                                                      |
>    | ----------                | -------    | -------------                                                                                                    |
>    | `certbot_domains`         | list       | List of domains to include as subject alternative names (the first will also be the certificate's *common name*) |
>    | `certbot_dns_provider`    | string     | Name of [Certbot DNS plugin][certbot-dns-plugins] to use                                                         |
>    | `certbot_dns_credentials` | dictionary | Plugin-specific credentials for performing dynamic DNS updates                                                   |
>    | `certbot_expand`          | boolean    | Whether to "expand" an existing certificate (add new domain names to it)                                         |
>
>    - Add a new item to the **existing** `certbot_domains` list so it matches:
>
>      {% raw %}
>      ```yaml
>      certbot_domains:
>        - "{{ inventory_hostname }}"
>        - "*.ep.interactivetool.{{ inventory_hostname }}"
>      ```
>      {% endraw %}
>
>    - Comment out the existing `certbot_auth_method` like so:
>
>      ```yaml
>      #certbot_auth_method: --webroot
>      ```
>
>        Although this is not explicitly required (setting `certbot_dns_provider` as we do overrides this setting), doing so is less confusing in the future, since it makes it clear that the "webroot" method for Let's Encrypt WEB-01 challenges is no longer in use for this server.
>
>    - Add the following lines to your `group_vars/galaxyservers.yml` file:
>
>      ```yaml
>      certbot_dns_provider: route53
>      aws_cli_credentials:
>        - access_key: "<SECRET PROVIDED BY INSTRUCTOR>"
>          secret_key: "<SECRET PROVIDED BY INSTRUCTOR>"
>          homedir: /root
>          owner: root
>          group: root
>      ```
>
> 5. Run the playbook **with `certbot_expand`**:
>
>    ```
>    ansible-playbook galaxy.yml -e certbot_expand=true
>    ```
>
>    > <question-title></question-title>
>    >
>    > What is the `-e` flag to `ansible-playbook` and why did we use it?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > As per `ansible-playbook --help`:
>    > >
>    > > ```
>    > >   -e EXTRA_VARS, --extra-vars EXTRA_VARS
>    > >                         set additional variables as key=value or YAML/JSON, if
>    > >                         filename prepend with @
>    > > ```
>    > >
>    > > We used this flag because `certbot_expand` only needs to be set *once*, when we are adding a new domain to the certificate. It should not be enabled on subsequent runs of the playbook, or else we would request a new certificate on each run! Thus, it does not make sense to add it to a vars file.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}

You can verify that your certificate has been expanded using your browser's developer tools:

![Wildcard Certificate Dialog](../../images/interactive-tools/wildcard-cert.png "Wildcard Certificate Dialog")

# Enabling Interactive Tools in Galaxy

A few Interactive Tool wrappers are provided with Galaxy, but they are [commented out in Galaxy's default tool config file][tool-conf-gxits]. As a result, we need to instruct the [galaxyproject.galaxy role][galaxy-role] to install a tool panel configuration file containing at least one of these tools in order to try them out. For the purposes of this tutorial, a good choice is the [EtherCalc][ethercalc] GxIT, because it has a relatively small [Docker image][ethercalc-docker-image].

[tool-conf-gxits]: https://github.com/galaxyproject/galaxy/blob/9d788cb144a18570e7e5e948ac1f7bc04fec4e75/lib/galaxy/config/sample/tool_conf.xml.sample#L126
[galaxy-role]: https://github.com/galaxyproject/ansible-galaxy
[ethercalc]: https://ethercalc.net/
[ethercalc-docker-image]: https://hub.docker.com/r/shiltemann/ethercalc-galaxy-ie

> <hands-on-title>Enabling Interactive Tools in Galaxy</hands-on-title>
>
> 1. Rather than modifying the default tool configuration file, we'll add a new one that only references the Interactive Tools. This way, the default set of tools will still load without us having to incorporate the entire default tool config into our playbook.
>
>    If the folder does not exist, create `templates/galaxy/config` next to your `galaxy.yml` (`mkdir -p templates/galaxy/config/`)
>
>    Create `templates/galaxy/config/tool_conf_interactive.xml.j2` with the following contents:
>
>    ```xml
>    <toolbox monitor="true">
>        <section id="interactivetools" name="Interactive Tools">
>            <tool file="interactive/interactivetool_ethercalc.xml" />
>            <tool file="interactive/interactivetool_jupyter_notebook_1.0.1.xml" />
>        </section>
>    </toolbox>
>    ```
>
> 2. We need to instruct Galaxy on how to run Interactive Tools (and specifically, how to run them in Docker). Let's begin with a basic job configuration:
>
>    ```yaml
>    # ... other settings ... #
>    # Galaxy Job Configuration
>    galaxy_job_config:
>      runners:
>        local_runner:
>          load: galaxy.jobs.runners.local:LocalJobRunner
>          workers: 4
>      execution:
>        default: local_env
>        environments:
>          local_env:
>            runner: local_runner
>    # ... other settings ... #
>    ```
>
>    > <comment-title>Note</comment-title>
>    > Depending on the order in which you are completing this tutorial in relation to other tutorials, you may have already defined the `galaxy_job_config` mapping in your `group_vars/galaxyservers.yml` file, as well as `galaxy_config_templates`. If this is the case, be sure to **merge the changes in this section with your existing playbook**.
>    {: .comment}
>
> 3. Next, we need to configure the interactive tools destination. Add a destination for submitting jobs as Docker containers using the [advanced sample job configuration][job-conf-docker] as a guide. Finally, use the [EtherCalc GxIT's][ethercalc-tool-wrapper] tool ID to route executions of the EtherCalc GxIT to the newly created destination:
>
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>         execution:
>           default: local_env
>           environments:
>             local_env:
>               runner: local_runner
>    +        local_interactive_env:
>    +          runner: local_runner
>    +          params:
>    +            docker_enabled: true
>    +            docker_volumes: $defaults
>    +            docker_sudo: false
>    +            docker_net: bridge
>    +            docker_auto_rm: true
>    +            docker_set_user: ""
>    +            require_container: true
>    +    tools:
>    +      - id: interactive_tool_ethercalc
>    +        environment: local_interactive_env
>    +      - id: interactive_tool_jupyter_notebook
>    +        environment: local_interactive_env
>    ```
>
>    Of considerable note is the `docker_volumes` param: the variable expansions are explained in the [advanced sample job configuration][job-conf-docker].  We'll use this volume configuration for now but it has some considerable data security problems. We'll discuss a better solution at the end of this tutorial.
>
> 4. Inform `galaxyproject.galaxy` of what tool configuration files to load in your group variables (`group_vars/galaxyservers.yml`):
>
>    {% raw %}
>    ```yaml
>    galaxy_tool_config_files:
>      - "{{ galaxy_server_dir }}/config/tool_conf.xml.sample"
>      - "{{ galaxy_config_dir }}/tool_conf_interactive.xml"
>    ```
>    {% endraw %}
>
>    Next, deploy the new config template using the `galaxy_config_templates` var in your group vars:
>
>    {% raw %}
>    ```yaml
>    galaxy_config_templates:
>      # ... possible existing config file definitions
>      - src: templates/galaxy/config/tool_conf_interactive.xml.j2
>        dest: "{{ galaxy_config_dir }}/tool_conf_interactive.xml"
>    ```
>    {% endraw %}
>
> 5. If you've configured container resolvers, make sure the explit docker container resolver is enabled, e.g in templates/galaxy/config/container_resolvers_conf.yml.j2
>    {% raw %}
>    ```diff
>    --- a/templates/galaxy/config/container_resolvers_conf.yml.j2
>    +++ b/templates/galaxy/config/container_resolvers_conf.yml.j2
>    @@ -9,3 +9,4 @@
>     - type: build_mulled_singularity
>       auto_install: False
>       cache_directory: "{{ galaxy_mutable_data_dir }}/cache/singularity/built/"
>     +- type: explicit
>    ```
>    {% endraw %}
>
> 6. Run the playbook:
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
>
> 7. Follow the Galaxy logs with `galaxyctl follow`.
>
{: .hands_on}

[job-conf-docker]: https://github.com/galaxyproject/galaxy/blob/a2c847439528df0becf9e382c323f1dfe34d9633/lib/galaxy/config/sample/job_conf.xml.sample_advanced#L495
[ethercalc-tool-wrapper]: https://github.com/galaxyproject/galaxy/blob/6622ad1acb91866febb3d2f229de7cfb8af3a9f6/tools/interactive/interactivetool_ethercalc.xml

# Run an Interactive Tool

You should now be ready to run an Interactive Tool in Galaxy!

> <hands-on-title>Running an Interactive Tool</hands-on-title>
>
> 1. Ensure that you are logged in to your Galaxy server by checking the **User** menu in the masthead.
> 2. We'll need an input for our test GxIT (EtherCalc). Any tabular file can be used, such as Galaxy's [1.tabular][1-tabular] test data. Copy this file's URL:
>
>    ```
>    https://raw.githubusercontent.com/galaxyproject/galaxy/release_20.01/test-data/1.tabular
>    ```
>
> 3. Click {% icon galaxy-upload %} **Upload** at the top of the tool panel (on the left side of the Galaxy UI).
> 4. In the resulting modal dialog, click the **Paste/Fetch data** button.
> 5. Paste the URL in the **text field** that has just appeared.
> 6. Give the new dataset a name such as `tabular`, if you like.
> 7. Click **Start** and then **Close**.
> 8. From the tool menu, click the **Interactive Tools** section, then click **EtherCalc** {% icon tool %}.
> 9. Ensure that your newly uploaded tabular dataset is selected as the input **Some tabular dataset**, then click **Execute**.
> 10. Monitor the blue **info box** on the next page, which will inform you when the Interactive Tool is accessible and provide you with a link to access it.
>
>     <br/>
>
>     If you navigate away from this page, you can view your running Interactive Tools from the **Active InteractiveTools** menu item in the **User** menu.
> 11. Click the **click here to display** link.
>
{: .hands_on}

If everything has worked correctly, your browser will load EtherCalc with your tabular data preloaded. Once you're done working with the data, return to Galaxy and stop EtherCalc by deleting its output dataset from your history, or stopping it via the interface from the **Active InteractiveTools** menu item in the **User** menu.

[1-tabular]: https://raw.githubusercontent.com/galaxyproject/galaxy/release_20.01/test-data/1.tabular

# Securing Interactive Tools

Inspecting the Docker container of a running Interactive Tool shows the volume configuration expanded from `$galaxy_root` in the job destination's `docker_volumes` param:

```console
$ docker inspect $(docker ps -q) | jq '.[0].HostConfig.Binds'
[
  "/srv/galaxy/server:/srv/galaxy/server:ro",                                     # Galaxy server dir
  "/srv/galaxy/server/tools/interactive:/srv/galaxy/server/tools/interactive:ro", # EtherCalc tool wrapper parent dir
  "/srv/galaxy/jobs/000/1:/srv/galaxy/jobs/000/1:ro",                             # Per-job root dir
  "/srv/galaxy/jobs/000/1/outputs:/srv/galaxy/jobs/000/1/outputs:rw",             # Job outputs dir
  "/srv/galaxy/jobs/000/1/configs:/srv/galaxy/jobs/000/1/configs:rw",             # Job config files dir
  "/srv/galaxy/jobs/000/1/working:/srv/galaxy/jobs/000/1/working:rw",             # Job working (cwd) dir
  "/data:/data:rw",                                                               # GALAXY USER DATASETS DIR (RW!)
  "/srv/galaxy/server/tool-data:/srv/galaxy/server/tool-data:ro"                  # Galaxy reference data dir
]
```

As hinted earlier, there is a concerning state here: The directory containing **all** of the the user-generated data in Galaxy (not just the data for this job) has been mounted **read-write** in to the container. **This configuration grants users running interactive tools full access to all the data in Galaxy, which is a very bad idea.** Unlike standard Galaxy tools, where the tool's design prevents users from writing to arbitrary paths, Interactive Tools are fully user controllable. Although EtherCalc does not provide a mechanism for writing to this path, other Interactive Tools (such as Jupyter Notebook) do.

Two solutions are discussed in the [advanced sample job configuration][job-conf-docker]:

1. Use the `outputs_to_working_directory` job configuration option, which allows you to mount datasets read-only: this prevents manipulation, but still allows GxIT users to read any dataset in your Galaxy server.
2. Use [Pulsar][pulsar], Galaxy's remote job execution engine, to provide full job isolation: this avoids all access to Galaxy data, with the performance penalty of copying input dataset(s) to the job directory.

Because we want to maintain dataset privacy, Pulsar is the better choice here. And in fact, we don't even need to set up a Pulsar server: because we only need Pulsar's input staging and isolation features, we can use [Embedded Pulsar][job-conf-pulsar-embedded], which runs a Pulsar server within the Galaxy application to perform these tasks. Embedded Pulsar can even interface with your distributed resource manager (aka cluster scheduler) of choice, as long as your Galaxy server and cluster both have access to a common filesystem (otherwise, you will need to use Pulsar in standalone mode; see the [Running Jobs on Remote Resources with Pulsar]({% link topics/admin/tutorials/pulsar/tutorial.md %}) tutorial).

[pulsar]: https://github.com/galaxyproject/pulsar
[job-conf-pulsar-embedded]: https://github.com/galaxyproject/galaxy/blob/6622ad1acb91866febb3d2f229de7cfb8af3a9f6/lib/galaxy/config/sample/job_conf.xml.sample_advanced#L106

> <hands-on-title>Running Interactive Tools with Embedded Pulsar</hands-on-title>
>
> 1. Create a configuration file *template* for the Pulsar application at `templates/galaxy/config/pulsar_app.yml.j2`.
>
>    <br/>
>
>    If the folder does not exist, create `templates/galaxy/config` next to your `galaxy.yml` (`mkdir -p templates/galaxy/config/`).
>
>    <br/>
>
>    Add the following contents to the template:
>
>    {% raw %}
>    ```yaml
>    ---
>
>    # The path where per-job directories will be created
>    staging_directory: "{{ galaxy_job_working_directory }}/_interactive"
>
>    # Where Pulsar state information will be stored (e.g. currently active jobs)
>    persistence_directory: "{{ galaxy_mutable_data_dir }}/pulsar"
>
>    # Where to find Galaxy tool dependencies
>    tool_dependency_dir: "{{ galaxy_tool_dependency_dir }}"
>
>    # How to run jobs (see https://pulsar.readthedocs.io/en/latest/job_managers.html)
>    managers:
>      _default_:
>        type: queued_python
>        num_concurrent_jobs: 1
>    ```
>    {% endraw %}
>
> 2. Modify the job configuration in `group_vars/galaxyservers.yml` to configure Interactive Tools to use the embedded Pulsar runner.
>
>    <br/>
>
>    **Add** the embedded Pulsar runner to the `galaxy_job_config` section of the file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>         runners:
>           local_runner:
>             load: galaxy.jobs.runners.local:LocalJobRunner
>             workers: 4
>    +      pulsar_embedded:
>    +        load: galaxy.jobs.runners.pulsar:PulsarEmbeddedJobRunner
>    +        pulsar_config: "{{ galaxy_config_dir }}/pulsar_app.yml"
>         execution:
>    ```
>    {% endraw %}
>
>    Next, **modify** the `interactive_local` destination to use the new runner and set the new parameter `container_monitor_result` to `callback` (explained in more detail in the next step):
>
>    > <warning-title>Untrusted SSL Certificates</warning-title>
>    > If you are completing this tutorial as part of a [Galaxy Admin Training][gat] course, you will also need the `env` setting shown below to prevent problems with the untrusted SSL certificates in use during the course. Galaxy servers with valid SSL certificates *do not need this option*.
>    {: .warning}
>
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>         execution:
>           default: local_env
>           environments:
>             local_env:
>               runner: local_runner
>             local_interactive_env:
>    -          runner: local_runner
>    +          runner: pulsar_embedded
>               params:
>                 docker_enabled: true
>                 docker_volumes: $defaults
>                 docker_sudo: false
>                 docker_net: bridge
>                 docker_auto_rm: true
>                 docker_set_user: ""
>                 require_container: true
>    +            container_monitor_result: callback
>    +          env:
>    +            - name: REQUESTS_CA_BUNDLE
>    +              value: /etc/ssl/certs/ca-certificates.crt
>         tools:
>    ```
>
> 3. Open your `galaxyservers` group variables file and instruct `galaxyproject.galaxy` to install the Pulsar configuration file:
>
>    > <comment-title>Note</comment-title>
>    > Depending on the order in which you are completing this tutorial in relation to other tutorials, you may have already defined `galaxy_config_templates`. If this is the case, be sure to **merge the changes in this step with your existing playbook**.
>    {: .comment}
>
>    {% raw %}
>    ```yaml
>    galaxy_config_templates:
>      - src: templates/galaxy/config/pulsar_app.yml.j2
>        dest: "{{ galaxy_config_dir }}/pulsar_app.yml"
>    ```
>    {% endraw %}
>
>    Additionally, you will need to set the `galaxy_infrastructure_url` config option:
>
>    {% raw %}
>    ```yaml
>    galaxy_config:
>      galaxy:
>        # ... existing configuration options in the `galaxy` section ...
>        galaxy_infrastructure_url: "https://{{ inventory_hostname }}/"
>    ```
>    {% endraw %}
>
>    > <details-title>Detail: Infrastructure URL/Callback</details-title>
>    > Galaxy must be made aware of the randomly selected port Docker has assigned after the GxIT begins operating, in order to update the proxy map. By default, this is done by writing a JSON file in the job directory. This method does not work with Pulsar since Pulsar uses a different job directory from the Galaxy job directory. As a result, Pulsar jobs use the `callback` method configured in the previous step to make a request to Galaxy's API, the URL for which is set in `galaxy_infrastructure_url`.
>    {: .details}
>
> 4. Run the playbook:
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
>
{: .hands_on}

Once the playbook run is complete and your Galaxy server has restarted, run the EtherCalc Interactive Tool again.

> <question-title></question-title>
>
> Once EtherCalc is running, check the mounts of its container. What do you observe?
>
> > <solution-title></solution-title>
> >
> > ```console
> > $ docker inspect $(docker ps -q) | jq '.[0].HostConfig.Binds'
> > [
> >   "/srv/galaxy/jobs/_interactive/2:/srv/galaxy/jobs/_interactive/2:ro",                       # Per-job root dir
> >   "/srv/galaxy/jobs/_interactive/2/tool_files:/srv/galaxy/jobs/_interactive/2/tool_files:ro", # EtherCalc tool wrapper parent dir
> >   "/srv/galaxy/jobs/_interactive/2/outputs:/srv/galaxy/jobs/_interactive/2/outputs:rw",       # Job outputs dir
> >   "/srv/galaxy/jobs/_interactive/2/working:/srv/galaxy/jobs/_interactive/2/working:rw",       # Job working (cwd) dir
> >   "/srv/galaxy/server/tool-data:/srv/galaxy/server/tool-data:ro"                              # Galaxy reference data dir
> > ]
> > ```
> >
> > Of note, the user data directory, `/data`, is no longer mounted in the container!
> >
> {: .solution }
>
{: .question }

# High availability setup with PostgresSQL (Optional)

> <comment-title></comment-title>
> This section is **only relevant if you are running a high-availability** setup, meaning that you have multiple copies of Galaxy running behind a load balancer.
>
> If you have installed Galaxy following the [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, or are completing this tutorial as part of a [Galaxy Admin Training][gat] course, please skip this section, as you are then _not_ running a high-availability setup.
{: .comment}

In a _high availability_ setup, multiple redundant copies of Galaxy run simultaneously behind a load balancer to minimize downtime and service interruptions.

As explained in [one of the previous sections](#enabling-the-interactive-tools-proxy), the Galaxy Interactive Tools Proxy redirects requests to each Interactive Tool's host and port. By default, the mapping of GxIT invocations to their corresponding host/port is kept in a SQLite database known as the _Interactive Tools Session Map_.

By design, [SQLite is the wrong choice for high availability setups][sqlite_situations_where_a_client_server_rdbms_may_work_better], the showstopper being that the SQLite database file would have to be shared over a network filesystem, which are usually associated with too high latencies for RDBMS use. For this reason, Galaxy and the Interactive Tools Proxy can also store the **Session Map in a PostgreSQL database**. 

[sqlite_situations_where_a_client_server_rdbms_may_work_better]: https://www.sqlite.org/whentouse.html#situations_where_a_client_server_rdbms_may_work_better

> <hands-on-title>Preparing the database</hands-on-title>
> 
> First, you need to create a database for the Interactive Tools Proxy. 
>
> > <warning-title></warning-title>
> > Do **not** use the Galaxy database for this purpose. The main Galaxy database is reserved for Galaxy's core functionality, and Interactive Tools have not yet reached this stage. Since Galaxy does not expect to find the Interactive Tools Session Map in this database, storing it there can lead to errors.
> {: .warning }
>
> <br>
>
> 1. **On your database server**, access PostgresSQL and create a `gxitproxy` database to store the Interactive Tools Session Map. For simplicity, the same user that operates on the Galaxy main database, typically named `galaxy`, is also going to operate on this one and will own the new database.
> 
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > # one-liner that connects to Postgres, creates the database and assigns ownership
>    > sudo -u postgres createdb -O galaxy gxitproxy
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 2. Connect to the `gxitproxy` database as `galaxy`.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > sudo -iu galaxy psql -d gxitproxy
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    > <code-out-title>SQL</code-out-title>
>    > ```
>    > psql (10.12 (Ubuntu 10.12-0ubuntu0.18.04.1))
>    > Type "help" for help.
>    >
>    > gxitproxy=#
>    > ```
>    {: .code-out}
>
> 3. Create a `gxitproxy` table in the new database.
>
>    > <code-in-title>SQL</code-in-title>
>    > ```sql
>    > CREATE TABLE IF NOT EXISTS gxitproxy (key TEXT, key_type TEXT, token TEXT, host TEXT, port INTEGER, info TEXT, PRIMARY KEY (key, key_type));
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    > <code-out-title>SQL</code-out-title>
>    > ```
>    > CREATE TABLE
>    > ```
>    {: .code-out}
>
> 4. This is enough to let Galaxy and the Interactive Tool Proxy store the Interactive Tools Session Map in PostgreSQL. But there is a catch: when the Interactive Tool Proxy uses SQLite, it knows the database has changed because it watches the file for changes. When using Postgres, this mechanism is not available. By default, the proxy simply polls the database at regular intervals. To let the user access interactive tools as fast as possible, the proxy can also be notified of updates via [PostgreSQL asynchronous notifications](https://www.postgresql.org/docs/16/libpq-notify.html). To enable them, you have to create a PostgreSQL trigger that sends a NOTIFY message to the channel `gxitproxy` every time the table `gxitproxy` changes.
>
>    Run the following commands to create to create a function that sends a NOTIFY message to the channel `gxitproxy` and a trigger that runs the function every time the table `gxitproxy` changes.
>
>    > <code-in-title>SQL</code-in-title>
>    > ```sql
>    > CREATE OR REPLACE FUNCTION notify_gxitproxy()
>    > RETURNS trigger AS $$
>    > BEGIN
>    >   PERFORM pg_notify('gxitproxy', 'Table "gxitproxy" changed');
>    >   RETURN NEW;
>    > END;
>    > $$ LANGUAGE plpgsql;
>    >
>    > CREATE TRIGGER gxitproxy_notify
>    > AFTER INSERT OR UPDATE OR DELETE ON gxitproxy
>    > FOR EACH ROW EXECUTE FUNCTION notify_gxitproxy();
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
>    > <code-out-title>SQL</code-out-title>
>    > ```
>    > CREATE FUNCTION
>    > CREATE TRIGGER
>    > ```
>    {: .code-out}
>
{: .hands_on}

The next step is configuring Galaxy and the Interactive Tool Proxy to use the new database.

> <hands-on-title>Configure Galaxy and the Interactive Tool Proxy</hands-on-title>
>
> 1. Adjust your `group_vars/galaxyservers.yml` file as follows.
>
>    {% raw %}
>    ```yaml
>    # ... existing configuration options ... #
>
>    galaxy_config:
>      galaxy:
>        # ... existing configuration options in the `galaxy` section ...
>        # interactivetools_map: "{{  gxit_proxy_sessions_path }}"  # comment, remove or leave this line in place (it will be overridden by the option below)
>        interactivetoolsproxy_map: "{{ gxit_proxy_sessions_path }}"
>        # ... other existing configuration options in the `galaxy` section ...
>
>    # ... other existing configurations ... #
>
>    gxit_proxy_sessions_path: "postgresql:///gxitproxy?host=/var/run/postgresql"
>    ```
>    {% endraw %}
>
> 2. Run the playbook:
>
>    ```
>    ansible-playbook galaxy.yml
>    ```
>
{: .hands_on}

That's it, once the playbook run is complete, both Galaxy and the Interactive Tools Proxy will be storing the Interactive Tools Session Map in PostgreSQL.
