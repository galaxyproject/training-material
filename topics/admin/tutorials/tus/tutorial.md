---
layout: tutorial_hands_on

title: "Performant Uploads with TUS"
zenodo_link: ""
questions:
objectives:
  - Setup TUSd
  - Configure Galaxy to use it to process uploads
time_estimation: "30M"
key_points:
  - Use TUS to make uploads more efficient, especially for large uploads over unstable connections.
contributors:
  - mvdbeek
  - hexylena
  - lldelisle
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

Here you'll learn to setup [TUS](https://tus.io/) an open source resumable file upload server to process uploads for Galaxy. We use an external process here to offload the main Galaxy processes for more important work and not impact the entire system during periods of heavy uploading.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# TUS and Galaxy

To allow your user to upload via TUS, you will need to:

- configure Galaxy to know where the files are uploaded.
- install TUSd
- configure Nginx to proxy TUSd

## Installing and Configuring

> ### {% icon hands_on %} Hands-on: Setting up ftp upload with Ansible
>
> 1. In your playbook directory, add the `galaxyproject.proftpd` role to your `requirements.yml`
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -12,3 +12,5 @@
>       version: 0.3.0
>     - src: usegalaxy_eu.certbot
>       version: 0.1.5
>    +- name: galaxyproject.tusd
>    +  version: 0.0.1
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement"}
>
> 2. Install the role with:
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Configure it in your group variables
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -122,3 +122,16 @@ nginx_conf_http:
>     nginx_ssl_role: usegalaxy_eu.certbot
>     nginx_conf_ssl_certificate: /etc/ssl/certs/fullchain.pem
>     nginx_conf_ssl_certificate_key: /etc/ssl/user/privkey-nginx.pem
>    +
>    +# TUS
>    +galaxy_tusd_port: 1080
>    +tusd_instances:
>    +  - name: main
>    +    user: "{{ galaxy_user.name }}"
>    +    group: "{{ galaxy_user.group }}"
>    +    args:
>    +      - "-host=localhost"
>    +      - "-port={{ galaxy_tusd_port }}"
>    +      - "-upload-dir=/data/tus/"
>    +      - "-hooks-http=https://{{ inventory_hostname }}/api/upload/hooks"
>    +      - "-hooks-http-forward-headers=X-Api-Key,Cookie"
>    {% endraw %}
>    ```
>    {: data-commit="Configure TUS in your group variables"}
>
> 4. We proxy the service next to Galaxy. As it resides "under" the Galaxy path, clients will send cookies and authentication headers to TUS, which it can use to process the uploads before telling Galaxy when they're done.
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -16,6 +16,22 @@ server {
>             include uwsgi_params;
>         }
>     
>    +    location /api/upload/resumable_upload {
>    +        # Disable request and response buffering
>    +        proxy_request_buffering     off;
>    +        proxy_buffering             off;
>    +        proxy_http_version          1.1;
>    +
>    +        # Add X-Forwarded-* headers
>    +        proxy_set_header X-Forwarded-Host   $host;
>    +        proxy_set_header X-Forwarded-Proto  $scheme;
>    +
>    +        proxy_set_header Upgrade            $http_upgrade;
>    +        proxy_set_header Connection         "upgrade";
>    +        client_max_body_size        0;
>    +        proxy_pass http://localhost:{{ galaxy_tusd_port }}/files;
>    +    }
>    +
>         # Static files can be more efficiently served by Nginx. Why send the
>         # request to uWSGI which should be spending its time doing more useful
>         # things like serving Galaxy!
>    {% endraw %}
>    ```
>    {: data-commit="Proxy it via NGINX"}
>
> 5. Run the playbook
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
>
{: .hands_on}

Congratulations, you've set up TUS for Galaxy.

## Check it works

> ### {% icon hands_on %} Hands-on: Check that it works.
>
> 1. SSH into your machine
>
> 2. Check the active status of proftpd by `systemctl status proftpd`.
>
> 3. Check the port has been correctly attributed by `sudo lsof -i -P -n`.
>
>    > ### {% icon question %} Question
>    >
>    > What do you see?
>    >
>    > > ### {% icon solution %} Solution
>    > > You should see all the ports used by the server. What interests us is the line with proftpd.
>    > > You should see TCP *:21 (LISTEN).
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 4. Check the directory `/data/uploads/` has been created and is empty.
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > sudo tree /data/uploads/
>    > ```
>    {: .code-in}
>
{: .hands_on}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}
