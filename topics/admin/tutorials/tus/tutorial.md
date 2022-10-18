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

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="tus" %}

# TUS and Galaxy

To allow your user to upload via TUS, you will need to:

- configure Galaxy to know where the files are uploaded.
- install TUSd
- configure Nginx to proxy TUSd

## Installing and Configuring

> <hands-on-title>Setting up ftp upload with Ansible</hands-on-title>
>
> 1. In your playbook directory, add the `galaxyproject.tusd` role to your `requirements.yml`
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
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Install the role with:
>
>    > <code-in-title>Bash</code-in-title>
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
>    @@ -60,6 +60,8 @@ galaxy_config:
>         allow_user_impersonation: true
>         # Tool security
>         outputs_to_working_directory: true
>    +    # TUS
>    +    tus_upload_store: /data/tus
>       gravity:
>         galaxy_root: "{{ galaxy_root }}/server"
>         app_server: gunicorn
>    @@ -139,3 +141,16 @@ nginx_conf_http:
>     nginx_ssl_role: usegalaxy_eu.certbot
>     nginx_conf_ssl_certificate: /etc/ssl/certs/fullchain.pem
>     nginx_conf_ssl_certificate_key: /etc/ssl/user/privkey-nginx.pem
>    +
>    +# TUS
>    +galaxy_tusd_port: 1080
>    +tusd_instances:
>    +  - name: main
>    +    user: "{{ galaxy_user.name }}"
>    +    group: "galaxy"
>    +    args:
>    +      - "-host=localhost"
>    +      - "-port={{ galaxy_tusd_port }}"
>    +      - "-upload-dir={{ galaxy_config.galaxy.tus_upload_store }}"
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
>    @@ -28,6 +28,22 @@ server {
>             proxy_set_header Upgrade $http_upgrade;
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
>         # request to Gunicorn which should be spending its time doing more useful
>         # things like serving Galaxy!
>    {% endraw %}
>    ```
>    {: data-commit="Proxy it via NGINX"}
>
> 5. Add to the end of your Galaxy playbook
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -19,3 +19,4 @@
>           become: true
>           become_user: "{{ galaxy_user.name }}"
>         - galaxyproject.nginx
>    +    - galaxyproject.tusd
>    {% endraw %}
>    ```
>    {: data-commit="Add the role to the playbook"}
>
> 6. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
>
{: .hands_on}

Congratulations, you've set up TUS for Galaxy.

## Check it works

> <hands-on-title>Check that it works.</hands-on-title>
>
> 1. SSH into your machine
>
> 2. Check the active status of tusd by `systemctl status tusd-main`.
>
> 3. Upload a small file! (Pasted text will not pass via TUS)
>
> 4. Check the directory `/data/tus/` has been created and it's contents
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > sudo tree /data/tus/
>    > ```
>    {: .code-in}
>
> 5. You'll see files in that directory, a file that's been uploaded and an 'info' file which contains metadata about the upload.
>
{: .hands_on}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}

{% snippet topics/admin/faqs/missed-something.md step=2 %}
