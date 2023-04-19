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
>    @@ -14,3 +14,6 @@
>     # gxadmin (used in cleanup, and later monitoring.)
>     - src: galaxyproject.gxadmin
>       version: 0.0.12
>    +# TUS (uploads)
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
>    @@ -67,6 +67,9 @@ galaxy_config:
>         # Tool security
>         outputs_to_working_directory: true
>         new_user_dataset_access_role_default_private: true # Make datasets private by default
>    +    # TUS
>    +    galaxy_infrastructure_url: "https://{{ inventory_hostname }}"
>    +    tus_upload_store: "{{ galaxy_tus_upload_store }}"
>       gravity:
>         process_manager: systemd
>         galaxy_root: "{{ galaxy_root }}/server"
>    @@ -87,6 +90,10 @@ galaxy_config:
>         celery:
>           concurrency: 2
>           loglevel: DEBUG
>    +    tusd:
>    +      enable: true
>    +      tusd_path: /usr/local/sbin/tusd
>    +      upload_dir: "{{ galaxy_tus_upload_store }}"
>         handlers:
>           handler:
>             processes: 2
>    @@ -156,3 +163,7 @@ nginx_conf_http:
>     nginx_ssl_role: usegalaxy_eu.certbot
>     nginx_conf_ssl_certificate: /etc/ssl/certs/fullchain.pem
>     nginx_conf_ssl_certificate_key: /etc/ssl/user/privkey-www-data.pem
>    +
>    +# TUS
>    +galaxy_tusd_port: 1080
>    +galaxy_tus_upload_store: /data/tus
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
>     		proxy_set_header Upgrade $http_upgrade;
>     	}
>     
>    +	location /api/upload/resumable_upload {
>    +		# Disable request and response buffering
>    +		proxy_request_buffering     off;
>    +		proxy_buffering             off;
>    +		proxy_http_version          1.1;
>    +
>    +		# Add X-Forwarded-* headers
>    +		proxy_set_header X-Forwarded-Host   $host;
>    +		proxy_set_header X-Forwarded-Proto  $scheme;
>    +
>    +		proxy_set_header Upgrade            $http_upgrade;
>    +		proxy_set_header Connection         "upgrade";
>    +		client_max_body_size        0;
>    +		proxy_pass http://localhost:{{ galaxy_tusd_port }}/files;
>    +	}
>    +
>     	# Static files can be more efficiently served by Nginx. Why send the
>     	# request to Gunicorn which should be spending its time doing more useful
>     	# things like serving Galaxy!
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
>    @@ -30,6 +30,7 @@
>             name: ['tmpreaper']
>           when: ansible_os_family == 'Debian'
>       roles:
>    +    - galaxyproject.tusd
>         - galaxyproject.galaxy
>         - role: galaxyproject.miniconda
>           become: true
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
> 2. Check the active status of tusd by `systemctl status galaxy-tusd`.
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

{% snippet topics/admin/faqs/git-commit.md page=page %}

{% snippet topics/admin/faqs/missed-something.md step=4 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="tus" %}
