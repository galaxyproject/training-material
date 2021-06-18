---
layout: tutorial_hands_on

title: "Enable upload via FTP"
zenodo_link: ""
questions:
objectives:
  - Configure galaxy and install a FTP server.
  - Use an Ansible playbook for this.
time_estimation: "TODO"
key_points:
contributors:
  - ldelisle
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

# Overview
{:.no_toc}

This tutorial will guide you to setup a ftp server so galaxy user can use it to upload large files. Indeed, as written on the [galaxy community hub](https://galaxyproject.org/ftp-upload/), uploading data directly from the browser can be unreliable and cumbersome. FTP will allow you to monitor the upload status as well as resume interrupted transfers.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# FTP
FTP stands for File Transfer Protocol. It is a communication protocol. It requires a server (here our galaxy server) and a client (user computers). The FTP server requires to have at least 2 ports accessible from outside (one for the commands and one for the transfer). Usually the port for the command is 21. 
% Possibly we will change it:
Here we will use 55000 because the VM we use in the galaxy admin training does not have the port 21 open.
You also need to know which are the other ports open so you can use them for the transfer (PassivePorts). In this training we consider that 56k to 60k are open.

# FTP and Galaxy

To allow your user to upload via FTP, you will need to:

- configure Galaxy to know where the files are uploaded.
- install a FTP server
- allow your FTP server to read Galaxy's database so users can use their credential and upload in the good directory.

For secure transmission we will use SSL/TLS (FTPS), not the SSH File Transfer Protocol (SFTP) as the Galaxy users don't correspond to users on the machine.

## Installing and Configuring

Luckily for us, there is an ansible role written by the Galaxy Project for this purpose. It will install proftpd. Firstly, we need to install the role and then update our playbook for using it.

If the terms "Ansible", "role" and "playbook" mean nothing to you, please checkout [the Ansible introduction slides]({% link topics/admin/tutorials/ansible/slides.html %}) and [the Ansible introduction tutorial]({% link topics/admin/tutorials/ansible/tutorial.md %})

{% snippet topics/admin/faqs/ansible_local.md %}

> ### {% icon hands_on %} Hands-on: Setting the ftp upload with Ansible
>
> 1. In your working directory, add the proftpd role to your `requirements.yml`
>
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -36,3 +36,5 @@
>       version: 0.12.0
>     - src: usegalaxy_eu.tiaas2
>       version: 0.0.6
>    +- src: galaxyproject.proftpd
>    +  version: 0.3.1
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
> 3. As in this training we are using certbot, we ask a private key for proftpd. Add the following line to your `group_vars/galaxyserver.yml` file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -159,6 +162,7 @@ certbot_well_known_root: /srv/nginx/_well-known_root
>     certbot_share_key_users:
>       - nginx
>       - rabbitmq
>    +  - proftpd
>     certbot_post_renewal: |
>         systemctl restart nginx || true
>         systemctl restart rabbitmq-server || true
>    {% endraw %}
>    ```
>    {: data-commit="Add proftpd in certbot"}
>
> 4. We will setup galaxy to specify that there are files uploaded by ftp. Add the following line to your `group_vars/galaxyserver.yml` file in the galaxy_config/galaxy section:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -93,6 +93,9 @@ galaxy_config:
>         allow_user_impersonation: true
>         # Tool security
>         outputs_to_working_directory: true
>    +    # FTP
>    +    ftp_upload_dir: /uploads
>    +    ftp_upload_site: "{{ inventory_hostname }}"
>       uwsgi:
>         socket: 127.0.0.1:5000
>         buffer-size: 16384
>    {% endraw %}
>    ```
>    {: data-commit="Add ftp vars in galaxy"}
>
> To check the other options for setting galaxy relative to ftp, please check the [galaxy options](https://docs.galaxyproject.org/en/master/admin/galaxy_options.html?highlight=ftp_upload_site#ftp-upload-dir).
>
> 5. Then we will set the different variable for proftpd. Add the following line to your `group_vars/galaxyserver.yml` file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -223,6 +227,27 @@ rabbitmq_users:
>         password: "{{ vault_rabbitmq_password_vhost }}"
>         vhost: /pulsar/galaxy_au
>     
>    +# Proftpd:
>    +proftpd_galaxy_auth: yes
>    +galaxy_ftp_upload_dir: /uploads
>    +proftpd_display_connect: |
>    +  example.org FTP server
>    +
>    +  Unauthorized access is prohibited
>    +proftpd_create_ftp_upload_dir: yes
>    +proftpd_options:
>    +  - User: galaxy
>    +  - Group: galaxy
>    +  - Port: 55000
>    +proftpd_sql_db: galaxy@/var/run/postgresql
>    +proftpd_sql_user: galaxy
>    +proftpd_conf_ssl_certificate: /etc/ssl/certs/cert.pem
>    +proftpd_conf_ssl_certificate_key: /etc/ssl/user/privkey-proftpd.pem
>    +proftpd_global_options:
>    +  - PassivePorts: 56000 60000
>    +proftpd_use_mod_tls_shmcache: false
>    +proftpd_tls_options: NoSessionReuseRequired
>    +
>     # Telegraf
>     telegraf_plugins_extra:
>       listen_galaxy_routes:
>    {% endraw %}
>    ```
>    {: data-commit="Add proftpd variables"}
>
> Here is a description of the set variables:
>    | Variable             | Description                                                                                                                                                                    |
>    | ----------           | -------------                                                                                                                                                                  |
>    | `proftpd_galaxy_auth`               | Attempt to authenticate users against a Galaxy database. |
>    | `galaxy_ftp_upload_dir`             | Path to the Galaxy FTP upload directory, should match `ftp_upload_dir` in your Galaxy config.  |
>    | `proftpd_display_connect`           |  Message to display when users connect to the FTP server. This should be the message, not the path to a file.   |
>    | `proftpd_create_ftp_upload_dir`     | Whether to allow the role to create this with owner `galaxy_user`.  |
>    | `proftpd_options`                   | Any option for proftpd, we will just set up the user and group of the `galaxy_user`.  |
>    | `proftpd_sql_db`                    |  Database name to connect to for authentication info.                             |
>    | `proftpd_sql_user`                  |  (default: the value of galaxy_user): Value of the username parameter to SQLConnectInfo. |
>    | `proftpd_conf_ssl_certificate`      | Path on the remote host where the SSL certificate file is. |
>    | `proftpd_conf_ssl_certificate_key`  | Path on the remote host where the SSL private key file is. |
>    | `proftpd_global_options`            | Set arbitrary options in the <Global> context. We set here the PassivePorts range. |
>    | `proftpd_use_mod_tls_shmcache`      | By default proftpd uses `mod_tls_shmcache` which is not installed on the server so we just disable it. |
>    | `proftpd_tls_options`               | Additional options for tls. We will use NoSessionReuseRequired |
>
>    > ### {% icon tip %} Why NoSessionReuseRequired?
>    > `mod_tls` only accepts SSL/TLS data connections that reuse the SSL session of the control connection, as a security measure. Unfortunately, there are some clients (e.g. curl/Filezilla) which do not reuse SSL sessions.
>    > To relax the requirement that the SSL session from the control connection be reused for data connections we set `NoSessionReuseRequired`.
>    {: .tip}
>
> 6. Add the new role to the list of roles under the `roles` key in your playbook, `galaxy.yml`:
>
>    {% raw %}
>    ```diff
>    --- a/galaxy.yml
>    +++ b/galaxy.yml
>    @@ -31,6 +31,7 @@
>           become_user: "{{ galaxy_user.name }}"
>         - usegalaxy_eu.rabbitmq
>         - galaxyproject.nginx
>    +    - galaxyproject.proftpd
>         - galaxyproject.cvmfs
>         - galaxyproject.gxadmin
>         - dj-wasabi.telegraf
>    {% endraw %}
>    ```
>    {: data-commit="Add role to playbook"}
>
> 5. Run the playbook
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
{: .hands_on}

Congratulations, you've set up FTP for Galaxy.

## Check it works

> ### {% icon hands_on %} Hands-on: Checking proftpd from the server
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
>    > > You should see TCP *:55000 (LISTEN).
>    > > If you set up the port 55000.
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 4. Check the directory `/uploads/` has been created and is empty.
>
>    > ### {% icon code-in %} Input: Bash
>    > ```
>    > sudo tree /uploads/
>    > ```
>    {: .code-in}
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Checking galaxy detected the ftp possibility
>
> 1. Open your galaxy in a browser.
>
> 2. Log in with a user (FTP is only possible for logged sessions).
>
> 3. Click on the upload button.
> You should now see on the bottom "Choose FTP files"
>
> 4. Click on the Choose FTP files button
> You should see a message "Your FTP directory does not contain any files."
>
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Upload your first file
>
> 1. Follow the [tutorial](https://galaxyproject.org/ftp-upload/) to upload a file.
> Don't forget to specify the port to 55000.
> You will have a message which ask you to approve the certificate, approve it.
>
>    > ### {% icon tip %} If you don't have a FTP client installed?
>    > You can use locally lftp to test the ftp.
>    > Install lftp with `sudo apt-get install lftp`.
>    > Add the public certificate to the list of known:
>    >    > ### {% icon code-in %} Input: Bash
>    >    > ```
>    >    > mkdir .lftp
>    >    > echo "set ssl:ca-file \"/etc/ssl/certs/cert.pem\"" > .lftp/rc
>    >    > ```
>    >    {: .code-in}
>    > Connect to the server with for example the admin account:
>    > `lftp -p 55000 admin@example.org@$HOSTNAME`
>    > Enter the password of the admin@example.org galaxy user.
>    > Put a random file:
>    > `put /srv/galaxy/server/CITATION`
>    > Check it is there with `ls`.
>    {: .tip}
>
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Use it in galaxy
>
> 1. Open your galaxy in a browser.
>
> 2. Log in with the user you used to upload the file.
>
> 3. Click on the upload button.
>
> 4. Click on the Choose FTP files button
> You should see your file.
>
> 5. Click on it and click on Start to launch the upload.
> It should go to your history as a new dataset.
>
> 6. Click again on Choose FTP files button.
> Your file has disappeared. By default, the files are removed from the ftp at import.
>
>    > ### {% icon tip %} You want to change this behaviour?
>    > You just need to add `ftp_upload_purge: false` to the galaxy_config/galaxy variables (next to `ftp_upload_dir`).
>    {: .tip}
{: .hands_on}
