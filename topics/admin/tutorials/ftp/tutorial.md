---
layout: tutorial_hands_on

title: "Enable upload via FTP"
zenodo_link: ""
questions:
  - How can I setup FTP to be easy for my users?
  - Can I authenticate ftp users with Galaxy credentials?
objectives:
  - Configure galaxy and install a FTP server.
  - Use an Ansible playbook for this.
time_estimation: "1h"
key_points:
  - FTP is easy to deploy thanks to the role
  - Users can be authenticated with their Galaxy credentials simplifying the user management process significantly
contributors:
  - lldelisle
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
abbreviations:
  FTP: File Transfer Protocol
  NAT: Network Address Translation
tags:
  - data
  - git-gat
---

This tutorial will guide you to setup an {FTP} server so galaxy users can use it to upload large files. Indeed, as written on the [galaxy community hub](https://galaxyproject.org/ftp-upload/), uploading data directly from the browser can be unreliable and cumbersome. FTP will allow users to monitor the upload status as well as resume interrupted transfers.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="ftp" %}

# FTP

{FTP} is a very old and reliable communication protocol that has been around since 1971 {% cite rfc114 %}. It requires a server (here our galaxy server) and a client (user's computer). The FTP server requires to have at least 2 ports accessible from outside (one for the commands and one for the transfer). Usually the port for the command is 21.

FTP supports two different modes: active, and passive. Active mode requires that the user's computer be reachable from the internet, which in the age of {NAT} and firewalls is usually unusable. So passive mode is the most commonly used. In passive mode, a client connects to the FTP server, and requests a channel for sending files. The server responds with an IP and port, from its range of "Passive Ports".

> <comment-title>Requirements for Running This Tutorial</comment-title>
>
> Your VM or wherever you are installing Galaxy needs to have the following ports available:
>
> - 21
> - Some high range of ports not used by another service, e.g. 56k-60k
>
> You need to know which ports are open so you can use them for the transfer (PassivePorts). In this training we assume that 56k to 60k are open.
>
> Which ports precisely is not important, and these numbers can differ between sites.
{: .comment}

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

> <hands-on-title>Setting up ftp upload with Ansible</hands-on-title>
>
> 1. In your playbook directory, add the `galaxyproject.proftpd` role to your `requirements.yml`
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -38,3 +38,5 @@
>       version: 0.12.0
>     - src: galaxyproject.tiaas2
>       version: 2.1.3
>    +- src: galaxyproject.proftpd
>    +  version: 0.3.1
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement"}
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
> 3. As in this training we are using certbot, we will ask for a private key for proftpd. Add the following line to your `group_vars/galaxyserver.yml` file:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -159,9 +159,11 @@ certbot_well_known_root: /srv/nginx/_well-known_root
>     certbot_share_key_users:
>       - nginx
>       - rabbitmq
>    +  - proftpd
>     certbot_post_renewal: |
>         systemctl restart nginx || true
>         systemctl restart rabbitmq-server || true
>    +    systemctl restart proftpd || true
>     certbot_domains:
>      - "{{ inventory_hostname }}"
>     certbot_agree_tos: --agree-tos
>    {% endraw %}
>    ```
>    {: data-commit="Add proftpd in certbot"}
>
>    This will make a copy of the current letsencrypt key available as `/etc/ssl/user/privkey-proftpd.pem`, and automatically restart proftpd every time the key is updated.
>
> 4. We will configure Galaxy to enable ftp file upload. Add the following line to your `group_vars/galaxyserver.yml` file in the galaxy_config/galaxy section:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -100,6 +100,9 @@ galaxy_config:
>         outputs_to_working_directory: true
>         # TUS
>         tus_upload_store: /data/tus
>    +    # FTP
>    +    ftp_upload_dir: /data/uploads
>    +    ftp_upload_site: "{{ inventory_hostname }}"
>       gravity:
>         galaxy_root: "{{ galaxy_root }}/server"
>         app_server: gunicorn
>    {% endraw %}
>    ```
>    {: data-commit="Add ftp vars in galaxy"}
>
> To check the other options for setting up ftp in Galaxy, please check the [Galaxy configuration documentation](https://docs.galaxyproject.org/en/master/admin/galaxy_options.html?highlight=ftp_upload_site#ftp-upload-dir).
>
> 5. Then we will set the different variables for proftpd. Add the following lines to your `group_vars/galaxyserver.yml` file. Please replace the PassivePorts below with the range of ports that are appropriate for your machine!
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -250,6 +250,27 @@ rabbitmq_users:
>         password: "{{ vault_rabbitmq_password_vhost }}"
>         vhost: /pulsar/galaxy_au
>     
>    +# Proftpd:
>    +proftpd_galaxy_auth: yes
>    +galaxy_ftp_upload_dir: "{{ galaxy_config.galaxy.ftp_upload_dir }}"
>    +proftpd_display_connect: |
>    +  {{ inventory_hostname }} FTP server
>    +
>    +  Unauthorized access is prohibited
>    +proftpd_create_ftp_upload_dir: yes
>    +proftpd_options:
>    +  - User: galaxy
>    +  - Group: galaxy
>    +  - Port: 21
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
>    Here is a description of the set variables:
>
>    | Variable                           | Description                                                                                                  |
>    | ----------                         | -------------                                                                                                |
>    | `proftpd_galaxy_auth`              | Attempt to authenticate users against a Galaxy database.                                                     |
>    | `galaxy_ftp_upload_dir`            | Path to the Galaxy FTP upload directory, should match `ftp_upload_dir` in your Galaxy config.                |
>    | `proftpd_display_connect`          | Message to display when users connect to the FTP server. This should be the message, not the path to a file. |
>    | `proftpd_create_ftp_upload_dir`    | Whether to allow the role to create this with owner `galaxy_user`.                                           |
>    | `proftpd_options`                  | Any option for proftpd, we will just set up the user and group of the `galaxy_user`.                         |
>    | `proftpd_sql_db`                   | Database name to connect to for authentication info.                                                         |
>    | `proftpd_sql_user`                 | (default: the value of galaxy_user): Value of the username parameter to SQLConnectInfo.                      |
>    | `proftpd_conf_ssl_certificate`     | Path on the remote host where the SSL certificate file is.                                                   |
>    | `proftpd_conf_ssl_certificate_key` | Path on the remote host where the SSL private key file is.                                                   |
>    | `proftpd_global_options`           | Set arbitrary options in the <Global> context. We set here the PassivePorts range.                           |
>    | `proftpd_use_mod_tls_shmcache`     | By default proftpd uses `mod_tls_shmcache` which is not installed on the server so we just disable it.       |
>    | `proftpd_tls_options`              | Additional options for tls. We will use `NoSessionReuseRequired`                                             |
>
>    > <tip-title>Why NoSessionReuseRequired?</tip-title>
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
>    @@ -32,6 +32,7 @@
>         - usegalaxy_eu.rabbitmq
>         - galaxyproject.tiaas2
>         - galaxyproject.nginx
>    +    - galaxyproject.proftpd
>         - galaxyproject.tusd
>         - galaxyproject.cvmfs
>         - galaxyproject.gxadmin
>    {% endraw %}
>    ```
>    {: data-commit="Add role to playbook"}
>
> 5. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
>
{: .hands_on}

Congratulations, you've set up FTP for Galaxy.

## Check it works

> <hands-on-title>Checking proftpd from the server</hands-on-title>
>
> 1. SSH into your machine
>
> 2. Check the active status of proftpd by `systemctl status proftpd`.
>
> 3. Check the port has been correctly attributed by `sudo lsof -i -P -n`.
>
>    > <question-title></question-title>
>    >
>    > What do you see?
>    >
>    > > <solution-title></solution-title>
>    > > You should see all the ports used by the server. What interests us is the line with proftpd.
>    > > You should see TCP *:21 (LISTEN).
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
> 4. Check the directory `/data/uploads/` has been created and is empty.
>
>    > <code-in-title>Bash</code-in-title>
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

> <hands-on-title>Checking galaxy detected the ftp possibility</hands-on-title>
>
> 1. Open your galaxy in a browser.
>
> 2. Log in with a user (FTP is only possible for authenticated sessions).
>
> 3. Click on the upload button. You should now see on the bottom "Choose FTP files"
>
> 4. Click on the Choose FTP files button. You should see a message "Your FTP directory does not contain any files."
>
{: .hands_on}

It's working!

> <hands-on-title>Upload your first file</hands-on-title>
>
> There are three options for uploading files, you can choose whichever is easiest for you.
>
> 1. **FileZilla**
>
>    1. Follow the [tutorial](https://galaxyproject.org/ftp-upload/) to upload a file.
>    2. You will have a message which ask you to approve the certificate, approve it.
>
> 2. **lftp**
>
>    You can use locally lftp to test the ftp.
>
>    1. Install lftp with `sudo apt-get install lftp`.
>    2. Add the public certificate to the list of known certificates (only for LetsEncrypt Staging Certificates!):
>       > <code-in-title>Bash</code-in-title>
>       > ```
>       > mkdir .lftp
>       > echo "set ssl:ca-file \"/etc/ssl/certs/cert.pem\"" > .lftp/rc
>       > ```
>       {: .code-in}
>
>    3. Connect to the server with for example the admin account:
>       > <code-in-title>Bash</code-in-title>
>       > ```
>       > lftp admin@example.org@$HOSTNAME
>       > ```
>       {: .code-in}
>
>    4. Enter the password of the admin@example.org galaxy user.
>    5. Put a random file:
>
>       `put /srv/galaxy/server/CITATION`
>
>    6. Check it is there with `ls`.
>    7. Leave lftp with `quit`.
>    {: .tip}
>
> 3. **Curl**
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > curl -T {"/srv/galaxy/server/CITATION"} ftp://localhost --user admin@example.org:password --ssl -k
>    > ```
>    > Here `-T` says to upload a file, `--ssl` ensures that the FTP connection is SSL/TLS encrypted, and `-k` ignores any certificate issues as the hostname `localhost` will not match the certificate we have.
>    {: .code-in}
>
{: .hands_on}

> <hands-on-title>Check where the file has been uploaded</hands-on-title>
>
> 1. SSH into your machine
>
> 4. Check the directory `/uploads/`.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > sudo tree /uploads/
>    > ```
>    {: .code-in}
>
>    > <question-title></question-title>
>    >
>    > What do you see?
>    >
>    > > <solution-title></solution-title>
>    > > As I uploaded a file called `CITATION` with the admin@example.org user I see:
>    > > ```
>    > > /uploads/
>    > > └── admin@example.org
>    > >     └── CITATION
>    > > ```
>    > >
>    > {: .solution }
>    >
>    {: .question}
>
{: .hands_on}

> <hands-on-title>Use it in galaxy</hands-on-title>
>
> 1. Open your galaxy in a browser.
>
> 2. Log in with the user you used to upload the file.
>
> 3. Click on the upload button.
>
> 4. Click on the Choose FTP files button. You should see your file.
>
> 5. Click on it and click on Start to launch the upload. It should go to your history as a new dataset.
>
> 6. Click again on Choose FTP files button. Your file has disappeared. By default, the files are removed from the FTP at import.
>
>    > <tip-title>You want to change this behaviour?</tip-title>
>    > You just need to add `ftp_upload_purge: false` to the galaxy_config/galaxy variables (next to `ftp_upload_dir`).
>    {: .tip}
>
{: .hands_on}

Congratulations! Let your users know this is an option, many of them will prefer to start large uploads from an FTP client.

{% snippet topics/admin/faqs/missed-something.md step=14 %}
