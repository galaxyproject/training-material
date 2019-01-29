---
layout: tutorial_hands_on

title: "Upstream Authentication"
questions:
- How can I connect Galaxy with CAS, SAML, etc.
objectives:
- be familiar with configuring Galaxy to use an upstream (proxy) authentication provider
- be able to log in to your Galaxy server with a file-configured user.
time_estimation: "2h"
key_points:
- Remote auth is not complex to set up and can help you meet institutional requirements
contributors:
  - natefoo
  - nsoranzo
  - erasche
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

For this exercise we will use a basic password file method for authenticating - this is probably not a very useful method in production, but it demonstrates how the proxy server can be configured to provide the correct header to Galaxy, and how Galaxy integrates with upstream authentication providers. This same method can be used with NGINX and Apache modules for CAS or SAML authentication.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Configuring Authentication

> ### {% icon hands_on %} Hands-on: Configuring everything
>
> 1. Edit your `galaxyservers` group variables file and update the main location block defined for serving galaxy. Add the parameters:
>      - `auth_basic galaxy;`
>      - `auth_basic_user_file /etc/nginx/passwd;`
>      - `uwsgi_param HTTP_REMOTE_USER $reomte_user;`
>
>    It should look like:
>
>    ```nginx
>        location / {
>            uwsgi_pass           127.0.0.1:8080;
>            uwsgi_param          UWSGI_SCHEME $scheme;
>            include              uwsgi_params;
>            auth_basic           galaxy;
>            auth_basic_user_file /etc/nginx/passwd;
>            uwsgi_param          HTTP_REMOTE_USER $remote_user;
>        }
>    ```
>
>    `auth_basic` enables validation of username and password using the "HTTP Basic Authentication" protocol. Its value `galaxy` is used as a realm name to be displayed to the user when prompting for credentials.
>
>    `auth_basic_user_file` specifies the file that keeps usernames and passwords, in the following format:
>
>    ```
>    # comment
>    name1:password1
>    name2:password2:comment
>    name3:password3
>    ```
>
>    `uwsgi_param` adds `HTTP_REMOTE_USER` to the special variables passed by nginx to uwsgi, with value `$remote_user`, which is a nginx embedded variable containing the username supplied with the Basic authentication.
>
> 2. Add a pre_task using the [`htpasswd`](https://docs.ansible.com/ansible/2.4/htpasswd_module.html) which sets up a password file in `/etc/nginx/passwd`, with owner and group set to root, and a name and password, and a mode of 0640.
>    > ### {% icon question %} Question
>    >
>    > How does your final configuration look?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > ```
>    > > - htpasswd:
>    > >     path: /etc/nginx/passwd
>    > >     name: helena
>    > >     password: 'squeamish ossifrage'
>    > >     owner: root
>    > >     group: root
>    > >     mode: 0640
>    > > ```
>    > {: .solution }
>    >
>    {: .question}
>
> 3. Galaxy needs to be instructed to expect authentication to come from the upstream proxy. In order to do this, set the following two options in your Galaxy group variables:
>
>    ```yaml
>    ...
>    galaxy_config:
>      galaxy:
>        ...
>        use_remote_user: true
>        remote_user_maildomain: "{{ hostname }}"
>    ```
>
>    Set the `remote_user_maildomain` option to the appropriate domain name for your site.
>
> 4. Run the playbook
>
{: .hands_on}

> ### {% icon tip %} Tip: Access denied
>
> If you see this message, it is because nginx is not correctly sending the `REMOTE_USER` variable
>
> ![access denied message](../../images/access_denied.png)
>
{: .tip}


# Testing

You should now be presented with a password dialog when attempting to load the Galaxy UI.

> ### {% icon hands_on %} Hands-on:
>
> 1. Log in using the username and password you provided when creating the `passwd` file. If your username and the value of `remote_user_maildomain` match an existing user, you will be logged in to that account. If not, a new account will be created for that user.
>
{: .hands_on}

Note that some user features are not available when remote user support is enabled.

Try logging out by selecting **User** -> **Logout**. You will discover that when returning to the user interface, you are still logged in. This is because Galaxy has no way of logging you out of the proxy's authentication system. Instead, you should set `remote_user_logout_href` in `galaxy.ini` to point to the URL of your authentication system's logout page.

# (Workshop Only) Reverting

We don't want to leave Galaxy this way for the rest of our workshop.

> ### {% icon hands_on %} Hands-on: Reverting the changes
>
> 1. Edit your group variables file and comment out:
>
>    - the NGINX changes
>    - `use_remote_user: true`
>
> 2. Run the playbook
>
{: .hands_on}
