---
title:  Getting your API key
area: user preferences
box_type: tip
layout: faq
contributors: [hexylena,shiltemann,bebatut]
---

{% if include.dev %}
A quick way to get an API key is to create a new user on our local Galaxy server. Navigate to the Galaxy base directory and execute the [run.sh](https://github.com/galaxyproject/galaxy/blob/dev/run.sh) script. This starts the Galaxy server.
{% elsif include.admin %}
Galaxy admin accounts are specified as a comma-separated email list in the `admin_users` directive of `galaxy.yml` . If you have set up your Galaxy server using the [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, this is set to `admin@example.org`.
{% endif %}

1. In your browser, open your Galaxy homepage{% if dev %}, [http://localhost:8080/](http://localhost:8080/){% endif %}
2. Log in, or register a new account, if it's the first time you're logging in
3. Go to `User -> Preferences` in the top menu bar, then click on `Manage API key`
4. If there is no current API key available, click on `Create a new key` to generate it
5. Copy your API key to somewhere convenient, you will need it throughout this tutorial
