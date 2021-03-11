---
title:  Getting your admin API key
area: user preferences
box_type: tip
layout: faq
---


Galaxy admin accounts are specified as a comma-separated email list in the `admin_users` directive of `galaxy.yml` . If you have set up your Galaxy server using the [Galaxy Installation with Ansible]({% link topics/admin/tutorials/ansible-galaxy/tutorial.md %}) tutorial, this is set to `admin@example.org`.

1. In your browser, open your Galaxy homepage
2. Log in using the admin email, or register a new account with it if it is the first time you use it
3. Go to `User -> Preferences` in the top menu bar, then click on `Manage API key`
4. If there is no current API key available, click on `Create a new key` to generate it
5. Copy your API key to somewhere convenient, you will need it throughout this tutorial
