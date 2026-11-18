---
redirect_from: [/faqs/galaxy/galaxy_creating_an_account, /faqs/galaxy/creating_account]
title: How do I create an account on a public Galaxy instance?
area: account
layout: faq
box_type: tip
contributions:
  authorship: [jennaj, bernandez, samanthaanjei, wm75]
---

1. To create an account at any public Galaxy instance, choose your server from the available [list of Galaxy Platforms](https://galaxyproject.org/use/).

    There are several UseGalaxy servers:
    {% assign servers = nil | list_usegalaxy_servers_shuffle %}
    {% for server in servers %}
    - [{{ server.name }}]({{ server.url }}) (**{{ server.id | upcase }}**)
    {% endfor %}

2. Click on **Register** in the top bar on the server.

   ![Register on the top panel]({% link faqs/galaxy/images/login_register.png %})

3. Fill in the registration form, then click on **Create**.

   Your account should now get created, but will remain inactive until you verify the email address you provided in the registration form.

4. Check for a **Confirmation Email** in the email you used for account creation.

    Missing? Check your **Trash** and **Spam** folders.

5. Click on the **Email confirmation link** to fully activate your account.

   {% icon galaxy-info %} Delivery of the confirmation email is blocked by your email provider or you mistyped the email address in the registration form?

   Please do **not** register again. Click **Resend Verification** in the banner to get a new confirmation email, or follow the [instructions to change the email address]({% link faqs/galaxy/user_preferences_manage_information.md %}) registered with your account! The confirmation email will be resent to your new address once you have changed it.

   Trouble logging in later? Account email addresses and public names are caSe-sensiTive. Check your activation email for formats.
