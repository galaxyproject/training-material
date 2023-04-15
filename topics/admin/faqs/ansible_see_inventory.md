---
title: How do I see what variables are set for a host?
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

If you are using a simple `group_vars` file only, per group, and no other variable sources, then it's relatively easy to tell what variables are getting set for your host! Just look at that one file.

But if you have graduated into using a more complex setup, perhaps with multiple sets of variables, like for example:

```
├── group_vars
│   ├── all
│   │   ├── all.yml
│   │   └── secret.yml
│   ├── galaxyservers.yml
│   └── pulsarservers.yml
├── hosts
├── host_vars
│   ├── galaxy.example.org
│   │   ├── all.yml
│   │   └── secret.yml
│   ├── pulsar.example.org
│   │   ├── all.yml
│   │   ├── pulsar.yml
│   │   └── secret.yml
...
```

Then it might be harder to figure out what variables are being set, in full. This is where [`ansible-inventory`](https://docs.ansible.com/ansible/latest/cli/ansible-inventory.html) command can be useful.

**Graph** shows you the structure of your host groups:

```
$ ansible-inventory --graph
@all:
  |--@cluster:
  |  |--allie.example.com
  |  |--bob.example.com
  |  |--charlie.example.com
[...]
```

Here is a relatively simple, flat example, but this can be more complicated if you nest sub-groups of hosts:

```
@all:
  |--@local:
  |  |--localhost
  |--@ungrouped:
  |--@workshop_instances:
  |  |--@workshop_eu:
  |  |  |--gat-0.eu.training.galaxyproject.eu
  |  |  |--gat-1.eu.training.galaxyproject.eu
  |  |--@workshop_oz:
  |  |--@workshop_us:
```

**List** shows you all defined variables:

{% raw %}
```
$ ansible-inventory --host galaxy.example.com | head
[WARNING]: While constructing a mapping from
/group_vars/galaxyservers.yml, line 3, column
1, found a duplicate dict key (tiaas_templates_dir). Using last defined value
only.
{
    "ansible_connection": "local",
    "ansible_user": "ubuntu",
    "certbot_agree_tos": "--agree-tos",
    "certbot_auth_method": "--webroot",
    "certbot_auto_renew": true,
    "certbot_auto_renew_hour": "{{ 23 |random(seed=inventory_hostname)  }}",
    "certbot_auto_renew_minute": "{{ 59 |random(seed=inventory_hostname)  }}",
{% endraw %}
```

And, helpfully, if variables are overridden in precedence you can see that as well with the above warnings.
