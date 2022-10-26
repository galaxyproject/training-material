---
title: 'Error: "skipping: no hosts matched"'
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

There can be multiple reasons this happens, so we'll step through all of them.
We'll start by assuming you're running the command

```
ansible-playbook galaxy.yml
```

The following things can cause issues:

1. Within your `galaxy.yml`, you've referred to a host group that doesn't exist or is misspelled. Check the `hosts: galaxyservers` to ensure it matches the host group defined in the `hosts` file.
2. Vice-versa, the group in your `hosts` file should match the hosts selected in the playbook, `galaxy.yml`.
3. If neither of these are the issue, it's possible Ansible doesn't know to check the `hosts` file for the inventory. Make sure you've specified `inventory = hosts` in your `ansible.cfg`.
