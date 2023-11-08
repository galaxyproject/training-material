---
title:  Running Ansible on your remote machine
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

It is possible to have ansible installed on the remote machine and run it there, not just from your local machine connecting to the remote machine.

Your hosts file will need to use `localhost`,  and whenever you run playbooks with `ansible-playbook -i hosts playbook.yml`, you will need to add `-c local` to your command.

Be **certain** that the playbook that you're writing on the remote machine is stored somewhere safe, like your user home directory, or backed up on your local machine. The cloud can be unreliable and things can disappear at any time.
