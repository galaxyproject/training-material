---
title: Using Git With Ansible Vaults
area: gat
box_type: hands_on
layout: faq
contributors: [hexylena]
---

When looking at `git log` to see what you changed, you cannot easily look into
Ansible Vault changes: you just see the changes in the encrypted versions which
is unpleasant to read.

Instead we can use [`.gitattributes`](https://www.git-scm.com/docs/gitattributes) to tell `git` that we want to use a
different program to visualise differences between two versions of a file,
namely `ansible-vault`.

1. Check your `git log -p` and see how the Vault changes look (you can type `/vault` to search). Notice that they're just changed encoded content.
1. Create the file `.gitattributes` in the same folder as your `galaxy.yml` playbook, with the following contents:

   ```
   group_vars/secret.yml diff=ansible-vault merge=binary
   ```

1. Try again to `git log -p` and look for the vault changes. Note that you can now see the decrypted content! Very useful.
