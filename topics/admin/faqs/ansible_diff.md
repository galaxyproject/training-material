---
title: What if you forget `--diff`?
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

If you forget to use `--diff`, it is not easy to see what has changed. Some modules like the `copy` and `template` modules have a `backup` option. If you set this option, then it will keep a backup copy next to the destination file.

However, most modules do not have such an option, so if you want to know what changes, always use `--diff`.
