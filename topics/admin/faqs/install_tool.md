---
title: Install tools via the Admin UI
area: galaxy admin interface
box_type: tip
layout: faq
contributors: [shiltemann]
---


1. Open Galaxy in your browser and type `{{ include.query }}` in the tool search box on the left. If "{{ include.name }}" is among the search results, you can skip the following steps.
2. Access the Admin menu from the top bar (you need to be logged-in with an email specified in the `admin_users` setting)
3. Click "Install and Uninstall", which can be found on the left, under "Tool Management"
4. Enter `{{ include.query }}` in the search interface
5. Click on the first hit, having `devteam` as owner
6. Click the "Install" button for the latest revision
7. Enter "{{ include.section }}" as the target section and click "OK".
