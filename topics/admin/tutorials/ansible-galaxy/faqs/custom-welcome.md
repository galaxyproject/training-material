---
title: Customising the welcome page
area: ansible-galaxy
box_type: tip
layout: faq
contributors: [hexylena]
---

Customising the `welcome.html` page is very easy!

1. Create a custom `welcome.html` in `templates/galaxy/config/` with whatever content you want, e.g.

   ```html
   <html>
     <body>
       <h1>Welcome to Galaxy</h1>
       <iframe src="https://galaxyproject.org/#" width="100%" height="500px" />
     </body>
   </html>
   ```

2. Add it to your deployed configuration

   {% raw %}
   ```bash
   --- a/group_vars/galaxyservers.yml
   +++ b/group_vars/galaxyservers.yml
   @@ -83,12 +83,7 @@ certbot_agree_tos: --agree-tos
    galaxy_config_templates:
    - src: templates/galaxy/config/job_conf.yml.j2
      dest: "{{ galaxy_config.galaxy.job_config_file }}"
    - src: templates/galaxy/config/reports.yml
      dest: "{{ galaxy_reports_path }}"
   +- src: templates/galaxy/config/welcome.html
   +  dest: "{{ galaxy_config_dir }}"
   ```
   {% endraw %}

3. Update the `templates/nginx/galaxy.j2` to point to the new location

   {% raw %}
   ```bash
   --- a/templates/nginx/galaxy.j2
   +++ b/templates/nginx/galaxy.j2
   @@ -83,12 +83,7 @@
    location config/welcome.html {
   -    alias {{ galaxy_server_dir }}/static/welcome.html.sample;
   +    alias {{ galaxy_config_dir }}/welcome.html;
        expires 24h;
    }
   ```
   {% endraw %}

4. Run the playbook
5. Done! You might need to refresh / clear your cache before the changes are visible.
