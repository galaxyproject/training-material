---
title: Updating from 22.01 to 23.0 with Ansible
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

Galaxy introduced a number of changes in 22.05 and 23.0 that are extremely important to be aware of during the upgrade process. Namely a new database migration system, and a new required running environment (gunicorn instead of uwsgi).

The scripts to migrate to the new database migration system are only compatible with release 22.05, and then were subsequently removed, so it is **mandatory** to upgrade to 22.05 if you want to go further.

Here is the recommended update procedure with ansible:

1. Update to 22.01 normally
2. Change the release to 22.05, and run the upgrade
   1. Galaxy will probably not start correctly here, ignore it.
   2. Run the database migration manually

      ```
      GALAXY_CONFIG_FILE=/srv/galaxy/config/galaxy.yml sh manage_db.sh -c /srv/galaxy/config/galaxy.yml upgrade
      ```

3. Update your system's ansible, you probably need something with a major version greater than 2.
3. Set the release to `23.0` and make other required changes. There are a lot of useful changes, but the easiest procedure is probably something like:

   1. git clone https://github.com/hexylena/git-gat/
   2. git checkout step-4
   3. Diff and sync (e.g. `vimdiff group_vars/galaxyservers.yml git-gat/group_vars/galaxyservers.yml`) for the main configuration files:

      - group_vars/all.yml
      - group_vars/dbservers.yml
      - galaxy.yml
      - requirements.yml
      - hosts
      - templates/nginx/galaxy.j2

   But the main change is the swap from uwsgi to gravity+gunicorn

   ```diff
   -  uwsgi:
   -    socket: 127.0.0.1:8080
   -    buffer-size: 16384
   -    processes: 1
   -    threads: 4
   -    offload-threads: 2
   -    static-map:
   -      - /static={{ galaxy_server_dir }}/static
   -      - /favicon.ico={{ galaxy_server_dir }}/static/favicon.ico
   -    static-safe: client/galaxy/images
   -    master: true
   -    virtualenv: "{{ galaxy_venv_dir }}"
   -    pythonpath: "{{ galaxy_server_dir }}/lib"
   -    module: galaxy.webapps.galaxy.buildapp:uwsgi_app()
   -    thunder-lock: true
   -    die-on-term: true
   -    hook-master-start:
   -      - unix_signal:2 gracefully_kill_them_all
   -      - unix_signal:15 gracefully_kill_them_all
   -    py-call-osafterfork: true
   -    enable-threads: true
   -    mule:
   -      - lib/galaxy/main.py
   -      - lib/galaxy/main.py
   -    farm: job-handlers:1,2
   +  gravity:
   +    process_manager: systemd
   +    galaxy_root: "{{ galaxy_root }}/server"
   +    galaxy_user: "{{ galaxy_user_name }}"
   +    virtualenv: "{{ galaxy_venv_dir }}"
   +    gunicorn:
   +      # listening options
   +      bind: "unix:{{ galaxy_mutable_config_dir }}/gunicorn.sock"
   +      # performance options
   +      workers: 2
   +      # Other options that will be passed to gunicorn
   +      # This permits setting of 'secure' headers like REMOTE_USER (and friends)
   +      # https://docs.gunicorn.org/en/stable/settings.html#forwarded-allow-ips
   +      extra_args: '--forwarded-allow-ips="*"'
   +      # This lets Gunicorn start Galaxy completely before forking which is faster.
   +      # https://docs.gunicorn.org/en/stable/settings.html#preload-app
   +      preload: true
   +    celery:
   +      concurrency: 2
   +      loglevel: DEBUG
   +    handlers:
   +      handler:
   +        processes: 2
   +        pools:
   +          - job-handlers
   +          - workflow-schedulers
   ```

   Some other important changes include:
   - uchida.miniconda is replaced with galaxyproject.conda
   - usegalaxy_eu.systemd is no longer needed
   - `galaxy_user_name` is defined in all.yml in the latest git-gat
   - git-gat also separates out the DB serving into a `dbservers.yml` host group

4. Backup your `venv`, `mv /srv/galaxy/venv/ /srv/galaxy/venv-old/`, as your NodeJS is probably out of date and Galaxy doesn't handle that gracefully
5. Do any local customs for luck (knocking on wood, etc.)
6. Run the playbook
7. Things might go wrong with systemd units
   - try running `galaxyctl -c /srv/galaxy/config/galaxy.yml update` as root
   - you may also need to `rm /etc/systemd/system/galaxy.service` which is then no longer needed
   - you'll have a `galaxy.target` and you can instead `systemctl daemon-reload` and `systemctl start galaxy.target`
