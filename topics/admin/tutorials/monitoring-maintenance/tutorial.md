---
layout: tutorial_hands_on
topic_name: admin
tutorial_name: monitoring-maintenance
---

# Monitoring and maintenance

## Runing the Reports Application

### Section 1 - Configure reports

Begin by making a copy of the reports config file to your config directory, and editing it:

```console
$ sudo -u galaxy cp /srv/galaxy/server/config/reports.ini.sample /srv/galaxy/config/reports.ini
$ sudo -u galaxy -e /srv/galaxy/config/reports.ini
```

Since we serve Galaxy at the root of our webserver, we'll need to serve Reports from a subdirectory: `/reports`. This is the default if we enable the `proxy-prefix` filter, all we need to do is uncomment the `proxy-prefix` setting. We also need to point the reports application at Galaxy's PostgreSQL database:

```ini
filter-with = proxy-prefix
cookie_path = /reports
database_connection = postgresql:///galaxy?host=/var/run/postgresql
file_path = /srv/galaxy/data
```

### Section 2 - Configure nginx

We now have to configure nginx to serve the reports app at `/reports`. This is done in `/etc/nginx/sites-available/galaxy`. Add the following upstream:

```nginx
upstream reports {
    server localhost:9001;
}
```

And in addition, add this new section to the existing `server { ... }` block:

```nginx
    location /reports {
        proxy_pass           http://reports;
        proxy_set_header     X-Forwarded-Host $host;
        proxy_set_header     X-Forwarded-For  $proxy_add_x_forwarded_for;
    }
```

Then, restart nginx with:

```console
$ sudo systemctl restart nginx
```

### Section 3 - Start reports

We need a way to start and stop the reports application. This can be done with supervisor. Add the following to `/etc/supervisor/conf.d/galaxy.conf`:

```ini
[program:reports]
command         = /srv/galaxy/venv/bin/python ./scripts/paster.py serve /srv/galaxy/config/reports.ini --log-file=/srv/galaxy/log/reports.log
directory       = /srv/galaxy/server
autostart       = true
autorestart     = true
startsecs       = 10
user            = galaxy
```
