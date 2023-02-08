---
layout: tutorial_hands_on

title: "How I learned to stop worrying and love the systemd"
questions:
- Unix is supposed to be about FILES™
- What is this systemd stuff?
- Why should I love it?
- I have so many worries!
objectives:
- Have an objective understanding of systemd allowing the user to obtain the benefits of this new system
- Realise the joys of journald, and how it makes logging easier and simpler
time_estimation: "30m"
key_points:
- systemd units are actually kinda nice
- In most versions you even can limit memory, cpu usage of process groups. In recent versions, disk usage!
- No more PID files! Cgroups mean the processes are actually cleaned up!
- journalctl makes accessing logs incredibly easy; view multiple logs simultaneously, properly interlaced.
- journalctl-vacuum replaces onerous cleaning processes and even can clean by relative timestamps

contributions:
  authorship:
  - hexylena
  editing:
  - natefoo
  testing:
  - natefoo
tags:
  - ansible
  - systemd
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible-galaxy
---

Many linux sysadmins with years and years of experience bemoan systemd ("it's infecting everything! Now it wants to mess with time? And DNS???") and journalctl ("unix was supposed to be about files!") and while those are fair complaints and make systemd and friends wildly more opaque than traditional SysV and logging to files, there are some benefits that can be obtained, and may be interesting even to the wise old admins. There is a lot of convenience in systemd that can make the tradeoffs worth it.

This tutorial assumes a working knowledge of unix systems (files, directories, services, logging, `/var/log`, etc.).

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# systemd

## Why systemd units, and not service shell scripts?

- CGroups
- Hard Memory Limits that include all children (also CPU, disk.)
- Dependencies (is your database up, before your service that needs it?)
- Fancy restarting (always/n times, delays, timeouts.)
- [Monitor it with Telegraf](https://github.com/influxdata/telegraf/tree/master/plugins/inputs/systemd_units)

## Unit Files

```
[Unit]
Description=Galaxy
After=network.target
After=time-sync.target

[Service]
UMask=022
Type=simple
User=galaxy
Group=galaxy
WorkingDirectory=/srv/galaxy/server
TimeoutStartSec=60
ExecStart=/srv/galaxy/venv/bin/galaxyctl start --foreground --quiet

Environment=HOME=/srv/galaxy
Environment=VIRTUAL_ENV=/srv/galaxy/venv

MemoryAccounting=yes
MemoryLimit=16G
CPUAccounting=yes
IOAccounting=yes

[Install]
WantedBy=multi-user.target
```

Here's a simple unit file. It has a name, and some targets that must be up before it should be started, then the service definition featuring all of the common fields like:

- User/group and umask restriction
- Working directory
- Process to start (and optionally reload/stop commands.)
- Some environment variables
- Enabling of CPU and Mmeory accounting
- As well as limiting the memory to 16GB, at which time the child will be killed and restarted.

Unit files are very easy to create, be sure to see the [examples section](https://man.archlinux.org/man/systemd.service.5#EXAMPLES) of the documentation, and stick with `Type=simple` services which stay in the foreground for maximum convenience.


## Status

Status is generally your first step looking at any unit file.

```bash
$ systemctl galaxy
● galaxy.service - Galaxy
     Loaded: loaded (/etc/systemd/system/galaxy.service; enabled; vendor preset: enabled)
     Active: active (running) since Tue 2022-06-28 15:17:04 CEST; 5 days ago
   Main PID: 217446
      Tasks: 73 (limit: 9521)
     Memory: 1.7G (limit: 16.0G)
        CPU: 18h 16min 20.685s
     CGroup: /system.slice/galaxy.service
             ├─217446 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/galaxyctl start --foreground --quiet
             ├─217467 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/supervisord -c /srv/galaxy/var/gravity/supervisor/supervisord.conf --nodaemon
             ├─217484 /bin/tail -f /srv/galaxy/var/gravity/supervisor/supervisord.log
             ├─217487 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/gunicorn galaxy.webapps.galaxy.fast_factory:factory() --timeout 300 --pythonpath lib -k galaxy.webapps.galaxy.workers.Worker -b unix:/srv/galaxy/var/config/gunicorn.sock --workers=2 --config python:galaxy.web_stack.gunicorn_config --preload --forwarded-allow-ips=*
             ├─245320 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/gunicorn galaxy.webapps.galaxy.fast_factory:factory() --timeout 300 --pythonpath lib -k galaxy.webapps.galaxy.workers.Worker -b unix:/srv/galaxy/var/config/gunicorn.sock --workers=2 --config python:galaxy.web_stack.gunicorn_config --preload --forwarded-allow-ips=*
             ├─245322 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/gunicorn galaxy.webapps.galaxy.fast_factory:factory() --timeout 300 --pythonpath lib -k galaxy.webapps.galaxy.workers.Worker -b unix:/srv/galaxy/var/config/gunicorn.sock --workers=2 --config python:galaxy.web_stack.gunicorn_config --preload --forwarded-allow-ips=*
             ├─245334 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/celery --app galaxy.celery worker --concurrency 2 --loglevel DEBUG --pool threads --queues celery,galaxy.internal,galaxy.external
             ├─245401 /srv/galaxy/venv/bin/python -c from multiprocessing.resource_tracker import main;main(7)
             ├─245582 /srv/galaxy/venv/bin/python /srv/galaxy/venv/bin/celery --app galaxy.celery beat --loglevel DEBUG --schedule /srv/galaxy/var/gravity/celery-beat-schedule
             ├─245806 python ./lib/galaxy/main.py -c /srv/galaxy/config/galaxy.yml --server-name=handler_0 --attach-to-pool=job-handler --attach-to-pool=workflow-scheduler --pid-file=/srv/galaxy/var/gravity/supervisor/handler_0.pid
             ├─246274 python ./lib/galaxy/main.py -c /srv/galaxy/config/galaxy.yml --server-name=handler_1 --attach-to-pool=job-handler --attach-to-pool=workflow-scheduler --pid-file=/srv/galaxy/var/gravity/supervisor/handler_1.pid
             └─246737 python ./lib/galaxy/main.py -c /srv/galaxy/config/galaxy.yml --server-name=handler_2 --attach-to-pool=job-handler --attach-to-pool=workflow-scheduler --pid-file=/srv/galaxy/var/gravity/supervisor/handler_2.pid
```

Helpful things to note:

- the file defining the service is listed after `Loaded`, if you want to see what it does or how it's configured, or to edit that configuration
- You can see how long it's been running
- How many tasks (children) there are
- Memory and the limit, if memory limiting is configured.
- CPU time used
- And then many of the children.

## Restarting

Here's one of the nice new things you can do, once you've moved away from SysV, and to the joyous world of systemd, restart multiple processes at once! If you have multiple Galaxy processes, use a wildcard to restart them all.

```bash
systemctl restart nginx galaxy-*
```

## Failed Units

Easily check if any service has failed with:

```bash
systemctl --failed
```

## Editing Units or Overriding

Sometimes one of the system units provided will have some weird behaviour that you need to override (or your ansible role doesn't expose it), then you can use `systemctl edit unit` to override some settings. 

Additional directives can be supplied, e.g. making your service start after another service is started, if you need to sequence their starts.

Note that if it was defined in the base unit file, then you need to set it to an empty value first, before overriding. E.g. `ExecStart=\nExecStart=command...`


## Masking Units

I never needed this until one day I did. Masking a unit makes it impossible to start, but if someone's gone and installed a GUI on your server, inclusive of power management systems that trigger hibernation when not in use, then this might be a useful trick. You can mask anything, service, or even targets like I had to in that case.

This command will prevent the device from reaching or activating any of those targets. Useful for servers when you're in a bind and don't know how to remove power management:

```console
systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target
```

## Unit Security Optimisation

Because systemd uses cgroups, it can also give us a nice overview of any security issues that might be worth looking into. Here we see the Galaxy unit has a lot of 

```console
ubuntu@gat-1:~$ systemd-analyze security galaxy
  NAME                                       DESCRIPTION                                         EXPOSURE
✗ PrivateNetwork=                            Service has access to the host's network            0.5
✓ User=/DynamicUser=                         Service runs under a static non-root user identity  
....
✗ CapabilityBoundingSet=~CAP_SYS_CHROOT      Service may issue chroot()                          0.1
✗ ProtectHostname=                           Service may change system host/domainname           0.1
✗ CapabilityBoundingSet=~CAP_BLOCK_SUSPEND   Service may establish wake locks                    0.1
✗ CapabilityBoundingSet=~CAP_LEASE           Service may create file leases                      0.1
✗ CapabilityBoundingSet=~CAP_SYS_PACCT       Service may use acct()                              0.1
✗ CapabilityBoundingSet=~CAP_SYS_TTY_CONFIG  Service may issue vhangup()                         0.1
✗ CapabilityBoundingSet=~CAP_WAKE_ALARM      Service may program timers that wake up the system  0.1
✗ RestrictAddressFamilies=~AF_UNIX           Service may allocate local sockets                  0.1

→ Overall exposure level for galaxy.service: 9.2 UNSAFE 😨
```

Here we see a lot of low hanging fruit for restricting the abilities of Galaxy's processes. This can also be applied to e.g. the job runner, or your webserver for additional security.

## Testing with systemd-run

Would you like to test some of these limits manually? Or just to run a one-off process in a very sandboxed environment? Enter `systemd-run`, which allows just that. Here you can set a lot of the capabilities or restrictions:

```bash
systemd-run -p IOAccounting=yes -p CPUAccounting=true -p MemoryAccounting=true -p TasksAccounting=true -p ProtectSystem=strict -p PrivateDevices=yes -t -S
```

And `-S` puts you in an interactive terminal, allowing you to play around with this environment. You get a cleaned out environment, working inside of a cgroup (visible in `systemd-cgtop`, and `systemctl status`) and when you exit it tells how you many resources you've used. This is both an easy way to test processes, and an option to sandbox rarely run processes, by setting an extremely restrictive environment.

## Optimising Slow Boots

systemd features an analyze command to check why your booting is slow

```console
$ systemd-analyze
Startup finished in 3.841s (firmware) + 3.925s (loader) + 11.104s (kernel) + 1min 5.331s (userspace) = 1min 24.203s
graphical.target reached after 11.547s in userspace
```

You can even have it generate a fancy plot for you!

```bash
systemd-analyze plot > plot.svg
```

![waterfall graphic of systemd boot time, cloud-init has the most, largest red bars](images/plot.svg "Notice how many units wait for network to finish being ready, and are stuck until that is ready. This image was generated from a training VM. Can you spot Galaxy? It doesn't show the red bars, because it does not currently notify systemd when it's done starting up.")

## Pro/Cons of systemd

> > <code-out-title>Pros</code-out-title>
> > - CGroups, to reap children!
> > - CGroups, for hard memory/cpu limits
> > - CGroups, for security benefits, and dropping capabilities!
> > - Service dependencies and simultaneously booting multiple units.
> > - Easily override tasks
> {: .code-out}
>
> > <code-in-title>Cons</code-in-title>
> > - It's new and requires time to learn
> >     - But all the distros have mostly switched so, might be time to learn.
> > - It is significantly more complicated
> > - Even if that complexity does have a lot of benefits
> {: .code-in}
{: .code-2col}

## Further Reading

- [Archi wiki page](https://wiki.archlinux.org/title/Systemd)
- [man 1 systemctl](https://man.archlinux.org/man/systemctl.1)
- [man 5 systemd.service](https://man.archlinux.org/man/systemd.service.5)

# timers

Do you use cron aggressively, like any other sysadmin? It's an amazing part of the system, no? Doing things periodically, and repeatably. But you've probably had to learn a lot of hard lessons along the way, right?

- Doing things idempotently
- Checking that you don't have a period cron job running, before starting this one, if it's a long one.
- What if there was a reboot?

Would you like to worry less about those issues? Try systemd timers today!

[The arch page](https://wiki.archlinux.org/title/Systemd/Timers) is an excellent reference on why they're interesting and useful. We'll reproduce their pro/con list below:

> > <code-out-title>Pros</code-out-title>
> > - Easier testing, you can trigger units any time
> > - Cgroups!
> > - dependencies (e.g. network!)
> > - automatic logging!!
> > - two flavours of timers: clock time based, and monotonic (e.g. time since last execution)
> {: .code-out}
>
> > <code-in-title>Cons</code-in-title>
> > - Two files, instead of a single line. It *is* annoying, you're not wrong. Maybe write an ansible role that translates it.
> > - No `MALITO` functionality built in.
> {: .code-in}
{: .code-2col}

## Listing Timers

We have a lovely timer listing, which tells you when it'll next run, how much time is left, when it last ran, how much time since that passed.

```
$ systemctl list-timers
NEXT                         LEFT          LAST                         PASSED       UNIT                         ACTIVATES
Mon 2022-07-04 12:32:45 CEST 6min left     Mon 2022-07-04 11:34:59 CEST 51min ago    anacron.timer                anacron.service
Mon 2022-07-04 14:07:08 CEST 1h 40min left Fri 2022-07-01 13:38:20 CEST 2 days ago   motd-news.timer              motd-news.service
Mon 2022-07-04 17:12:05 CEST 4h 45min left Mon 2022-07-04 10:20:49 CEST 2h 5min ago  ua-timer.timer               ua-timer.service
Mon 2022-07-04 17:14:33 CEST 4h 48min left Fri 2022-07-01 19:22:14 CEST 2 days ago   apt-daily.timer              apt-daily.service
Tue 2022-07-05 00:00:00 CEST 11h left      Mon 2022-07-04 09:32:16 CEST 2h 54min ago logrotate.timer              logrotate.service
Tue 2022-07-05 00:00:00 CEST 11h left      Mon 2022-07-04 09:32:16 CEST 2h 54min ago man-db.timer                 man-db.service
Tue 2022-07-05 01:44:11 CEST 13h left      Mon 2022-07-04 11:07:49 CEST 1h 18min ago fwupd-refresh.timer          fwupd-refresh.service
Tue 2022-07-05 06:28:44 CEST 18h left      Mon 2022-07-04 09:33:59 CEST 2h 52min ago apt-daily-upgrade.timer      apt-daily-upgrade.service
Tue 2022-07-05 10:12:36 CEST 21h left      Mon 2022-07-04 09:47:01 CEST 2h 39min ago systemd-tmpfiles-clean.timer systemd-tmpfiles-clean.service
Sun 2022-07-10 03:10:47 CEST 5 days left   Mon 2022-07-04 09:32:52 CEST 2h 53min ago e2scrub_all.timer            e2scrub_all.service
Mon 2022-07-11 00:00:00 CEST 6 days left   Mon 2022-07-04 09:32:16 CEST 2h 54min ago fstrim.timer                 fstrim.service
```

In the last column we can see the timer units, and the services which they activate. In systemd everything is a unit, and hopefully by this point in the tutorial you've seen some of the benefits of that.

Let's look at the `motd-news` timer:

```console
$ systemctl cat motd-news.timer
[Unit]
Description=Message of the Day

[Timer]
OnCalendar=00,12:00:00
RandomizedDelaySec=12h
Persistent=true
OnStartupSec=1min

[Install]
WantedBy=timers.target
```

It has a description, an `OnCalendar` time since it's a clock based timer, potentially a random delay as well (despite being `Sec`, can be expressed in any unit). The associated `motd-news.service` does not need to be mentioned, it will be started by default based on the name of the timer.

What about logs? They're in journalctl now! Let's take advantage of some of the things we learned earlier, like, wildcards!

```console
journalctl -u motd-news.* --since '5 days ago'
```

You'll see the timer unit start mostly around reboots, usually it doesn't have much useful in it. The service unit has what we're more interested in, any logged output.

## Fixing those Cron Issues

Before we mentioned a few issues that admins dealt with with cron

- Doing things idempotently
- Checking that you don't have a period cron job running, before starting this one, if it's a long one.
- What if there was a reboot?

While we can't help with the first, the second is taken care of by having systemd services be state machines, if there's an already running copy of the `.service`, then the `.timer` won't do anything. No more `pgrep` to check for specific processes. And as for reboots, `Persistent=true` can cause a timer to trigger immediately, if it missed the last start time.

## Notifying you

As mentioned in the Arch wiki, you can setup a "notification" service. This is a service, typically a simple bash script, that, when run, sends a notification to a channel of your choice (Slack, XMPP, Email, etc.) This can be a really simple service:

```
/etc/systemd/system/failure-notification@.service

[Unit]
Description=Send a notification
After=network.target

[Service]
Type=simple
ExecStart=/path/to/failure-notification.sh %i
```

And then for any unit (timer's associated service, normal units, etc) where you'd like to be informed of the failure, simply update or override that to do the notification:

```
[Unit]
OnFailure=failure-notification@%n
```

Here it will be templated out with `%n` meaning the unit name, which will replace the `%i` in the notification unit, and be included in the notification to you.


## Further Reading 

- [The arch page](https://wiki.archlinux.org/title/Systemd/Timers)
- [man 5 systemd.timer](https://man.archlinux.org/man/systemd.timer.5)
- [man 7 systemd.time](https://man.archlinux.org/man/systemd.time.7)

# journald

> Where did my files go?
{:.quote}

In a number of systems they're still there, to be honest, you've still got all of your old log files.

But now instead of search N different log files, or memorising which files have what contents, there's one firehose for you to watch instead: journalctl. So we need to make this manageable with filters!

## Per Unit Logs

> <code-out-title>The journalctl Way</code-out-title>
> ```bash
> journalctl -u galaxy
> ```
{: .code-out}

returns all of the logs for Galaxy.  The first line even tells us how far back our logs go:

```console
$ journalctl -u galaxy | head -n1
-- Logs begin at Tue 2022-06-28 08:31:06 UTC, end at Mon 2022-07-04 09:30:01 UTC. --
```

> > <code-out-title>Pros</code-out-title>
> > You did not need to know or remember:
> > 
> > - Where the log files were (`/var/log`? `/srv/galaxy/log`? somewhere else?)
> > - If any of the old logs were compressed (or use `zless`/`zcat`/`zgrep`)
> {: .code-out}
>
> > <code-in-title>Cons</code-in-title>
> > - Can't just `cat` it.
> {: .code-in}
{: .code-2col}

## Logs for a Group of Units

Sometimes you want to tail the logs of multiple proceses as once. Have you ever been debugging an issue, and wanted to see the logs of `nginx`, **all** `galaxy` processes, and `postgres`? Well you can't because `nginx` needs special configuration to write to journalctl (see, your files are still there, and it makes things extra confusing.)

Adding these two lines to your `nginx` configuration will place logs in journald

```nginx
error_log syslog:server=unix:/dev/log;
access_log syslog:server=unix:/dev/log;
```

And then you can simply

> <code-out-title>The journalctl Way</code-out-title>
> ```bash
> journalctl -f -u galaxy* -u nginx -u postgresql
> ```
{: .code-out}

If you have multiple similarly named units, the wildcard feature is incredibly helpful. UseGalaxy.eu had multiple named processes for: web handlers, workflow schedulers, job handlers. We could see logs simultaneously across all of these units.

> > <code-out-title>Pros</code-out-title>
> > - Tail multiple units at once
> > - Again don't have to know where the files are
> > - Or switch tools if you want to look at older logs
> {: .code-out}
>
> > <code-in-title>Cons</code-in-title>
> > ? 
> {: .code-in}
{: .code-2col}

## Dropping Logs

journald has the habit of throwing away logs if it gets too many messages simultaneously, to prevent {DOS} and denial-of-storage attacks. To disable that behaviour, you can tune

```
RateLimitBurst=0
```

(see `man 5 journald.conf`)

## Logs during a specific time period

This is my absolute favourite feature of journalctl, and I've regularly used it to narrow down logs from dozens of processes, especially when I don't know the root cause of an issue and need to watch all of the logs.

> > <code-in-title>The Old Way</code-in-title>
> > ```bash
> > zcat /var/log/* | egrep "(Jul\s*1 11:4.:..|2022-07-01 11:4.:..)"
> > ```
> > Hope everything logs in the same timezone! And no one uses UTC when there's a different system time configured.
> {: .code-in}
>
> > <code-out-title>The journalctl Way</code-out-title>
> > ```bash
> > journalctl --since "2022-07-01 11:40 CEST" --until "2022-07-01 11:50 CEST"
> > ```
> {: .code-out}
{: .code-2col}

This is very useful if you know when an error is, easily pull out entire system logs from around a specific timepoint. Everything is forced through journald's logging format, so you can be sure the timestamps are all correct, and `--since` and `--when` can be made timezone aware if you have reason to use that (e.g. a user reporting an issue in their current local time.) You can mix and match explicit time zone locales for even more fun:

```bash
journalctl --since "2022-06-21 14:24 Pacific/Auckland" --until "2022-06-21 14:30 Europe/Amsterdam"
```

The logs will be displayed in the last mentioned TZ.

## Logs during a specific time period, relatively

This is another flag that I use literally every day:

```bash
journalctl --since '1 hour ago'
```

If you know relatively when an issue started happening, it's a great way to filter on the logs that are interesting, without having to do clock math.

## Previous Boot Logs

Want to see what happened during a previous run of the system, during a previous boot? Use `--list-boots` to list all of the previous system restarts

```console
$ journalctl --list-boots
-14 143fdb7ac7a54370b41196c02affc3e8 Wed 2022-04-06 20:03:38 CEST—Mon 2022-04-11 08:26:00 CEST
-13 45801d31b0654fe3a4a7d5541e6ff71a Mon 2022-04-11 09:06:07 CEST—Wed 2022-04-20 13:37:31 CEST
-12 661e959635954aa0bbd6b800433e5097 Wed 2022-04-20 13:37:53 CEST—Thu 2022-04-28 17:26:32 CEST
-11 cedaef23e95849c39cfcdb1611b616a2 Thu 2022-04-28 17:40:33 CEST—Thu 2022-04-28 17:43:23 CEST
-10 6ac084d5256d4134bd6e395b490e99f2 Thu 2022-04-28 18:44:57 CEST—Mon 2022-05-09 11:21:33 CEST
 -9 680c9440c75d4699b6de48a56891d835 Mon 2022-05-09 11:21:54 CEST—Tue 2022-05-17 15:21:36 CEST
 -8 b8581fd6394144a9abc59b2832056037 Wed 2022-05-18 09:08:54 CEST—Fri 2022-05-20 17:33:36 CEST
 -7 8bb5c9d7d4ef4117ad8181d0314de7e1 Sun 2022-05-22 17:45:23 CEST—Thu 2022-06-02 16:35:45 CEST
 -6 5325407857054a828dbba6bb104efc80 Fri 2022-06-03 12:04:28 CEST—Thu 2022-06-09 10:15:40 CEST
 -5 0428cc6d11bc4aa9bd8e3999f9bd2c7f Thu 2022-06-09 10:20:05 CEST—Fri 2022-06-17 12:26:06 CEST
 -4 990910c70e79458394f5d9a97ae724ab Sat 2022-06-18 11:21:16 CEST—Tue 2022-06-21 12:14:47 CEST
 -3 d9a9088afacd4c6c9f85db49e9b321ce Tue 2022-06-21 14:29:20 CEST—Tue 2022-06-21 15:07:26 CEST
 -2 c3f0ec4fbc93410280f1e086d6533286 Tue 2022-06-21 15:07:45 CEST—Thu 2022-06-23 12:58:11 CEST
 -1 5d4d9c039daf4e5795a13211de4d6fd4 Thu 2022-06-23 12:58:31 CEST—Fri 2022-07-01 20:10:15 CEST
  0 490a1609fb5d422cad1e7135db88efe7 Mon 2022-07-04 09:32:12 CEST—Mon 2022-07-04 11:56:09 CEST
```

You can then see logs for those specific timeperiods with 

```
$ journalctl -b -3 | head
-- Logs begin at Wed 2022-04-06 20:03:38 CEST, end at Mon 2022-07-04 11:59:13 CEST. --
jun 21 14:29:20 cosima kernel: microcode: microcode updated early to revision 0x88, date = 2021-03-31
jun 21 14:29:20 cosima kernel: Linux version 5.13.0-44-generic (buildd@lcy02-amd64-107) (gcc (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0, GNU ld (GNU Binutils for Ubuntu) 2.34) #49~20.04.1-Ubuntu SMP Wed May 18 18:44:28 UTC 2022 (Ubuntu 5.13.0-44.49~20.04.1-generic 5.13.19)
```

All of the other filtering commands (`--since`, `--when`) natively search across boots so there's generally no need to supply this flag unless you want to restrict to a specific time period based on restarting.


## Advanced Log Filtering

There's a lot of additional information journald tracks with your logs, you can read them in the `verbose` format to see

```console
$ journalctl -o verbose
Wed 2022-04-06 20:03:38.281794 CEST [s=1b57c151fb2e49ad98a84ac45ed2e245;i=a310796;b=143fdb7ac7a54370b41196c02affc3e8;m=51917bf226;t=5dc002e4a3f42;x=4ae28fb4218f3691]
    _TRANSPORT=stdout
    _STREAM_ID=45b0adde513f46348564e6f9ce85e4b2
    PRIORITY=6
    SYSLOG_FACILITY=3
    SYSLOG_IDENTIFIER=influxd-systemd-start.sh
    MESSAGE=[httpd] 127.0.0.1 - - [06/Apr/2022:20:03:38 +0200] "POST /write?db=telegraf HTTP/1.1 " 204 0 "-" "Telegraf/1.21.4 Go/1.17.7" e265a273-b5d3-11ec-88d5-80fa5b94d68a>
    _PID=2491
    _UID=1001
    _GID=1001
    _COMM=influxd
    _EXE=/usr/bin/influxd
    _CMDLINE=/usr/bin/influxd -config /etc/influxdb/influxdb.conf
    _CAP_EFFECTIVE=0
    _SELINUX_CONTEXT=unconfined
    _SYSTEMD_CGROUP=/system.slice/influxdb.service
    _SYSTEMD_UNIT=influxdb.service
    _SYSTEMD_SLICE=system.slice
    _SYSTEMD_INVOCATION_ID=0ba7986947854acc938a93f020f96cee
    _BOOT_ID=143fdb7ac7a54370b41196c02affc3e8
    _MACHINE_ID=08854338f5f046d5a0d73f81ca5c2d79
    _HOSTNAME=cosima
```

All of these attributes can be used to filter your logs, just add them to your `journalctl` command like so:

```bash
journalctl _GID=1001
```

Or, for example, to get every log, by any process that logged, running under a specific user, like the Galaxy user:

```bash
journalctl _UID=$(getent passwd | grep galaxy | cut -f 3 -d : )
```

Or to see every log from a python program?

```bash
journalctl _COMM=python3
```

## Cleaning Logs

Journald comes with a nice tool for cleaning up logs, `journalctl --vacuum-...`. There are a couple different ways to invoke this. It's very useful if you left Galaxy poorly configured, reading bad XML files and crashing every minute, logging endlessly until you checked it again. That sort of behaviour generates an incredible number of logs.

First we can check how much disk space our logs are using:

```
$ journalctl --disk-usage
Archived and active journals take up 120.0M in the file system.
``` 

And then we can clean it!

Based on **size**:

```bash
journalctl --vacuum-size=4GB
```

Based on **age**

```bash
journalctl --vacuum-time="2 days"
```

Or based on number of **files**

```bash
journalctl --vacuum-files=3 # I've never used this one.
```

Note that this only applies to "archived" logs, not "active" ones, but you can use `journalctl --rotate` to force them to archived, before cleaning..


## Log Verification

Want to check if they've been tampered with (and someone forgot to recalculate checksums?). The `--verify` flag can be used to do that.

```
$ journalctl --verify
PASS: /var/log/journal/08854338f5f046d5a0d73f81ca5c2d79/user-1001@0005df9b9d17faa5-c4276c13359e16e1.journal~
File corruption detected at /var/log/journal/08854338f5f046d5a0d73f81ca5c2d79/system@0005df9b9c8ebfda-d41419e76ede7f2c.journal~:3f070f8 (of 67108864 bytes, 98%).
FAIL: /var/log/journal/08854338f5f046d5a0d73f81ca5c2d79/system@0005df9b9c8ebfda-d41419e76ede7f2c.journal~ (Bad message)
PASS: /var/log/journal/08854338f5f046d5a0d73f81ca5c2d79/system@ede9cabbb3874bd98bcb0b16ada49d8f-000000000a6ce684-0005de680b29a077.journal
7e5abf8: Unused data (entry_offset==0)
..  94%
```

## Further Reading

- [man 5 journald.conf](https://manpages.ubuntu.com/manpages/bionic/man5/journald.conf.5.html)
- [man 1 journalctl](https://manpages.ubuntu.com/manpages/bionic/man1/journalctl.1.html)
- [man 7 systemd.journal-fields](https://manpages.ubuntu.com/manpages/bionic/man7/systemd.journal-fields.7.html)
