---
title: Debugging Memory Leaks
area: ansible
box_type: tip
layout: faq
contributors: [mvdbeek, hexylena]
---

[memray](https://github.com/bloomberg/memray) is a great memory profiler for debugging memory issues.

In the context of Galaxy, this is significantly easier for job handlers. Install it in your `virtualenv` and 

```console
memray run  --trace-python-allocators -o the_dump <your_handler_startup_command_here>
```

Once you've collected enough data,

```console
memray flamegraph --leaks --temporal the_dump -o the_dump.html
```

would then produce a report that shows allocation made but not freed over time.

It might also be useful to just check what the process is doing with py-spy dump.

You can follow web workers in gunicorn with

```console
memray run --follow-fork -o the_dump gunicorn 'galaxy.webapps.galaxy.fast_factory:factory()' --timeout 600 --pythonpath lib -k galaxy.webapps.galaxy.workers.Worker -b localhost:8082 --config python:galaxy.web_stack.gunicorn_config -w 1 --preload
```

the traced app will run on port 8082, you can then for instance in an upstream nginx section direct a portion of the traffic to your profiled app.
