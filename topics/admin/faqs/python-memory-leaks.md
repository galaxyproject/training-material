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

If you're debugging a Galaxy process, and a user is splitting a dataset into a million element collection that could use some memory, however we have some pretty stringent limits for this now and haven't had a problem in a while.
