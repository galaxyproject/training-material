---
title: Tool missing from Galaxy
area: deployment
box_type: tip
layout: faq
contributors: [natefoo]
---

First, restart Galaxy and watch the log for lines like:

```console
Loaded tool id: toolshed.g2.bx.psu.edu/repos/iuc/sickle/sickle/1.33, version: 1.33 into tool panel....
```

After startup, check `integrated_tool_panel.xml` for a line like the following to be sure it was loaded properly and added to the toolbox (if not, check the logs further)

```xml
<tool id="toolshed.g2.bx.psu.edu/repos/iuc/sickle/sickle/1.33" />
```

If it is a toolshed tool, check `shed_tool_conf.xml` for

```xml
<tool file="toolshed.g2.bx.psu.edu/repos/iuc/sickle/43e081d32f90/sickle/sickle.xml" guid="toolshed.g2.bx.psu.edu/repos/iuc/sickle/sickle/1.33">
...
</tool>
```

Additionally if you have multiple job handlers, *sometimes*, rarely they don't all get the update. Just restart them if that's the case. Alternatively you can send an (authenticated) API requested:

```
curl -X PUT https://galaxy.example.org/api/configuration
```
