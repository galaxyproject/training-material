---
layout: tutorial_hands_on

title: Provide your custom JupyterLab via Galaxy
subtopic: practises
draft: true
time_estimation: 15m
questions:
  - As a Galaxy admin, how can I provide my users with a custom Notebook Server?
objectives:
  - Learn how to add additional interactive tools for each Jupyter customization.
key_points:
  - Each user-base has different needs. It's easy to serve them providing customized Jupyter Notebook Servers.
contributors:
  - mittler-works
---

# Introduction

(WIP)

# Make your Image available

E.g. on Docker Hub, Quay, GitLab or any other public accessible Registry

# Create a new interactive tool definition (xml)

--> You may take the available Jupyter Notebook tool definition as starting point: https://github.com/galaxyproject/galaxy/blob/dev/tools/interactive/interactivetool_jupyter_notebook.xml

--> Tweak it to your needs

--> Adjust the container image to point to your container image on the registry.

# Include your new interactive tool definition into your galaxy installation

Official docs: https://docs.galaxyproject.org/en/master/admin/special_topics/interactivetools.html

--> Include your newly created definition in your `config/tool_conf.xml` like described in the docs:

```xml
...
<toolbox monitor="true">
    ...
    <section id="interactivetools" name="Interactive tools">
      ...
      <tool file="PATH_TO_YOUR_NEWLY_CREATED_TOOL_DEFINITION" />
      ...
    </section>
    ...
</toolbox>
```
