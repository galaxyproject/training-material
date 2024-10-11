---
title: Annotating Pre-requisites
area: gtn
box_type: tip
layout: faq
contributors: [hexylena]
---

If you are adding a tutorial, annotating the pre-requisites is an important task! It will help ensure learners know what they need to know before starting the tutorial. They also let instructors plan a schedule optimally.

Internal requirements often include specific features of Galaxy you plan to use in your training material, and let learners know which tutorials to follow first, before starting your tutorial. 

```yaml
requirements:
  - type: "internal"
    topic_name: galaxy-interface
    tutorials:
      - collections
      - upload-rules
```

Or you can have external requirements, which link to another site.

```yaml
requirements:
  -
    type: "external"
    title: "Trackster"
    link: "https://wiki.galaxyproject.org/Learn/Visualization"
```

Least commonly needed are software requirements. These are usually used in e.g. Galaxy Admin Training tutorials, but if you have specific software requirements, you can list them here:

```yaml
requirements:
- type: none
  title: "Web browser"
- type: none
  title: "A linux-based machine or linux emulator"
- type: none
```
