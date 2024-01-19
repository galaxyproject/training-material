---
title: Using data source tools with Pulsar
area: deployment
box_type: tip
layout: faq
contributors: [martinwolst]
---

Data source tools such as UCSC Main will fail if Pulsar is the default destination.

To prevent the error from occureing you can force individual tools to run on a different destination or handler by adding to the job_conf file like:


```xml
<tools>
    <tool id="longbar" destination="my-local" />
</tools>
```

```yml
tools:
- id: longbar
  handler: my-local
```