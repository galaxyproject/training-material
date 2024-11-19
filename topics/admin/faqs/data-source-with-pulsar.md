---
title: Using data source tools with Pulsar
area: deployment
box_type: tip
layout: faq
contributors: [martinwolst]
---

Data source tools such as UCSC Main will fail if Pulsar is the default destination.

To fix this issue you can force individual tools to run on a specific destination or handler by adding to your job_conf file:

For job_conf.xml
```xml
<tools>
    <tool id="ucsc_table_direct1" destination="my-local" />
</tools>
```
For job_conf.yml
```yml
tools:
- id: ucsc_table_direct1
  handler: my-local
```