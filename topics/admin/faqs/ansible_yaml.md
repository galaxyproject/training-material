---
title: Is YAML sensitive to True/true/False/false
area: ansible
box_type: tip
layout: faq
contributors: [hexylena]
---

By [this reference](https://yaml.org/refcard.html), YAML doesn't really care:
```
{ Y, true, Yes, ON   }    : Boolean true
{ n, FALSE, No, off  }    : Boolean false
```
