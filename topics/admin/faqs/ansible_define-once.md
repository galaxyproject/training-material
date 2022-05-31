---
title: Define once, reference many times
area: ansible
box_type: tip
layout: faq
contributors: [natefoo, shiltemann]
---

Using practices like those shown above helps to avoid problems caused when paths are defined differently in multiple places. The datatypes config file will be copied to the same path as Galaxy is configured to find it in, because that path is only defined in one place. Everything else is a reference to the original definition! If you ever need to update that definition, everything else will be updated accordingly.
