---
title: Define once, reference many times
area: ansible
box_type: tip
layout: faq
contributors: [natefoo, shiltemann]
---

Using variables, either by defining them ahead of time, or simply accessing them via existing data structures that have been defined, e.g.:

```yaml
# defining a variable that gets reused is great!
galaxy_user: galaxy

galaxy_config:
  galaxy:
    # Re-using the galaxy_config_dir variable saves time and ensures everything
    # is in sync!
    datatypes_config_file: "{{ galaxy_config_dir  }}/datatypes_conf.xml"

# and now we can re-use "{{ galaxy_config.galaxy.datatypes_config_file }}"
# in other places!


galaxy_config_templates:
  - src: templates/galaxy/config/datatypes_conf.xml
    dest: "{{ galaxy_config.galaxy.datatypes_config_file }}"
```

Practices like those shown above help to avoid problems caused when paths are defined differently in multiple places. The datatypes config file will be copied to the same path as Galaxy is configured to find it in, because that path is only defined in one place. Everything else is a reference to the original definition! If you ever need to update that definition, everything else will be updated accordingly.
