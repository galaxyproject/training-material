---
title: Mapping Jobs to Specific Storage By User
area: ansible
box_type: tip
layout: faq
contributors: [hexylena, bgruening, natefoo]
---

It is possible to map your jobs to use specific storage backends based on user! If you have e.g. specific user groups that need their data stored separately from other users, for whatever political reasons, then in your dynamic destination you can do something like:

```python
job_destination = app.job_config.get_destination(destination_id)
if user == "alice":
    job_destination.params['object_store_id'] = 'foo' # Maybe lookup the ID from a mapping somewhere
```

If you manage to do this in production, please let us know and we can update this FAQ with any information you encounter.
