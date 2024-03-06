---
title: Debugging tool memory errors
area: deployment
box_type: tip
layout: faq
contributors: [natefoo]
---

Often the tool output contains one of:

```console
MemoryError                 # Python
what():  std::bad_alloc     # C++
Segmentation Fault          # C - but could be other problems too
Killed                      # Linux OOM Killer
```

Solutions:
- Change input sizes or params
  - Map/reduce?
- Decrease the amount of memory the tool needs
- Increase the amount of memory available to the job
  - Request more memory from cluster scheduler
  - Use job resubmission to automatically rerun with a larger memory allocation
- Cross your fingers and rerun the job
