---
title: Debugging tool errors
area: deployment
box_type: tip
layout: faq
contributors: [natefoo]
---

Tool stdout/stderr is available in UI under "i" icon on history dataset

1. Set `cleanup_job` to `onsuccess`
2. Cause a job failure
3. Go to job working directory (find in logs or `/data/jobs/<hash>/<job_id>`)
4. Poke around, try running things (`srun --pty bash` considered useful)

Familiarize yourself with the places Galaxy keeps things
