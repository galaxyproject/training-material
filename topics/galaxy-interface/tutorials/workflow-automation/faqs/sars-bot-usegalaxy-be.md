---
title: What should I do special if on usegalaxy.be?
box_type: tip
layout: faq
contributors: [bedroesb]
---

**Note for anyone trying to follow the tutorial on usegalaxy.be:**

In step 3 of the hands-on section of setting up the sars-cov-2 analysis bot, when suggested to run

```
planemo run vcf2lineage.ga vcf2lineage-job.yml --profile planemo-tutorial --history_name "vcf2lineage test"
```

please use directly the workflow ID `814dd8d1c056bc54` instead of `vcf2lineage.ga`. This ID points to a public workflow that's using the version of the pangolin tool installed on usegalaxy.be`.

