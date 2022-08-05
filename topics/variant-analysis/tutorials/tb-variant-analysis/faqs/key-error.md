---
title: "TB Variant Report crashes (with an error about KeyError: 'protein')"
box_type: question
layout: faq
contributors: [pvanheus, hexylena]
---

This is a bug present in TB Variant Report (aka tbvcfreport) version 0.1.8 and earlier. In this case it is triggered by the presence of variants in Rv3798. You only see this bug, however, if you forget to run `tb_variant_filter` (TB Variant Filter). Rv3798 is a suspected transposase and any variants in this gene region would be filtered out by `tb_variant_filter`, so if you see this crash, make sure you have run the filter step before the TB Variant Report step.
