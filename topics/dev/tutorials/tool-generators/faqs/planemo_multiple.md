---
title: Only one Planemo test runs at a time. Why doesn't the server allow more than one at once?
box_type: question
layout: faq
contributors: [fubar2]
---

- When a new dependency is being installed in the Planemo Conda repository, there is no locking to prevent a second process from overwriting or otherwise interfering with it's own independent repository update.
- The result is not pretty.
- Allowing two tests to run at once has proven to be unstable so the Appliance is currently limited to one.


