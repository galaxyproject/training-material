---
title: First job I submitted remains grey or running for a long time - is it broken?
box_type: question
layout: faq
contributors: [fubar2]
---

- Check with `top` or your system monitor - if Conda is running, things are working but it's slow the first time a dependency is installed.
- The first run generally takes a while to install all the needed dependencies.
- Subsequent runs should start immediately with all dependencies already in place.
- Installing new Conda dependencies just takes time so tools that have new Conda packages will take longer to run the first time if they must be installed.
- In general, a `planemo_test` job usually takes around a minute - planemo has to build and tear down a new Galaxy for generating test results and then again for testing properly. Longer if the tool has Conda dependencies.
- The very first test in a fresh appliance may take 6 minutes so be patient.



