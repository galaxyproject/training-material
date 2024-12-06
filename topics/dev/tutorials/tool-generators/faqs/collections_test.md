---
title: I want to use a collection for outputs but it always passes the test even when the script fails. Why?
box_type: question
layout: faq
contributors: [fubar2]
---

- Collections are tricky for generating tests.
  - The contents appear only after the tool has been run and even then may vary with settings.
- A manual test override is currently the only way to test collections properly.
- Automation is hard. If you can help, pull requests are welcomed.
- Until it's automated, please take a look at the plotter sample.
- It is recommended that you modify the test over-ride that appears in that sample form. Substitute one or more of the file names you expect to see after the collection is filled by your new tool for the `<element.../>` used in the plotter sample's tool test.

