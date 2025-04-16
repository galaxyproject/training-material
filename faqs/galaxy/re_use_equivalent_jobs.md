---
title: How do I re-use equivalent jobs in Galaxy (aka Job Cache)?
area: features
layout: faq
box_type: tip
contributors:
  - dadrasarmin
  - teresa-m
---

We can reuse the reproducibility of Galaxy to detect if a tool has been run with the exact same parameters and inputs before. In this case, we can simply skip the computational step and just reuse the
data we have previously computed. We call this feature the `job cache`. Part of the `job cache` is all your personal data and all data in public histories. 
This can be highly helpful, e.g., for training events, if the instructor makes a respective training history public before the event. 
If the trainee activates this option in their account and uses the same input and parameters, they will immediately receive the results. This feature reduces the waiting time in the training sessions,
saves energy and computational resources, and therefore reduces environmental impact.

To activate this feature, take the following steps:

1. To activate this option for your account, click on your username at the top right of the page.
- Select `Preferences` and navigate to your user-references.
- In your middle panel search for `Manage Information` and select them. You can also navigate to "https://<your-galaxy-server.domain>/user" — for example, https://usegalaxy.eu/user.
- Find the grey box: `Do you want to be able to re-use equivalent jobs?` 
- Within the box, change the slider from `no` to `yes`.
- Scroll down to the bottom of the page and click the `Save` button.

2. For every tool you want to run now, you will notice the option `Attempt to re-use jobs with identical parameters?`. To test this:
- Click on any tool you would like to run
- If you scroll down to the end of the  `Tool Parameters` section until you see the `Run tool` button, you will notice the new option `Attempt to re-use jobs with identical parameters?` above the `Run tool` button.
- You can enable this option by sliding the `No` to `Yes`
- Once you click on the `Run tool`, Galaxy will check if this tool was run before with the exact same parameters and inputs. If so, the results will be retrieved from the `job cache` and **not** be calculated.

⚠️ At the moment, this feature only works with data shared/reused inside Galaxy. If you upload the same file twice, we can not detect that it is the same file.
