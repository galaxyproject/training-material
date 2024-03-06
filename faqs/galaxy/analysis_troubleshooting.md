---
title: Troubleshooting errors
area: analysis
description: When you get a red dataset in your history, it means something went wrong. But how can you find out what it was? And how can you report errors?
box_type: tip
layout: faq
contributors: [shiltemann]
---

When something goes wrong in Galaxy, there are a number of things you can do to find out what it was. Error messages can help you figure out whether it was a problem with one of the settings of the tool, or with the input data, or maybe there is a bug in the tool itself and the problem should be reported. Below are the steps you can follow to troubleshoot your Galaxy errors.

1. **Expand** the red history dataset by clicking on it.
   - Sometimes you can already see an error message here


2. **View the error message** by clicking on the **bug icon** {% icon galaxy-bug %}


3. **Check the logs.** Output (stdout) and error logs (stderr) of the tool are available:
   - Expand the history item
   - Click on the {% icon details %} icon
   - Scroll down to the **Job Information** section to view the 2 logs:
     - Tool Standard Output
     - Tool Standard Error
   - For more information about specific tool errors, please see the [Troubleshooting section]({% link faqs/galaxy/index.md %}#troubleshooting)


4. **Submit a bug report!** If you are still unsure what the problem is.
   - Click on the bug icon {% icon galaxy-bug %}
   - Write down any information you think might help solve the problem
     - See [this FAQ]({% link faqs/galaxy/analysis_troubleshooting_reporting.md %}) on how to write good bug reports
   - Click **{% icon galaxy-bug %} Report** button

5. **Ask for help!**
   - Where?
     - In the [GTN Gitter Channel](https://gitter.im/Galaxy-Training-Network/Lobby)
     - In the [Galaxy Gitter Channel](https://gitter.im/galaxyproject/Lobby)
     - Browse the [Galaxy Help Forum](https://help.galaxyproject.org/) to see if others have encountered the same problem before (or post your question).
   - When asking for help, it is useful to [share a link to your history]({% link faqs/galaxy/histories_sharing.md %})


