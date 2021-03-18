---
layout: tutorial_hands_on

title: UseGalaxy.eu's Training Infrastructure as a Service
time_estimation: 10m
questions:
  - Is this service appropriate for my event?
objectives:
  - Identify if it is appropriate
  - Interact with the UseGalaxy.eu admins to arrange for infrastructure
key_points:
  - Infrastructure is available for running Galaxy trainings for free from UseGalaxy.eu
  - This can be easier than setting up a local Galaxy and may have more resources available
  - But have a backup plan
contributors:
  - hexylena
---

# Introduction
{:.no_toc}

![TIaaS Logo](../../images/tiaas-logo.png){: width="500px"}

UseGalaxy.eu has developed Training Infrastructure as a Service (TIaaS for short) which allows you to use UseGalaxy.eu with a private queue for your training event. Your trainees' jobs won't wait in the main queue, and can be processed much more quickly than they might be otherwise. This can provide the experience of a local, private Galaxy combined with a public Galaxy that you are not responsible for maintaining. Additionally if something goes wrong, you can conveniently blame the UseGalaxy.eu admins, rather than feeling the stress of debugging and fixing your private Galaxy.

> ### Agenda
>
> In this tutorial, we will see:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Identify if TIaaS is Appropriate For Your Training

First consider the requirements for your training to see if TIaaS from UseGalaxy.eu is a good fit for you:

- Do you need really special tools that are not already available on UseGalaxy.eu?
  - The UseGalaxy.eu server supports most training workflows and has one of the largest toolboxes
  - But we only install publicly available and open tools: [usegalaxy.eu tool repository](https://github.com/usegalaxy-eu/usegalaxy-eu-tools)
  - If you need something special you might want to [contact the admins](mailto:contact@usegalaxy.eu), and they can potentially accomodate your needs.
- Do you need extra guarantees that the server will be online?
  - The UseGalaxy.eu server has good uptime for a server run by an academic group, but it cannot make promises regarding availability.
  - Additionally this server can experience occasional slow downs due to usage by other groups and users. If you need more guarantees, please find an alternative.

# How TIaaS Works

We have several groups of virtual machines (VMs) attached to UseGalaxy.eu that run user jobs. For trainings we attach a new group of VMs that is specially labelled for that training. When normal users run tools on our server, these jobs are instructed to avoid the training pools by default.

When your users join a training, using a special URL provided to you, they then are placed in a special training group. Their jobs will then preferentially run on a training machine, and, in the event there is no more capacity, they will run on the main queue. If a spot on a training VM opens up first, they will run there rather than continuing to wait in the main queue.

# The Application Process

> ### {% icon hands_on %} Hands-on: Apply for TIaaS Training
>
> 1. [Fill out the form](https://usegalaxy.eu/request-tiaas)
>
> 2. The UseGalaxy.eu admin team will review the request and get in touch with you as needed to identify the compute resources you need for your training
>
> 3. The UseGalaxy.eu admin team will inform you of the URL you should provide to your participants during the training. They can open this URL and they will be added to the special training group
>
{: .hands_on}

# The Student's Process

We have a "test" TIaaS training setup which will never have compute resources associated with it, but you can use it to test the process of signing up, if you wish:

> ### {% icon hands_on %} Hands-on: Join a TIaaS Training
>
> 1. Click on this link: [https://usegalaxy.eu/join-training/test](https://usegalaxy.eu/join-training/test)
>
> 2. In the background you are added to a group in Galaxy. If you were to run jobs they would be tagged with this special queue.
>
{: .hands_on}

That's it, you're now in the "test" TIaaS group. It's really that easy for students.

# The Status Dashboard

Once your students are registered, and you're running your training, a common question instructors ask are "Are you all done?" and students are often not as vocal or repsonsive as we would like. So the TIaaS service now has a dashboard you can view which shows you the queue status for everyone in your training group. This works by finding all of the members of that training group, and checking all jobs that were created in the last 3 hours. For all of those jobs, this is displayed as a simple dashboard with the status of these jobs:

It shows:

- Overview of queue (how many are in state new/queued/ok/error)
- Overview split by tools (how many people are done running Fastqc?)
- A full listing of the queue

![TIaaS Queue Status](../../images/tiaas-status.png "The TIaaS Status dashboard gives you an overview of all jobs states (are they ok or not), as well as a breakdown by tool. This is useful for finding out if everyone is finished running FastQC this morning and if they mostly worked OK. Finally it gives you a detailed breakdown, shown in the order they were submitted. This can give you a more detailed feeling for how the students are progressing through the tutorial.")


> ### {% icon hands_on %} Hands-on: View the dashboard
>
> 1. The status dashboard is just the same URL as to join the group, with `/status` at the end: [https://usegalaxy.eu/join-training/test/status](https://usegalaxy.eu/join-training/test/status)
>
{: .hands_on}

# Before Your Training

1. Run through your planned trainings on the UseGalaxy.eu server to ensure everything is available including data libraries and tools.
2. If anything is missing, contact the admins via [Gitter](https://gitter.im/usegalaxy-eu/Lobby) or [Email](mailto:contact@usegalaxy.eu) and they will work to resolve the missing training material resources.
3. Follow the rest of the training material's [guide to preparing a workshop]({% link topics/instructors/tutorials/organize-workshop/tutorial.md %})

# During Your Training

1. Check out your status dashboard
2. Watch for problems and contact the admins via [Gitter](https://gitter.im/usegalaxy-eu/Lobby)
   - Jobs spending abnormally long in the queue
   - Jobs failing unexpectedly
2. Remind your participants to submit bug reports if they experience any tool errors
3. Follow the rest of the training material's [guide to running a workshop]({% link topics/instructors/tutorials/running-workshop/tutorial.md %})

# Conclusion
{:.no_toc}
