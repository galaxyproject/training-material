---
layout: tutorial_hands_on
redirect_from:
  - /topics/instructors/tutorials/setup-tiaas-for-training/tutorial

title: Training Infrastructure as a Service
subtopic: prepare
time_estimation: 10m
questions:
  - Is this service appropriate for my event?
objectives:
  - Identify if it is appropriate
  - Interact with the service administrators to arrange for infrastructure
key_points:
  - Infrastructure is available for running Galaxy trainings for free from UseGalaxy.eu, UseGalaxy.org, and UseGalaxy.org.au
  - This can be easier than setting up a local Galaxy and may have more resources available
  - But have a backup plan
contributors:
  - hexylena
tags:
  - cyoa
---

UseGalaxy.eu has developed Training Infrastructure as a Service (TIaaS for short) which allows you to use Galaxy with a private queue for your training event. Your trainees' jobs won't wait in the main queue, and can be processed much more quickly than they might be otherwise. This can provide the experience of a local, private Galaxy combined with a public Galaxy that you are not responsible for maintaining. Additionally if something goes wrong, you can conveniently blame the Galaxy admins, rather than feeling the stress of debugging and fixing your private Galaxy.

![TIaaS Logo](../../images/tiaas-logo.png){: width="500px"}

UseGalaxy.eu has developed Training Infrastructure as a Service (TIaaS for short) which allows you to use Galaxy with a private queue for your training event. Your trainees' jobs won't wait in the main queue, and can be processed much more quickly than they might be otherwise. This can provide the experience of a local, private Galaxy combined with a public Galaxy that you are not responsible for maintaining. Additionally if something goes wrong, you can conveniently blame the Galaxy admins, rather than feeling the stress of debugging and fixing your private Galaxy.

> <agenda-title></agenda-title>
>
> In this tutorial, we will see:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Where are you

{% include _includes/cyoa-choices.html option1="Europe, the Middle East, and Africa (EMEA)" option2="Asia/Pacific (APAC)" option3="the Americas (AMER)" text="Galaxy has servers across the world, but it's best to choose one close to you. Here we customise this tutorial based on where you are, so, pick the location closest to where you'll be running your training." default="Europe, the Middle East, and Africa (EMEA)" %}

# Identify if TIaaS is Appropriate For Your Training

First consider the requirements for your training to see if TIaaS is a good fit for you:

- Do you need really special tools that are not already available on <span class="Europe-the-Middle-East-and-Africa-EMEA"><a href="https://usegalaxy.eu">UseGalaxy.eu</a></span><span class="AsiaPacific-APAC"><a href="https://usegalaxy.org.au">UseGalaxy.org.au</a></span><span class="the-Americas-AMER"><a href="https://usegalaxy.org">UseGalaxy.org</a></span>?
  - The UseGalaxy.* servers support most training workflows.
  - But we only install publicly available and open tools: <span class="Europe-the-Middle-East-and-Africa-EMEA"><a href="https://github.com/usegalaxy-eu/usegalaxy-eu-tools">EU Tool Repository</a></span><span class="AsiaPacific-APAC"><a href="https://github.com/usegalaxy-au/usegalaxy-au-tools">Australian Tool Repository</a></span><span class="the-Americas-AMER"><a href="https://github.com/galaxyproject/usegalaxy-tools/">American Tool Repository</a></span>
- Do you need extra guarantees that the server will be online?
  - The big UseGalaxy.* servers have good uptime, for services run by academic groups, but it cannot make promises regarding availability.
  - Additionally this server can experience occasional slow downs due to usage by other groups and users. If you need more guarantees, please find an alternative.

# How TIaaS Works

<p class="Europe-the-Middle-East-and-Africa-EMEA">
We have several groups of virtual machines (VMs) attached to UseGalaxy.eu that run user jobs. For trainings we attach a new group of VMs that is specially labelled for that training. When normal users run tools on our server, these jobs are instructed to avoid the training pools by default.
</p>
<p class="AsiaPacific-APAC">
We have a pool of virtual machines (VMs) dedicated to training, that run training specific jobs only.
</p>
<p class="the-Americas-AMER">
Jobs submitted while in a TIaaS group have an artificially capped memory limited, and walltime. By lowering these parameters we can be sure the jobs will execute more quickly, but this only works for training data. If you run a training with "real data" it may require too much memory or computation time, and those jobs will be killed.
</p>

When your users join a training, using a special URL provided to you, they then are placed in a special training group. Their jobs will then preferentially run on a training machine, and, in the event there is no more capacity, they will run on the main queue. If a spot on a training VM opens up first, they will run there rather than continuing to wait in the main queue.

# The Application Process

> <hands-on-title>Apply for TIaaS Training</hands-on-title>
>
> 1. <span class="Europe-the-Middle-East-and-Africa-EMEA"><a href="https://usegalaxy.eu/tiaas/new/">Fill out the form</a></span><span class="AsiaPacific-APAC"><a href="https://usegalaxy.org.au/tiaas/new/">Fill out the form</a></span><span class="the-Americas-AMER"><a href="https://usegalaxy.org/tiaas/new/">Fill out the form</a></span>
>
> 2. The admin team will review the request and get in touch with you as needed to identify the compute resources you need for your training
>
> 3. The admin team will inform you of the URL you should provide to your participants during the training. They can open this URL and they will be added to the special training group
>
{: .hands_on}

# The Student's Process

We have a "test" TIaaS training setup on UseGalaxy.eu, which will never have compute resources associated with it, but you can use it to test the process of signing up, if you wish:

> <hands-on-title>Join a TIaaS Training</hands-on-title>
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


> <hands-on-title>View the dashboard</hands-on-title>
>
> 1. The status dashboard is just the same URL as to join the group, with `/status` at the end: [https://usegalaxy.eu/join-training/test/status](https://usegalaxy.eu/join-training/test/status)
>
{: .hands_on}

# Before Your Training

**At least 2 weeks** before your training:

1. Run through your planned trainings on the UseGalaxy.eu server to ensure everything is available including data libraries and tools.
2. If anything is missing, contact the admins and they will work to resolve the missing training material resources.
3. Follow the rest of the training material's [guide to preparing a workshop]({% link topics/teaching/tutorials/organize-workshop/tutorial.md %})

# During Your Training

1. Check out your status dashboard
2. Watch for problems and contact the admins if
   - Jobs are spending abnormally long in the queue
   - Jobs are failing unexpectedly
2. Remind your participants to submit bug reports if they experience any tool errors
3. Follow the rest of the training material's [guide to running a workshop]({% link topics/teaching/tutorials/running-workshop/tutorial.md %})
