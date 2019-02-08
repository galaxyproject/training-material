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
contributors:
  - erasche
---

# Introduction
{:.no_toc}

UseGalaxy.eu has developed Training Infrastructure as a Service (TIaaS for short) which allows you to use UseGalaxy.eu with a private queue for your training event. Your trainees' jobs won't wait in the main queue, and can be processed much more quickly than they might be otherwise. This can provide the experience of a local, private Galaxy combined with a public Galaxy that you are not responsible for maintaining. Additionally if something goes wrong, you can conveniently blame the EU admins, rather than feeling the stress of debugging and fixing your private Galaxy.

> ### Agenda
>
> In this tutorial, we will see:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Identify if TIaaS is Appropriate For Your Training

Consider the requirements for your training:

- Do you need really special tools that are not already available on EU? The EU server supports [most training workflows]({{ site.url }}/badges/#galaxy-europe), but if you need something special you might want to [contact the admins](mailto:contact@usegalaxy.eu) to see if they can accomodate your needs.
- Do you have special private training data that should be made public? If so the EU server is not appropriate
- Do you need extra guarantees that the server will be online? The EU server has very good uptime for a server run by an academic group, but it cannot make promises regarding availability.

# How TIaaS Works

We have several "pools" of VMs attached to UseGalaxy.eu that run user jobs. For trainings we attach a new pool of VMs that is specially labelled for that training. When normal users run tools on our server, these jobs are instructed to avoid the training pools by default.

When your users join a training, using a special URL provided to you, they then are placed in a special training group. Their jobs will then preferentially run on a training machine, and, in the event there is no more capacity, they will run on the main queue. If a spot on a training VM opens up first, they will run there rather than continuing to wait in the main queue.

# The Application Process

> ### {% icon hands_on %} Hands-on: Apply for TIaaS Training
>
> 1. [Fill out the form](https://usegalaxy.eu/request-tiaas)
>
> 2. The EU admin team will review the request and get in touch with you as needed to identify the compute resources you need for your training
>
> 3. The EU admin team will inform you of the URL you should provide to your participants during the training. They can open this URL and they will be added to the special training group
>
{: .hands_on}

# Things to do: Before Your Training

1. Run through your planned trainings on the EU server to ensure everything is available including data libraries and tools.
2. If anything is missing, contact the admins via [Gitter](https://gitter.im/usegalaxy-eu/Lobby) or [Email](mailto:contact@usegalaxy.eu) and they will work to resolve the missing training material resources.

# Things to do: During Your Training

1. Watch for any problems and contact the admins via [Gitter](https://gitter.im/usegalaxy-eu/Lobby)

   - Jobs spending abnormally long in the queue
   - Jobs failing

2. Remind your participants to submit bug reports if they experience any tool errors

# Conclusion
{:.no_toc}
