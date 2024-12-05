---
title: What Galaxy instance should I use for my training?
area: instructors
layout: faq
box_type: tip
contributors: [bebatut,shiltemann]
---


To teach the hands-on tutorials you need a Galaxy server to run the examples on.

Each tutorial is annotated with the information on which [public Galaxy servers](https://galaxyproject.org/public-galaxy-servers/) it can be run. These servers are available to anyone on the world wide web and some may have all the tools that are needed by a specific tutorial. If you choose this option then you should work with that server's admins to confirm that the server can handle the workload for a workshop. For example, the [usegalaxy.eu](https://usegalaxy.eu/)

If your organization/consortia/community has its own Galaxy server, then you may want to run tutorials on that. This can be ideal because then the instance you are teaching on is the same as your participants will be using after the training. They'll also be able to revisit any analysis they did during the training. If you pursue this option you'll need to work with your organization's Galaxy Admins to confirm that

- the server can support a room full of people all doing the same analysis at the same time.
- all tools and reference datasets needed in the tutorial are locally installed.  To learn how to setup a Galaxy instance for a tutorial, you can follow our [dedicated tutorial]({% link topics/teaching/tutorials/setup-galaxy-for-training/tutorial.md %}).
- all participants will be able to create/use accounts on the system.

Some training topics have a Docker image that can be installed and run on all participants' laptops.  These images contain Galaxy instances that include all tools and datasets used in a tutorial, as well as saved analyses and repeatable workflows that are relevant.

Finally, you can also run your tutorials on cloud-based infrastructures.  Galaxy is [available on many national research infrastructures](https://galaxyproject.org/galaxy-services/) such as [Jetstream](https://galaxyproject.org/cloud/jetstream/) (United States), [GenAP](https://www.genap.ca/) (Canada), [GVL](https://launch.genome.edu.au/launch) (Australia), [CLIMB](http://www.climb.ac.uk/) (United Kingdom), and more.


