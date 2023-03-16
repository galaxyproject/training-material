---
title: Where can I run the hands-on tutorials?
area: learners
layout: faq
box_type: tip
contributors: [bebatut,shiltemann]
---

To run the hands-on tutorials you need a Galaxy server to run them on.

Each tutorial is annotated with information about which [public Galaxy servers](https://galaxyproject.org/public-galaxy-servers/) it can be run on. These servers are available to anyone on the world wide web and some may have all the tools that are needed by a specific tutorial.

If your organization/consortia/community has its own Galaxy server, then you may  want to run tutorials on that. You will need to confirm that all necessary tools and reference genomes are available on your server and possible install missing tools and data. To learn how to do that, you can follow our [dedicated tutorial]({% link topics/teaching/tutorials/setup-galaxy-for-training/tutorial.md %}).

Some topics have a [Docker](https://www.docker.com/) image that can be installed and run on participants' laptops.  These Docker images contain Galaxy instances that include all tools and datasets used in a tutorial, as well as saved analyses and repeatable workflows that are relevant. You will need to [install Docker](https://docs.docker.com/install/).

Finally, you can also run your tutorials on cloud-based infrastructures.  Galaxy is [available on many national research infrastructures](https://galaxyproject.org/galaxy-services/) such as [Jetstream](https://galaxyproject.org/cloud/jetstream/) (United States), [GenAP](https://www.genap.ca/) (Canada), [GVL](https://www.gvl.org.au/) (Australia), [CLIMB](http://www.climb.ac.uk/) (United Kingdom), and more.  These instances are typically easy to launch, and easy to shut down when you are done.

If you are already familiar with, and have an account on [Amazon Web Services](https://aws.amazon.com/) then you can also launch a Galaxy server there using [CloudLaunch](https://launch.usegalaxy.org/).

