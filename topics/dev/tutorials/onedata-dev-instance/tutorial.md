---
layout: tutorial_hands_on

title: "Setting up a dev Onedata instance"
subtopic: advanced
tags:
  - storage
contributors:
  - lopiola

time_estimation: "20m"
level: Introductory

requirements:
 - type: "internal"
   topic_name: galaxy-interface
   tutorials:
     - onedata-getting-started

questions:
- How to set up a minimal Onedata dev instance?
- How can I test the Galaxy & Onedata integration in a sandbox environment?
objectives:
- Quickly set up a zero-config dockerized Onedata deployment.
- Run Onedata locally for Galaxy development and testing purposes.
key_points:
- A local Onedata deployment can be easily set up using the dockerized "demo mode".
- You may set up Onedata and use the sandbox credentials to connect it to Galaxy.
- A local Onedata deployment may be useful to test your changes to the data-related logic of Galaxy.
---


# Prerequisites

It's recommended that you have basic knowledge about Onedata. If needed, follow 
[this tutorial]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html %})
first!


# Introduction

If you are familiar with **docker**, you can easily set up a local Onedata
sandbox deployment. It may come in handy if you wish to test Galaxy & Onedata
integration or develop features related to data management in Galaxy.


# Demo mode

Use the official Onedata tutorial for running the system in
[demo mode](https://onedata.org/#/home/documentation/topic/stable/demo-mode).

If you are new to this, run a Onezone and Oneprovider instance in foreground.
When the setup is finished, log in to the Onezone interface, using the IP
address and credentials displayed in Onezone docker logs (as explained in 
[the guide](https://onedata.org/#/home/documentation/topic/stable/demo-mode)).

Use the Onezone **docker IP** as the **Onezone domain** in all Onedata related
configs in Galaxy. Use the predefined Space called **demo-space** for Storage
Location configs (global or BYOS).

> <warning-title>The demo mode runs on self-signed certificates</warning-title>
> You will need to accept the untrusted connection in your browser to enter the UI. 
> 
> Remember to disable TLS certificate validation in the Galaxy configs!
{: .warning}


Follow the 
[getting started tutorial]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html#access-tokens %})
to acquire an access token, or use the convenience scripts of the demo mode
(as explained in 
[the guide](https://onedata.org/#/home/documentation/topic/stable/demo-mode)).


# Related materials

The demo mode is used for setting up Onedata in Galaxy's Object Store
integration tests. 
[Take a look at the code](https://github.com/galaxyproject/galaxy/blob/dev/test/integration/objectstore/_base.py)
(search for "onedata").
