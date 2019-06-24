---
layout: tutorial_hands_on

title: "Galaxy Installation on Kubernetes"
questions:
- How do I deploy Galaxy on Kubernetes using Helm?
- How can I simple create a replica of usegalaxy.org?
objectives:
- Have an understanding of how to use Galaxy's Helm chart
- Be able to use Helm to install different flavors of Galaxy for different purposes
time_estimation: "30m"
key_points:
- Stock deployment of production Galaxy components on Kubernetes is simple
- Helm chart allows for a lot of deployment flexibility
contributors:
  - pcm32
  - afgane
  - nuwang
  - almahmoud
  - ic4f
tags:
  - kubernetes
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - kubernetes
---

# Overview
{:.no_toc}

TODO
- K8s
- Helm
- Galaxy Helm chart goals

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Galaxy Helm Chart

## Prerequisites

We'll be using the [Galaxy Helm chart](https://github.com/CloudVE/galaxy-kubernetes)
to install and manage a Galaxy deployment. To be able to use this chart, we'll
need access to a Kubernetes cluster, with Helm installed. For development and
testing purposes this can be easily achieved by installing
[Docker Desktop locally and enabling Kubernetes](https://www.docker.com/products/kubernetes). Afterwards, also install [Helm](https://helm.sh).
