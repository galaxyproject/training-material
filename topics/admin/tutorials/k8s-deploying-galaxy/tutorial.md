---
layout: tutorial_hands_on

title: "Galaxy Installation on Kubernetes"
level: "Intermediate"
questions:
- How do I deploy Galaxy on Kubernetes using Helm?
- How can I create a simple replica of usegalaxy.org?
objectives:
- Have an understanding of how to use Galaxy's Helm chart
- Be able to use Helm to install different flavors of Galaxy for different purposes
time_estimation: "30m"
key_points:
- Stock deployment of production Galaxy components on Kubernetes is simple
- Helm chart allows easy configuration changes
contributors:
  - pcm32
  - afgane
  - nuwang
  - almahmoud
  - jdavcs
tags:
  - kubernetes
subtopic: cloud
follow_up_training:
  - type: "internal"
    topic_name: admin
    tutorials:
      - k8s-managing-galaxy
priority: 2
---

# Galaxy Helm Chart

## Overview

This tutorial describes how to use the Galaxy Helm Chart to deploy a production
grade instance of Galaxy on a Kubernetes cluster. The Helm Chart has been designed
to follow best practices adopted by the community, including the usegalaxy.* federation,
and will install a Galaxy with the following features by default:
- Zero-downtime configuration changes and upgrades
- Scalable web and job handlers
- Automatic failure recovery based on liveness and readiness probes
- A built-in nginx for efficiently serving large files
- TUSD for resumable uploads
- Celery for background jobs
- Access to CVMFS reference data
- A toolset matching the usegalaxy.* federation (also served off CVMFS)
- Interactive tools (wildcard DNS mapping required)
- Minimal privileges, with jobs running as non-root and only having access to datasets they need
- Automatic maintenance scripts to cleanup the galaxy database and tmp directories

Optionally, the chart can be configured with
- High-availability components - this includes trivial scaling of clustered Postgres, Rabbit MQ etc.
- Replacement components - You can replace the built-in operators with a managed or existing Postgres database (e.g. Amazon RDS), RabbitMQ cluster etc.
- Use S3 as an alternative to CVMFS
- Automatic scraping of metrics which can be sent to Influxdb

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Prerequisites
Some familarity with Kubernetes is assumed. This includes general administrative
familarity and how to install and configure Helm Charts.

A running Kubernetes cluster is also required (1.27 or higher), with Helm
(3.5 or higher) configured to access it.
For development and testing purposes this can be easily achieved by installing
[Docker Desktop locally and enabling Kubernetes][DfD]. Afterwards, install
[Helm](https://helm.sh).

For production deployments, we'll also need some storage resources for data
persistence. This can be done by either defining a [storage class][SC] or
creating a [Persistent Volume and a corresponding Persistent Volume Claim][PV].
Once created, just keep a note of the resources Persistent Volume Claim ID and
to use later.

## Deploying the Default Configuration
The default set of values for the Galaxy chart configures only a minimal set
of Galaxy options necessary. The configured options are required for suitable
operation of the system. Setting other options will depend on the environment
and it's best to refer to the general [Galaxy documentation][GxyDocs]; we'll
also take a look at how to make configuration changes in the context of the
chart later in this tutorial.

> <hands-on-title>Deploying the Galaxy Helm Chart</hands-on-title>
>
> 1. First, we need to add the helm repository for the chart. The chart is
> automatically packaged, versioned and uploaded to a helm repository on github
> with each accepted PR. Therefore, the latest version of the chart can be directly
> installed from that repository.
>
>    {% raw %}
>    ```bash
>    helm repo add galaxyproject https://raw.githubusercontent.com/galaxyproject/helm-charts/master/
>    helm repo update
>    ```
>    {% endraw %}
>
> 2. We can now deploy Galaxy via the Chart. Running this command will create a new
>    Helm release (i.e., chart installation) called `mygalaxy`.
>
>    {% raw %}
>    ```bash
>    helm install mygalaxy galaxyproject/galaxy
>    ```
>    {% endraw %}
>
> 3. It will take about a minute or two for the necessary containers to download,
>    the database to initialize, and Galaxy processes to start. Ultimately,
>    while this may depend on the Kubernetes cluster setup you are using,
>    Galaxy should be available at https://<root URI>/galaxy for the given machine. We can
>    always check the status of our release by typing `helm status galaxy`.
>
{: .hands_on}

## Setting the admin user and changing the brand
The chart is designed to follow standard Kubernetes and Helm idioms, and therefore,
it should be intuitively similar to the steps required to change configuration in
any other Helm chart. For example, ingress paths, resource allocations, container
images etc. can be changed following standard helm conventions. The list of
available configuration options are also documented in the Galaxy Helm
_[Chart repository](https://github.com/galaxyproject/galaxy-helm/tree/master?tab=readme-ov-file#configuration)_

To change Galaxy specific configuration, such as setting the admin user or change the brand in `galaxy.yml`,
we can follow the following steps. Once done, we will also rollback our change to demonstrate how Helm manages
configuration.

> <hands-on-title>Setting admin user and changing the brand</hands-on-title>
>
> 1. Modify the following entries in your `mygalaxy.yml`. Make sure to add these
>    keys under the `configs:` section of the file.
>
>    {% raw %}
>    ```
>    configs:
>      galaxy.yml:
>        galaxy:
>          brand: "Hello World"
>          admin_users: "admin@mydomain.com"
>    ```
>    {% endraw %}
>
> 2. Now, let’s upgrade the chart to apply the new configuration.
>
>    {% raw %}
>    ```bash
>    helm upgrade --reuse-values -f mygalaxy.yml mygalaxy galaxyproject/galaxy
>    ```
>    {% endraw %}
>
> 3. Inspect the currently set Helm values by:
>
>    {% raw %}
>    ```bash
>    helm get values mygalaxy
>    ```
>    {% endraw %}
>
> 4. List the installed Helm charts again and note that the revision of the chart has changed as expected.
>
>    {% raw %}
>    ```bash
>    helm list
>    NAME  	  REVISION	UPDATED                 	STATUS  	CHART                 	APP VERSION	NAMESPACE
>    mygalaxy	2       	Wed Jun 26 14:51:17 2023	DEPLOYED	galaxy-5.14.2          	v24.0.2    	default
>    ```
>    {% endraw %}
>
> 5. Revisit the Galaxy Application in your browser to check whether the settings have changed. This will
>    take a short while (< 1 minute) for the new container to come up. You should experience no downtime.
>
> 6. Let’s now roll back to the previous revision.
>
>    {% raw %}
>    ```bash
>    helm rollback mygalaxy 1
>    ```
>    {% endraw %}
>
>    Use `helm get values` again to observe that the values have reverted to
>    the previous revision. After a short while, once the new container is up
>    and running, Kubernetes will automatically switch over to it and you can
>    see that the previous configuration has been restored.
>
{: .hands_on}

## Deleting a Deployed Helm Release
By default, the Helm chart is designed to install all required dependencies, so that it's easy
to get an instance up and running quickly for experimentation. However, in production, we
recommend installing the dependency charts separately, once per cluster, by installing
Galaxy with helm options
`--set postgresql.deploy=false --set s3csi.deploy=false --set cvmfs.deploy=false --set rabbitmq.deploy=false`.

This is particularly important during uninstallation, where orderly destruction of dependencies is often required
For example, if the rabbitmq operator is uninstalled before the rest of the Galaxy helm chart is deleted, there will be
no operator left to cleanup rabbitmq resources. Installing the aforementioned operators separately sidesteps this problem.

> <hands-on-title>Deleting a Deployed Helm Release</hands-on-title>
>
>    {% raw %}
>    ```bash
>    helm delete mygalaxy
>    helm delete mycvmfs # and any other operators
>    ```
>    {% endraw %}
>
{: .hands_on}

# Next Steps
This tutorial covers the basics of getting Galaxy deployed on Kubernetes using
Helm. There is a lot more to understanding all the configuration options for
the chart and the available deployment models. For more info on some of these
topics, take a look at the [Galaxy Helm chart] repository as well as
[other tutorials](../..) tagged with _kubernetes_. Also, feel free to reach out
on Gitter: [https://gitter.im/galaxyproject/FederatedGalaxy][fedG].


[Galaxy Helm chart]: https://github.com/galaxyproject/galaxy-helm
[DfD]: https://www.docker.com/products/kubernetes
[SC]: https://kubernetes.io/docs/concepts/storage/storage-classes/
[PV]: https://kubernetes.io/docs/concepts/storage/persistent-volumes/
[CSI]: https://kubernetes.io/blog/2019/01/15/container-storage-interface-ga/
[GxyDocs]: https://galaxyproject.org/admin/#configuration
[CVMFS]: https://training.galaxyproject.org/training-material/topics/admin/tutorials/cvmfs/tutorial.html
[BioContainers]: https://biocontainers.pro/#/
[Namespace]: https://kubernetes.io/docs/concepts/overview/working-with-objects/namespaces/
[reclaim policy]: https://kubernetes.io/docs/tasks/administer-cluster/change-pv-reclaim-policy/
[fedG]: https://gitter.im/galaxyproject/FederatedGalaxy
