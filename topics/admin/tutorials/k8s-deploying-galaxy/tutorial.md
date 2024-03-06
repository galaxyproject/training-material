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


Galaxy has a minimal number of required dependencies, which makes its basic
installation quick for both users and developers. However, configuring a
multi-user production instance is a complex undertaking due to Galaxy’s many
interacting and dependent systems, components, and configurations. Software
containerization has become the preferred method of addressing deployment
challenges across operating environments. Containerization also requires
orchestration, so that multiple containers can work together to deliver a
complex application. [Kubernetes](https://kubernetes.io/) has emerged as the
primary container orchestration technology, as it is both container agnostic and
widely adopted. Kubernetes allows managing, scaling, and deploying different
pieces of an application–in a standardized way–while providing excellent tooling
for doing so.

In this tutorial, we'll take a look at Kubernetes and [Helm](https://helm.sh/)
as tools for deploying containerized Galaxy. The goals for this model of
deploying Galaxy is to use best-practices from the Galaxy community on how to
deploy the Galaxy application in a well-defined package. This model can simplify
deployment and management burden of running Galaxy. While it is possible to
follow this tutorial by simply copying and pasting supplied commands, and a
production-grade Galaxy will be installed, it is desirable to have a basic
understanding of the container concepts and Kubernetes and Helm technologies.

Some of the goals for deploying and running Galaxy in this mode include:
- Design a mostly stateless model for running Galaxy where processes can be
  horizontally scaled as needed
- Integrate components from the Galaxy project ecosystem to leverage existing
  resources
- Provide a unified handling of Galaxy configurations
- Minimize customized dependencies
- Minimize the need to build custom components

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Prerequisites
We'll be using the [Galaxy Helm chart] to install and manage a Galaxy
deployment. To be able to use this chart, we'll need access to a Kubernetes
cluster, with Helm installed. For development and testing purposes this can be
easily achieved by installing
[Docker Desktop locally and enabling Kubernetes][DfD]. Afterwards, also install
[Helm](https://helm.sh).

For production deployments, we'll also need some storage resources for data
persistence. This can be done by either defining a [storage class][SC] or
creating a [Persistent Volume and a corresponding Persistent Volume Claim][PV].
Once created, just keep a note of the resources Persistent Volume Claim ID and
to use later.

For the CVMFS-enabled version of the chart (more on this below), it is also
necessary to run Kubernetes version 1.13 or newer because we'll be using the
[Container Storage Interface (CSI)][CSI].

# Downloading the Galaxy Helm Chart
The Galaxy Helm Chart is currently under active development with enhancements
continuously trickling in. As a result, there are no regular releases yet and
instead we recommend just cloning the GitHub repository with the chart
implementation. This will be the easiest method to keep up with chart changes
for the time being.

> <hands-on-title>Download the chart</hands-on-title>
>
>   Clone the chart repository from the machine where you would like to deploy
>   Galaxy and change into the chart directory.
>
>   {% raw %}
>   ```bash
>   git clone https://github.com/galaxyproject/galaxy-helm
>   cd galaxy-helm/galaxy
>   ```
>   {% endraw %}
{: .hands_on}

# Deploying Galaxy
The Galaxy Helm chart packages best-practice solutions for deploying Galaxy
into a single package that can be readily deployed as a unit. Behind the
scenes, all the supporting services are started and configured into an
interoperable system. Specifically, this involves starting a database service
based on Postgres, using Nginx as a web proxy, and running an independently
scalable set of web and job handler processes for Galaxy. This follows
the production-quality deployment recommendation setup for Galaxy and leverages
some of the Kubernetes features to help with running long-term services (e.g.,
liveness probes that automatically restart stuck processes).

## Deploying the Default Configuration
The default set of values for the Galaxy chart configures only a minimal set
of Galaxy options necessary. The configured options are required for suitable
operation of the system. Setting other options will depend on the environment
and it's best to refer to the general [Galaxy documentation][GxyDocs]; we'll
also take a look at how to make configuration changes in the context of the
chart later in this tutorial.

> <hands-on-title>Deploying the Galaxy Helm Chart</hands-on-title>
>
> 1. First, we need to fetch any dependencies for the chart. One of the
>    advantages of using Helm is that we can reuse best-practice deployment
>    methods for other software right out of the box by relying on published
>    charts and integrating them into the Galaxy chart.
>
>    {% raw %}
>    ```bash
>    helm dependency update
>    ```
>    {% endraw %}
>
> 2. We can now deploy Galaxy via the Chart. Before running this command make
>    sure you are in the chart source code directory (where `values.yaml` file
>    resides) and note the trailing dot. Running this command will create a new
>    Helm release (i.e., chart installation) called `galaxy`.
>
>    {% raw %}
>    ```bash
>    helm install --name galaxy .
>    ```
>    {% endraw %}
>
> 3. It will take about a minute or two for the database to be initialized,
>    necessary containers downloaded, and Galaxy processes started. Ultimately,
>    while this may depend on the Kubernetes cluster setup you are using,
>    Galaxy should be available at the root URI for the given machine. We can
>    always check the status of our release by typing `helm status galaxy`.
>
{: .hands_on}

## Deploying a CVMFS-enabled Configuration
The Galaxy Helm chart also comes with a more comprehensive set of configuration
options that leverage more of the Galaxy project ecosystem. In practice this
means deploying Galaxy with the same toolset as that of
_[usegalaxy.org](https://usegalaxy.org/)_ right out of the box. It's important to note
that this deployment configuration leverages all the same chart components but
just defines more configuration options. Namely, we attach to the
[Galaxy CVMFS][CVMFS] ready-only file system for retrieving the tool
configurations while leveraging [BioContainers] for resolving tool dependencies.

> <hands-on-title>Deploying the CVMFS-enabled Configuration</hands-on-title>
>
> 1. If you are following this tutorial sequentially and have a release of
>    Galaxy already running, let's delete it (assuming that's fine and you have
>    no data to keep). More details about the deletion process are available in
>    the [Deleting a Deployment section](#deleting-a-deployed-helm-release). If
>    you're just playing around, run `helm delete --purge galaxy`.
>
> 2. The CVMFS variant of the Galaxy chart has an additional dependency on the
>    [Galaxy CVMFS chart](https://github.com/CloudVE/galaxy-cvmfs-csi-chart).
>    We'll deploy this chart into its own [Namespace] to keep its resources
>    nicely grouped. We'll also fetch the chart from a packaged chart
>    repository instead of its GitHub repo.
>
>    {% raw %}
>    ```bash
>    kubectl create namespace cvmfs
>    helm repo add galaxy https://raw.githubusercontent.com/CloudVE/helm-charts/master/
>    helm repo update
>    helm install --name cvmfs --namespace cvmfs galaxy/galaxy-cvmfs-csi
>    ```
>    {% endraw %}
>
> 3. We can now install the CVMFS-enabled set of values.
>
>    {% raw %}
>    ```bash
>    helm install --name galaxy galaxy/galaxy
>    ```
>    {% endraw %}
>
> 4. Again, it will take a few minutes for Galaxy to start up. This time most of
>    the waiting is due to the tool definition files to be cached on CVMFS and
>    loaded into the tool panel. We can check the status of the deployment by
>    running `helm status galaxy`. We can also watch the boot process by tailing
>    the logs of the relevant container with a command similar to
>    `kubectl logs -f galaxy-web-7568c58b94-hjl9w` where the last argument is
>    the name of the desired pod, as printed following the `helm install`
>    command. Once the boot process has completed, we can access Galaxy at
>    `/galaxy/` URI (note the trailing `/`; it's significant).
>
{: .hands_on}

## Deleting a Deployed Helm Release
After we're done experimenting with an installation of the chart, we can just
as easily delete all the resources as we've created them. However, that may
not be desirable so make sure you understand the system you're working on to
avoid undesired surprises. Namely, deleting and recreating a Helm release is
generally not a problem where the processes will just respawn and everything will go back to operational; however, underlying storage configuration may
interfere here with all the application data being potentially lost. This
predominantly depends on how the relevant storage class was configured.

> <hands-on-title>Deleting a Deployed Helm Release</hands-on-title>
>
> 1. Before we delete a deployment, let's ensure we understand what will happen
>    with the underlying storage used by Galaxy.
>
>    {% raw %}
>    ```bash
>    $ kubectl get pv
>    NAME                                       CAPACITY   ACCESS MODES   RECLAIM POLICY   STATUS   CLAIM                                   STORAGECLASS      REASON   AGE
>    cvmfs-cache-pv                             1000Mi     RWX            Retain           Bound    cvmfs/cvmfs-cache-pvc                   manual                     31m
>    pvc-55806281-96c6-11e9-8e96-0251cc6c62f4   1Gi        ROX            Delete           Bound    default/galaxy-cvmfs-gxy-data-pvc       cvmfs-gxy-data             28m
>    pvc-5580c830-96c6-11e9-8e96-0251cc6c62f4   1Gi        ROX            Delete           Bound    default/galaxy-cvmfs-gxy-main-pvc       cvmfs-gxy-main             28m
>    pvc-55814757-96c6-11e9-8e96-0251cc6c62f4   10Gi       RWX            Delete           Bound    default/galaxy-galaxy-pvc               nfs-provisioner            28m
>    pvc-70d4cc48-96be-11e9-8e96-0251cc6c62f4   8Gi        RWO            Delete           Bound    default/data-galaxy-galaxy-postgres-0   nfs-provisioner            84m
>    pvc-8cb27bc9-9679-11e9-8e96-0251cc6c62f4   100Gi      RWO            Delete           Bound    cloudman/data-nfs-provisioner-0         ebs-provisioner            9h
>    ```
>    {% endraw %}
>
>    As we can see in the command output, the storage resources associated with
>    the current deployment have the reclaim policy set to `Delete`, which will
>    happen once no resources are using the given resource. If what you see is
>    the not the intended behavior, you can change the [reclaim policy].
>
> 2. Once we're ok with the state of the resources and are ready to delete a
>    a deployment, we can do so with the following commands:
>
>
>    {% raw %}
>    ```bash
>    helm delete --purge galaxy
>    helm delete --purge cvmfs
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
