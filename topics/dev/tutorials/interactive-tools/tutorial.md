---
layout: tutorial_hands_on

title: "Galaxy Interactive Tools"
questions:
  - "What is an Interactive Tool on Galaxy (GxIT)"
  - "How do GxIT work?"
  - "How do I build an GxIT?"
  - "How do I debug an GxIT?"
  - "How do I deploy an GxIT?"
objectives:
  - "Discover what Galaxy Interactive Tools (GxIT) are"
  - "Understand when a GxIT is appropriate"
  - "Understand how GxITs are structured"
  - "Be able to wrap a dockerised application as a GxIT"
  - "Be able to install the new GxIT to an existing Galaxy instance"
requirements:
  -
    type: "internal"
    topic_name: Galaxy tool development
    tutorials:
      - tool-from-scratch
  -
    title: "Docker basics"
    type: "none"
  -
    title: "A Galaxy instance configured to serve interactive tools"
    type: "None"
time_estimation: '1h'
subtopic: tooldev
key_points:
  - "Galaxy Interactive Tools (GxIT) provides an easy way to integrate external interactive tools into Galaxy"
  - "Example GxITs are Jupyter notebooks, RStudio or R Shiny apps"
  - "GxITs rely heavily on container technologies like Docker"
  - "If a tool is already containerized, it can be integrated rapidly into Galaxy as a GxIT"
  - "In theory, any containerized web application can be wrapped as a GxIT"
contributors:
  - eancelet
  - yvanlebras
  - neoformit

---

## Introduction

In this tutorial we will demonstrate how to build and deploy a Galaxy Interactive Tool (GxIT). Installed GxITs are accessible through the Galaxy tool panel like any Galaxy tool. To see some example GxITs, take a look at the "live" subdomain of Galaxy EU: https://live.usegalaxy.eu/.

In this tutorial, we will look at a concrete example, an app named `HeatMap`. This Galaxy interactive tool is composed of 4 elements:
  - a [dockerfile](https://raw.githubusercontent.com/workflow4metabolomics/gie-shiny-heatmap/master/Dockerfile)
  - a xml
  - a pinch of config ;)
  - The original [R script](https://raw.githubusercontent.com/workflow4metabolomics/gie-shiny-heatmap/master/app.R) responsible to launch the R shiny app inside the container.


## How do interactive tools work?

Interactive tools are a special breed of Galaxy tool which have seen relatively
little air time.
They enable the user to run an entire web application through Galaxy, which
opens as a new tab in the browser. This allows users to explore and manipulate
data in a rich interface, such as Jupyter notebooks and RStudio.
Interactive tool development builds on the canonical tool-wrapping process.
Instead of running a command, the tool feeds user input to a docker container
which can then be accessed through a unique subdomain. The user opens their
interactive, makes the observations or manipulations that they require and then
terminates the tool. On termination, the docker container is stopped and
removed.

## Technical development Stack

> Not sure what goes here? Tech requirements?

### Dockerfile

At first let's discover the Dockerfile of this amazing `HeatMap` Galaxy Interactive Tool (GxIT)!

source : [dockerfile](https://raw.githubusercontent.com/workflow4metabolomics/gie-shiny-heatmap/master/Dockerfile)

Github: [heatmap](https://github.com/workflow4metabolomics/gie-shiny-heatmap/)

Mais d'abord, qu'est-ce qu'un dockerfile? Et une image docker? Et un container? Bref, c'est quoi docker?

Commentaire de Lain sur la structure du dockerfile :

```txt
Docker, why?
Isolate a process by building a file system with process isolation based on namespaces.
A docker image is not a VM.
A container on a linux machine can only run native programs.
A docker file is used to build a filesystem with the necessary resources to perform the defined actions

```

```txt

### our based Docker image
FROM quay.io/workflow4metabolomics/gie-shiny:latest

# Installing packages needed for check traffic on the container and kill if none
RUN apt-get update && \
    apt-get install --no-install-recommends -y r-cran-car

# Install R packages
COPY ./packages.R /tmp/packages.R
RUN Rscript /tmp/packages.R

# Build the app
RUN rm -rf /srv/shiny-server/sample-apps && \
    rm /srv/shiny-server/index.html && \
    mkdir -p /srv/shiny-server/samples/heatmap && \
    chmod -R 755 /srv/shiny-server/samples/heatmap && \
    chown shiny.shiny /srv/shiny-server/samples/heatmap

COPY ./app.R /srv/shiny-server/samples/heatmap/app.R
COPY ./static/css/styles.css /srv/shiny-server/samples/heatmap/styles.css

```


**Notes :**

Recommendations for version numbering

There are 2 docker repositories:
https://github.com/abretaud/geoc_gxit_ansible/
hub.docker.com and http://quay.io/

The default port of dockerized RShiny app is 3838

The environment_variables would be needed to retrieve data from Galaxy History. In the example, this is therefore not necessary.

### Tool script


### XML Wrapper + a bit of configuration

Create your first Galaxy Interactive Tool using an existing Docker image!

> ### {% icon hands_on %} Hands-on
>
> 1. Create a Galaxy tool xml file named `interactive_heatmap.xml` :
>
>    > ### {% icon tip %} Tip: some Galaxy and XML related tips
>    >
>    > * You can take inspiration from [the askomics tool XML](https://github.com/galaxyproject/galaxy/blob/dev/tools/interactive/interactivetool_askomics.xml)
>    > * To check the xml syntax https://www.xmlvalidation.com/ or https://www.w3schools.com/xml/xml_validator.asp
>    >
>    {: .tip}
>    > > ### {% icon solution %} Solution
>    > >
>    > >```xml=
>    > >
>    > ><tool id="interactive_tool_heatmap" tool_type="interactive" name="heatmap" version="0.1">
>    > >    <description>An interactive tool for heatmap</description>
>    > >    <requirements>
>    > >        <container type="docker">emetabohub/heatmap</container>
>    > >    </requirements>
>    > >    <entry_points>
>    > >        <entry_point name="heatmap" requires_domain="True">
>    > >            <port>8765</port>
>    > >            <url>/</url>
>    > >        </entry_point>
>    > >    </entry_points>
>    > >    <environment_variables>
>    > >        <environment_variable name="HISTORY_ID">$__history_id__</environment_variable>
>    > >        <environment_variable name="REMOTE_HOST">$__galaxy_url__</environment_variable>
>    > >        <environment_variable name="GALAXY_WEB_PORT">8080</environment_variable>
>    > >        <environment_variable name="GALAXY_URL">$__galaxy_url__</environment_variable>
>    > >    </environment_variables>
>    > >    <command><![CDATA[
>    > >
>    > >        mkdir -p /srv/shiny-server/samples/heatmap &&
>    > >        cp '$infile' /srv/shiny-server/samples/heatmap/inputdata.tsv &&
>    > >        cd /srv/shiny-server/samples/heatmap/ &&
>    > >        Rscript app.R > '$outfile'
>    > >
>    > >    ]]>
>    > >    </command>
>    > >    <inputs>
>    > >        <param name="infile" type="data" format="tabular,csv" label="tsv file"/>
>    > >    </inputs>
>    > >    <outputs>
>    > >        <data name="outfile" format="txt"/>
>    > >    </outputs>
>    > >    <tests>
>    > >    </tests>
>    > >    <help>
>    > ><![CDATA[
>    > >
>    > >Some help is always of interest ;).
>    > >
>    > >]]>
>    > >    </help>
>    > >
>    > ></tool>
>    > >```
>    > {: .solution}

{: .hands_on}

![chronometre](https://i.imgur.com/pXWdClv.png)


## Testing locally

## Installing on a Galaxy production instance (set up with Ansible)

## Distributing your GxIT

Two types of materials have to be distributed:
- the Galaxy tool xml
- the Docker image

Here we are focusing on the Galaxy tool but the Docker image must be available on Dockerhub.
If planemo and toolshed supported it, it looks like you need to deposit on iuc. In the meantime we can say to put it in https://github.com/galaxyproject/galaxy/? https://github.com/galaxyproject/galaxy/tree/dev/tools/interactive

## Other examples

## Final Notes
Is it the same as Conclusion?

## Conclusion
