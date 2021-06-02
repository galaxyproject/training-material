---
layout: tutorial_hands_on

title: "Galaxy Interactive Tools"
questions:
  - "What are Galaxy Interactive Tools (GxIT)"
  - "How to create them?"
objectives:
  - "Discover what Galaxy Interactive Tools (GxIT) are"
  - "Be able to create Galaxy Interactive Tools (GxIT)"
  - "Be able to add a Galaxy Interactive Tools (GxIT) in a Galaxy instance"
time_estimation: '1h'
contributors:
 - eancelet
 - yvanlebras
key_points:
  - "A Galaxy Interactive Tools (GxIT) provides an easy way to integer external interactive tools into Galaxy"
  - "Galaxy Interactive Tools (GxIT) can be for example Jupyter notebooks, RStudio or R Shiny apps"
  - "Galaxy Interactive Tools (GxIT) are intensively using containers technology, notamment Docker"
  - "With a minimal amount of code you can integrate external tools into Galaxy if already containerized"
---

## Introduction

In this tutorial we are going to demonstrate how to integer an extrenal interactive environment as a Galaxy Interactive Tool (GxIT). Galaxy Interactive Tools (GxIT) are accessible on the Galaxy tool panel as any Galaxy tool. This functionality has been used to integer notably JupyterLab or RStudio server as some R Shiny apps. To see Galaxy Interactive Tools example, you can have a look at the "live subdomain" European Galaxy instance https://live.usegalaxy.eu/.

In this tutorial, we will look at a concrete example, an app named `HeatMap`. This Galaxy interactive tool is composed of 4 elements:
  - a [dockerfile](https://raw.githubusercontent.com/workflow4metabolomics/gie-shiny-heatmap/master/Dockerfile)
  - a xml
  - a pinch of config ;)
  - The original [R script](https://raw.githubusercontent.com/workflow4metabolomics/gie-shiny-heatmap/master/app.R) responsible to launch the R shiny app inside the container.




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


### "A first bunch of code!" (XML + conf)

Create your first Galaxy Interactive Tool using an existing Docker image!

> ### {% icon hands_on %} Hands-on
>
> 1. Create a Galaxy tool xml file named `interactive_heatmap.xml` :
>
>    > ### {% icon tip %} Tip: some Galaxy and XML related tips
>    >
>    > * You can take inspiration from [xml d'askomics](https://github.com/galaxyproject/galaxy/blob/dev/tools/interactive/interactivetool_askomics.xml)
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



### Distributing your GxIT

Two types of materials have to be distributed:
- the Galaxy tool xml
- the Docker image

Here we are focusing on the Galaxy tool but the Docker image must be available on Dockerhub. 
If planemo and toolshed supported it, it looks like you need to deposit on iuc. In the meantime we can say to put it in https://github.com/galaxyproject/galaxy/? https://github.com/galaxyproject/galaxy/tree/dev/tools/interactive



## Conclusion
