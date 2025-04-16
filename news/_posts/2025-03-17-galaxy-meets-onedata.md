---
title: Galaxy meets Onedata — distributed storage for science
contributions:
  authorship:
    - lopiola
  funding:
    - eurosciencegateway
    - egi
tags: [esg-wp4, esg, contributing]
cover: news/images/2025-03-17-galaxy-and-onedata.png
coveralt: Galaxy now features comprehensive integration with Onedata
layout: news
---

Thanks to the efforts undertaken in the
[EuroScienceGateway](https://galaxyproject.org/projects/esg/) project, Galaxy
now offers integration with Onedata. It can be used as 
a **remote source for data import/export** (a.k.a. Files Source Plugin) and as 
a **storage location for Galaxy datasets** (an Object Store). 
The integration includes 
[BYOD]({% link topics/galaxy-interface/tutorials/onedata-remote-import/tutorial.html %})
(Bring Your Own Data) and
[BYOS]({% link topics/galaxy-interface/tutorials/onedata-byos/tutorial.html %})
(Bring Your Own Storage) templates.

All the good stuff has been integrated in Galaxy version **24.2** and is now
**live** at the EU server: [https://usegalaxy.eu](https://usegalaxy.eu).


## Is for me?

You can give it a go even if you have never used Onedata before. Take a look at
the [tutorial]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html#accessing-onedata-services %})
that will take you through creating your first free account in the Onedata
sandbox (demo) environment. If you like it and want more, 
[let us know](https://onedata.org/#/home/contact) — we'll see what can be done to
organize a fully-fledged Onedata ecosystem for your use cases.

If your organization is already a part of a Onedata ecosystem or you have access to
Onedata services for science, like the [EGI DataHub](https://datahub.egi.eu),
then it's definitely worth looking into.


## What's the point?

Good question! It's 
[all explained in this tutorial]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html#whats-the-point %}).


## Getting started

There is a range of tutorials concerning Onedata on the 
[Galaxy training portal]({% link search2?query=onedata %}). We recommend to 
[start with this one]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html %}),
which will give you a basic understanding about Galaxy & Onedata integration and
guide you through the next steps.

The below image shows the Onedata file browser UI and the Galaxy **Upload** tool
making it possible to import data stored in Onedata Spaces.
![Onedata & Galaxy interfaces](../../../../topics/galaxy-interface/images/onedata-getting-started/galaxy-onedata-side-by-side.png)


## A quick overview of Onedata

Onedata ([https://onedata.org](https://onedata.org)) is a data management platform that
provides easy and unified access to globally distributed storage resources. It's an 
[open-source project](https://github.com/onedata), developed by the team from the
Cyfronet Computing Center in Krakow (Poland), since 2013.

Onedata creates a POSIX-compatible, virtual file system layer spanning
geographically distributed data providers hosting heterogeneous storage
resources. The virtualized data can be access using multiple interfaces: Web
GUI, fuse-based POSIX mount, REST API, CDMI API, Python libraries, or S3.
Regardless of the interface, the user gets the same, unified view of all his
data.

![Onedata virtual FS](../../../../topics/galaxy-interface/images/onedata-getting-started/onedata-virtual-fs.png)

The data management capabilities of Onedata and its ability to facilitate
cross-organizational collaborative data sharing may be a great asset for your
data processing pipelines in Galaxy. 

Read on [in the tutorial]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html %})!
