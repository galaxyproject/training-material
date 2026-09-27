---
layout: tutorial_hands_on
title: Display marine Antarctic species occurences & environemental information on QGIS
questions:
- How to use QGIS through Galaxy to create maps showing species occurences on environmental context
objectives:
- Deal with QGIS execution management
- Access, filter, and import GIS data using WFS webservice
time_estimation: 0H30M
key_points:
- From WFS webservice, you can now access GIS remote data, build a query to select a subpart of the data, and import it to Galaxy through QGIS Galaxy interactive tool.
tags:
  - GIS
  - Geographical Information System
  - spatial data
  - maps
  - OGC
  - biodiversity
  - ecology
  - data visualization
contributions:
  authorship:
  - Meyane
  - C0c0-C00k13
  - TuturBaba
  - yvanlebras
  funding:
  - mnhn
  - pndb
  - sorbonneuniv
  - ISYEB

subtopic: ecologyviz
---

Dealing with Biodiveristy data, notably species occurences, scientists are intensively using QGIS solution to display such biological information on a geographical context. To illustrate this kind of task, we will use Antarctic geographical context to create a map presenting marine species occurences obtained from OBIS, environmental information from IBSCO, and basemap from NPI/Quantarctica.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Managing QGIS Galaxy interactive tool

QGIS is integrated in Galaxy as an interactive tool. This kind of tool works differently than classical tools as it allows the user to interact with a dedicated graphical interface.

To use QGIS, you need to use the {% tool [dedicated form](interactive_tool_qgis) %}, you can specify input datasets from your history you want to use in QGIS, or not ;), then press the execute button to launch a dedicated QGIS instance. When the graphical user interface of QGIS is ready to be used, a URL will be displayed at the top of the Galaxy center panel. If you don't see it, you can see and access it through the "Interactive Tools" space of the left panel or you can click on {% icon galaxy-eye %} on the tool in the history.

Once you finish your work on QGIS, if you want to reuse data and/or the entire project, you need to save files in the "output" folder (which you can find in the "working" directory). 
Then, quit QGIS properly through the "Project" Menu tab top left and click on "Exit QGIS".

# Importing data to Galaxy

The starting point is to gather occurences and environmental data with an Antarctic base map in your Galaxy history

## Antarctic species occurences

> <hands-on-title>Seraching and downloading marine occurences data from OBIS</hands-on-title>
>
> 1. 
{: .hands_on}

## Antarctic bathymetry


> <hands-on-title>Downloading bathymetry from GEBCO IBSCO </hands-on-title>
>
> 1. 
{: .hands_on}

> <question-title></question-title>
>
> 1. Are you seeing how QGIS shows the fact that a layer is resulting from a query made on a larger layer?
>
> > <solution-title></solution-title>
> >
> > 1. There is a dedicated "filter" icon next to the name of the layer.
> >
> {: .solution}
>
{: .question}

## Antarctic basemap

> <hands-on-title>Downloading bathymetry from GEBCO IBSCO </hands-on-title>
>
> 1. 
{: .hands_on}

# Launching a QGIS instance on your data

# Using QGIS to create a Antarctic biodiversity map


# Conclusion

You just did a classical GIS operation, importing species occurences, environemtal and geographical context basemap data into QGIS, creating a QGIS project, displaying features of interest saved in a QGIS project archive you can share broadly through Galaxy!
