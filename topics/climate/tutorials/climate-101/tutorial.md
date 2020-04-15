---
layout: tutorial_hands_on
enable: true
title: Getting your hands-on climate data
zenodo_link: ''
questions:
- What is climate?
- What type of data is available?
objectives:
- Learn about the terminology
- Learn about the different source of climate data
- Learn about climate observations, reanalysis, climate predictions and climate projections
time_estimation: 1H
key_points:
- Weather versus Climate
- Essential Climate Variables
- Observations, reanalysis, predictions and projections.
contributors:
- annefou

---


# Introduction
{:.no_toc}

> ### {% icon comment %} Comment
>
> This tutorial is significantly based on [Getting your hands-on Climate data](https://nordicesmhub.github.io/climate-data-tutorial/).
>
{: .comment}

The practical aims at familiarzing you with Climate Science and the terminology used by climate scientists. The target audience is not a climate scientist but 
anyone interested in learning about climate. 

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### {% icon comment %} Background
>
> [European Copernicus Climate Change Service (C3S)](https://climate.copernicus.eu/) provide authoritative information about the past, present
> and future climate. C3S is one of the many services provided by Copernicus, the European Union's Earth Observation Programme, looking
> at our planet and its environment for the ultimate benefit of all European citizens.
> The C3S [Climate Data Store (CDS)](https://cds.climate.copernicus.eu/#!/home) provides a single point of access to a wide range of 
> quality-assured climate datasets distributed in the cloud.
> Access to the CDS data is open, free and unrestricted.
> We will be using freely available datasets from the CDS, including 
> observations, historical climate data records, estimates of Essential Climate Variables (ECVs) derived from Earth observations, 
> global and regional climate reanalyses of past observations, seasonal forecasts and climate projections. 
{:  .comment}

In this tutorial, we will be using data downloaded from [C3S](https://climate.copernicus.eu/):
- Data collected from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/#!/home)
- Data prepared and issued by [Copernicus Climate bulletins](https://climate.copernicus.eu/climate-bulletins).

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *climate101* for example...
>    {% include snippets/create_new_history.md %}
> 2. Import the file from [Zenodo]() or from the shared data library
>
>    ```
>    https://zenodo.org/record/3697454/files/ecv_1979.nc
>    https://zenodo.org/record/3697454/files/ecv_2018.nc
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Check that the datatype is **netcdf**
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 4. Add a tag to the dataset corresponding to `copernicus`
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# What is climate?

According to [wikipedia](https://en.wikipedia.org/wiki/Climate), 
Climate is defined as the average state of everyday's weather condition **over a period of 30 years**. It is measured by assessing
the patterns of variation in temperature, humidity, atmospheric pressure, wind, precipitation, atmospheric particle count and
other meteorological variables in a given region over long periods of time.
Climate differs from weather, in that weather only describes the short-term conditions of these variables in a given region. 


# Climate versus Weather

Quantities that climate scientists are interested in are similar to those used to assess the weather (temperature, precipitation, etc.). 
But there is a big difference between climate and weather: **weather** varies from hour to hour and from day to day whereas **climate** 
is defined as the average of weather over several decades or longer.

The figure below shows a woman walking her dog and we can use it to make an analogy to illustrate the difference between weather and climate. 
if you focus your attention on the dog, you can see that it is all over the place, sometimes upwards, sometimes downwards: this can represent the weather and its 
variability. The dog (weather) is not following a fully random pattern and varies around a main direction (trend) that is given by the woman: the woman is representing
the climate and gives us an indication of where both the woman and dog are likely to be in the future.
  
 ![Illustrate the difference beteen weather and climate](../../images/weather_versus_climate.png "Weather versus Climate")

*Source: [Animated short on statistics](https://youtu.be/e0vj-0imOLw) from Norwegian infotainment program Siffer. Produced by TeddyTV for NRK. Animation by Ole Christoffer Haga*

You can also watch this [Video](https://youtu.be/e0vj-0imOLw) to get an animated illustration of the difference between climate and weather.

Climate questions are different from weather questions:
- Will it rain tomowrow?                  *<- weather*
- By how much will global temperature at the end of the century be warmer than the beginning of the century?                 *<- climate*
- What could happen if CO2 emissions double within the next century?                 *<- climate*

# Climate variables

*Temperature* is often the first variable that comes to mind when we talk about climate. However, it is insufficient to fully characterize the climate, and scientists have agreed on a number of variables to systematically observe Earth`s changing climate. 

That is what we call *Essential Climate Variables*.

## Essential Climate Variables

The [Global Climate Observing System](https://gcos.wmo.int/) (GCOS) and its GCOS expert panels maintain definitions of [Essential Climate Variables](https://gcos.wmo.int/en/essential-climate-variables) (ECVs). 

GCOS is co-sponsored by the [World Meteorological Organization](https://public.wmo.int/en) (WMO), the [Intergovernmental Oceanographic Commission of the United Nations Educational, Scientific and Cultural Organization](http://www.ioc-unesco.org/) (IOC-UNESCO), the [United Nations Environment Programme](https://www.unenvironment.org/) (UN Environment), and the [International Science Council](https://council.science/) (ISC). It regularly assesses the status of global climate observations of the atmosphere, land and ocean and produces guidance for its improvement.

At the moment, there are [54 ECVs](https://gcos.wmo.int/en/essential-climate-variables/ecv-factsheets) divided in 3 categories:
- Atmosphere
- Land
- Ocean

<img src="../fig/ECVs_GCOS.png" width="70%"/> 
*Source: [https://gcos.wmo.int/en/essential-climate-variables/ecv-factsheets](https://gcos.wmo.int/en/essential-climate-variables/ecv-factsheets)*

> ## Which Climate variables do you plan to use?
> Write in the workshop etherpad the list of variables you think would be useful for your study and at the same time write what you would like to do with climate data (area of application). 
> This is a preliminary list of variables and we will re-discuss it later after learning what is generally on offer.
>
{: .challenge}

# Types of climate data resources

When we talk about climate data, the type of data can vary significantly. We have very little actual observations at the scale of climate and usually not covering a large area. 

The type of climate data you will be using greatly depends on the period of time you are interested in:

- Observations
- Re-analyses
- Climate models


## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *Panoply* for example...
>    {% include snippets/create_new_history.md %}
> 2. Import the file from [Zenodo]() or from the shared data library
>
>    ```
>    https://zenodo.org/record/3697454/files/ecv_1979.nc
>    https://zenodo.org/record/3697454/files/ecv_2018.nc
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Check that the datatype is **netcdf**
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 4. Add a tag to the dataset corresponding to `copernicus`
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Conclusion
{:.no_toc}

We have learnt to differentiate climate from weather and got an overview of the terminology used by climate scientists to identify the 
various source of climate data. 
