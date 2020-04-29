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

For the purpose of this tutorial, sample datasets have been created from data downloaded from [C3S](https://climate.copernicus.eu/) through
[Copernicus Climate Data Store](https://cds.climate.copernicus.eu/#!/home):
- [E-OBS daily gridded meteorological data for Europe from 1950 to present derived from in-situ observations](https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-gridded-observations-europe?tab=overview)
- [Fifth generation ECMWF reanalysis for the global climate and weather (ERA5)](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means)
- [Essential climate variables for assessment of climate variability from 1979 to present](https://cds.climate.copernicus.eu/cdsapp#!/dataset/ecv-for-climate-change?tab=overview)

To reduce the volume of data, the data resolution (in space and/or time) has been significantly reduced and/or data has been selected on sample locations (Paris, Oslo and 
Freiburg). The data format may also have been changed (for instance to tabular) to ease processing. 

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *climate101*.
>    {% include snippets/create_new_history.md %}
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    https://zenodo.org/record/3697454/files/tg_ens_mean_0.1deg_reg_v20.0e_Paris_daily.csv
>    https://zenodo.org/record/3697454/files/ts_cities.csv
>    https://zenodo.org/record/3697454/files/era5_Paris.csv
>    https://zenodo.org/record/3697454/files/ecv_1979.nc
>    https://zenodo.org/record/3697454/files/ecv_2018.nc
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Check that the datatype is **csv** or **netcdf** (files ending with `.nc`)
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


## Climate versus Weather

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

# What is the weather like in Paris?

In order to answer this question, we are going to inspect and visualize the dataset `tg_ens_mean_0.1deg_reg_v20.0e_Paris_daily.csv` using simple galaxy tools.

> ### {% icon hands_on %} Hands-on: Daily temperature time series
>    > ### {% icon comment %} Tip: search for the tool
>    >
>    > Many different tools can be used to answer to the questions. Here we give you some guidelines to help you to choose.
>    > Use the **tools search box** at the top of the tool panel to find **Select lines that match an expression** {% icon tool %} and **Datamash** {% icon tool %}.
>    {: .tip}
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What was the average temperature in Paris on the 14th of July 2003?
>    > 2. What is the minimum and maximum temperatures in Paris?
>    > 3. On which dates did the minimum and maximum temperatures occured?
>    > 4. Can we observe a trend (cooling/warming)?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. The average temperature in Paris on the 14th of July 2003 was 26.73 degrees Celcius. It can be found by using **Select lines that match an expression** {% icon tool %} with parameter **"the pattern"** set to 2003-07-14.
>    > > 2. The minimum temperature in Paris is -11.6799995 degrees celcius and the maximum temperature in Paris is 33.579998 degrees celcius. To find out, you can use **Datamash** {% icon tool %} with the following parameters:
>    > >      - {% icon param-file %} *"Input tabular dataset"*: `tg_ens_mean_0.1deg_reg_v20.0e_Paris_daily.csv`
>    > >      - *"Input file has a header line"*: `Yes`
>    > >      - *"Print header line"*: `Yes`
>    > >      - "Print all fields from input file": `No`
>    > >      - In *"Operation to perform on each group"*:
>    > >          - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>    > >              - *"Type"*: `minimum`
>    > >              - *"On column"*: `c2`
>    > >          - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>    > >              - *"Type"*: `maximum`
>    > >              - *"On column"*: `c2`
>    > > 
>    > > 3. The minimum temperature (-11.6799995 degrees celcius) was observed on January 16 1985. The maximum temperature (33.579998 degrees celcius) was observed on July 25 2019. You can use different Galaxy tools to find out the solution and here we show you how to use **Datamash**  {% icon tool %} with the following parameters:
>    > >      - {% icon param-file %} *"Input tabular dataset"*: `tg_ens_mean_0.1deg_reg_v20.0e_Paris_daily.csv`
>    > >      - *"Input file has a header line"*: `Yes`
>    > >      - *"Print header line"*: `Yes`
>    > >      - "Print all fields from input file": `Yes`
>    > >      - In *"Operation to perform on each group"*:
>    > >          - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>    > >              - *"Type"*: `minimum`
>    > >              - *"On column"*: `c2`
>    > > 
>    > > For the **maximum**, repeat **Datamash**  {% icon tool %} with the following parameters:
>    > >      - {% icon param-file %} *"Input tabular dataset"*: `tg_ens_mean_0.1deg_reg_v20.0e_Paris_daily.csv`
>    > >      - *"Input file has a header line"*: `Yes`
>    > >      - *"Print header line"*: `Yes`
>    > >      - "Print all fields from input file": `Yes`
>    > >      - In *"Operation to perform on each group"*:
>    > >          - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>    > >              - *"Type"*: `maximum`
>    > >              - *"On column"*: `c2`
>    > > 
>    > {: .solution}
>    {: .question}
> 
{: .hands_on}

# What is the climate in Paris?

To get some information about the (past and current) climate in Paris, we will first look at monthly averages.

## Seasonality

> ### {% icon hands_on %} Hands-on: What is the monthly climatological temperature in Paris?
>
>   To answer to this question, we will compute the global average temperatures over the entire period 1950 and 2019 for each month (January, February, etc.). Indeed, 
>   this period of time is sufficiently long for computing monthly climatological temperature (more than 30 years).
>    > ### {% icon question %} Questions
>    > 
>    > 1. What is the warmest summer month e.g. between June, July and August (JJA) in Paris?
>    > 2. What is the coolest winter month e.g. between December, January and February (DJF) in Paris?
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. The warmest summer month in Paris is July (19.921018171429 degrees celcius). However, it is interesting to remark that on our dataset we see very little difference in the mean temperature between July and August.
>    > > 2. The coolest winter month in Paris is January (4.4669169722484 degrees celcius).
>    > >
>    > > Below, we show you how we found these results.
>    > > We will first split all the dates (first column) from YYYY-MM-DD (where YYYY is the year, MM the month and DD the day) to three column to get 3 columns: one for the year, one for the month and one for the day:
>    > > 
>    > > Please note that you may use other Galaxy tools to reach the same results.
>    > > Results can be slightly different when using different source of climate information. However, you will always observe the same pattern e.g. cool month in winter and warm month on summer. We can also clearly see that Paris has a mild climate with on average no extreme temperatures.
>    > {: .solution}
>    {: .question}
> 
{: .hands_on}

> ### {% icon tip %} Tip: Using existing climatologies
>
> In this tutorial, we compute manually the monthly climatological temperatures to explain you the algorithm used behing.
>  However, many data providers have pre-computed climatologies and can be directly downloaded. For instance, on the [CDS](https://cds.climate.copernicus.eu/cdsapp#!/search?type=dataset), climatologies are provided for [Essential climate variables for assessment of climate variability from 1979 to present](https://cds.climate.copernicus.eu/cdsapp#!/dataset/ecv-for-climate-change?tab=overview).
{: .tip}


## Yearly average

> ### {% icon hands_on %} Hands-on: What is the trend (cooling/warming) in the climate for Paris between 1950 and 2019?
>
>  To answer to this question, we will compute yearly mean of the temperature in Paris and visualize it.
>
>    > ### {% icon question %} Questions
>    > 
>    > 1. Can we easily observe a trend?
>    > > ### {% icon solution %} Solution
>    > >
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

## Anomalies


Climate and weather data can have **biases** (explain why we use anomalies rather than absolute values in climate); discuss uncertainties around observations and numerical models.


> ### {% icon hands_on %} Hands-on: Climate stripes for Paris
>
>    > ### {% icon question %} Questions
>    > 
>    > 1. Do you observe a warming or cooling between 1950 and 2019?
>    > > ### {% icon solution %} Solution
>    > >
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

> ### {% icon tip %} Tip: Copernicus Climate Bulletin
>
> [Copernicus Climate Bulletins](https://climate.copernicus.eu/climate-bulletins) presents the current condition of the climate using key climate change indicators. 
> They also provide data, analysis of the maps and guidance on how they are produced. Datasets for temperature anomalies can be found and are 
> regularly updated (with recent dates). For instance, in March 2020, the corresponding dataset can be found [here](https://climate.copernicus.eu/sites/default/files/2020-04/ts_1month_anomaly_Global_ea_2t_202003_v01.csv).
>
{: .tip}


## Ensemble

Explain why we use ensembles (more than one source of information and/or perturbations) in climate.


> ### {% icon hands_on %} Hands-on: Heatmap for different models (use ERA5 data)
>
>    > ### {% icon question %} Questions
>    > 
>    > 1. Do you observe a warming or cooling between 1979 and 2019?
>    > > ### {% icon solution %} Solution
>    > >
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


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
 
*Source: [https://gcos.wmo.int/en/essential-climate-variables/ecv-factsheets](https://gcos.wmo.int/en/essential-climate-variables/ecv-factsheets)*

> ### {% icon hands_on %} Hands-on: Essential Climate Variables
>
>    > ### {% icon question %} Questions
>    > 
>    > 1. XXX
>    > > ### {% icon solution %} Solution
>    > >
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}


## Past, present and future climate?


When we talk about climate data, the type of data can vary significantly. We have very little actual observations at the scale of climate and usually not covering a large area. In addition to observations, we can make use of:
- Re-analyses where observations and numerical modelling are combined together. 
- Climate models.

Observations and re-analyses provide information about the past and current climate while climate models can provide past, current and future climate information.
When it comes to future climate, we usually need to make some assumptions (such as how much CO2 emissions, etc.) and make different scenarios e.g. we run climate models using different assumptions and look at future trends under each of these scenarios: this is what we call **climate projections**. Climate projections will be discussed in a separate Galaxy tutorial.


# Conclusion

{:.no_toc}

We have learnt to differentiate climate from weather and got an overview of the terminology used by climate scientists to identify the 
various source of climate data. 
