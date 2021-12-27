---
layout: tutorial_hands_on

title: Pangeo ecosystem 101 for everyone - Introduction to Xarray Galaxy Tools
zenodo_link: 'https://doi.org/10.5281/zenodo.5805953'
questions:
- What Xarray Galaxy Tools can I use in Galaxy and what for?
- What is an Xarray?
- How do I use Xarray in Galaxy?
- How to get metadata information?
- How to make a selection?
- How to visualize?
- How to filter? 
- How to make reduction operations (mean, max, min)?
- How to resample my data?
objectives:
- Understand what Pangeo and Xarray are
- Learn to get metadata information using Xarray Galaxy Tools
- Learn to select data
- Learn to visualize geographical data on a map
- Learn to filter, make reduction operations (mean, max, min)
- Learn to resample my data
time_estimation: 1H
key_points:
- Xarray Tools in Galaxy
contributors:
- annefou

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

[Pangeo](https://pangeo.io/) is a project that effectively began in 2016 with a workshop at Columbia University. The mission for Pangeo developed at that workshop is still valid nowadays:

*Our mission is to cultivate an ecosystem in which the next generation of open-source analysis tools for ocean, atmosphere and climate science can be developed, distributed, and sustained. These tools must be scalable in order to meet the current and future challenges of big data, and these solutions should leverage the existing expertise outside of the geoscience community.*

In this tutorial, you will learn how to manipulate [netCDF](https://en.wikipedia.org/wiki/NetCDF) data file using [Xarray](https://xarray.pydata.org/en/stable/) Galaxy Tools. 


> ### {% icon comment %} Xarray and Earth Science
> Xarray works with labelled multi-dimensional arrays and can be used for a very wide range of data and data formats. However, in this training material, we focus on the usage of Xarray for Earth Science data following [CF-Convention](https://cfconventions.org/). However, some Galaxy Tools also work for non Earth Science datasets and if needed current Xarray Galaxy Tools could be extended to accomodate new usage.
{: .comment}

In this tutorial, we will be using data from [Copernicus Atmosphere Monitoring Service](https://ads.atmosphere.copernicus.eu/)
and more precisely PM2.5 ([Particle Matter < 2.5 μm](https://en.wikipedia.org/wiki/Particulates#Size,_shape_and_solubility_matter)) 4 days forecast from December, 22 2021. This dataset is very small and there is no need to parallelize our data analysis. Parallel data analysis with Pangeo is not covered in this tutorial. 


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Create a history

> ### {% icon hands_on %} Hands-on: Create history
>
> 1. Make sure you start from an empty analysis history.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. **Rename your history** to be meaningful and easy to find. For instance, you can choose **Pangeo 101 for everyone - Xarray** as the name of your new history.
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}


##  Upload CAMS PM2.5 data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/record/5805953/files/CAMS-PM2_5-20211222.netcdf
>    ```
>
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add a tag corresponding to ads (for Atmosphere Data Service)
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}


# Understanding our dataset


> ### {% icon details %} More information about CAMS PM2.5 forecast datasetss
>
> Our CAMS PM2.5 forecast dataset is in [netCDF](https://en.wikipedia.org/wiki/NetCDF) format. You could find the same dataset in different formats such as [GRIdded Binary or General Regularly-distributed Information in Binary form (GRIB)](https://en.wikipedia.org/wiki/GRIB) or [geoTIFF](https://en.wikipedia.org/wiki/GeoTIFF). The same Xarray Tools can be used with these other data formats.
{: .details}


To understand what is contained in our dataset, we will first use Xarray metadata Galaxy Tool. That will give us all the metadata information about the dataset.


## Get metadata

### Global metadata information

> ### {% icon hands_on %} Hands-on: netCDF dataset with Xarray metadata Galaxy Tool
>
>
{: .hands_on}


We can identify 4 different sections:
1. **Dimensions**: name of dimensions and corresponding number of elements;
2. **Coordinates**: contains coordinate arrays (longitude, latitude, level and time) with their values.
3. **Data variables**: contains all the variables available in the dataset. Here, we only have one variable. For each variable, we get information on its shape and values.
4. **Attributes**: at this level, we get all the attributes of the dataset. 


> ### {% icon question %} Questions CAM PM2.5 Dataset
>
> What is the name of the variable for Particle matter < 2.5 μm and its physical units?
>
> > ### {% icon solution %} Solution
> > 1. Information about variable names and units can be found in **info file** that was generated by Xarray metadata Galaxy Tool. 
> >      - Variable name: `mass_concentration_of_pm2p5_ambient_aerosol_in_air`
> >      - Units: `µg/m3`
> >
> > 
> > > ### {% icon code-out %} Output
> > > ```bash
> > > xarray.Dataset {
> > > dimensions:
> > > 	latitude = 400 ;
> > > 	level = 1 ;
> > > 	longitude = 700 ;
> > > 	time = 97 ;
> > > 
> > > variables:
> > > 	float32 longitude(longitude) ;
> > > 		longitude:long_name = longitude ;
> > > 		longitude:units = degrees_east ;
> > > 	float32 latitude(latitude) ;
> > > 		latitude:long_name = latitude ;
> > > 		latitude:units = degrees_north ;
> > > 	float32 level(level) ;
> > > 		level:long_name = level ;
> > > 		level:units = m ;
> > > 	timedelta64[ns] time(time) ;
> > > 		time:long_name = FORECAST time from 20211222 ;
> > > 	float32 pm2p5_conc(time, level, latitude, longitude) ;
> > > 		pm2p5_conc:species = PM2.5 Aerosol ;
> > > 		pm2p5_conc:units = µg/m3 ;
> > > 		pm2p5_conc:value = hourly values ;
> > > 		pm2p5_conc:standard_name = mass_concentration_of_pm2p5_ambient_aerosol_in_air ;
> > > 
> > > // global attributes:
> > > 	:title = PM25 Air Pollutant FORECAST at the Surface ;
> > > 	:institution = Data produced by Meteo France ;
> > > 	:source = Data from ENSEMBLE model ;
> > > 	:history = Model ENSEMBLE FORECAST ;
> > > 	:FORECAST = Europe, 20211222+[0H_96H] ;
> > > 	:summary = ENSEMBLE model hourly FORECAST of PM25 concentration at the Surface from 20211222+[0H_96H] on Europe ;
> > > 	:project = MACC-RAQ (http://macc-raq.gmes-atmosphere.eu) ;
> > > }
> > > ```
> > {: .code-out}
> {: .solution }
{: .question }

### Coordinates information

> ### {% icon hands_on %} Hands-on: Get Coordinate information with Xarray Coordinate
>
>
{: .hands_on}


> ### {% icon question %} Understanding PM2.5 forecast coordinates
>
>  1. What is the units of the `time` coordinate?
>  2. What is the frequency of PM2.5 forecasts? 
>  3. What is the range of values for latitudes and longitudes?
>
> > ### {% icon solution %} Solution
> > 
> {: .solution }
{: .question }

# Plotting our dataset on a geographical map

> ### {% icon hands_on %} Hands-on: Map plot
>  We will use Xarray mapplot Galaxy Tool to plot PM2.5 on December 22, 2021.
>
{: .hands_on}


> ### {% icon question %} Visualize and compare
>
> Make a plot to Visualize the forecast for December, 24th 2021 at 12:00 UTC. Do you see any obvious differences with the plot from December 22, 2021 at 00:00 UTC?
>
> > ### {% icon solution %} Solution
> > 
> > Data starts on December, 22nd 2021 at 00:00 UTC so we need to add 2 days and 12 hours to select the correct time index.
> >
> {: .solution }
{: .question }

# Select / Subset from coordinates


> ### {% icon hands_on %} Hands-on: NetCDF xarray operations manipulate xarray from netCDF and save back to netCDF
>  
>
{: .hands_on}


> ### {% icon question %} PM2.5 over Italy
>
> Using a Multi-plot between Rome and Naples, can you tell us if the forecasted PM2.5 will increase or decrease during the first 24 hours?
>
> > ### {% icon solution %} Solution
> > 
> >
> {: .solution }
{: .question }

> ### {% icon comment %}  `latitude=slice(47.3, 36.5)` and not `latitude=slice(36.5, 47.3)`
> Why did we slice latitudes with `latitude=slice(47.3, 36.5)` and not `latitude=slice(36.5, 47.3)`?
> - because when using slice, you need to specify values using the same order than in the coordinates. Latitudes are specified in 
> decreasing order for CAMS.
>
{: .comment}

# Masking with Where statement

- Sometimes we may want to make more complex selections with criteria on the values of a given variable and not only on its coordinates. For this we use `where`.
- For instance, we may want to only keep PM2.5 if values are greater than 25 μm.m-3 (or any threshold you would like to choose)

> ### {% icon hands_on %} Hands-on: Plot where PM2.5 is greater than 0
> 
{: .hands_on}



> ### {% icon question %} PM2.5 over Italy over 35 μm.m-3
>
> Over the first 24 hour forecast, will PM2.5 exceed 35 μm.m-3?
>
> > ### {% icon solution %} Solution00:00 UTC so we need to add 2 days and 12 hours to select the correct time index.`
> >
> {: .solution }
{: .question }

# From Xarray to Tabular Data

> ### {% icon hands_on %} Hands-on: Xarray selection
> 
{: .hands_on}



> ### {% icon question %} PM2.5 at Naples over the 4 forecasted days
>
> From a qualitative point of view, can you say if PM2.5 may increase or decrease over the 4 forecasted days?
>
> > ### {% icon solution %} Solution
> >  - Select a single location with Xarray selection, use a scatter plot with ggplot2 and/or climate stripes
> >
> {: .solution }
{: .question }


# Conclusion
{:.no_toc}

{% icon trophy %} Well done! In this tutorial, Xarray Galaxy Tools have been introduced and we learned to use these tools on a real dataset from Copernicus Atmosphere Monitoring Service. We encourage you to try with your own datasets. 
