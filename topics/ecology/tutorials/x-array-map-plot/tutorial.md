---
layout: tutorial_hands_on
title: Visualization of Climate Data using  NetCDF Xarray Map Plotting
zenodo_link: '10.5281/zenodo.6600820'
questions:
- What is xarray map plotting?
- How to plot NetCDF data using xarray map plotting tool?
- What are the different types of projections and colormap options available in the tool?
objectives:
- Learn about plotting Climate Data
- Learn about NetCDF xarray map plotting tool
- Learn about visualzing the climate data using NetCDF xarray map plotting using the different kinds of projections aand colormap options  
time_estimation: 1H
key_points:
- NetCDF xarray climate data visualization
- NetCDF xarray map plotting projections
- NetCDF xarray map plotting colormaps
contributors:
- Quickbeasts51429

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->
>### {% icon comment %}Comment
>
> The tutorial aims for establishing a good knowledge about meaningful visualization of climate data. It is beginner friendly and requires not much knowledge about the tool.
>
{: .comment}

It is beginner friendly and requires not much knowledge about the tool.

> ### Agenda
>
> In this tutorial we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

>### {% icon comment %} Background
>
>According to [UN](https://www.un.org/en/climatechange/what-is-climate-change) , Climate is the long term shift in temperature and weather patterns which may be due to natural or artificial causes. To learn more about climate, refer this [tutorial](https://training.galaxyproject.org/training-material/topics/climate/tutorials/climate-101/tutorial.html) from the GTN. Due to the frequently changing nature of the weather patterns, the data collected is huge in size.
The climate data is mainly represented in these three categories : NetCDF (Network Common Data Form), HDF (Hierarchical Data Format) , GRIB (GRIdded Binary or General Regularly-distributed Information in Binary form).
>
>NetCDF datasets are basically used for storing multidimensional data which generally consits variables such as temperature, precipitation , direction , etc. The variation of climate variables over a period of time is suitably plotted using this dataset. The entire earth is divided into both horizontal as well as verticle coordinates which makes plotting of the variables such as the oocean temperatures possible.
>
>The coordinate system, types of projections and colormaps are some of the very important considerations in achieving the most suitable visualization option.



# Plotting air temperature at 2 metres using the ECMWF Reanalysis Data

The data used in this tutorial is ECMWF Reanalysis . We will be concerned with the following variables : air temperature at 2 metres, latitude, longitude and time. Our main objective is to plot the global air temperature at 2 metres with respect to time. For this we will be using the netCDF xarray tool exsisting in the Galaxy Europe ( or your favourable Galaxy Instance ) server.

It will be a fun learning experience for anyone who loves visualization ! 

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. Name it as per your choice. My suggestions : *ECMWF_Reanalysis*.
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 2. Import the files from
>    the remote files
>    ```
>    (`Upload Data` -> `Choose remote files` -> `ECMWF ERA5 Reanalysis` >   ->`2022` -> `05` -> `data` -> `air_temperature_at_2_metres.nc`)
>    ```
>
> 3. Check that the datatype of uploaded data is **netCDF**. 
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
>    If it is not `netCDF` make sure to convert it using the Galaxy built-in format converters.
>
>    {% snippet faqs/galaxy/datasets_convert_datatype.md conversion="Convert h5 to netCDF" %}
>
> 4. Rename Datasets if you find the name to be too long or you have somthing more meaningful in your mind.
>
> 5. Add to each database a tag corresponding to `ecmwf`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# What is Xarray?
Xarray , formerly known as xray, is a python package which enable us to play with gridded data. This package shares most of its features from numpy , but in a more convinient manner by keeping track of labels in arrays. The gridded data is mainly available in netCDF data format. Thus **xarray comes very handy while dealing with netCDF files.**


## What is Visualisation in Xarray?
Xarray uses **Cartopy** and **Matplotlib** as the two main tools for creating detailed and informative plots.
Cartopy is a python package used for geospatial data analysis. In the Python library, matplotlib is the most used scientific plotting library.
For a multidimensional data consisting of latitudes and longitudes along with the other variables, xarray has the capability of appling cartopy map projections.

# Splitting the dataset  using  seltimestep and splithour 

After loading the required data, it comes to obtaining the meta info or meta data of the file. The very purpose of these steps are to obtain the information about dimensions, variables, global attributes, etc. The coordinate info helps to know about the actual data entries present under the various variables.  

Follow the below steps :

## Finding the  **NetCDF xarray Metadata Info**

> ### {% icon hands_on %} Hands-on: netCDF dataset with Xarray metadata Galaxy Tool
>
> 1. {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `air_temperature_at_2_metres.nc` (Input dataset)
>
>
> 2. View {% icon galaxy-eye%} the two generated outputs:
>    - `Metadata infos` is a `tabular` providing the list of variables, their dimension names and number of elements per dimension. This file is used by other Xarray Tools. 
>    - The second file `info file` provide a summary of the **Xarray Dataset** contained in your netCDF file.
{: .hands_on}

In `info file` output file, we can identify 4 different sections:
1. **Dimensions**: name of dimensions and corresponding number of elements;
2. **Coordinates**: contains coordinate arrays (longitude, latitude, level and time) with their values.
3. **Data variables**: contains all the variables available in the dataset. Here, we only have one variable. For each variable, we get information on its shape and values.
4. **Global Attributes**: at this level, we get the global attributes of the dataset. Each attribute has a name and a value.  


> ### {% icon question %} Questions 
>
> 1. What is the name of the variable for air temperature at 2 metres and its physical units?
>
>
> > ### {% icon solution %} Solution
> > 1. Information about variable names and units can be found in **info file** that was generated by Xarray metadata Galaxy Tool. 
> >      - Variable name: `air_temperature_at_2_metres`
> >      - Units: `K`
> >
> > 
> > > ### {% icon code-out %} Output
> > > ```bash
> > > xarray.Dataset {
> > > dimensions:
> > > lat = 721 ;
> > > lon = 1440 ;
> > > time0 = 595 ;
> > > 
> > > 
> > > variables:
> > >   float32 lon(lon) ;
> > >     lon:standard_name = longitude ;
> > >     lon:long_name = longitude ;
> > >     lon:units = degrees_east ;
> > >  
> > >   float32 lat(lat) ;
> > >     lat:standard_name = latitude ;
> > >     lat:long_name = latitude ;
> > >     lat:units = degrees_north ;
> > >  
> > >   datetime64[ns] time0(time0) ;
> > >  		time0:standard_name = time ;
> > >  
> > >  	float32 air_temperature_at_2_metres(time0, lat, lon) ;
> > >  		air_temperature_at_2_metres:standard_name = air_temperature ;
> > >  		air_temperature_at_2_metres:units = K ;
> > >  		air_temperature_at_2_metres:long_name = 2 metre temperature ;
> > >  		air_temperature_at_2_metres:nameECMWF = 2 metre temperature ;
> > >  		air_temperature_at_2_metres:shortNameECMWF = 2t ;
> > >  		air_temperature_at_2_metres:nameCDM = 2_metre_temperature_surface;
> > >  		air_temperature_at_2_metres:product_type = analysis ;
> > >  
> > >  
> > > // global attributes:
> > >  
> > >   :source = Reanalysis ;
> > >   :institution = ECMWF ;
> > > 	:title = ERA5 forecasts ;
> > > }
> > > ```
> > {: .code-out}
> {: .solution }
{: .question }

## Sub-step with **NetCDF xarray Coordinate Info**

> ### {% icon hands_on %} Hands-on:  Get Coordinate information with Xarray Coordinate
>
> 1. {% tool [NetCDF xarray Coordinate Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_coords_info/xarray_coords_info/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `output` (Input dataset)
>2. View {% icon galaxy-eye%} the 5 generated outputs:
>    - `lat`: a tabular file containing all the latitude values of our Xarray dataset;
>    - `lon`: a tabular file containing all the longitudes values;
>    - `time0`: this tabular file contains all the forecast times. In our case, these are relative to 25/05/2022, 18:00:00 UTC.;
>    - `version`: this is a text file returning the Xarray package version. It is useful when publishing your Galaxy workflow.

>
>    > ### {% icon comment %} Comment
>    >
>    > This tool returns as many tabular files as the number of coordinate variables present in your input file. The values are decoded from the netCDF input file and no further processing is done. So units for instance for latitudes, longitudes, level and time may vary from one file to another depending on how it was coded in the original input file.
>    {: .comment}
>
{: .hands_on}


> ### {% icon question %} Understanding air temperature coordinates at 2 metres
>
>
> 1. What is the format of time coordinate?
> 2. What is the range of values for latitude and longitude?
>
> > ### {% icon solution %} Solution
> > 1. The `info file` tells us that `time0` is coded as `timedelta64[ns]` e.g. as differences in times (here in nanoseconds). Here the reference time is ,May, 25, 2022. If we look at the tabular file named `time0` (generated by `NetCDF xarray Coordinate Info`), we see that these times are automatically converted to human readable time format when printed:
> >
> > > ### {% icon code-out %} Output
> > > ```bash
> > >    0 2022-05-01 00:00:00
> > > 1	2022-05-01 01:00:00
> > > 2	2022-05-01 02:00:00
> > > 3	2022-05-01 03:00:00
> > > 4	2022-05-01 04:00:00
> > > 
> > > ```
> > {: .code-out}
> > This tells us that we have hourly forecast data. The last forecast time is `4 days 00:00:00` which means that the last forecast is in 4 days at 00:00 UTC (from May 25, 2022).
> {: .solution }
{: .question }

## Operations on Climate data using **CDO Operations**

> ### {% icon hands_on %} Hands-on: Defining a particular time range using seltimestep 
>
> 1. {% tool [CDO Operations](toolshed.g2.bx.psu.edu/repos/climate/cdo_operations/cdo_operations/2.0.0+galaxy0) %} with the following parameters:
>    - In *"CDO Operators"*:
>        - {% icon param-repeat %} *"Insert CDO Operators"*
>            - *"Select cdo operator"*: `seltimestep (Select timesteps)`
>                - *"Timesteps for selection"*: `594/595`
>                - {% icon param-file %} *"Additional input file"*: `output` (Input dataset)
>
>   Our main aim is to plot the last two hour data on the last day present in the dataset. The data present is in the form of daily hourly time frame. Thus we need to split the dataset into smaller part upto which we want to plot. 
>
>
>
>    > ### {% icon comment %} Comment
>    >
>    > The syntax of using the `seltimestep` is `(initial data number / final data entry)`. An important thing to pay attention is how your data entries are number: are they numbered starting from 0 or 1. Accordingly we can add or skip adding 1 to the data number to attain the desired result.
>    {: .comment}
>
{: .hands_on}



> ### {% icon question %} Questions
>
> 1. How will you plot the last two hours separately?
> 
>
> > ### {% icon solution %} Solution
> >
> > 1. {% tool [CDO Operations](toolshed.g2.bx.psu.edu/repos/climate/cdo_operations/cdo_operations/2.0.0+galaxy0) %} with the following parameters:
> >- In *"CDO Operators"*:
>        - {% icon param-repeat %} *"Insert CDO Operators"*
> >  - *"Select cdo operator"*: `splithour (Split hours)`
>                - *"Timesteps for selection"*: `594/595`
>                - {% icon param-file %} *"Additional input file"*: `outfile.netcdf` (Input dataset) generated from the previous step.
>
> >
> {: .solution}
>
{: .question}


## Finding the  **NetCDF xarray Metadata Info**

> ### {% icon hands_on %} Hands-on: netCDF dataset with Xarray metadata Galaxy Tool for the hourly plots
>
> 1. {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `outfile_17.nc` (Input dataset). The meta data generated from this step can be used for second file too due to same meta-data of both the datasets(which is quite obvious).
>
>
> 2. View {% icon galaxy-eye%} the two generated outputs:
>    - `Metadata infos` is a `tabular` providing the list of variables, their dimension names and number of elements per dimension. This file is used by other Xarray Tools. 
>    - The second file `info file` provide a summary of the **Xarray Dataset** contained in your netCDF file.
{: .hands_on}
## Map Plotting using **NetCDF xarray map plotting**

> ### {% icon hands_on %} Hands-on: Map Plotting
> The air temperatures at the last two hours of the day are plotted here : 
>
> 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %}
*"Netcdf file"*: `outfile_17.nc`
>    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from outfile_17.nc`
 *"Tabular of variables"*: `output` (output of **NetCDF xarray Metadata Info** {% icon tool %})
>    - *"Choose the variable to plot"*: `air_temperature_at_2_metres`
>    - *"Name of latitude coordinate"*: `lat`
>    - *"Name of longitude coordinate"*: `lon`
>    - *"Datetime selection"*: `No`
>    - *"Range of values for plotting e.g. minimum value and maximum value (minval,maxval) (optional)"*: `220,320`
>    - *"Add country borders with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add coastline with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add ocean with alpha value [0-1] (optional)"*: `1.0`
>    - *"Specify plot title (optional)"*: `Projection :  Mercator 17:00 UTC `
>    - *"Specify which colormap to use for plotting (optional)"*: `lajolla`
>    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator', 'central_longitude': 12.0}`
>
>    > ### {% icon comment %} Why shifting longitudes?
>    >
>    > Longitudes are coded from 0 to 360 degrees. As we do not have global data but only covering Europe, we need to shift longitudes so that `NetCDF xarray map plotting` can plot properly our dataset. 
>    {: .comment}
>
> 
{: .hands_on}



> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **NetCDF xarray map plotting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %}
*"Netcdf file"*: `outfile_18.nc`
>    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from outfile_17.nc`
 *"Tabular of variables"*: `output` (output of **NetCDF xarray Metadata Info** {% icon tool %})
>    - *"Choose the variable to plot"*: `air_temperature_at_2_metres`
>    - *"Name of latitude coordinate"*: `lat`
>    - *"Name of longitude coordinate"*: `lon`
>    - *"Datetime selection"*: `No`
>    - *"Range of values for plotting e.g. minimum value and maximum value (minval,maxval) (optional)"*: `220,320`
>    - *"Add country borders with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add coastline with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add ocean with alpha value [0-1] (optional)"*: `1.0`
>    - *"Specify plot title (optional)"*: `Projection :  Mercator 18:00 UTC `
>    - *"Specify which colormap to use for plotting (optional)"*: `lajolla`
>    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator', 'central_longitude': 12.0}`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}




# Conclusion
{:.no_toc}

We have learnt about the xarray map plotting tool dealing with the netcdf data set. The tutorial also discussed about the types of climate datasets. One of the tutorial is info about usage of different [colormaps](https://github.com/Quickbeasts51429/Xarray_ColorMaps/blob/main/index.md#color-maps) and [projections]() in xarray. 