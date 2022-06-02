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
> {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from
>    the remote files
>    ```
>    (`Upload Data` -> `Choose remote files` -> `ECMWF ERA5 Reanalysis` >->      `2022` -> `05` -> `data` -> `air_temperature_at_2_metres.nc`)
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

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
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

## Sub-step with **NetCDF xarray Coordinate Info**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [NetCDF xarray Coordinate Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_coords_info/xarray_coords_info/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
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

## Sub-step with **CDO Operations**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [CDO Operations](toolshed.g2.bx.psu.edu/repos/climate/cdo_operations/cdo_operations/2.0.0+galaxy0) %} with the following parameters:
>    - In *"CDO Operators"*:
>        - {% icon param-repeat %} *"Insert CDO Operators"*
>            - *"Select cdo operator"*: `seltimestep (Select timesteps)`
>                - *"Timesteps for selection"*: `594/595`
>                - {% icon param-file %} *"Additional input file"*: `output` (Input dataset)
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
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

## Sub-step with **NetCDF xarray map plotting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular of variables"*: `output` (output of **NetCDF xarray Metadata Info** {% icon tool %})
>    - *"Choose the variable to plot"*: ``
>    - *"Name of latitude coordinate"*: ``
>    - *"Name of longitude coordinate"*: ``
>    - *"Datetime selection"*: `No`
>    - *"Range of values for plotting e.g. minimum value and maximum value (minval,maxval) (optional)"*: `220,320`
>    - *"Add country borders with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add coastline with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add ocean with alpha value [0-1] (optional)"*: `1.0`
>    - *"Specify plot title (optional)"*: `Projection :  Mercator 17:00 UTC `
>    - *"Specify which colormap to use for plotting (optional)"*: `lajolla`
>    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator', 'central_longitude': 12.0}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
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

## Sub-step with **NetCDF xarray map plotting**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Tabular of variables"*: `output` (output of **NetCDF xarray Metadata Info** {% icon tool %})
>    - *"Choose the variable to plot"*: ``
>    - *"Name of latitude coordinate"*: ``
>    - *"Name of longitude coordinate"*: ``
>    - *"Datetime selection"*: `No`
>    - *"Range of values for plotting e.g. minimum value and maximum value (minval,maxval) (optional)"*: `220,320`
>    - *"Add country borders with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add coastline with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add ocean with alpha value [0-1] (optional)"*: `1.0`
>    - *"Specify plot title (optional)"*: `Projection :  Mercator 18:00 UTC `
>    - *"Specify which colormap to use for plotting (optional)"*: `lajolla`
>    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator', 'central_longitude': 12.0}`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
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