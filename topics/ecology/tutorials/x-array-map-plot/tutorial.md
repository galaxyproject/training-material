---
layout: tutorial_hands_on
title: Visualization of Climate Data using  NetCDF xarray Map Plotting
zenodo_link: "https://doi.org/10.5281/zenodo.6621460"
requirements:
  -
    type: "internal"
    topic_name: climate
    tutorials:
        - climate-101
questions:
- What is xarray map plotting?
- How to plot NetCDF data using xarray map plotting tool?
- What are the different types of projections and colormap options available in the tool?
objectives:
- Learn about plotting Climate Data
- Learn about NetCDF xarray map plotting tool
- Learn about visualizing the climate data using NetCDF xarray map plotting using the different kinds of projections and colormap options
time_estimation: 1H
key_points:
- NetCDF xarray climate data visualization
- NetCDF xarray map plotting projections
- NetCDF xarray map plotting colormaps
tags:
  - pangeo
contributors:
- Quickbeasts51429

---



><comment-title></comment-title>
>
> The tutorial aims for establishing a good knowledge about meaningful visualization of climate data.
>
{: .comment}

 It is beginner friendly and does not require much knowledge about the tool.

> <agenda-title></agenda-title>
>
> In this tutorial we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

><comment-title>Background</comment-title>
>
>According to [UN](https://www.un.org/en/climatechange/what-is-climate-change) , Climate is the long term shift in temperature and weather patterns which may be due to natural or artificial causes. To learn more about climate, refer this [tutorial]({% link topics/climate/tutorials/climate-101/tutorial.md %}) from the GTN. Due to the frequently changing nature of the weather patterns, the size of the collected data is huge.
The climate data is mainly represented in these three categories : [NetCDF](https://en.wikipedia.org/wiki/NetCDF) (Network Common Data Form), [HDF](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) (Hierarchical Data Format) , [GRIB](https://en.wikipedia.org/wiki/GRIB) (GRIdded Binary or General Regularly-distributed Information in Binary form).
>
>The NetCDF file format is basically used for storing multidimensional data which generally consists of variables such as temperature, precipitation, wind direction, etc. The variation of climate variables over a period of time is suitably plotted using this dataset. The entire earth is divided into both horizontal as well as vertical coordinates which makes plotting of the variables such as the ocean temperatures possible.
>
>The coordinate system, types of projections and colormaps are some of the very important considerations in achieving the most suitable visualization option.
{: .comment}


# Introduction

The data used in this tutorial is [ECMWF Reanalysis 5 hourly data on single levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview). We are interested in the following variables: air temperature at 2 metres, latitude, longitude and time. Our main objective is to plot the global air temperature at 2 metres with respect to time. For this we will be using the netCDF xarray tool available in the Galaxy Europe (or your favourite Galaxy Instance) server.

It will be a fun learning experience for anyone who loves visualization !

# Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial. Name it as per your choice. My suggestions : *ECMWF_Reanalysis*.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. **Import remote files**
>    - location: `ECMWF ERA5 Reanalysis` ->`2022` -> `05` -> `data` -> `air_temperature_at_2_metres.nc`
>
>    {% snippet faqs/galaxy/datasets_import_from_remote_files.md location="ECMWF ERA4 Reanalysis -> 2022 -> 05 -> data -> air_temperature_at_2_metres.nc" %}
>
>    **or**, if this is not available, you can import from [Zenodo]({{page.zenodo_link}}):
>
>     ```
>     https://zenodo.org/record/6621460/files/air_temperature_at_2_metres.netcdf
>     ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>
> 3. **Check the datatype** of uploaded data, it should be `netCDF`.
>    -  If it is not `netCDF`, make sure to convert it using the Galaxy built-in format converters.
>
>    {% snippet faqs/galaxy/datasets_convert_datatype.md conversion="Convert h5 to netCDF" %}
>
> 4. **Rename Datasets** {% icon galaxy-pencil %} if you find the name to be too long or you have something more meaningful in your mind.
>
> 5. **Tag Datasets** {% icon galaxy-tags %} with the tag `#ecmwf`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# What is xarray?

[xarray](https://docs.xarray.dev/en/stable/), formerly known as xray, is a python package which enables us to play with gridded data. This package shares most of its features from [numpy](https://numpy.org/), but in a more convenient manner by keeping track of labels in arrays. The gridded data is mainly available in netCDF data format. Thus **xarray comes very handy while dealing with netCDF files.**


## What is Visualisation in xarray?

xarray uses **[Cartopy](https://scitools.org.uk/cartopy/docs/latest/)** and **[Matplotlib](https://matplotlib.org/)** as the two main tools for creating detailed and informative plots.
Cartopy is a python package used for geospatial data analysis. In the Python library, matplotlib is the most used scientific plotting library.
For a multidimensional data consisting of latitudes and longitudes along with the other variables, xarray has the capability of appling cartopy map projections.

# Splitting the dataset  using  seltimestep, splithour and plotting

After loading the required data, the following stage is to obtain the meta info or meta data of the file. The very purpose of these steps are to obtain the information about dimensions, variables, global attributes, etc. The coordinate info helps to know about the actual data entries present under the various variables.

Follow the below steps:

## **NetCDF xarray Metadata Info**

> <hands-on-title>netCDF dataset with xarray metadata Galaxy Tool</hands-on-title>
>
> 1. In the tools search for {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %} with the following parameter:
>    - {% icon param-file %} *"Netcdf file"*: `air_temperature_at_2_metres.nc`
>  and click on `Execute`.
>
> 2. View {% icon galaxy-eye%} the two generated outputs:
>    - `Metadata infos` is a `tabular` providing the list of variables, their dimension names and number of elements per dimension. This file is used by other xarray Tools.
>    - The second file `info file` provides a summary of the **xarray Dataset** contained in your netCDF file.
{: .hands_on}




> <question-title></question-title>
>
> 1. What is the name of the variable for air temperature at 2 metres. What are its physical units?
>
>
> > <solution-title></solution-title>
> > 1. Information about variable names and units can be found in **info file** that was generated by xarray metadata Galaxy Tool.
> >      - Variable name: `air_temperature_at_2_metres`
> >      - Units: `K`
> >
> >
> > > <code-out-title></code-out-title>
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

> <hands-on-title> Get Coordinate information with xarray Coordinate</hands-on-title>
>
> 1. {% tool [NetCDF xarray Coordinate Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_coords_info/xarray_coords_info/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `air_temperature_at_2_metres.nc`.
>2. View {% icon galaxy-eye%} the 4 generated outputs:
>    - `lat`: a tabular file containing all the latitude values of our Xarray dataset;
>    - `lon`: a tabular file containing all the longitudes values;
>    - `time0`: this tabular file contains all the times extracted from the netCDF file. In our case, these are relative to 01/05/2022, 00:00:00 UTC;
>    - `version`: this is a text file returning the xarray package version. It is useful when publishing your Galaxy workflow.
>
{: .hands_on}

> <comment-title></comment-title>
>
>The number of tabular files returned by this programme is proportional to the number of coordinate variables in your input file. No further processing is done after decoding the values from the netCDF input file. As a result, depending on how the original input file was coded, units for latitudes, longitudes, (level, if this dimension was present, which is not the case in this particular case) and time may differ from one file to the next.
>
{: .comment}




> <question-title>Understanding air temperature coordinates at 2 metres</question-title>
>
>
> 1. What is the format of time coordinate?
> 2. What is the range of values for latitude and longitude?
>
> > <solution-title></solution-title>
> > 1. The `info file` tells us that `time0` is coded as `timedelta64[ns]` i.e. the time differences are in ns (here in nanoseconds). The format of time is `yy-mm-dd  hh:mm:ss`. If we look at the tabular file named `time0` (generated by `NetCDF xarray Coordinate Info`), we see that these times are automatically converted to human readable time format when printed:
> >
> > > <code-out-title></code-out-title>
> > > ```
> > > 0  2022-05-01 00:00:00
> > > 1  2022-05-01 01:00:00
> > > 2  2022-05-01 02:00:00
> > > 3  2022-05-01 03:00:00
> > > 4  2022-05-01 04:00:00
> > > ```
> > >
> > {: .code-out}
> >
> > This tells us that we have hourly forecast data.
> >
> > 2.If we look at the tabular file named `lat` and `lon` (generated by `NetCDF xarray Coordinate Info`), we see that these files are displayed with two columns e.g. index numbers and corresponding values. From these, we infer that the range of values for `lat` is  `90.0`  to `-90.0`and for `lon` is `0` to `359.75`.
> >
> {: .solution }
{: .question }


## Operations on Climate data using **CDO Operations**

We have hourly data. In order to plot it, we must first extract the hours from the bigger dataset. This is done using the `seltimestep` and the `splithour` options available in the `CDO Operations` tool. `splithour` is used when we wish to plot more than an hour. Our main aim is to plot the last hour data on the last day present in the dataset.

> <hands-on-title>Defining a particular time range using seltimestep </hands-on-title>
>
> 1. {% tool [CDO Operations](toolshed.g2.bx.psu.edu/repos/climate/cdo_operations/cdo_operations/2.0.0+galaxy0) %} with the following parameters:
>    - In *"CDO Operators"*:
>        - {% icon param-repeat %} *"Insert CDO Operators"*
>            - *"Select cdo operator"*: `seltimestep (Select timesteps)`
>            - *"Timesteps for selection"*: `744/744`
>            - {% icon param-file %} *"Additional input file"*: `air_temperatures_at_2_metres.nc`
>
>
>
> > <comment-title></comment-title>
> >
> > The syntax of using the `seltimestep` is `(initial data number / final data entry)`. An important thing to pay attention to is how the data entries are numbered: are they numbered starting from 0 or 1. Accordingly we can add or skip adding 1 to the data number to achieve the desired result.
> >
Although we are not using `splithour` here, you can find below the syntax for future uses.
> >
> >1. {% tool [CDO Operations](toolshed.g2.bx.psu.edu/repos/climate/cdo_operations/cdo_operations/2.0.0+galaxy0) %} with the following parameters:
> >- In *"CDO Operators"*:
> >- {% icon param-repeat %} *"Insert CDO Operators"*
> >- *"Select cdo operator"*: `splithour (Split hours)`
> >- {% icon param-file %} *"Additional input file"*: `outfile.netcdf` generated from the previous step.
> >
> > This step generates that `N` number of `outfiles.netcdf` files where `N` is the range of selection.
> > Suppose your selected range was `744/744` for the `seltimestep` , then it will generate `2` files which can be plotted further.
> >
> {: .comment}
>
{: .hands_on}






## Finding the  **NetCDF xarray Metadata Info**

> <hands-on-title>netCDF dataset with xarray metadata Galaxy Tool for the hourly plots</hands-on-title>
>
> 1. {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `outfile.netcdf`. .
>
> 2. View {% icon galaxy-eye%} the two generated outputs:
>    - `Metadata infos` is a `tabular` providing the list of variables, their dimension names and number of elements per dimension. This file is used by other xarray Tools.
>    - The second file `info file` provide a summary of the **xarray Dataset** contained in your netCDF file.
{: .hands_on}

# Map Plotting using **NetCDF xarray map plotting**

> <hands-on-title>Plotting the data of the last hour of the day</hands-on-title>
> The air temperatures corresponding to the 744th time step from the original netCDF file, namely `23:00:00` for 31st May 2022 is plotted here :
>
> 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.20.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `outfile.netcdf`
>    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from outfile.netcdf` (output of **NetCDF xarray Metadata Info** >{% icon tool %})
>    - *"Choose the variable to plot"*: `air_temperature_at_2_metres`
>    - *"Name of latitude coordinate"*: `lat`
>    - *"Name of longitude coordinate"*: `lon`
>    - *"Datetime selection"*: `No`
>    - *"Range of values for plotting e.g. minimum value and maximum value (minval,maxval) (optional)"*: `220,320`
>    - *"Add country borders with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add coastline with alpha value [0-1] (optional)"*: `1.0`
>    - *"Add ocean with alpha value [0-1] (optional)"*: `1.0`
>    - *"Specify plot title (optional)"*: `Air_temperature_at_2_metres___23:00:00_UTC`
>    - *"Specify which colormap to use for plotting (optional)"*: `lajolla`
>    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator'}`
>
> ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 UTC](../../images/x-array-map-plot/plot1.png)
>
>
{: .hands_on}

> <details-title>Syntax for Projections</details-title>
> **Xarray Map Plotting** comes with a variety of projections. There are basically four types of projections available namely **Azimuthal, Conic, Cylindrical and Pseudo-cylindrical.**
>
>The **Azimuthal projection** depicts the surface of the earth with the help of a flat plane. It is also known as '**Plane Projection**.' This simple projection type forms a family of projections by considering the poles in "normal aspect." Another form of this projection is **Azimuthal Equidistant Projection**. It preserves both directions and distance from the central point of the earth. The **stereographic projection** which is a type of azimuthal projection was created before 150AD and deployed while mapping areas over the poles.
> The **conic type** is based on the concept of projecting the surface of the Earth on a conical surface which is unrolled into a plane  surfaced map. Some of the examples of conic type of projections are : **Albers Equal Area Conic, Equidistant Conic, Lambert Conformal Conic, and Polyconic**. This type of projection basically has a big implementation in aviation navigation.
>
>
> The **cylindrical type** is based on the concept of plotting the geographical features on the surface of the cylinder and then it is unrolled to present as a **flat projection**. Some examples of cylindrical projections are : **Cylindrical Equal Area, Behrmann Cylindrical Equal-Area , Stereographic Cylindrical, Peters, Mercator, and Transverse Mercator**. **Mercator** is one of the most famous and preferred type of projection. The straight lines (Rhumb Lines) are best for the navigation process.
>
> In the **pseudo cylindrical** type of projection , the meridians are straight instead of being curved. These kinds of projections are characterised by straight horizontal lines for parallels of latitude and (usually) equally-spaced curved meridians of longitude. **Oval projections** pinched or flattened at the poles are its identifying features.
>
>
> **Projections**:
>
> There are about **30** different kinds of projections available in the proj option. The basic knowledge of projections is very essential in plotting out the best fit map. Below is the list of plots deployed using the various kinds of projections supported by the NetCDF xarray mapplotting tool.
>
>
> **PlateCarree**
>**Syntax: {"proj":"PlateCarree"}**
>![PlateCarree](../../images/x-array-map-plot/projections/PlateCarree.png)
>
> **EquidistantConic**
>**Syntax: {"proj":"EquidistantConic", "central_longitude": 20.0, "central_latitude": 70.0 }**
>![EquidistantConic](../../images/x-array-map-plot/projections/EquidistantConic.png)
>
> **AlbersEqualArea**
> **Syntax: {"proj":"AlbersEqualArea", "central_longitude": 20.0, "central_latitude": 70.0 }**
> ![AlbersEqualArea](../../images/x-array-map-plot/projections/AlbersEqualArea.png)
>
> **EuroPP**
> **Syntax: {"proj":"EuroPP"}**
> ![EuroPP](../../images/x-array-map-plot/projections/EuroPP.png)
>
> **LambertConformal**
> **Syntax:{"proj":"LambertConformal"  }**
> ![LambertConformal](../../images/x-array-map-plot/projections/LambertConformal.png)
>
> **AzimuthalEquidistant**
> **Syntax: {"proj":"AzimuthalEquidistant" }**
> ![AzimuthalEquidistant](../../images/x-array-map-plot/projections/AzimuthalEquidistant.png)
>
> **LambertCylindrical**
> **Syntax:{"proj":"LambertCylindrical"}**
> ![LambertCylindrical](../../images/x-array-map-plot/projections/LambertCylindrical.png)
>
> **Mercator**
> **Syntax: {"proj":"Mercator" }**
> ![Mercator](../../images/x-array-map-plot/projections/Mercator.png)
>
> **Miller**
> **Syntax: {"proj":"Miller"  }**
> ![Miller](../../images/x-array-map-plot/projections/Miller.png)
>
> **Mollweide**
> **Syntax: {"proj":"Mollweide" }**
> ![Mollweide](../../images/x-array-map-plot/projections/Mollweide.png)
>
> **Orthographic**
> **Syntax:{"proj":"Orthographic"  }**
> ![Orthographic](../../images/x-array-map-plot/projections/Orthographic.png)
>
> **Robinson**
> **Syntax: {"proj":"Robinson"  }**
> ![Robinson](../../images/x-array-map-plot/projections/Robinson.png)
>
> **Sinusoidal**
> **Syntax: {"proj":"Sinusoidal"  }**
> ![Sinusoidal](../../images/x-array-map-plot/projections/sinusoidal.png)
>
> **Stereographic**
> **Syntax: {"proj":"Stereographic"  }**
> ![Stereographic](../../images/x-array-map-plot/projections/Stereographic.png)
>
> **TransverseMercator**
> **Syntax:{"proj":"TransverseMercator"  }**
> ![TransverseMercator](../../images/x-array-map-plot/projections/TransverseMercator.png)
>
> **InterruptedGoodeHomolosine**
> **Syntax: {"proj":"InterruptedGoodeHomolosine"  }**
> ![InterruptedGoodeHomolosine](../../images/x-array-map-plot/projections/InterruptedGoodeHomolosine.png)
>
> **RotatedPole**
> **Syntax: {"proj":"RotatedPole"  }**
> ![RotatedPole](../../images/x-array-map-plot/projections/RotatedPole.png)
>
> **OSGB**
> **Syntax:{"proj":"OSGB"  }**
> ![OSGB](../../images/x-array-map-plot/projections/OSGB.png)
>
> **Geostationary**
> **Syntax: {"proj":"Geostationary"  }**
> ![Geostationary](../../images/x-array-map-plot/projections/Geostationary.png)
>
> **NearsidePerspective**
> **Syntax:{"proj":"NearsidePerspective" }**
> ![NearsidePerspective](../../images/x-array-map-plot/projections/NearsidePerspective.png)
>
> **EckertI**
> **Syntax: {"proj":"EckertI" }**
> ![EckertI](../../images/x-array-map-plot/projections/EckertI.png)
>
> **EckertII**
> **Syntax: {"proj":"EckertII" }**
> ![EckertII](../../images/x-array-map-plot/projections/EckertII.png)
>
> **EckertIII**
> **Syntax: {"proj":"EckertIII" }**
> ![EckertIII](../../images/x-array-map-plot/projections/EckertIII.png)
>
> **EckertIV**
> **Syntax: {"proj":"EckertIV" }**
> ![EckertIV](../../images/x-array-map-plot/projections/EckertIV.png)
>
> **EckertV**
> **Syntax: {"proj":"EckertV" }**
> ![EckertV](../../images/x-array-map-plot/projections/EckertV.png)
>
> **EckertVI**
> **Syntax: {"proj":"EckertVI" }**
> ![EckertVI](../../images/x-array-map-plot/projections/EckertVI.png)
>
> **EqualEarth**
> **Syntax:{"proj":"EqualEarth" }**
> ![EqualEarth](../../images/x-array-map-plot/projections/EqualEarth.png)
>
> **Gnomonic**
> **Syntax: {"proj":"Gnomonic" }**
> ![Gnomonic](../../images/x-array-map-plot/projections/Gnomonic.png)
>
> **LambertAzimuthalEqualArea**
> **Syntax: {"proj":"LambertAzimuthalEqualArea" }**
> ![LambertAzimuthalEqualArea](../../images/x-array-map-plot/projections/LambertAzimuthalEqualArea.png)
>
> **NorthPolarStereo**
> **Syntax: {"proj":"NorthPolarStereo" }**
> ![NorthPolarStereo](../../images/x-array-map-plot/projections/NorthPolarStereo.png)
>
> **OSNI**
> **Syntax: {"proj":"OSNI" }**
> ![OSNI](../../images/x-array-map-plot/projections/OSNI.png)
>
> **SouthPolarStereo**
> **Syntax: {"proj":"SouthPolarStereo" }**
> ![SouthPolarStereo](../../images/x-array-map-plot/projections/SouthPolarStereo.png)
>
{: .details}


> <question-title></question-title>
>
> What are the different kinds of projections that can be used?
>
> > <solution-title></solution-title>
> >
> > There are many projections which can be used in the `NetCDF xarray map plotting` tool. Different projections have different purposes and need to be carefully chosen.
> >
> > Plotting different projections using **NetCDF xarray map plotting**:
> > > <hands-on-title>Plotting the major projections using the xarray tool.</hands-on-title>
> > >
> > > 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.20.2+galaxy0) %} with the following parameters:
> > >    - {% icon param-file %} *"Netcdf file"*: `outfile.netcdf`
> > >    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from outfile.netcdf` (output of **NetCDF xarray Metadata Info** {% icon tool %})
> > >    - *"Choose the variable to plot"*: `air_temperature_at_2_metres`
> > >    - *"Name of latitude coordinate"*: `lat`
> > >    - *"Name of longitude coordinate"*: `lon`
> > >    - *"Datetime selection"*: `No`
> > >    - *"Range of values for plotting e.g. minimum value and maximum value (minval,maxval) (optional)"*: `220,320`
> > >    - *"Add country borders with alpha value [0-1] (optional)"*: `1.0`
> > >    - *"Add coastline with alpha value [0-1] (optional)"*: `1.0`
> > >    - *"Add ocean with alpha value [0-1] (optional)"*: `1.0`
> > >    - *"Specify plot title (optional)"*: `Air_temperature_at_2_metres___23:00:00_UTC_InterruptedGoodeHomolosine`
> > >    - *"Specify which colormap to use for plotting (optional)"*: `lajolla`
> > >    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{"proj":"InterruptedGoodeHomolosine" }`
> > >
> > >
> > >
> > > The final plot is shown below:
> > >
> > > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 : Interrupted_GoodeHomolosine](../../images/x-array-map-plot/Interrupted_GoodeHomolosine.png)
> > >
> > >
> > >
> > > Some other potentially interesting types of projections can be found below :
> > >
> > > *{"proj":"LambertCylindrical" }*
> > > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 : LambertCylindrical](../../images/x-array-map-plot/Lam.png)
> > > *{"proj":"Orthographic"  }*
> > > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 : Orthographic](../../images/x-array-map-plot/Ortho.png)
> > > *{"proj":"Sinusoidal" }*
> > > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 : Sinusoidal](../../images/x-array-map-plot/sinu.png)
> > > *{"proj":"EquidistantConic"}*
> > > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 : EquidistantConic](../../images/x-array-map-plot/Equidi.png)
> > > *{"proj":"LambertConformal" }*
> > > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 : LambertConformal](../../images/x-array-map-plot/Lamcon.png)
> > > *{"proj":"AzimuthalEquidistant" }*
> > > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 at 23:00:00 : AzimuthalEquidistant](../../images/x-array-map-plot/Azi.png)
> > {: .hands_on}
> {: .solution}
{: .question}


> <details-title>Available ColorMaps</details-title>
> There are about **150 different kinds of maps available in the  colormaps** options. This document contains a list of all the types of projections and colormaps available to the users while using the Xarray Map Plotting. These depict the global air temperature at 2 metres on a Kelvin Scale of range 220 K TO 320 K. Let's get started!
> Below is the complete list for reference.
>
> **ACCENT**
> ![Accent](../../images/x-array-map-plot/colormaps/Accent.png)
>  **BLUES**
> ![Blues]( ../../images/x-array-map-plot/colormaps/Blues.png )
>  **BrBg**
> ![BrBg]( ../../images/x-array-map-plot/colormaps/BrBG.png )
>  **BuGn**
> ![BuGn]( ../../images/x-array-map-plot/colormaps/BuGn.png )
>  **BuPu**
> ![BuPu](../../images/x-array-map-plot/colormaps/BuPu.png  )
>  **CMRmap**
> ![CMRmap]( ../../images/x-array-map-plot/colormaps/CMRmap.png )
>  **Dark2**
> ![Dark2](../../images/x-array-map-plot/colormaps/Dark2.png  )
>  **Acton**
> ![Acton](../../images/x-array-map-plot/colormaps/acton.png  )
>  **Acton_r**
> ![Acton_r]( ../../images/x-array-map-plot/colormaps/acton_r.png )
>  **Afmhot**
> ![Afmhot]( ../../images/x-array-map-plot/colormaps/afmhot.png )
>  **Autumn**
> ![Autumn]( ../../images/x-array-map-plot/colormaps/autumn.png )
>  **Bam**
> ![BamO]( ../../images/x-array-map-plot/colormaps/bam.png  )
>  **BamO**
> ![BamO]( ../../images/x-array-map-plot/colormaps/bam0.png )
>  **BamO_r**
> ![bamO_r]( ../../images/x-array-map-plot/colormaps/bamO_r.png )
>  **Bamako**
> ![Accent]( ../../images/x-array-map-plot/colormaps/bamako.png )
>  **Bamako_r**
> ![bamako_r]( ../../images/x-array-map-plot/colormaps/bamako_r.png )
>  **Batlow**
> ![Batlow]( ../../images/x-array-map-plot/colormaps/batlow.png  )
>  **BatlowK**
> ![BatlowK]( ../../images/x-array-map-plot/colormaps/batlowK.png )
>  **BatlowK_r**
> ![BatlowK_r]( ../../images/x-array-map-plot/colormaps/batlowK_r.png )
>  **BatlowW**
> ![batlowW]( ../../images/x-array-map-plot/colormaps/batlowW.png )
>  **BatlowW_r**
> ![batlowW_r]( ../../images/x-array-map-plot/colormaps/batlowW_r.png )
>  **Batlow_r**
> ![Batlow_r]( ../../images/x-array-map-plot/colormaps/batlow_r.png )
>  **Berlin**
> ![Berlin]( ../../images/x-array-map-plot/colormaps/berlin.png )
>  **berlin_r**
> ![Berlin_r]( ../../images/x-array-map-plot/colormaps/berlin_r.png )
>  **Bilbao**
> ![bilbao](  ../../images/x-array-map-plot/colormaps/bilbao.png)
>  **Bilbao_r**
> ![Bilbao_r]( ../../images/x-array-map-plot/colormaps/bilbao_r.png )
>  **Binary**
> ![Binary]( ../../images/x-array-map-plot/colormaps/binary.png )
>  **Bone**
> ![Bone]( ../../images/x-array-map-plot/colormaps/bone.png )
>  **Brg**
> ![Brg]( ../../images/x-array-map-plot/colormaps/brg.png )
>  **Broc**
> ![Broc]( ../../images/x-array-map-plot/colormaps/broc.png )
>  **BrocO**
> ![BrocO]( ../../images/x-array-map-plot/colormaps/brocO.png )
>  **BrocO_r**
> ![BrocO_r]( ../../images/x-array-map-plot/colormaps/brocO_r.png )
>  **Broc_r**
> ![Broc_r]( ../../images/x-array-map-plot/colormaps/broc_r.png )
>  **Buda**
> ![Buda]( ../../images/x-array-map-plot/colormaps/buda.png )
>  **Buda_r**
> ![Buda_r]( ../../images/x-array-map-plot/colormaps/buda_r.png )
>  **Bukavu**
> ![Bukavu]( ../../images/x-array-map-plot/colormaps/bukavu.png )
>  **Bukavu_r**
> ![Bukavu_r]( ../../images/x-array-map-plot/colormaps/bukavu_r.png )
>  **Bwr**
> ![Bwr]( ../../images/x-array-map-plot/colormaps/bwr.png )
>  **Cool**
> ![Cool]( ../../images/x-array-map-plot/colormaps/cool.png )
>  **Coolwarm**
> ![Coolwarm]( ../../images/x-array-map-plot/colormaps/coolwarm.png )
>  **Copper**
> ![Copper]( ../../images/x-array-map-plot/colormaps/copper.png )
>  **Cork**
> ![Cork]( ../../images/x-array-map-plot/colormaps/cork.png )
>  **CorkO**
> ![CorkO]( ../../images/x-array-map-plot/colormaps/corkO.png )
>  **CorkO_r**
> ![CorkO_r]( ../../images/x-array-map-plot/colormaps/corkO_r.png )
>  **Cork_r**
> ![Cork_r]( ../../images/x-array-map-plot/colormaps/cork_r.png )
>  **Cubehelix**
> ![Cubehelix]( ../../images/x-array-map-plot/colormaps/cubehelix.png )
>  **Davos**
> ![Davos]( ../../images/x-array-map-plot/colormaps/davos.png )
>  **GnBu**
> ![GnBu]( ../../images/x-array-map-plot/colormaps/GnBu.png )
>  **Greens**
> ![Greens]( ../../images/x-array-map-plot/colormaps/Greens.png )
>  **Greys**
> ![Greys]( ../../images/x-array-map-plot/colormaps/Greys.png )
>  **Imola**
> ![Imola]( ../../images/x-array-map-plot/colormaps/Imola.png )
>  **OrRd**
> ![OrRd]( ../../images/x-array-map-plot/colormaps/OrRd.png )
>  **Oranges**
> ![Oranges]( ../../images/x-array-map-plot/colormaps/Oranges.png )
>  **Paired**
> ![Paired]( ../../images/x-array-map-plot/colormaps/Paired.png )
>  **Pastel1**
> ![Pastel1]( ../../images/x-array-map-plot/colormaps/Pastel1.png )
>  **Pastel2**
> ![Pastel2]( ../../images/x-array-map-plot/colormaps/Pastel2.png )
>  **Davos_r**
> ![Davos_r]( ../../images/x-array-map-plot/colormaps/davos_r.png )
>  **Devon_r**
> ![Devon_r]( ../../images/x-array-map-plot/colormaps/devon_r.png )
>  **Devon**
> ![Devon]( ../../images/x-array-map-plot/colormaps/devon.png )
>  **Fes**
> ![Fes]( ../../images/x-array-map-plot/colormaps/fes.png )
>  **Fes_r**
> ![Fes_r]( ../../images/x-array-map-plot/colormaps/fes_r.png )
>  **Flag**
> ![Flag]( ../../images/x-array-map-plot/colormaps/flag.png )
>  **Gist_earth**
> ![Gist_earth]( ../../images/x-array-map-plot/colormaps/gist_earth.png )
>  **Gist_gray**
> ![Gist_gray]( ../../images/x-array-map-plot/colormaps/gist_gray.png )
>  **Gist_heat**
> ![Gist_heat]( ../../images/x-array-map-plot/colormaps/gist_heat.png )
>  **Gist_ncar**
> ![Gist_ncar]( ../../images/x-array-map-plot/colormaps/gist_ncar.png )
>  **Gist_rainbow**
> ![Gist_rainbow]( ../../images/x-array-map-plot/colormaps/gist_rainbow.png)
>  **Gist_stern**
> ![Gist_stern]( ../../images/x-array-map-plot/colormaps/gist_stern.png)
>  **Gist_yarg**
> ![Gist_yarg]( ../../images/x-array-map-plot/colormaps/gist_yarg.png)
>  **Gnuplot**
> ![Gnuplot]( ../../images/x-array-map-plot/colormaps/gnuplot.png)
>  **Gnuplot2**
> ![Gnuplot2]( ../../images/x-array-map-plot/colormaps/gnuplot2.png)
>  **Gray**
> ![Gray]( ../../images/x-array-map-plot/colormaps/gray.png)
>  **GrayC**
> ![GrayC](../../images/x-array-map-plot/colormaps/grayC.png )
>  **GrayC_r**
> ![GrayC_r]( ../../images/x-array-map-plot/colormaps/grayC_r.png)
>  **Hawaii**
> ![Hawaii]( ../../images/x-array-map-plot/colormaps/hawaii.png)
>  **Hawaii_r**
> ![Hawaii_r]( ../../images/x-array-map-plot/colormaps/hawaii_r.png)
>  **Hot**
> ![Hot]( ../../images/x-array-map-plot/colormaps/hot.png )
>  **Hsv**
> ![Hsv]( ../../images/x-array-map-plot/colormaps/hsv.png )
>  **Imola_r**
> ![Imola_r]( ../../images/x-array-map-plot/colormaps/imola_r.png )
>  **Jet**
> ![Jet]( ../../images/x-array-map-plot/colormaps/jet.png )
>  **Lajolla**
> ![Lajolla]( ../../images/x-array-map-plot/colormaps/lajolla.png )
>  **Lajolla_r**
> ![Lajolla_r]( ../../images/x-array-map-plot/colormaps/lajolla_r.png )
>  **Lapaz**
> ![Lapaz]( ../../images/x-array-map-plot/colormaps/lapaz.png )
>  **Lapaz_r**
> ![Lapaz_r]( ../../images/x-array-map-plot/colormaps/lapaz_r.png )
>  **Lisbon**
> ![Lisbon](  ../../images/x-array-map-plot/colormaps/lisbon.png)
>  **Lisbon_r**
> ![Lisbon_r]( ../../images/x-array-map-plot/colormaps/lisbon_r.png )
>  **Nipy_spectral**
> ![Nipy_spectral]( ../../images/x-array-map-plot/colormaps/nipy_spectral.png )
>  **Nuuk**
> ![Nuuk]( ../../images/x-array-map-plot/colormaps/nuuk.png )
>  **Nuuk_r**
> ![Nuuk_r]( ../../images/x-array-map-plot/colormaps/nuuk_r.png )
>  **Ocean**
> ![Ocean]( ../../images/x-array-map-plot/colormaps/ocean.png )
>  **Oleron**
> ![Oleron]( ../../images/x-array-map-plot/colormaps/oleron.png )
>  **Oslo**
> ![Oslo]( ../../images/x-array-map-plot/colormaps/oslo.png )
>  **Oslo_r**
> ![Oslo_r]( ../../images/x-array-map-plot/colormaps/oslo_r.png )
>  **Pink**
> ![Pink]( ../../images/x-array-map-plot/colormaps/pink.png )
>  **PRGn**
> ![PRGn]( ../../images/x-array-map-plot/colormaps/PRGn.png )
>  **PiYG**
> ![PiYG]( ../../images/x-array-map-plot/colormaps/PiYG.png )
>  **PuBu**
> ![PuBu]( ../../images/x-array-map-plot/colormaps/PuBu.png )
>  **PuBuGn**
> ![PuBuGn](  ../../images/x-array-map-plot/colormaps/PuBuGn.png)
>  **PuOr**
> ![PuOr]( ../../images/x-array-map-plot/colormaps/PuOr.png )
>  **PuRd**
> ![PuRd]( ../../images/x-array-map-plot/colormaps/PuRd.png )
>  **Purples**
> ![Purples]( ../../images/x-array-map-plot/colormaps/Purples.png )
>  **RdBu**
> ![RdBu]( ../../images/x-array-map-plot/colormaps/RdBu.png )
>  **RdBu_r**
> ![RdBu_r]( ../../images/x-array-map-plot/colormaps/RdBu_r.png )
>  **RdBu_r**
> ![RdBu_r]( ../../images/x-array-map-plot/colormaps/RdBu_r.png)
>  **RdGy_r**
> ![RdGy_r]( ../../images/x-array-map-plot/colormaps/RdGy_r.png )
>  **RdPu**
> ![RdPu]( ../../images/x-array-map-plot/colormaps/RdPu.png)
> **RdPu_r**
> ![RdPu_r](../../images/x-array-map-plot/colormaps/RdPu_r.png )
> **RdYIBu**
> ![RdYIBu]( ../../images/x-array-map-plot/colormaps/RdYIBu.png)
> **RdYIGn**
> ![RdYIGn]( ../../images/x-array-map-plot/colormaps/RdYIGn.png)
> **Reds**
> ![Reds](../../images/x-array-map-plot/colormaps/Reds.png )
> **Set1**
> ![Set1]( ../../images/x-array-map-plot/colormaps/Set1.png)
> **Set2**
> ![Set2]( ../../images/x-array-map-plot/colormaps/Set2.png )
> **Set3**
> ![Set3]( ../../images/x-array-map-plot/colormaps/Set3.png )
> **Spectral**
> ![Spectral]( ../../images/x-array-map-plot/colormaps/Spectral.png )
> **Wistia**
> ![Wistia]( ../../images/x-array-map-plot/colormaps/Wistia.png )
> **YIGn**
> ![YIGn]( ../../images/x-array-map-plot/colormaps/YIGn.png )
> **YIGnBu**
> ![YIGnBu]( ../../images/x-array-map-plot/colormaps/YIGnBu.png )
> **YIOrBr**
> ![YIOrBr]( ../../images/x-array-map-plot/colormaps/YIOrBr.png )
> **YIOrRd**
> ![YIOrRd]( ../../images/x-array-map-plot/colormaps/YIOrRd.png )
> **Prism**
> ![Prism]( ../../images/x-array-map-plot/colormaps/prism.png )
> **Rainbow**
> ![Rainbow]( ../../images/x-array-map-plot/colormaps/rainbow.png )
> **Roma**
> ![Roma]( ../../images/x-array-map-plot/colormaps/roma.png )
> **RomaO**
> ![RomaO]( ../../images/x-array-map-plot/colormaps/romaO.png )
> **RomaO_r**
> ![RomaO_r]( ../../images/x-array-map-plot/colormaps/romaO_r.png )
> **Roma_r**
> ![Roma_r]( ../../images/x-array-map-plot/colormaps/roma_r.png )
> **Seismic**
> ![Seismic]( ../../images/x-array-map-plot/colormaps/seismic.png )
> **Spring**
> ![Spring]( ../../images/x-array-map-plot/colormaps/spring.png  )
> **Summer**
> ![Summer]( ../../images/x-array-map-plot/colormaps/summer.png )
> **Tab10**
> ![Tab10]( ../../images/x-array-map-plot/colormaps/tab10.png )
> **Tab20**
> ![Tab20]( ../../images/x-array-map-plot/colormaps/tab20.png )
> **Tab20b**
> ![Tab20b]( ../../images/x-array-map-plot/colormaps/tab20b.png )
> **Tab20c**
> ![Tab20c]( ../../images/x-array-map-plot/colormaps/tab20c.png )
> **Terrain**
> ![Terrain]( ../../images/x-array-map-plot/colormaps/terrain.png )
> **Tofino**
> ![Tofino]( ../../images/x-array-map-plot/colormaps/tofino.png)
> **Tofino_r**
> ![Tofino_r]( ../../images/x-array-map-plot/colormaps/tofino_r.png )
> **Tokyo**
> ![Tokyo]( ../../images/x-array-map-plot/colormaps/tokyo.png )
> **Tokyo_r**
> ![Tokyo_r]( ../../images/x-array-map-plot/colormaps/tokyo_r.png )
> **Turku**
> ![Turku]( ../../images/x-array-map-plot/colormaps/turku.png )
> **Turku_r**
> ![Turku_r]( ../../images/x-array-map-plot/colormaps/turku_r.png )
> **Vanimo**
> ![Vanimo]( ../../images/x-array-map-plot/colormaps/vanimo.png )
> **Vanimo_r**
> ![Vanimo_r]( ../../images/x-array-map-plot/colormaps/vanimo_r.png )
> **Vik**
> ![Vik]( ../../images/x-array-map-plot/colormaps/vik.png )
> **VikO**
> ![VikO]( ../../images/x-array-map-plot/colormaps/vikO.png )
> **VikO_r**
> ![VikO_r]( ../../images/x-array-map-plot/colormaps/vikO_r.png )
> **Vik_r**
> ![Vik_r]( ../../images/x-array-map-plot/colormaps/vik_r.png )
> **Winter**
> ![Winter]( ../../images/x-array-map-plot/colormaps/winter.png )
>
{: .details}

> <question-title></question-title>
>
> 1. What are the different kinds of colormaps that can be used?
>
> > <solution-title></solution-title>
> >
> >When it comes to conveying the correct information, through visualisation, colors play a major role. Suppose you want to display a cold region, its an obvious practice of using cooler tones such as a `blue`. Thus it is important to understand the choices we have.
> >
> >
> >## Plotting different colormaps using **NetCDF xarray map plotting**
 > > > <hands-on-title>Plotting the major colormaps using the xarray tool.</hands-on-title>
> > >
> > > 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.20.2+galaxy0) %} with the following parameters:
> > >    - {% icon param-file %} *"Netcdf file"*: `outfile.netcdf`
> > >    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from outfile.netcdf` (output of **NetCDF xarray Metadata Info** {% icon tool %})
> > >    - *"Choose the variable to plot"*: `air_temperature_at_2_metres`
> > >    - *"Name of latitude coordinate"*: `lat`
> > >    - *"Name of longitude coordinate"*: `lon`
> > >    - *"Datetime selection"*: `No`
> > >    - *"Range of values for plotting e.g. minimum value and maximum value (minval,maxval) (optional)"*: `220,320`
> > >    - *"Add country borders with alpha value [0-1] (optional)"*: `1.0`
> > >    - *"Add coastline with alpha value [0-1] (optional)"*: `1.0`
> > >    - *"Add ocean with alpha value [0-1] (optional)"*: `1.0`
> > >    - *"Specify plot title (optional)"*: `Air_temperature_at_2_metres___23:00:00_UTC_batlow`
> > >    - *"Specify which colormap to use for plotting (optional)"*: `batlow`
> > >
> > >    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator'}`
> > >
> > >
> > > The final plot is shown below:
> > >    ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-25 >at 18:00:00](../../images/x-array-map-plot/bat.png)
> > >Some other important color variants of the same map can be found below :
> > > - *colormap : oslo*
> > >![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-25 >at 18:00:00](../../images/x-array-map-plot/os.png)
> > >  - *colormap : bamako*
> > >![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-25 >at 18:00:00](../../images/x-array-map-plot/bama.png)
> > >  - *colormap : vik*
> > >![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-25 >at 18:00:00](../../images/x-array-map-plot/vi.png)
> > > - *colormap : tokyo*
> > >![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-25 >at 18:00:00](../../images/x-array-map-plot/tok.png)
> > > - *colormap : hawaii*
> > >![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-25 >at 18:00:00](../../images/x-array-map-plot/hawa.png)
> > {: .hands_on}
> >
> >
> >
> {: .solution}
>
{: .question}



> <question-title></question-title>
>
> 1. What insights can be driven from the data and how is the data important on a local level?
> 2. What are some giveaways from the tutorial ?
>
> > <solution-title></solution-title>
> >
> > 1. Every piece of data recites a story. The air temperature at a certain height has a lot of significance in major commercial and day to day activities. Read this data-blog on the above analysis. [Click Here](https://quickbeasts51429.github.io/Outreachy_Galaxy_Community_contributor/).
> > 2. The tutorial has summed up a proper way of plotting data from a netcdf file. It has discussed everything from loading of data to its final display. Some other key points to keep in mind are :
> > > 1. It may take some time while plotting the maps. It depends on traffic / load on the Galaxy server. It is suggested to have a 64-bit processor with 8GB RAM storage. Be patient.
> > > 2.  You can view as well as download the generated plots to use further.
> > > 3. Plotting over global maps is very convinient as you saw above. But many a times,  you want to plot a specific region, it becomes very easy using CDO tool. Refer to [this tutorial]({% link topics/climate/tutorials/pangeo/tutorial.md %}) for more info.
> >
> > 3. If you wish to present all the plotted maps at one place for comparision or analysis. It is a short and simple step and can be doe using the tool name `Image Montage`.
> >
> > {% tool [Image Montage](toolshed.g2.bx.psu.edu/repos/bgruening/graphicsmagick_image_montage/graphicsmagick_image_montage/1.3.31+galaxy1) %} with the following parameters:
> >   - {% icon param-files %} *"Images"*: `Map plots`
> >   - {% icon param-text %} *"# of images wide"*: `4`
> > Here you can see that all the projections we plotted above have been shown in a single image using `Image Montage`.
> > ![ECMWF Reanalysis Air temperature a 2 metres on 2022-05-31 Monatge for 23:00:00 ](../../images/x-array-map-plot/montage.png)
> {: .solution}
>
{: .question}




# Conclusion

We have learnt about the `xarray map plotting tool` dealing with the `netcdf data set`. The tutorial also discussed about the types of climate datasets. One of the tutorial is info about usage of different [colormaps](https://github.com/Quickbeasts51429/Xarray_ColorMaps/blob/main/index.md#color-maps) and [projections](https://github.com/Quickbeasts51429/xarray_projection) in xarray.
