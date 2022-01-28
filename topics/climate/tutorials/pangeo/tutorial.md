---
layout: tutorial_hands_on

title: Pangeo ecosystem 101 for everyone - Introduction to Xarray Galaxy Tools
zenodo_link: https://doi.org/10.5281/zenodo.5805953
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
time_estimation: '1H'
key_points:
- Xarray Tools in Galaxy
tags:
  - pangeo
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
> Xarray works with labelled multi-dimensional arrays and can be used for a very wide range of data and data formats. In this training material, we focus on the usage of Xarray for Earth Science data following [CF-Convention](https://cfconventions.org/). However, some Galaxy Tools also work for non Earth Science datasets and if needed current Xarray Galaxy Tools could be extended to accomodate new usage.
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
> 1. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
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
> 2. If needed rename the datasets to `CAMS-PM2_5-20211222.netcdf`
> 3. Check that the datatype is `netcdf`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 4. Add a tag corresponding to ads (for Atmosphere Data Service)
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}


# Understanding our dataset


> ### {% icon details %} More information about CAMS PM2.5 forecast datasetss
>
> Our CAMS PM2.5 forecast dataset is in [netCDF](https://en.wikipedia.org/wiki/NetCDF) format. You could find the same dataset in different formats such as [GRIdded Binary or General Regularly-distributed Information in Binary form (GRIB)](https://en.wikipedia.org/wiki/GRIB) or [geoTIFF](https://en.wikipedia.org/wiki/GeoTIFF). The same Xarray Tools can be used with these other data formats. More information about this particular data set can be found on the [CAMS European air quality forecast webpage](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-europe-air-quality-forecasts?tab=overview).
{: .details}


To understand what is contained in our dataset, we will first use Xarray metadata Galaxy Tool. That will give us all the metadata information about the dataset.


## Get metadata

### Global metadata information

> ### {% icon hands_on %} Hands-on: netCDF dataset with Xarray metadata Galaxy Tool
>
> 1. {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %} with the following parameters:
>
>    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222.netcdf`
> 2. View {% icon galaxy-eye%} the two generated outputs:
>    - `Metadata infos` is a `tabular` providing the list of variables, their dimension names and number of elements per dimension. This file is used by other Xarray Tools. 
>    - The second file `info file` provide a summary of the **Xarray Dataset** contained in your netCDF file.
{: .hands_on}

In `info file` output file, we can identify 4 different sections:
1. **Dimensions**: name of dimensions and corresponding number of elements;
2. **Coordinates**: contains coordinate arrays (longitude, latitude, level and time) with their values.
3. **Data variables**: contains all the variables available in the dataset. Here, we only have one variable. For each variable, we get information on its shape and values.
4. **Global Attributes**: at this level, we get the global attributes of the dataset. Each attribute has a name and a value. 


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
> 1. {% tool [NetCDF xarray Coordinate Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_coords_info/xarray_coords_info/0.18.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222.netcdf`
>
> 2. View {% icon galaxy-eye%} the 5 generated outputs:
>    - `latitude`: a tabular file containing all the latitude values of our Xarray dataset;
>    - `longitude`: a tabular file containing all the longitudes values;
>    - `level`: this file contains information on all the different levels (here, we have surface data so level=0 meter);
>    - `time`: this tabular file contains all the forecast times. In our case, these are relative to December 22, 2021;
>    - `version`: this is a text file returning the Xarray package version. It is useful when publishing your Galaxy workflow.
>
>    > ### {% icon comment %} Comment
>    >
>    > This tool returns as many tabular files as the number of coordinate variables present in your input file. The values are decoded from the netCDF input file and no further processing is done. So units for instance for latitudes, longitudes, level and time may vary from one file to another depending on how it was coded in the original input file.
>    {: .comment}
>
{: .hands_on}


> ### {% icon question %} Understanding PM2.5 forecast coordinates
>
>  1. What is the units of the `time` coordinate?
>  2. What is the frequency of PM2.5 forecasts? 
>  3. What is the range of values for latitudes and longitudes?
>
> > ### {% icon solution %} Solution
> > 1. `info file` tells us that `time` is coded as `timedelta64[ns]` e.g. as differences in times (here in nanoseconds). Here the reference time is December 22, 2021. If we look at the tabular file named `time` (generated by `NetCDF xarray Coordinate Info`), we see that these times are automatically converted to human readable time format when printed:
> >
> > > ### {% icon code-out %} Output
> > > ```bash
> > > 0	0 days 00:00:00
> > > 1	0 days 01:00:00
> > > 2	0 days 02:00:00
> > > 3	0 days 03:00:00
> > > 4	0 days 04:00:00
> > > ```
> > {: .code-out}
> > This tells us that we have hourly forecast data. The last forecast time is `4 days 00:00:00` which means that the last forecast is in 4 days at 00:00 UTC (from December 22, 2021).
> {: .solution }
{: .question }

# Plotting our dataset on a geographical map

> ### {% icon hands_on %} Hands-on: Map plot
>  We will use Xarray mapplot Galaxy Tool to plot PM2.5 on December 22, 2021.
>
> 1. {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.18.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222.netcdf`
>    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222.netcdf`
>    - *"Choose the variable to plot"*: `pm2p5_conc`
>    - *"Name of latitude coordinate"*: `latitude`
>    - *"Name of longitude coordinate"*: `longitude`
>    - *"Datetime selection"*: `Yes`
>        - {% icon param-file %} *"Tabular of time values"*: `time` 
>        - *"Choose the times to plot"*: `0 days 00:00:00`
>    - *"Shift longitudes [0,360] --> [-180,180]"*: `Yes`
>    - *"Range of values for plotting e.g. minimum value abd maximum value (minval,maxval) (optional)"*: `0,35`
>    - *"Add country borders with alpha value [0-1] (optional)"*: `0.2`
>    - *"Add coastline with alpha value [0-1] (optional)"*: `0.5`
>    - *"Specify which colormap to use for plotting (optional)"*: `roma_r`
>    - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator', 'central_longitude': 12.0}`
>
>   ![CAMS PM2.5 December, 22th 2021 at 00:00 UTC](../../images/PM2_5_galaxy_20211222-00UTC.png)
>
>    > ### {% icon comment %} Why shifting longitudes?
>    >
>    > Longitudes are coded from 0 to 360 degrees. As we do not have global data but only covering Europe, we need to shift longitudes so that `NetCDF xarray map plotting` can plot properly our dataset. 
>    {: .comment}
>
> 
{: .hands_on}


> ### {% icon question %} Visualize and Compare
>
> Make a plot to Visualize the forecast for December, 24th 2021 at 12:00 UTC. Do you see any obvious differences with the plot from December 22, 2021 at 00:00 UTC?
>
> > ### {% icon solution %} Solution
> > 
> > Data starts on December, 22nd 2021 at 00:00 UTC so we need to add 2 days and 12 hours to select the correct time index. We reuse the same `NetCDF xarray map plotting` with a different selection for time:
> >
> > {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.18.2+galaxy0) %} with the following parameters:
> >    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222.netcdf`
> >   - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222.netcdf`
> >   - *"Choose the variable to plot"*: `pm2p5_conc`
> >   - *"Name of latitude coordinate"*: `latitude`
> >   - *"Name of longitude coordinate"*: `longitude`
> >   - *"Datetime selection"*: `Yes`
> >       - {% icon param-file %} *"Tabular of time values"*: `time` 
> >       - *"Choose the times to plot"*: `2 days 12:00:00`
> >   - *"Shift longitudes [0,360] --> [-180,180]"*: `Yes`
> >   - *"Range of values for plotting e.g. minimum value abd maximum value (minval,maxval) (optional)"*: `0,35`
> >   - *"Add country borders with alpha value [0-1] (optional)"*: `0.2`
> >   - *"Add coastline with alpha value [0-1] (optional)"*: `0.5`
> >   - *"Specify which colormap to use for plotting (optional)"*: `roma_r`
> >   - *"Specify the projection (proj4) on which we draw e.g. {"proj":"PlateCarree"} with double quote (optional)"*: `{'proj': 'Mercator', 'central_longitude': 12.0}`
> >   ![CAMS PM2.5 December, 24th 2021 at 12:00 UTC](../../images/PM2_5_galaxy_20211224-12UTC.png)
> {: .solution }
{: .question }

# Select / Subset from coordinates


> ### {% icon hands_on %} Hands-on: NetCDF xarray operations manipulate xarray from netCDF and save back to netCDF
>
> 1. {% tool [NetCDF xarray operations](toolshed.g2.bx.psu.edu/repos/ecology/xarray_netcdf2netcdf/xarray_netcdf2netcdf/0.18.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222.netcdf`
>    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222.netcdf`
>    - *"Choose the variable to extract"*: `pm2p5_conc`
>    - In *"additional filter"*:
>        - {% icon param-repeat %} *"Insert additional filter"*
>            - *"Dimensions"*: `time`
>            - *"Comparator"*: `slice(threshold1,threshold2)`
>                - *"Choose the start value for slice"*: `0 days 00:00:00`
>                - *"Choose the end value for slice"*: `1 days 00:00:00`
> 2. Rename the output dataset to `CAMS-PM2_5-20211222_fc0-23h.netcdf`
> 3.  Add a tag corresponding to `0-23h`
> 4. {% tool [NetCDF xarray Coordinate Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_coords_info/xarray_coords_info/0.18.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222_fc0-23h.netcdf`
> 5. Check the generated outputs and in particular `time`. We see that the tabular file `time` only contains **24** lines with times from `0 days 00:00:00` to `0 days 23:00:00`
>
>    > ### {% icon comment %} slice *threshold2* not included in selection
>    >
>    > You may have noticed already but when selecting a range with `slice` the upper limit (here `1 days 00:00:00`) is not included.
>    {: .comment}
>  
>
{: .hands_on}


> ### {% icon question %} PM2.5 over Italy region
>
> Using a selection and making plots of PM2.5 over Italy (latitudes: 43.N, 40.N and longitudes: 11.E,15.E), can you tell us if the forecasted PM2.5 will increase or decrease during the next 24 hours between 10:00 UTC and 17:00 UTC? Over which town in Italy do you see high values?
>
> > ### {% icon solution %} Solution
> > 1. {% tool [NetCDF xarray operations](toolshed.g2.bx.psu.edu/repos/ecology/xarray_netcdf2netcdf/xarray_netcdf2netcdf/0.18.2+galaxy0) %} with the following parameters:
> >   - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222.netcdf`
> >   - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222.netcdf`
> >   - *"Choose the variable to extract"*: `pm2p5_conc`
> >   - In *"additional filter"*:
> >       - {% icon param-repeat %} *"Insert additional filter"*
> >           - *"Dimensions"*: `time`
> >           - *"Comparator"*: `slice(threshold1,threshold2)`
> >               - *"Choose the start value for slice"*: `0 days 10:00:00`
> >               - *"Choose the end value for slice"*: `0 days 18:00:00`
> >       - {% icon param-repeat %} *"Insert additional filter"*
> >           - *"Dimensions"*: `latitude`
> >           - *"Comparator"*: `slice(threshold1,threshold2)`
> >               - *"Choose the start value for slice"*: `43.05`
> >               - *"Choose the end value for slice"*: `40.05`
> >       - {% icon param-repeat %} *"Insert additional filter"*
> >           - *"Dimensions"*: `longitude`
> >           - *"Comparator"*: `slice(threshold1,threshold2)`
> >               - *"Choose the start value for slice"*: `11.05`
> >               - *"Choose the end value for slice"*: `15.05`
> > 2. Rename the output dataset to `CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
> > 3.  Add a tag corresponding to `0-23h-Italy`
> > 4. {% tool [NetCDF xarray Metadata Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_metadata_info/xarray_metadata_info/0.15.1) %} with the following parameters:
> >   - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
> > 5. {% tool [NetCDF xarray Coordinate Info](toolshed.g2.bx.psu.edu/repos/ecology/xarray_coords_info/xarray_coords_info/0.18.2+galaxy0) %} with the following parameters:
> >   - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
> > 6.  {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.18.2+galaxy0) %} with the following parameters:
> >   - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222_fc10-18h_Italy.netcdf`
> >   - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
> >   - *"Choose the variable to plot"*: `pm2p5_conc`
> >   - *"Name of latitude coordinate"*: `latitude`
> >   - *"Name of longitude coordinate"*: `longitude`
> >   - *"Datetime selection"*: `Yes`
> >   - {% icon param-file %} *"Tabular of time values"*: `time`
> >   - *"Choose the times to plot"*: **Tick Select all**
> >   - *"Shift longitudes [0,360] --> [-180,180]"*: `No`
> >   - *"Range of values for plotting e.g. minimum value abd maximum value (minval,maxval) (optional)"*: `0,35`
> >   - *"Add country borders with alpha value [0-1] (optional)"*: `0.2`
> >   - *"Add coastline with alpha value [0-1] (optional)"*: `0.5`
> >   - *"Specify which colormap to use for plotting (optional)"*: `roma_r` 
> > 8. {% tool [Image Montage](toolshed.g2.bx.psu.edu/repos/bgruening/graphicsmagick_image_montage/graphicsmagick_image_montage/1.3.31+galaxy1) %} with the following parameters:
> >   - {% icon param-files %} *"Images"*: `Map plots`
> >   - {% icon param-text %} *"# of images wide"*: `4`
> > 
> > ![CAMS PM2.5 Italy 10:00 - 17:00 December 22, 2021](../../images/PM2_5_galaxy_20211222_Italy-10-17UTC.png)
> >
> > From the plot there is no obvious trend over this entire region of Italy. Howver, we clearly see that PM2.5 is always higher over Naples and tend to increase during the day (on that particular date).
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
- For instance, we may want to only keep PM2.5 if values are greater than a chosen threshold.

> ### {% icon hands_on %} Hands-on: Plot where PM2.5 is greater than 30 μm.m-3
> 1.  {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.18.2+galaxy0) %} with the following parameters:
>   - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
>   - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
>   - *"Choose the variable to plot"*: `pm2p5_conc`
>   - *"Name of latitude coordinate"*: `latitude`
>   - *"Name of longitude coordinate"*: `longitude`
>   - *"Datetime selection"*: `Yes`
>   - {% icon param-file %} *"Tabular of time values"*: `time` 
>   - *"Choose the times to plot"*: `0 days 10:00:00` 
>   - *"Shift longitudes [0,360] --> [-180,180]"*: `No`
>   - *"Range of values for plotting e.g. minimum value abd maximum value (minval,maxval) (optional)"*: `0,35`
>   - *"Do not plot values below this threshold (optional)"*: `30`
>   - *"Add country borders with alpha value [0-1] (optional)"*: `0.2`
>   - *"Add coastline with alpha value [0-1] (optional)"*: `0.5`
>   - *"Specify which colormap to use for plotting (optional)"*: `roma_r`
>    ![CAMS PM2.5 Italy 10:00  December 22, 2021](../../images/PM2_5_galaxy_20211222_Italy-10UTC-30.png)
> Now we clearly see that values of PM2.5 >  30 μm.m-3 are only found over Naples on December 22, 10:00 UTC.
>
{: .hands_on}



> ### {% icon question %} PM2.5 over Italy over 30 μm.m-3
>
>Using the same geographical region over Italy, can you tell us if the forecasted PM2.5 will exceed 30 μm.m-3 between 10:00 UTC and 17:00 UTC on december 22, 2021? 
>
> > ### {% icon solution %} Solution
> > 1.  {% tool [NetCDF xarray map plotting](toolshed.g2.bx.psu.edu/repos/ecology/xarray_mapplot/xarray_mapplot/0.18.2+galaxy0) %} with the following parameters:
> >   - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
> >   - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222_fc10-17h_Italy.netcdf`
> >   - *"Choose the variable to plot"*: `pm2p5_conc`
> >   - *"Name of latitude coordinate"*: `latitude`
> >   - *"Name of longitude coordinate"*: `longitude`
> >   - *"Datetime selection"*: `Yes`
> >         - {% icon param-file %} *"Tabular of time values"*: `time` 
> >         - *"Choose the times to plot"*: **Tick Select all** 
> >   - *"Shift longitudes [0,360] --> [-180,180]"*: `No`
> >   - *"Range of values for plotting e.g. minimum value abd maximum value (minval,maxval) (optional)"*: `0,35`
> >   - *"Do not plot values below this threshold (optional)"*: `30`
> >   - *"Add country borders with alpha value [0-1] (optional)"*: `0.2`
> >   - *"Add coastline with alpha value [0-1] (optional)"*: `0.5`
> >   - *"Specify which colormap to use for plotting (optional)"*: `roma_r`
> >   
> > 7. {% tool [Image Montage](toolshed.g2.bx.psu.edu/repos/bgruening/graphicsmagick_image_montage/graphicsmagick_image_montage/1.3.31+galaxy1) %} with the following parameters:
> >   - {% icon param-files %} *"Images"*: `Map plots`
> >   - {% icon param-text %} *"# of images wide"*: `4`
> > 
> > ![CAMS PM2.5 Italy 10:00 - 17:00 December 22, 2021 with PM2.5 > 30 μm.m-3](../../images/PM2_5_galaxy_20211222_Italy-10-17UTC_30.png)
> > Using thresholds, we can clearly identify areas if any where we have "high" values of PM2.5. On that particular day, PM2.5 is not very high but there are a few pixels where we have values greater than 30 μm.m-3.
> {: .solution }
{: .question }

# From Xarray to Tabular Data

> ### {% icon hands_on %} Hands-on: Xarray selection
> 
> We will select a single location: Naples (40.8518° N, 14.2681° E) and select the grid point that is closest to Naples.
> 1. {% tool [NetCDF xarray operations](toolshed.g2.bx.psu.edu/repos/ecology/xarray_netcdf2netcdf/xarray_netcdf2netcdf/0.18.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222.netcdf`
>    - {% icon param-file %} *"Tabular of variables"*: `Metadata infos from CAMS-PM2_5-20211222.netcdf`
>    - *"Choose the variable to extract"*: `pm2p5_conc`
>    - In *"additional filter"*:
>        - {% icon param-repeat %} *"Insert additional filter"*
>            - *"Dimensions"*: `latitude`
>            - *"Comparator"*: `slice(threshold1,threshold2)`
>                - *"Choose the start value for slice"*: `40.95`
>                - *"Choose the end value for slice"*: `40.85`
>        - {% icon param-repeat %} *"Insert additional filter"*
>            - *"Dimensions"*: `longitude`
>            - *"Comparator"*: `slice(threshold1,threshold2)`
>                - *"Choose the start value for slice"*: `14.25`
>                - *"Choose the end value for slice"*: `14.35`
> 2. Rename the output dataset to `CAMS-PM2_5-20211222_Naples.netcdf`
> 3.  Add a tag corresponding to `Naples`
> 4. {% tool [NetCDF xarray Selection](toolshed.g2.bx.psu.edu/repos/ecology/xarray_select/xarray_select/0.15.1) %} with the following parameters:
>    - {% icon param-file %} *"Netcdf file"*: `CAMS-PM2_5-20211222_Naples.netcdf`
>    - *"Choose the variable to extract"*: `pm2p5_conc`
>    - *"Source of coordinates"*: `Manually enter coordinates`
>        - *"Geographical area"*: `Whole available region`
>    - In *"Select Time series"*:
>        - *"Datetime selection"*: `No`
>
> 5. Rename the output dataset to `CAMS-PM2_5-20211222_Naples.tabular`
> 6. View {% icon galaxy-eye%} the generated file. It is a tabular with timeseries of PM2.4 concentrations over Naples. The total number of lines is 97 but we only print the first 5 lines.
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > 	time	level	latitude	longitude	pm2p5_conc
>    > 	0	0 days 00:00:00.000000000	0.0	40.95000076293945	14.25	24.00212
>    > 	1	0 days 01:00:00.000000000	0.0	40.95000076293945	14.25	23.5767
>    > 	2	0 days 02:00:00.000000000	0.0	40.95000076293945	14.25	21.383186
>    > 	3	0 days 03:00:00.000000000	0.0	40.95000076293945	14.25	20.04839
>    > 	4	0 days 04:00:00.000000000	0.0	40.95000076293945	14.25	18.347801
>    > ```
>    {: .code-out}
>
{: .hands_on}



> ### {% icon question %} PM2.5 at Naples over the 4 forecasted days
>
> From a qualitative point of view, can you say if PM2.5 may increase or decrease over the 4 forecasted days?
>
> > ### {% icon solution %} Solution
> >  We can make a simple plot using Scatterplot with ggplot2 ot climate stripes:
> > 1. {% tool [Scatterplot with ggplot2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_point/ggplot2_point/2.2.1+galaxy2) %} with the following parameters:
> >   - {% icon param-file %} *Input in tabular format"*: `CAMS-PM2_5-20211222_Naples.tabular`
> >   - *"Column to plot on x-axis"*: `1`
> >   - *"Column to plot on y-axis"*: `6`
> >   - *"Label for x axis"*: `Forecast time (hour) from December, 22 2021`
> >   - *"Label for y axis"*: `Particule Matter < 2.5 um-3`
> >   - In *"Advanced Options"*:
> >       - *"Type of plot"*: `Lines only `
> >       - *"Data point options"*: `User defined point options`
> >           - *"Transparency of points (On a scale of 0-1; 0=transparent, 1=default)"*: `0.7`
> >       - *"Plotting multiple groups"*: `No thanks - just plot the data as one group`
> >       - *"Axis title options"*: `Default`
> >       - *"Axis text options"*: `Default`
> >       - *"Plot title options"*: `Default`
> >       - *"Axis scaling"*: `Automatic axis scaling`
> >     ![CAMS PM2.5 Naples](../../images/PM2_5_galaxy_20211222_Naples_ggplot2.png)
> > 2. {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.1) %} with the following parameters:
> >   - {% icon param-file %} *"Select cells from"*: `CAMS-PM2_5-20211222_Naples.tabular`
> >   - *"using column"*: `c2`
> >   - In *"Check"*:
> >       - {% icon param-repeat %} *"Insert Check"*
> >           - *"Find Regex"*: `0 days`
> >           - *"Replacement"*: `20211222`
> >       - {% icon param-repeat %} *"Insert Check"*
> >           - *"Find Regex"*: `1 days`
> >           - *"Replacement"*: `20211223`
> >       - {% icon param-repeat %} *"Insert Check"*
> >           - *"Find Regex"*: `2 days`
> >           - *"Replacement"*: `20211224`
> >       - {% icon param-repeat %} *"Insert Check"*
> >           - *"Find Regex"*: `3 days`
> >           - *"Replacement"*: `20211225`
> >       - {% icon param-repeat %} *"Insert Check"*
> >           - *"Find Regex"*: `4 days`
> >           - *"Replacement"*: `20211226`
> > 3. Rename your dataset to `CAMS-PM2_5-20211222_Naples_with_dates.tabular`
> > 4. {% tool [climate stripes](toolshed.g2.bx.psu.edu/repos/climate/climate_stripes/climate_stripes/1.0.1) %} with the following parameters:
> >   - *"column name to use for plotting"*: `pm2p5_conc`
> >   - *"plot title"*: `PM2.5 4 days forecast from December 22 2021 over Naples`
> >   - In *"Advanced Options"*:
> >       - *"column name to use for x-axis"*: `time`
> >       - *"format for input date/time column"*: `%Y%m%d %H:%M:%S.%f`
> >       - *"format for plotting dates on the x-axis"*: `%d %b %H hours`
> >       - *""*: `winter`
> >  ![CAMS PM2.5 Naples Stripes](../../images/PM2_5_galaxy_20211222_Naples_stripes.png)
> > From December 24, 2021 at 00:00 UTC and onwards, PM2-5 concentration are much lower than at the beginning. This is bothe visible on the 1D plot and stripes.
> {: .solution }
{: .question }


# Conclusion
{:.no_toc}

{% icon trophy %} Well done! In this tutorial, Xarray Galaxy Tools have been introduced and we learned to use these tools on a real dataset from Copernicus Atmosphere Monitoring Service. We encourage you to try with your own datasets. 
