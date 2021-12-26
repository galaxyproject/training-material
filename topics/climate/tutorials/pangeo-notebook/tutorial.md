---
layout: tutorial_hands_on

title: Pangeo Notebook in Galaxy - Introduction to Xarray
zenodo_link: ''
requirements:
  - type: "external"
    title: Programming with Python
    link: https://swcarpentry.github.io/python-novice-inflammation/
  - type: "external"
    title: Data Carpentry Geospatial Workshop
    link: https://datacarpentry.org/geospatial-workshop/
questions:
- What is Pangeo notebook?
- How to start Pangeo Notebook in Galaxy?
- What are the main software components of the Pangeo ecosystem?
- What is Xarray?
- How to manipulate Xarrays?
- How to print metadata information?
- How to make a selection?
- How to visualize?
- How to filter? 
- How to make reduction operations (mean, max, min)?
- How to resample my data?
- Where to go next?
objectives:
- Learn to get metadata information using Xarray Galaxy Tools
- Learn to select data
- Learn to visualize data
- Learn to filter, make reduction operations (mean, max, min)
- Learn to resample my data
- Learn to cite and contribute to Pangeo
time_estimation: 1H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
tags:
  - interactive-tools
contributors:
- annefou

---


# Introduction
{:.no_toc}

In this tutorial, we will learn about [Xarray](https://xarray.pydata.org/), one of the most used Python library from the [Pangeo](https://pangeo.io/) ecosystem. 

We will be using data from [Copernicus Ar Monitoring Service](https://ads.atmosphere.copernicus.eu/)
and more precisely PM2.5 ([Particle Matter < 2.5 μm](https://en.wikipedia.org/wiki/Particulates#Size,_shape_and_solubility_matter)) 4 days forecast from December, 22 2021.

> ### {% icon comment %} Remark
>
> This tutorial uses data on a regular latitude, longitude grid. More complex and irregular grids are not discussed in this tutorial. In addition,
> this tutorial is not meant to cover all the different possibilities offered by Xarrays but shows functionalities we find useful for day to day
> analysis.
>
{: .comment}

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
> 2. **Rename your history** to be meaningful and easy to find. For instance, you can choose **Xarray with Pangeo notebook** as the name of your new history.
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}


## Upload CAMS PM2.5 data

Data can be retrieved directly from [Copernicus Atmosphere Monitoring Service](https://ads.atmosphere.copernicus.eu/) but to make it easier, you can
download the tutorial data from Zenodo.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the dataset to **CAMS-PM2_5-20211222.netcdf**
> 4. Check that the datatype is **netcdf**
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

# Starting Galaxy Pangeo JupyterLab

> ### {% icon hands_on %} Hands-on: Launch Pangeo notebook JupyterLab in Galaxy
>
> Currently Pangeo notebook in Galaxy is available on [useGalaxy.eu](https://usegalaxy.eu) only. Make sure you login to Galaxy and search for Pangeo notebook and not the default JupyterLab in Galaxy to make sure you ahve all the Pangeo Software stack available. The default JupyterLab in Galaxy would not be sufficient for executing all the tasks in this tutorial.
>
> 1. Open the {% tool [Pangeo Notebook](interactive_tool_pangeo_notebook) %} by clicking [here](https://usegalaxy.eu/?tool_id=interactive_tool_pangeo_notebook){:target="_blank"}
> 2. *Include data into the environment*: `CAMS-PM2_5-20211222.netcdf`
> 3. Click on **Execute**: the tool will start running and will stay running 
> 4. Click on the **User** menu at the top and go to **Active Interactive Tools** and locate the Pangeo JupyterLab instance you started.
> 5. Click on your Pangeo JupyterLab instance
{: .hands_on}

# Analysis

## Import Python packages

> ### {% icon hands_on %} Hands-on: Import Python packages
>    > ### {% icon code-in %} Input: Python
>    > ```python
>    > import numpy as np
>    > import xarray as xr
>    > import cartopy.crs as ccrs
>    > import matplotlib.pyplot as plt
>    > import cmcrameri.cm as cmc
>    > import pandas as pd
>    >  ```
>    {: .code-in}
{: .hands_on}

## Open and read metadata

> ### {% icon hands_on %} Hands-on: netCDF dataset with Xarray
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    > dset = xr.open_dataset("CAMS-PM2_5-20211222.netcdf")
>    >  ```
>    {: .code-in}
{: .hands_on}


Once opened, we can get metadata using `print` statement.

> ### {% icon hands_on %} Hands-on: Get metadata
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset)
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    >   <xarray.Dataset>
>    >   Dimensions:     (longitude: 700, latitude: 400, level: 1, time: 97)
>    >   Coordinates:
>    >     * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
>    >     * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
>    >     * level       (level) float32 0.0
>    >     * time        (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
>    >   Data variables:
>    >       pm2p5_conc  (time, level, latitude, longitude) float32 0.4202 ... 7.501
>    >   Attributes:
>    >       title:        PM25 Air Pollutant FORECAST at the Surface
>    >       institution:  Data produced by Meteo France
>    >       source:       Data from ENSEMBLE model
>    >       history:      Model ENSEMBLE FORECAST
>    >       FORECAST:     Europe, 20211222+[0H_96H]
>    >       summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
>    >       project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
>    > ```
>    {: .code-out}
{: .hands_on}

We can identify 4 different sections:
1. **Dimensions**: name of dimensions and corresponding number of elements;
2. **Coordinates**: contains coordinate arrays (longitude, latitude, level and time) with their values.
3. **Data variables**: contains all the variables available in the dataset. Here, we only have one variable. For each variable, we get information on its shape and values.
4. **Attributes**: at this level, we get all the attributes of the dataset. 

We can also get metadata information for each coordinate and data variables using "." followed by the coordinate or data variable name.


> ### {% icon hands_on %} Hands-on: get metadata of a coordinate or data variable
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.time)
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.DataArray 'time' (time: 97)>
>    > array([              0,   3600000000000,   7200000000000,  10800000000000,
>    >         14400000000000,  18000000000000,  21600000000000,  25200000000000,
>    >         28800000000000,  32400000000000,  36000000000000,  39600000000000,
>    >         43200000000000,  46800000000000,  50400000000000,  54000000000000,
>    >         57600000000000,  61200000000000,  64800000000000,  68400000000000,
>    >         72000000000000,  75600000000000,  79200000000000,  82800000000000,
>    >         86400000000000,  90000000000000,  93600000000000,  97200000000000,
>    >        100800000000000, 104400000000000, 108000000000000, 111600000000000,
>    >        115200000000000, 118800000000000, 122400000000000, 126000000000000,
>    >        129600000000000, 133200000000000, 136800000000000, 140400000000000,
>    >        144000000000000, 147600000000000, 151200000000000, 154800000000000,
>    >        158400000000000, 162000000000000, 165600000000000, 169200000000000,
>    >        172800000000000, 176400000000000, 180000000000000, 183600000000000,
>    >        187200000000000, 190800000000000, 194400000000000, 198000000000000,
>    >        201600000000000, 205200000000000, 208800000000000, 212400000000000,
>    >        216000000000000, 219600000000000, 223200000000000, 226800000000000,
>    >        230400000000000, 234000000000000, 237600000000000, 241200000000000,
>    >        244800000000000, 248400000000000, 252000000000000, 255600000000000,
>    >        259200000000000, 262800000000000, 266400000000000, 270000000000000,
>    >        273600000000000, 277200000000000, 280800000000000, 284400000000000,
>    >        288000000000000, 291600000000000, 295200000000000, 298800000000000,
>    >        302400000000000, 306000000000000, 309600000000000, 313200000000000,
>    >        316800000000000, 320400000000000, 324000000000000, 327600000000000,
>    >        331200000000000, 334800000000000, 338400000000000, 342000000000000,
>    >        345600000000000], dtype='timedelta64[ns]')
>    > Coordinates:
>    >   * time     (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
>    > Attributes:
>    >     long_name:  FORECAST time from 20211222
>    > ```
>    {: .code-out}
{: .hands_on}



> ### {% icon question %} Questions CAM PM2.5 Dataset
>
> What is the name of the variable for Particle matter < 2.5 μm and its physical units?
>
> > ### {% icon solution %} Solution
> > 1. To get metadata information from `pm2p5_conc`Data variable, we use its variable name and print it. Printing it will only print metadata, not the values.
> >      - Variable name: `mass_concentration_of_pm2p5_ambient_aerosol_in_air`
> >      - Units: `µg/m3`
> >
> > > ### {% icon code-in %} Input: Python
> > > ```python
> > > print(dset.pm2p5_conc)
> > > ```
> > {: .code-in}
> > 
> > > ### {% icon code-out %} Output
> > > ```bash
> > >    <xarray.DataArray 'pm2p5_conc' (time: 97, level: 1, latitude: 400, longitude: 700)>
> > >    [27160000 values with dtype=float32]
> > >    Coordinates:
> > >      * longitude  (longitude) float32 335.0 335.1 335.2 335.4 ... 44.75 44.85 44.95
> > >      * latitude   (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
> > >      * level      (level) float32 0.0
> > >      * time       (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
> > >    Attributes:
> > >        species:        PM2.5 Aerosol
> > >        units:          µg/m3
> > >        value:          hourly values
> > >        standard_name:  mass_concentration_of_pm2p5_ambient_aerosol_in_air
> > > ```
> > {: .code-out}
> {: .solution }
{: .question }

> ### {% icon comment %} Different ways to access Data variables
>
> To access a variable or coordinate, we can use "." or specify its name as a string between squared brackets "[" "]". For example:
>  > ### {% icon code-in %} Input: Python
>  > ```python
>  > print(dset['pm2p5_conc'])
>  > ```
> {: .code-in}
>
> > ### {% icon code-out %} Output
>  > ```bash
>  >   <xarray.DataArray 'pm2p5_conc' (time: 97, level: 1, latitude: 400, longitude: 700)>
>  >   [27160000 values with dtype=float32]
>  >   Coordinates:
>  >     * longitude  (longitude) float32 335.0 335.1 335.2 335.4 ... 44.75 44.85 44.95
>  >     * latitude   (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
>  >     * level      (level) float32 0.0
>  >     * time       (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
>  >   Attributes:
>  >       species:        PM2.5 Aerosol
>  >       units:          µg/m3
>  >       value:          hourly values
>  >       standard_name:  mass_concentration_of_pm2p5_ambient_aerosol_in_air
>  > ```
> {: .code-out}
>
{: .comment}


## Select / Subset from coordinates

> ### {% icon hands_on %} Hands-on: Select elements from coordinates by index
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.isel(time=0))
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.Dataset>
>    > Dimensions:     (longitude: 700, latitude: 400, level: 1)
>    > Coordinates:
>    >   * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
>    >   * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
>    >   * level       (level) float32 0.0
>    >     time        timedelta64[ns] 00:00:00
>    > Data variables:
>    >     pm2p5_conc  (level, latitude, longitude) float32 0.4202 0.4331 ... 14.22
>    > Attributes:
>    >     title:        PM25 Air Pollutant FORECAST at the Surface
>    >     institution:  Data produced by Meteo France
>    >     source:       Data from ENSEMBLE model
>    >     history:      Model ENSEMBLE FORECAST
>    >     FORECAST:     Europe, 20211222+[0H_96H]
>    >     summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
>    >     project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
>    > ```
>    {: .code-out}
{: .hands_on}



> ### {% icon hands_on %} Hands-on: Select elements from coordinates by value
> When selecting elements by the value of the coordinate, we need to use the same datatype. For instance, to select and element from
> `time`, we need to use `timedelta64`. The code below will give the same result than `isel(time=0)`.
>
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.sel(time=np.timedelta64(0)))
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.Dataset>
>    > Dimensions:     (longitude: 700, latitude: 400, level: 1)
>    > Coordinates:
>    >   * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
>    >   * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
>    >   * level       (level) float32 0.0
>    >     time        timedelta64[ns] 00:00:00
>    > Data variables:
>    >     pm2p5_conc  (level, latitude, longitude) float32 0.4202 0.4331 ... 14.22
>    > Attributes:
>    >     title:        PM25 Air Pollutant FORECAST at the Surface
>    >     institution:  Data produced by Meteo France
>    >     source:       Data from ENSEMBLE model
>    >     history:      Model ENSEMBLE FORECAST
>    >     FORECAST:     Europe, 20211222+[0H_96H]
>    >     summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
>    >     project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
>    > ```
>    {: .code-out}
{: .hands_on}


> ### {% icon question %} Select a single time for PM2.5
>
> How to select the forecast for December, 24th 2021 at 12:00 UTC?
>
> > ### {% icon solution %} Solution
> > Data starts on December, 22nd 2021 at 00:00 UTC so we need to add 2 days and 12 hours to select the correct time index.`
> >
> > > ### {% icon code-in %} Input: Python
> > > ```python
> > > print(dset.sel(time=(np.timedelta64(2,'D')+ np.timedelta64(12,'h'))))
> > > ```
> > {: .code-in}
> > 
> > > ### {% icon code-out %} Output
> > > ```bash
> > > <xarray.Dataset>
> > > Dimensions:     (longitude: 700, latitude: 400, level: 1)
> > > Coordinates:
> > >   * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
> > >   * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
> > >   * level       (level) float32 0.0
> > >     time        timedelta64[ns] 2 days 12:00:00
> > > Data variables:
> > >     pm2p5_conc  (level, latitude, longitude) float32 0.4499 0.4421 ... 10.71
> > > Attributes:
> > >     title:        PM25 Air Pollutant FORECAST at the Surface
> > >     institution:  Data produced by Meteo France
> > >     source:       Data from ENSEMBLE model
> > >     history:      Model ENSEMBLE FORECAST
> > >     FORECAST:     Europe, 20211222+[0H_96H]
> > >     summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
> > >     project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
> > > ```
> > {: .code-out}
> {: .solution }
{: .question }


## Plotting

- To plot a map, you need to select a variable with data on geographical coordinates (latitude, longitude).
- In addition coordinates need to be sorted (preferably in increasing order). This is not the case for "longitude" because they are expressed as 0 to 360. Let's shift them to -180 to 180.

> ### {% icon hands_on %} Hands-on: Shift longitude from (0, 360) to (-180, 180)
> 
> We print the longitudes before and after shifting them so we can see what is happening.
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.longitude)
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.DataArray 'longitude' (longitude: 700)>
>    > array([335.05, 335.15, 335.25, ...,  44.75,  44.85,  44.95], dtype=float32)
>    > Coordinates:
>    >   * longitude  (longitude) float32 335.0 335.1 335.2 335.4 ... 44.75 44.85 44.95
>    > ```
>    {: .code-out}
>
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    > dset.coords['longitude'] = (dset['longitude'] + 180) % 360 - 180
>    > print(dset.longitude)
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.DataArray 'longitude' (longitude: 700)>
>    > array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
>    >         44.949997], dtype=float32)
>    > Coordinates:
>    >   * longitude  (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
>    > ```
>    {: .code-out}
> 
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Visualize on a map PM2.5 for December, 24th 2021 at 12:00 UTC
> 
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  dset.sel(time=(np.timedelta64(2,'D')+ np.timedelta64(12,'h'))).pm2p5_conc.plot()
>    >  ```
>    {: .code-in}
>
>  ![CAMS PM2.5 December, 24th 2021 at 12:00 UTC](../../images/PM2_5_default.png)
{: .hands_on}


> ### {% icon comment %} What about `level`
> Note that in the previous plot, we did not need to select `level` because there is one value only. However, if we had more than one level, we would need to add a selection on the level before plotting
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Customize your plot
> There are many ways to customize your plots and here we only give you what we think is important for creating publication ready figures:
> - Define the size of the figure
> - Choose to project data on a different projection.
> - Add coastline
> - Set the min and max values for plotting
> - Add a title, change colorbar title
> - Save figure into png
>
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  fig = plt.figure(1, figsize=[15,10])
>    >  
>    >  # We're using cartopy to project our data.
>    >  # (see documentation on cartopy)
>    >  ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator())
>    >  ax.coastlines(resolution='10m')
>    >  
>    >  # We need to project our data to the new projection and for this we use `transform`.
>    >  # we set the original data projection in transform (here PlateCarree)
>    >  dset.sel(time=(np.timedelta64(2,'D') + np.timedelta64(12,'h')))['pm2p5_conc'].plot(ax=ax, 
>    >                                                                                    transform=ccrs.PlateCarree(),
>    >                                                                                    vmin = 0, vmax = 35,
>    >                                                                                    cmap=cmc.roma_r)
>    >  # One way to customize your title
>    >  plt.title("Copernicus Monitoring Service PM2.5, 2 day forecasts\n 24th December 2021 at 12:00 UTC", fontsize=18)
>    >  plt.savefig("CAMS-PM2_5-fc-20211224.png")
>    >  ```
>    {: .code-in}
>
>  ![Customized plot for CAMS PM2.5 December, 24th 2021 at 12:00 UTC](../../images/CAMS-PM2_5-fc-20211224.png)
{: .hands_on}


> ### {% icon hands_on %} Hands-on: Multi-plots
> - Here we would like to plot several times on the same figure in different sub-plots
> - We will not plots all the times (too many) but the first 24 forecasted values.
>
>  1. Make a list of times and convert to pandas datetime to make it easier to format times when plotting
>
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  list_times = np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))
>    >  print(pd.to_datetime(list_times).strftime("%d %b %H:%S UTC"))
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > Index(['22 Dec 00:00 UTC', '22 Dec 01:00 UTC', '22 Dec 02:00 UTC',
>    >        '22 Dec 03:00 UTC', '22 Dec 04:00 UTC', '22 Dec 05:00 UTC',
>    >        '22 Dec 06:00 UTC', '22 Dec 07:00 UTC', '22 Dec 08:00 UTC',
>    >        '22 Dec 09:00 UTC', '22 Dec 10:00 UTC', '22 Dec 11:00 UTC',
>    >        '22 Dec 12:00 UTC', '22 Dec 13:00 UTC', '22 Dec 14:00 UTC',
>    >        '22 Dec 15:00 UTC', '22 Dec 16:00 UTC', '22 Dec 17:00 UTC',
>    >        '22 Dec 18:00 UTC', '22 Dec 19:00 UTC', '22 Dec 20:00 UTC',
>    >        '22 Dec 21:00 UTC', '22 Dec 22:00 UTC', '22 Dec 23:00 UTC',
>    >        '23 Dec 00:00 UTC'],
>    >       dtype='object')
>    > ```
>    {: .code-out}
>
>   2. We use the same plotting method than earlier but we pass additional parameters:
>        - `vmin = 0`and `vmax = 35` to set the minimum and maximum values when plotting (this is useful to highlight features in your plot);
>        - `subplot_kws={"projection": proj_plot}` to project data on a non-default projection. See [cartopy projection](https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html) for more information about projections;
>        - `col='time'` because we will plot several `time`;
>        -  `col_wrap=4` to have a maximum of 4 plots per row. If we have more times to plot, then the next figures will be on another row;
>        - `robust=True` and `aspect=dset.dims["longitude"] / dset.dims["latitude"]` are additional parameters to make each subplot with a "sensible" figsize;
>        - `cmap=cmc.roma_r` to select a non-default and color-blind driendly colormap (see [scientific colormaps](https://www.fabiocrameri.ch/colourmaps/)).
>
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  fig = plt.figure(1, figsize=[10,10])
>    >  
>    >  # We're using cartopy to project our data.
>    >  # (see documentation on cartopy)
>    >  proj_plot = ccrs.Mercator()
>    >  
>    >  # We need to project our data to the new projection and for this we use `transform`.
>    >  # we set the original data projection in transform (here PlateCarree)
>    >  p = dset.sel(time=slice(np.timedelta64(1,'h'),np.timedelta64(1,'D')))['pm2p5_conc'].plot(transform=ccrs.PlateCarree(),
>    >                                                                                       vmin = 0, vmax = 35,
>    >                                                                                       subplot_kws={"projection": proj_plot},
>    >                                                                                       col='time', col_wrap=4,
>    >                                                                                       robust=True,
>    >                                                                                       aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
>    >                                                                                       cmap=cmc.roma_r)
>    >  # We have to set the map's options on our axes
>    >  for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
>    >      ax.coastlines('10m')
>    >      ax.set_title("CAMS PM2.5 " + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
>    >  # Save your figure
>    >  plt.savefig("CAMS-PM2_5-fc-multi.png")
>    >  ```
>    {: .code-in}
>
>  In the second part of our plot, we customize each subplot (this is why we loop for each of them and get their axes) by adding:
>         -  `coastlines`: we pass a parameter `10m` to get coastlines with a high resolution (non-default);
>         - `set_title` to set a title for each subplot.
>    
>  ![Customized multi-plot](../../images/CAMS-PM2_5-fc-multi.png)
> 
{: .hands_on}



> ### {% icon question %} PM2.5 over Italy
>
> Using a Multi-plot between Rome and Naples, can you tell us if the forecasted PM2.5 will increase or decrease during the first 24 hours?
>
> > ### {% icon solution %} Solution
> > We will select a sub-area: 11. East to 15.0 East and 40. N to 43. N. PM2.5 will increase and reach values close to 35 μm.m-3
> >
> > > ### {% icon code-in %} Input: Python
> > > ```python
> > > fig = plt.figure(1, figsize=[10,10])
> > > 
> > > # We're using cartopy to project our data.
> > > # (see documentation on cartopy)
> > > proj_plot = ccrs.Mercator()
> > > 
> > > # We need to project our data to the new projection and for this we use `transform`.
> > > # we set the original data projection in transform (here PlateCarree)
> > > p = dset.sel(time=slice(np.timedelta64(1,'h'),np.timedelta64(1,'D'))).sel(latitude=slice(43., 40.),
> > >                                                                           longitude=slice(11.,15.))['pm2p5_conc'].plot(transform=ccrs.PlateCarree(),
> > >                                                                                      vmin = 0, vmax = 35,
> > >                                                                                      subplot_kws={"projection": proj_plot},
> > >                                                                                      col='time', col_wrap=4,
> > >                                                                                      robust=True,
> > >                                                                                      aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
> > >                                                                                      cmap=cmc.roma_r)
> > > # We have to set the map's options on all axes
> > > for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
> > >     ax.coastlines('10m')
> > >     ax.set_title("CAMS PM2.5 " + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
> > > # Save your figure
> > > plt.savefig("CAMS-PM2_5-fc-multi-Italy.png")
> > > ```
> > {: .code-in}
> > 
> >  ![PM2.5 over Italy](../../images/CAMS-PM2_5-fc-multi-Italy.png)
> {: .solution }
{: .question }


> ### {% icon comment %}  `latitude=slice(47.3, 36.5)` and not `latitude=slice(36.5, 47.3)`
> Why did we slice latitudes with `latitude=slice(47.3, 36.5)` and not `latitude=slice(36.5, 47.3)`?
> - because when using slice, you need to specify values using the same order than in the coordinates. Latitudes are specified in 
> decreasing order for CAMS.
>
{: .comment}


## Where 

- Sometimes we may want to make more complex selections with criteria on the values of a given variable and not only on its coordinates. For this we use `where`.
- For instance, we may want to only keep PM2.5 if values are greater than 25 μm.m-3 (or any threshold you would like to choose)

> ### {% icon hands_on %} Hands-on: Mask values that do not meet a criteria with `Where`
> 
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.where(dset['pm2p5_conc'] > 25))
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.Dataset>
>    > Dimensions:     (time: 97, level: 1, latitude: 400, longitude: 700)
>    > Coordinates:
>    >   * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
>    >   * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
>    >   * level       (level) float32 0.0
>    >   * time        (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
>    > Data variables:
>    >     pm2p5_conc  (time, level, latitude, longitude) float32 nan nan ... nan nan
>    > Attributes:
>    >     title:        PM25 Air Pollutant FORECAST at the Surface
>    >     institution:  Data produced by Meteo France
>    >     source:       Data from ENSEMBLE model
>    >     history:      Model ENSEMBLE FORECAST
>    >     FORECAST:     Europe, 20211222+[0H_96H]
>    >     summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
>    >     project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
>    > ```
>    {: .code-out}
{: .hands_on}


> ### {% icon comment %} What happened?
> You may not see any changes but if you look carefuly to `pm2p5_conc` values, you will see that many `nan`. In fact, we now have 
> `nan` values Whenever the criteria within the `where` statement is not met e.g. when PM2.5 <= 25
>
{: .comment}

Let's plot one time to better see what happened:

> ### {% icon hands_on %} Hands-on: Plotting with mask
> 
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  fig = plt.figure(1, figsize=[15,10])
>    >  
>    >  # We're using cartopy to project our data.
>    >  # (see documentation on cartopy)
>    >  ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator())
>    >  ax.coastlines(resolution='10m')
>    >  
>    >  # We need to project our data to the new projection and for this we use `transform`.
>    >  # we set the original data projection in transform (here PlateCarree)
>    >  dset.where(dset['pm2p5_conc'] > 25).isel(time=0)['pm2p5_conc'].plot(ax=ax,
>    >                                                                      transform=ccrs.PlateCarree(),
>    >                                                                      vmin = 0, vmax = 35,
>    >                                                                      cmap=cmc.roma_r)
>    >  # One way to customize your title
>    >  plt.title("Copernicus Monitoring Service PM2.5, 2 day forecasts\n 24th December 2021 at 12:00 UTC\n only values > 25", fontsize=18)
>    >  plt.savefig("CAMS-PM2_5-fc-20211224-25.png")
>    >  ```
>    {: .code-in}
>
>  ![PM2.5 over Italy with threshold at 25](../../images/CAMS-PM2_5-fc-20211224-25.png)
{: .hands_on}

We can then make the same multi-plot as earlier (over Italy) but with a `where` statement to mask values lower than 25 μm.m-3:

> ### {% icon hands_on %} Hands-on: Multi-plot over Italy using a mask
> 
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  fig = plt.figure(1, figsize=[10,10])
>    >  
>    >  # We're using cartopy to project our data.
>    >  # (see documentation on cartopy)
>    >  proj_plot = ccrs.Mercator()
>    >  
>    >  # We need to project our data to the new projection and for this we use `transform`.
>    >  # we set the original data projection in transform (here PlateCarree)
>    >  p = dset.where(dset['pm2p5_conc'] > 25).sel(time=slice(np.timedelta64(1,'h'),np.timedelta64(1,'D'))).sel(latitude=slice(43., 40.),
>    >                                                                            longitude=slice(11.,15.))['pm2p5_conc'].plot(transform=ccrs.PlateCarree(),
>    >                                                                                       vmin = 0, vmax = 35,
>    >                                                                                       subplot_kws={"projection": proj_plot},
>    >                                                                                       col='time', col_wrap=4,
>    >                                                                                       robust=True,
>    >                                                                                       aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
>    >                                                                                       cmap=cmc.roma_r)
>    >  # We have to set the map's options on all four axes
>    >  for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
>    >      ax.coastlines('10m')
>    >      ax.set_title("PM2.5 > 25 μm.m-3" + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
>    >  # Save your figure
>    >  plt.savefig("CAMS-PM2_5-fc-multi-Italy-25.png")
>    >  ```
>    {: .code-in}
>
>  ![Multi-plot of PM2.5 over Italy with threshold at 25](../../images/CAMS-PM2_5-fc-multi-Italy-25.png)
>
{: .hands_on}


## Reduction operations
- We often want to compute the mean of all our dataset, or along a dimension (for instance time).
- If you do not pass any argument to the operation then it is done over all dimension.


> ### {% icon hands_on %} Hands-on: Mean
>  When we do not specify any paramters, we get a single value
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).mean())
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.Dataset>
>    > Dimensions:     ()
>    > Data variables:
>    >     pm2p5_conc  float32 9.118
>    > ```
>    {: .code-out}
{: .hands_on}



> ### {% icon question %} Maximum PM2.5 over Italy
>
> What is the maximum forecasted PM2.5 value over the Rome-Naples region?
>
> > ### {% icon solution %} Solution
> >  We  select the same sub-area: 11. East to 15.0 East and 40. N to 43. N and compute the maximum with `max`. The maximum PM2.5 value is **59.13694382** μm.m-3 (that is rounded to **59.14**).
> >
> > > ### {% icon code-in %} Input: Python
> > > ```python
> > > dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).max()
> > > ```
> > {: .code-in}
> > 
> > > ### {% icon code-out %} Output
> > > ```bash
> > > xarray.Dataset
> > > Dimensions:
> > > Coordinates: (0)
> > > Data variables:
> > > pm2p5_conc
> > > ()
> > > float64
> > > 59.14
> > > array(59.13694382)
> > > Attributes: (0)
> > > ```
> > {: .code-out}
> {: .solution }
{: .question }





> ### {% icon question %} Find when the maximum PM2.5 is forecasted
>
> When is forecasted the maximum PM2.5 value?
>
> > ### {% icon solution %} Solution
> > We will select a sub-area: 11. East to 15.0 East and 40. N to 43. N and average over the entire selected area and search where the maximum PM2.5 value of 59.13694382 μm.m-3 is found. The maximum PM2.5 value occurs on 2021-12-22 at 20:00 UTC.
> >
> > > ### {% icon code-in %} Input: Python
> > > ```python
> > > dset_tmean = dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).max(dim=('latitude', 'longitude'))
> > > dset_tmean_max = dset_tmean.where(dset_tmean['pm2p5_conc'] == 59.13694382, drop=True)
> > > print(dset_tmean_max)
> > > ```
> > {: .code-in}
> > 
> > > ### {% icon code-out %} Output
> > > ```bash
> > > <xarray.Dataset>
> > > Dimensions:     (time: 1, level: 1)
> > > Coordinates:
> > >   * level       (level) float32 0.0
> > >   * time        (time) timedelta64[ns] 20:00:00
> > > Data variables:
> > >     pm2p5_conc  (time, level) float32 59.14
> > > ```
> > {: .code-out}
> {: .solution }
{: .question }

> ### {% icon comment %} Pixel size when averaging
> We average over a relatively small so we do not make a weighted average. Use weighted averages when averaging over the entire globe or over a large area where the pixel sizes may vary (depending on the latitude).
>
{: .comment}


## Resample
- We often want to resample (on time dimension) our data:
    - if your resampling frequency is lower than your original data, you would need to apply an operation on the data you group together, such as mean, min, max;
    - if your resampling frequency is higher than your original data, you would need to indicate how to fill the gaps, for instance interpolate and indicate which interpolation method to apply or select nearest values, etc.

> ### {% icon hands_on %} Hands-on: 1 day Resampling
> 
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.resample(time='1D').mean())
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.Dataset>
>    > Dimensions:     (time: 5, longitude: 700, latitude: 400, level: 1)
>    > Coordinates:
>    >   * time        (time) timedelta64[ns] 0 days 1 days 2 days 3 days 4 days
>    >   * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
>    >   * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
>    >   * level       (level) float32 0.0
>    > Data variables:
>    >     pm2p5_conc  (time, level, latitude, longitude) float32 0.4298 ... 7.501
>    > ```
>    {: .code-out}
{: .hands_on}

> ### {% icon hands_on %} Hands-on: 30 minute resampling
> 
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  print(dset.resample(time='30min').interpolate('linear'))
>    >  ```
>    {: .code-in}
>
>    > ### {% icon code-out %} Output
>    > ```bash
>    > <xarray.Dataset>
>    > Dimensions:     (longitude: 700, latitude: 400, level: 1, time: 193)
>    > Coordinates:
>    >   * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
>    >   * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
>    >   * level       (level) float32 0.0
>    >   * time        (time) timedelta64[ns] 00:00:00 00:30:00 ... 4 days 00:00:00
>    > Data variables:
>    >     pm2p5_conc  (time, level, latitude, longitude) float64 0.4202 ... 7.501
>    > Attributes:
>    >     title:        PM25 Air Pollutant FORECAST at the Surface
>    >     institution:  Data produced by Meteo France
>    >     source:       Data from ENSEMBLE model
>    >     history:      Model ENSEMBLE FORECAST
>    >     FORECAST:     Europe, 20211222+[0H_96H]
>    >     summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
>    >     project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)
>    > ```
>    {: .code-out}
{: .hands_on}



> ### {% icon comment %} Be careful when sub-sampling!
> Increasing the frequency of your data e.g. artificialy creating data may not be scientifically relevant. Please use carefuly! Interpolating is not always scientifically relevant and sometimes you may prefer to choose a different method, for instance to take the nearest value:
> 
>    > ### {% icon code-in %} Input: Python
>    >  ```python
>    >  dset.resample(time='30min').nearest()
>    >  ```
>    {: .code-in}
>
{: .comment}



> ### {% icon question %} PM2.5 over Italy in the next 4 days
>
> Using a Multi-plot between Rome and Naples, and making averages per day, can you tell us if forecasted PM2.5 will increase or decrease?
>
> > ### {% icon solution %} Solution
> > PM2.5 over Italy is overall decreasing over the next 4 forecasted days.
> >
> > > ### {% icon code-in %} Input: Python
> > > ```python
> > > fig = plt.figure(1, figsize=[10,10])
> > > 
> > > # We're using cartopy to project our data.
> > > # (see documentation on cartopy)
> > > proj_plot = ccrs.Mercator()
> > > 
> > > sub_dset = dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).resample(time='1D').mean()
> > > # We need to project our data to the new projection and for this we use `transform`.
> > > # we set the original data projection in transform (here PlateCarree)
> > > p = sub_dset['pm2p5_conc'].plot(transform=ccrs.PlateCarree(), 
> > >                                 vmin = 0, vmax = 35,
> > >                                 subplot_kws={"projection": proj_plot},
> > >                                 col='time', col_wrap=5,
> > >                                 robust=True,
> > >                                 aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
> > >                                 cmap=cmc.roma_r)
> > > # We have to set the map's options on all axes
> > > for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
> > >     ax.coastlines('10m')
> > >     ax.set_title("CAMS PM2.5 " + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
> > > # Save your figure
> > > plt.savefig("CAMS-PM2_5-fc-multi-Italy-mean-per-day.png")
> > > ```
> > {: .code-in}
> > 
> > >  ![Daily mean for PM2.5 over Italy](../../images/CAMS-PM2_5-fc-multi-Italy-mean-per-day.png)
> {: .solution }
{: .question }




> ### {% icon comment %} `Grouby` versus `resample`
>
> Use `groupby` instead of `resample` when you wish to group over a dimension that is not `time`. `groupby` is very similar to resample but can be applied to any coordinates and not only to time. 
>
{: .comment}
