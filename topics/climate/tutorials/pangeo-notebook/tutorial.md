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

## Import Python packages

> ### {% icon code-in %} Input: Python
>    ```python
>    import numpy as np
>    import xarray as xr
>    import cartopy.crs as ccrs
>    import matplotlib.pyplot as plt
>    import cmcrameri.cm as cmc
>    import pandas as pd
>    ```
{: .hands_on}

## Open and read metadata


```python
dset = xr.open_dataset("CAMS-PM2_5-20211222.netcdf")
```

### Getting metadata

#### Global metadata


```python
dset
```


#### Metadata for a particular variable or coordinate variable


```python
print(dset.time)
```

    <xarray.DataArray 'time' (time: 97)>
    array([              0,   3600000000000,   7200000000000,  10800000000000,
            14400000000000,  18000000000000,  21600000000000,  25200000000000,
            28800000000000,  32400000000000,  36000000000000,  39600000000000,
            43200000000000,  46800000000000,  50400000000000,  54000000000000,
            57600000000000,  61200000000000,  64800000000000,  68400000000000,
            72000000000000,  75600000000000,  79200000000000,  82800000000000,
            86400000000000,  90000000000000,  93600000000000,  97200000000000,
           100800000000000, 104400000000000, 108000000000000, 111600000000000,
           115200000000000, 118800000000000, 122400000000000, 126000000000000,
           129600000000000, 133200000000000, 136800000000000, 140400000000000,
           144000000000000, 147600000000000, 151200000000000, 154800000000000,
           158400000000000, 162000000000000, 165600000000000, 169200000000000,
           172800000000000, 176400000000000, 180000000000000, 183600000000000,
           187200000000000, 190800000000000, 194400000000000, 198000000000000,
           201600000000000, 205200000000000, 208800000000000, 212400000000000,
           216000000000000, 219600000000000, 223200000000000, 226800000000000,
           230400000000000, 234000000000000, 237600000000000, 241200000000000,
           244800000000000, 248400000000000, 252000000000000, 255600000000000,
           259200000000000, 262800000000000, 266400000000000, 270000000000000,
           273600000000000, 277200000000000, 280800000000000, 284400000000000,
           288000000000000, 291600000000000, 295200000000000, 298800000000000,
           302400000000000, 306000000000000, 309600000000000, 313200000000000,
           316800000000000, 320400000000000, 324000000000000, 327600000000000,
           331200000000000, 334800000000000, 338400000000000, 342000000000000,
           345600000000000], dtype='timedelta64[ns]')
    Coordinates:
      * time     (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
    Attributes:
        long_name:  FORECAST time from 20211222


### Question 1: What is the name of the variable for Particle matter < 2.5 μm and its physical units?

- **Answer**: to get metadata, we use its variable name and print it. Printing it will only print metadata, not the values.
     - Variable name: mass_concentration_of_pm2p5_ambient_aerosol_in_air
     - Units: µg/m3


```python
print(dset.pm2p5_conc)
```

    <xarray.DataArray 'pm2p5_conc' (time: 97, level: 1, latitude: 400, longitude: 700)>
    [27160000 values with dtype=float32]
    Coordinates:
      * longitude  (longitude) float32 335.0 335.1 335.2 335.4 ... 44.75 44.85 44.95
      * latitude   (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
      * level      (level) float32 0.0
      * time       (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
    Attributes:
        species:        PM2.5 Aerosol
        units:          µg/m3
        value:          hourly values
        standard_name:  mass_concentration_of_pm2p5_ambient_aerosol_in_air


**Remark**: to access a variable or coordinate, we can use "." or specify its name as a string between squared brackets "[" "]". For example:


```python
print(dset['pm2p5_conc'])
```

    <xarray.DataArray 'pm2p5_conc' (time: 97, level: 1, latitude: 400, longitude: 700)>
    [27160000 values with dtype=float32]
    Coordinates:
      * longitude  (longitude) float32 335.0 335.1 335.2 335.4 ... 44.75 44.85 44.95
      * latitude   (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
      * level      (level) float32 0.0
      * time       (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
    Attributes:
        species:        PM2.5 Aerosol
        units:          µg/m3
        value:          hourly values
        standard_name:  mass_concentration_of_pm2p5_ambient_aerosol_in_air


## Select

### Select elements from coordinates by index


```python
dset.isel(time=0)
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (longitude: 700, latitude: 400, level: 1)
Coordinates:
  * longitude   (longitude) float32 335.0 335.1 335.2 ... 44.75 44.85 44.95
  * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
  * level       (level) float32 0.0
    time        timedelta64[ns] 00:00:00
Data variables:
    pm2p5_conc  (level, latitude, longitude) float32 ...
Attributes:
    title:        PM25 Air Pollutant FORECAST at the Surface
    institution:  Data produced by Meteo France
    source:       Data from ENSEMBLE model
    history:      Model ENSEMBLE FORECAST
    FORECAST:     Europe, 20211222+[0H_96H]
    summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
    project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-698fc86a-6f7f-4728-80ff-1957c5aa3f20' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-698fc86a-6f7f-4728-80ff-1957c5aa3f20' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>longitude</span>: 700</li><li><span class='xr-has-index'>latitude</span>: 400</li><li><span class='xr-has-index'>level</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-0ee9125c-43fe-40cc-beb8-35efe5b63ac4' class='xr-section-summary-in' type='checkbox'  checked><label for='section-0ee9125c-43fe-40cc-beb8-35efe5b63ac4' class='xr-section-summary' >Coordinates: <span>(4)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>335.0 335.1 335.2 ... 44.85 44.95</div><input id='attrs-005be0e7-161f-4bbe-9767-dd7ed519b4b2' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-005be0e7-161f-4bbe-9767-dd7ed519b4b2' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-2b009679-2bc4-47a1-8ec2-02f36e0badf2' class='xr-var-data-in' type='checkbox'><label for='data-2b009679-2bc4-47a1-8ec2-02f36e0badf2' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>longitude</dd><dt><span>units :</span></dt><dd>degrees_east</dd></dl></div><div class='xr-var-data'><pre>array([335.05, 335.15, 335.25, ...,  44.75,  44.85,  44.95], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>69.95 69.85 69.75 ... 30.15 30.05</div><input id='attrs-c0ad34f6-99f0-455a-9289-cb8643f07d92' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-c0ad34f6-99f0-455a-9289-cb8643f07d92' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-af0ca814-8622-48c3-a3ab-981b3b71674b' class='xr-var-data-in' type='checkbox'><label for='data-af0ca814-8622-48c3-a3ab-981b3b71674b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd></dl></div><div class='xr-var-data'><pre>array([69.95, 69.85, 69.75, ..., 30.25, 30.15, 30.05], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-9f1862fa-a31d-4fe3-b1f7-6d0c06cd1fd9' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-9f1862fa-a31d-4fe3-b1f7-6d0c06cd1fd9' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-f4f994b3-24e2-4278-86b0-963206bd700e' class='xr-var-data-in' type='checkbox'><label for='data-f4f994b3-24e2-4278-86b0-963206bd700e' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>time</span></div><div class='xr-var-dims'>()</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>00:00:00</div><input id='attrs-0d93c8c6-b3cd-4b5b-8304-19b5562a4e89' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-0d93c8c6-b3cd-4b5b-8304-19b5562a4e89' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-c84741bb-fa54-4490-a1ce-7732533585fa' class='xr-var-data-in' type='checkbox'><label for='data-c84741bb-fa54-4490-a1ce-7732533585fa' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array(0, dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-930daa26-fbac-444d-bc74-9f98d3ecbd24' class='xr-section-summary-in' type='checkbox'  checked><label for='section-930daa26-fbac-444d-bc74-9f98d3ecbd24' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(level, latitude, longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>...</div><input id='attrs-1493a990-9479-4eea-a1a9-20f10c1814e7' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-1493a990-9479-4eea-a1a9-20f10c1814e7' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-f1dd635f-4a38-43f9-917e-f87e32208380' class='xr-var-data-in' type='checkbox'><label for='data-f1dd635f-4a38-43f9-917e-f87e32208380' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>species :</span></dt><dd>PM2.5 Aerosol</dd><dt><span>units :</span></dt><dd>µg/m3</dd><dt><span>value :</span></dt><dd>hourly values</dd><dt><span>standard_name :</span></dt><dd>mass_concentration_of_pm2p5_ambient_aerosol_in_air</dd></dl></div><div class='xr-var-data'><pre>[280000 values with dtype=float32]</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-186a1ef0-48c6-4cfe-8e64-408183cb8340' class='xr-section-summary-in' type='checkbox'  checked><label for='section-186a1ef0-48c6-4cfe-8e64-408183cb8340' class='xr-section-summary' >Attributes: <span>(7)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>title :</span></dt><dd>PM25 Air Pollutant FORECAST at the Surface</dd><dt><span>institution :</span></dt><dd>Data produced by Meteo France</dd><dt><span>source :</span></dt><dd>Data from ENSEMBLE model</dd><dt><span>history :</span></dt><dd>Model ENSEMBLE FORECAST</dd><dt><span>FORECAST :</span></dt><dd>Europe, 20211222+[0H_96H]</dd><dt><span>summary :</span></dt><dd>ENSEMBLE model hourly FORECAST of PM25 concentration at the Surface from 20211222+[0H_96H] on Europe</dd><dt><span>project :</span></dt><dd>MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</dd></dl></div></li></ul></div></div>



### Select elements from coordinates by value
- We need to use the same type of time. Here `timedelta64`.


```python
dset.sel(time=np.timedelta64(0))
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (longitude: 700, latitude: 400, level: 1)
Coordinates:
  * longitude   (longitude) float32 335.0 335.1 335.2 ... 44.75 44.85 44.95
  * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
  * level       (level) float32 0.0
    time        timedelta64[ns] 00:00:00
Data variables:
    pm2p5_conc  (level, latitude, longitude) float32 ...
Attributes:
    title:        PM25 Air Pollutant FORECAST at the Surface
    institution:  Data produced by Meteo France
    source:       Data from ENSEMBLE model
    history:      Model ENSEMBLE FORECAST
    FORECAST:     Europe, 20211222+[0H_96H]
    summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
    project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-1303c599-7334-49f4-adc6-3d87206a2f50' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-1303c599-7334-49f4-adc6-3d87206a2f50' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>longitude</span>: 700</li><li><span class='xr-has-index'>latitude</span>: 400</li><li><span class='xr-has-index'>level</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-b4deb462-27de-48d5-9969-c2b1fe450072' class='xr-section-summary-in' type='checkbox'  checked><label for='section-b4deb462-27de-48d5-9969-c2b1fe450072' class='xr-section-summary' >Coordinates: <span>(4)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>335.0 335.1 335.2 ... 44.85 44.95</div><input id='attrs-7e73c332-9bfb-414d-b173-1b4c031073ec' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-7e73c332-9bfb-414d-b173-1b4c031073ec' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-6a41ac8b-6df5-492c-acbf-2e0fa67e2582' class='xr-var-data-in' type='checkbox'><label for='data-6a41ac8b-6df5-492c-acbf-2e0fa67e2582' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>longitude</dd><dt><span>units :</span></dt><dd>degrees_east</dd></dl></div><div class='xr-var-data'><pre>array([335.05, 335.15, 335.25, ...,  44.75,  44.85,  44.95], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>69.95 69.85 69.75 ... 30.15 30.05</div><input id='attrs-48064df2-b083-434f-9608-93c5e4413c92' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-48064df2-b083-434f-9608-93c5e4413c92' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-53bc8caa-1cca-421b-898b-1262872410f3' class='xr-var-data-in' type='checkbox'><label for='data-53bc8caa-1cca-421b-898b-1262872410f3' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd></dl></div><div class='xr-var-data'><pre>array([69.95, 69.85, 69.75, ..., 30.25, 30.15, 30.05], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-d250a0e9-89db-4f8b-98aa-a3d01f89deac' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-d250a0e9-89db-4f8b-98aa-a3d01f89deac' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-c080cd54-0788-4a51-833a-c0afc86674ec' class='xr-var-data-in' type='checkbox'><label for='data-c080cd54-0788-4a51-833a-c0afc86674ec' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>time</span></div><div class='xr-var-dims'>()</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>00:00:00</div><input id='attrs-a7ac7fb1-67a0-4456-9249-d60a301f5f98' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-a7ac7fb1-67a0-4456-9249-d60a301f5f98' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-9802c41c-74fe-446f-8065-0eb7b8362fd8' class='xr-var-data-in' type='checkbox'><label for='data-9802c41c-74fe-446f-8065-0eb7b8362fd8' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array(0, dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-97778c6e-425e-4439-9116-5d5d75c4c853' class='xr-section-summary-in' type='checkbox'  checked><label for='section-97778c6e-425e-4439-9116-5d5d75c4c853' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(level, latitude, longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>...</div><input id='attrs-66d49731-6812-4d0c-a6df-78885ff445ee' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-66d49731-6812-4d0c-a6df-78885ff445ee' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-65d18cab-266a-4bb0-857f-0f8e852f4f0b' class='xr-var-data-in' type='checkbox'><label for='data-65d18cab-266a-4bb0-857f-0f8e852f4f0b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>species :</span></dt><dd>PM2.5 Aerosol</dd><dt><span>units :</span></dt><dd>µg/m3</dd><dt><span>value :</span></dt><dd>hourly values</dd><dt><span>standard_name :</span></dt><dd>mass_concentration_of_pm2p5_ambient_aerosol_in_air</dd></dl></div><div class='xr-var-data'><pre>[280000 values with dtype=float32]</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-11fbbc2b-64e9-4b13-9e78-f5edff831e20' class='xr-section-summary-in' type='checkbox'  checked><label for='section-11fbbc2b-64e9-4b13-9e78-f5edff831e20' class='xr-section-summary' >Attributes: <span>(7)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>title :</span></dt><dd>PM25 Air Pollutant FORECAST at the Surface</dd><dt><span>institution :</span></dt><dd>Data produced by Meteo France</dd><dt><span>source :</span></dt><dd>Data from ENSEMBLE model</dd><dt><span>history :</span></dt><dd>Model ENSEMBLE FORECAST</dd><dt><span>FORECAST :</span></dt><dd>Europe, 20211222+[0H_96H]</dd><dt><span>summary :</span></dt><dd>ENSEMBLE model hourly FORECAST of PM25 concentration at the Surface from 20211222+[0H_96H] on Europe</dd><dt><span>project :</span></dt><dd>MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</dd></dl></div></li></ul></div></div>



## Question 2: Select forecast for December, 24th 2021 at 12:00 UTC

-  **Answer**: . The initial date/time is December, 22 2021 at 00:00 UTC. To select December, 24th 2021, you can use either isel or sel:


```python
dset.sel(time=(np.timedelta64(2,'D')+ np.timedelta64(12,'h')))
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (longitude: 700, latitude: 400, level: 1)
Coordinates:
  * longitude   (longitude) float32 335.0 335.1 335.2 ... 44.75 44.85 44.95
  * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
  * level       (level) float32 0.0
    time        timedelta64[ns] 2 days 12:00:00
Data variables:
    pm2p5_conc  (level, latitude, longitude) float32 ...
Attributes:
    title:        PM25 Air Pollutant FORECAST at the Surface
    institution:  Data produced by Meteo France
    source:       Data from ENSEMBLE model
    history:      Model ENSEMBLE FORECAST
    FORECAST:     Europe, 20211222+[0H_96H]
    summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
    project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-7cd13b21-e914-4feb-940b-955052bc396b' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-7cd13b21-e914-4feb-940b-955052bc396b' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>longitude</span>: 700</li><li><span class='xr-has-index'>latitude</span>: 400</li><li><span class='xr-has-index'>level</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-6d9cad50-89ed-469e-84b0-eb3eafe33e95' class='xr-section-summary-in' type='checkbox'  checked><label for='section-6d9cad50-89ed-469e-84b0-eb3eafe33e95' class='xr-section-summary' >Coordinates: <span>(4)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>335.0 335.1 335.2 ... 44.85 44.95</div><input id='attrs-a102bce6-0d76-45a9-b820-766b8a962910' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-a102bce6-0d76-45a9-b820-766b8a962910' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-c8e599f5-4773-4ba0-9186-382113b9fc2f' class='xr-var-data-in' type='checkbox'><label for='data-c8e599f5-4773-4ba0-9186-382113b9fc2f' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>longitude</dd><dt><span>units :</span></dt><dd>degrees_east</dd></dl></div><div class='xr-var-data'><pre>array([335.05, 335.15, 335.25, ...,  44.75,  44.85,  44.95], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>69.95 69.85 69.75 ... 30.15 30.05</div><input id='attrs-2aeac11d-8f95-4d57-8d0f-871e963dcda0' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-2aeac11d-8f95-4d57-8d0f-871e963dcda0' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-25c431c7-2543-4cb6-aeac-a751b713cfb3' class='xr-var-data-in' type='checkbox'><label for='data-25c431c7-2543-4cb6-aeac-a751b713cfb3' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd></dl></div><div class='xr-var-data'><pre>array([69.95, 69.85, 69.75, ..., 30.25, 30.15, 30.05], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-f894477f-0756-4384-b22d-a7ce3fcd75d2' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-f894477f-0756-4384-b22d-a7ce3fcd75d2' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-2c615a6a-a134-4112-acc7-f61975fb44ff' class='xr-var-data-in' type='checkbox'><label for='data-2c615a6a-a134-4112-acc7-f61975fb44ff' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>time</span></div><div class='xr-var-dims'>()</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>2 days 12:00:00</div><input id='attrs-096db6ee-c7af-4273-b13e-ab681bbf5918' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-096db6ee-c7af-4273-b13e-ab681bbf5918' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-acb1bbab-fcf1-4efb-953f-5e4050a94f19' class='xr-var-data-in' type='checkbox'><label for='data-acb1bbab-fcf1-4efb-953f-5e4050a94f19' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array(216000000000000, dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-4011faf5-7b45-476d-8d76-ea3338bac924' class='xr-section-summary-in' type='checkbox'  checked><label for='section-4011faf5-7b45-476d-8d76-ea3338bac924' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(level, latitude, longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>...</div><input id='attrs-a4a0d1b1-68d4-40ad-b7a7-bb87d0dca65c' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-a4a0d1b1-68d4-40ad-b7a7-bb87d0dca65c' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-31515a08-2865-43bc-8c78-510cd633d17c' class='xr-var-data-in' type='checkbox'><label for='data-31515a08-2865-43bc-8c78-510cd633d17c' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>species :</span></dt><dd>PM2.5 Aerosol</dd><dt><span>units :</span></dt><dd>µg/m3</dd><dt><span>value :</span></dt><dd>hourly values</dd><dt><span>standard_name :</span></dt><dd>mass_concentration_of_pm2p5_ambient_aerosol_in_air</dd></dl></div><div class='xr-var-data'><pre>[280000 values with dtype=float32]</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-fe7f2cb0-960b-4c0f-8834-45e5fc0cbd2c' class='xr-section-summary-in' type='checkbox'  checked><label for='section-fe7f2cb0-960b-4c0f-8834-45e5fc0cbd2c' class='xr-section-summary' >Attributes: <span>(7)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>title :</span></dt><dd>PM25 Air Pollutant FORECAST at the Surface</dd><dt><span>institution :</span></dt><dd>Data produced by Meteo France</dd><dt><span>source :</span></dt><dd>Data from ENSEMBLE model</dd><dt><span>history :</span></dt><dd>Model ENSEMBLE FORECAST</dd><dt><span>FORECAST :</span></dt><dd>Europe, 20211222+[0H_96H]</dd><dt><span>summary :</span></dt><dd>ENSEMBLE model hourly FORECAST of PM25 concentration at the Surface from 20211222+[0H_96H] on Europe</dd><dt><span>project :</span></dt><dd>MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</dd></dl></div></li></ul></div></div>



## Plotting
- To plot a map, you need to select a variable with data on geographical coordinates (latitude, longitude).
- In addition coordinates need to be sorted (preferably in increasing order). This is not the case of "longitude" because they are expressed as 0 to 360. Let's shift them to -180 to 180.


```python
dset.longitude
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.DataArray &#x27;longitude&#x27; (longitude: 700)&gt;
array([335.05, 335.15, 335.25, ...,  44.75,  44.85,  44.95], dtype=float32)
Coordinates:
  * longitude  (longitude) float32 335.0 335.1 335.2 335.4 ... 44.75 44.85 44.95
Attributes:
    long_name:  longitude
    units:      degrees_east</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.DataArray</div><div class='xr-array-name'>'longitude'</div><ul class='xr-dim-list'><li><span class='xr-has-index'>longitude</span>: 700</li></ul></div><ul class='xr-sections'><li class='xr-section-item'><div class='xr-array-wrap'><input id='section-e87e9d41-11bc-46da-9e8a-608b924bf223' class='xr-array-in' type='checkbox' checked><label for='section-e87e9d41-11bc-46da-9e8a-608b924bf223' title='Show/hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-array-preview xr-preview'><span>335.0 335.1 335.2 335.4 335.5 335.5 ... 44.55 44.65 44.75 44.85 44.95</span></div><div class='xr-array-data'><pre>array([335.05, 335.15, 335.25, ...,  44.75,  44.85,  44.95], dtype=float32)</pre></div></div></li><li class='xr-section-item'><input id='section-4bd450d5-e928-4bf1-9802-42c5dda49bc8' class='xr-section-summary-in' type='checkbox'  checked><label for='section-4bd450d5-e928-4bf1-9802-42c5dda49bc8' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>335.0 335.1 335.2 ... 44.85 44.95</div><input id='attrs-0a2dc9a6-f847-4585-9823-385994645506' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-0a2dc9a6-f847-4585-9823-385994645506' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-2ed6dec1-0f4b-412f-acf1-723f229dfe27' class='xr-var-data-in' type='checkbox'><label for='data-2ed6dec1-0f4b-412f-acf1-723f229dfe27' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>longitude</dd><dt><span>units :</span></dt><dd>degrees_east</dd></dl></div><div class='xr-var-data'><pre>array([335.05, 335.15, 335.25, ...,  44.75,  44.85,  44.95], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-e499416d-39fb-444f-a16c-e7a00b087fb0' class='xr-section-summary-in' type='checkbox'  checked><label for='section-e499416d-39fb-444f-a16c-e7a00b087fb0' class='xr-section-summary' >Attributes: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>longitude</dd><dt><span>units :</span></dt><dd>degrees_east</dd></dl></div></li></ul></div></div>




```python
dset.coords['longitude'] = (dset['longitude'] + 180) % 360 - 180
```


```python
dset.longitude
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.DataArray &#x27;longitude&#x27; (longitude: 700)&gt;
array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
        44.949997], dtype=float32)
Coordinates:
  * longitude  (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.DataArray</div><div class='xr-array-name'>'longitude'</div><ul class='xr-dim-list'><li><span class='xr-has-index'>longitude</span>: 700</li></ul></div><ul class='xr-sections'><li class='xr-section-item'><div class='xr-array-wrap'><input id='section-215a93e4-c852-448c-b7f4-886cc5e9ef86' class='xr-array-in' type='checkbox' checked><label for='section-215a93e4-c852-448c-b7f4-886cc5e9ef86' title='Show/hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-array-preview xr-preview'><span>-24.95 -24.85 -24.75 -24.65 -24.55 ... 44.55 44.65 44.75 44.85 44.95</span></div><div class='xr-array-data'><pre>array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
        44.949997], dtype=float32)</pre></div></div></li><li class='xr-section-item'><input id='section-c133c0cd-c468-411e-b305-8a0fb57694d4' class='xr-section-summary-in' type='checkbox'  checked><label for='section-c133c0cd-c468-411e-b305-8a0fb57694d4' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>-24.95 -24.85 ... 44.85 44.95</div><input id='attrs-298f32c3-7b83-4565-b697-ae9e73a70206' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-298f32c3-7b83-4565-b697-ae9e73a70206' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-7b0ba0aa-5b2d-4e7d-bba3-cc437b9af57b' class='xr-var-data-in' type='checkbox'><label for='data-7b0ba0aa-5b2d-4e7d-bba3-cc437b9af57b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
        44.949997], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-0aec563d-141f-49da-91cc-56b8f88417a6' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-0aec563d-141f-49da-91cc-56b8f88417a6' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>




```python
dset.sel(time=(np.timedelta64(2,'D')+ np.timedelta64(12,'h'))).pm2p5_conc.plot()
```




    <matplotlib.collections.QuadMesh at 0x7f62c58a1bb0>




    
![png](xarray-Galaxy_files/xarray-Galaxy_24_1.png)
    


**Remark**: Noe that in the previous plot, we did not need to select `level` because there is one value only. However, if we had more than one level, we would need to add a selection on the level before plotting

### Customize your plot
There are many ways to customize your plots and here we only give you what we think is important for creating publication ready figures:
- Define the size of the figure
- Choose to project data on a different projection.
- Add coastline
- Set the min and max values for plotting
- Add a title, change colorbar title
- Save figure into png or svg


```python
fig = plt.figure(1, figsize=[15,10])

# We're using cartopy to project our data.
# (see documentation on cartopy)
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator())
ax.coastlines(resolution='10m')

# We need to project our data to the new projection and for this we use `transform`.
# we set the original data projection in transform (here PlateCarree)
dset.sel(time=(np.timedelta64(2,'D') + np.timedelta64(12,'h')))['pm2p5_conc'].plot(ax=ax, 
                                                                                  transform=ccrs.PlateCarree(),
                                                                                  vmin = 0, vmax = 35,
                                                                                  cmap=cmc.roma_r)
# One way to customize your title
plt.title("Copernicus Monitoring Service PM2.5, 2 day forecasts\n 24th December 2021 at 12:00 UTC", fontsize=18)
plt.savefig("CAMS-PM2_5-fc-20211224.png")
```

    /srv/conda/envs/notebook/lib/python3.9/site-packages/cartopy/crs.py:825: ShapelyDeprecationWarning: __len__ for multi-part geometries is deprecated and will be removed in Shapely 2.0. Check the length of the `geoms` property instead to get the  number of parts of a multi-part geometry.
      if len(multi_line_string) > 1:
    /srv/conda/envs/notebook/lib/python3.9/site-packages/cartopy/crs.py:877: ShapelyDeprecationWarning: Iteration over multi-part geometries is deprecated and will be removed in Shapely 2.0. Use the `geoms` property to access the constituent parts of a multi-part geometry.
      for line in multi_line_string:
    /srv/conda/envs/notebook/lib/python3.9/site-packages/cartopy/crs.py:944: ShapelyDeprecationWarning: __len__ for multi-part geometries is deprecated and will be removed in Shapely 2.0. Check the length of the `geoms` property instead to get the  number of parts of a multi-part geometry.
      if len(p_mline) > 0:



    
![png](xarray-Galaxy_files/xarray-Galaxy_27_1.png)
    


### Multi-plots
- Here we would like to plot several times on the same figure in different sub-plots
- We will not plots all the times (too many) but the first 24 forecasted values.


```python
dset.time.shape
```




    (97,)



#### Make a list of times and convert to pandas datetime to make it easier to format times when plotting


```python
list_times = np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))
```


```python
pd.to_datetime(list_times).strftime("%d %b %H:%S UTC")
```




    Index(['22 Dec 00:00 UTC', '22 Dec 01:00 UTC', '22 Dec 02:00 UTC',
           '22 Dec 03:00 UTC', '22 Dec 04:00 UTC', '22 Dec 05:00 UTC',
           '22 Dec 06:00 UTC', '22 Dec 07:00 UTC', '22 Dec 08:00 UTC',
           '22 Dec 09:00 UTC', '22 Dec 10:00 UTC', '22 Dec 11:00 UTC',
           '22 Dec 12:00 UTC', '22 Dec 13:00 UTC', '22 Dec 14:00 UTC',
           '22 Dec 15:00 UTC', '22 Dec 16:00 UTC', '22 Dec 17:00 UTC',
           '22 Dec 18:00 UTC', '22 Dec 19:00 UTC', '22 Dec 20:00 UTC',
           '22 Dec 21:00 UTC', '22 Dec 22:00 UTC', '22 Dec 23:00 UTC',
           '23 Dec 00:00 UTC'],
          dtype='object')




```python
fig = plt.figure(1, figsize=[10,10])

# We're using cartopy to project our data.
# (see documentation on cartopy)
proj_plot = ccrs.Mercator()

# We need to project our data to the new projection and for this we use `transform`.
# we set the original data projection in transform (here PlateCarree)
p = dset.sel(time=slice(np.timedelta64(1,'h'),np.timedelta64(1,'D')))['pm2p5_conc'].plot(transform=ccrs.PlateCarree(),
                                                                                     vmin = 0, vmax = 35,
                                                                                     subplot_kws={"projection": proj_plot},
                                                                                     col='time', col_wrap=4,
                                                                                     robust=True,
                                                                                     aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
                                                                                     cmap=cmc.roma_r)
# We have to set the map's options on all four axes
for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
    ax.coastlines('10m')
    ax.set_title("CAMS PM2.5 " + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
# Save your figure
plt.savefig("CAMS-PM2_5-fc-multi.png")
```


    <Figure size 720x720 with 0 Axes>



    
![png](xarray-Galaxy_files/xarray-Galaxy_33_1.png)
    


## Question 3: Using a Multi-plot between Rome and Naples, can you tell us if the forecasted PM2.5 will increase or decrease during the first 24 hours?
- You can use https://geojson.io/ to find the coordinates for Italy. We suggest you use: sub-area 11. East to 15.0 East and 40. N to 43. N.

-  **Answer**: We will select a sub-area: 11. East to 15.0 East and 40. N to 43. N. PM2.5 will increase and reach values close to 35 μm.m-3


```python
fig = plt.figure(1, figsize=[10,10])

# We're using cartopy to project our data.
# (see documentation on cartopy)
proj_plot = ccrs.Mercator()

# We need to project our data to the new projection and for this we use `transform`.
# we set the original data projection in transform (here PlateCarree)
p = dset.sel(time=slice(np.timedelta64(1,'h'),np.timedelta64(1,'D'))).sel(latitude=slice(43., 40.),
                                                                          longitude=slice(11.,15.))['pm2p5_conc'].plot(transform=ccrs.PlateCarree(),
                                                                                     vmin = 0, vmax = 35,
                                                                                     subplot_kws={"projection": proj_plot},
                                                                                     col='time', col_wrap=4,
                                                                                     robust=True,
                                                                                     aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
                                                                                     cmap=cmc.roma_r)
# We have to set the map's options on all four axes
for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
    ax.coastlines('10m')
    ax.set_title("CAMS PM2.5 " + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
# Save your figure
plt.savefig("CAMS-PM2_5-fc-multi-Italy.png")
```


    <Figure size 720x720 with 0 Axes>



    
![png](xarray-Galaxy_files/xarray-Galaxy_35_1.png)
    


**Remark**: Why did we slice latitudes with `latitude=slice(47.3, 36.5)` and not `latitude=slice(36.5, 47.3)`?
- because when using slice, you need to specify values using the same order than in the coordinates. Latitudes are specified in decreasing order for CAMS.

## Where
- Sometimes we may want to make more complex selections with criteria on the values of a given variable and not only on its coordinates. For this we use `where`.
- For instance, we may want to only keep PM2.5 if they are greater than 25  μm.m-3


```python
dset.where(dset['pm2p5_conc'] > 25)
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (time: 97, level: 1, latitude: 400, longitude: 700)
Coordinates:
  * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
  * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
  * level       (level) float32 0.0
  * time        (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
Data variables:
    pm2p5_conc  (time, level, latitude, longitude) float32 nan nan ... nan nan
Attributes:
    title:        PM25 Air Pollutant FORECAST at the Surface
    institution:  Data produced by Meteo France
    source:       Data from ENSEMBLE model
    history:      Model ENSEMBLE FORECAST
    FORECAST:     Europe, 20211222+[0H_96H]
    summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
    project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-7d1e582c-5faf-486e-80b9-633d08b603a8' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-7d1e582c-5faf-486e-80b9-633d08b603a8' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 97</li><li><span class='xr-has-index'>level</span>: 1</li><li><span class='xr-has-index'>latitude</span>: 400</li><li><span class='xr-has-index'>longitude</span>: 700</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-deba531c-1ce4-411e-a745-325de4eea7f2' class='xr-section-summary-in' type='checkbox'  checked><label for='section-deba531c-1ce4-411e-a745-325de4eea7f2' class='xr-section-summary' >Coordinates: <span>(4)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>-24.95 -24.85 ... 44.85 44.95</div><input id='attrs-9411d03c-8c80-4d76-9340-d90c1732875e' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-9411d03c-8c80-4d76-9340-d90c1732875e' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-3abab2ba-acbb-40b5-9021-a18121e872d6' class='xr-var-data-in' type='checkbox'><label for='data-3abab2ba-acbb-40b5-9021-a18121e872d6' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
        44.949997], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>69.95 69.85 69.75 ... 30.15 30.05</div><input id='attrs-3cbab486-0966-466a-aec9-ee63c137e6f1' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-3cbab486-0966-466a-aec9-ee63c137e6f1' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-12f7f26e-a79f-49e8-96ae-f53258230a1c' class='xr-var-data-in' type='checkbox'><label for='data-12f7f26e-a79f-49e8-96ae-f53258230a1c' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd></dl></div><div class='xr-var-data'><pre>array([69.95, 69.85, 69.75, ..., 30.25, 30.15, 30.05], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-71dc2b93-978b-469e-a9f9-9c62e615d350' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-71dc2b93-978b-469e-a9f9-9c62e615d350' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-d14e5cde-32b0-4f6d-8980-7e3869518dcb' class='xr-var-data-in' type='checkbox'><label for='data-d14e5cde-32b0-4f6d-8980-7e3869518dcb' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>00:00:00 ... 4 days 00:00:00</div><input id='attrs-6fc71eec-ccfe-4ec0-9fe5-5931ae4fd7d5' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-6fc71eec-ccfe-4ec0-9fe5-5931ae4fd7d5' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-52fdda26-2b5c-4054-8a25-a7d28706d370' class='xr-var-data-in' type='checkbox'><label for='data-52fdda26-2b5c-4054-8a25-a7d28706d370' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array([              0,   3600000000000,   7200000000000,  10800000000000,
        14400000000000,  18000000000000,  21600000000000,  25200000000000,
        28800000000000,  32400000000000,  36000000000000,  39600000000000,
        43200000000000,  46800000000000,  50400000000000,  54000000000000,
        57600000000000,  61200000000000,  64800000000000,  68400000000000,
        72000000000000,  75600000000000,  79200000000000,  82800000000000,
        86400000000000,  90000000000000,  93600000000000,  97200000000000,
       100800000000000, 104400000000000, 108000000000000, 111600000000000,
       115200000000000, 118800000000000, 122400000000000, 126000000000000,
       129600000000000, 133200000000000, 136800000000000, 140400000000000,
       144000000000000, 147600000000000, 151200000000000, 154800000000000,
       158400000000000, 162000000000000, 165600000000000, 169200000000000,
       172800000000000, 176400000000000, 180000000000000, 183600000000000,
       187200000000000, 190800000000000, 194400000000000, 198000000000000,
       201600000000000, 205200000000000, 208800000000000, 212400000000000,
       216000000000000, 219600000000000, 223200000000000, 226800000000000,
       230400000000000, 234000000000000, 237600000000000, 241200000000000,
       244800000000000, 248400000000000, 252000000000000, 255600000000000,
       259200000000000, 262800000000000, 266400000000000, 270000000000000,
       273600000000000, 277200000000000, 280800000000000, 284400000000000,
       288000000000000, 291600000000000, 295200000000000, 298800000000000,
       302400000000000, 306000000000000, 309600000000000, 313200000000000,
       316800000000000, 320400000000000, 324000000000000, 327600000000000,
       331200000000000, 334800000000000, 338400000000000, 342000000000000,
       345600000000000], dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-37bd24d7-7d6d-4e34-aba2-9bf4438d02c3' class='xr-section-summary-in' type='checkbox'  checked><label for='section-37bd24d7-7d6d-4e34-aba2-9bf4438d02c3' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(time, level, latitude, longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>nan nan nan nan ... nan nan nan nan</div><input id='attrs-bfd494e6-7fe5-421b-b701-c6c5a35dfe63' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-bfd494e6-7fe5-421b-b701-c6c5a35dfe63' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-110ecbe9-1b92-4561-90eb-cfa4981cdba6' class='xr-var-data-in' type='checkbox'><label for='data-110ecbe9-1b92-4561-90eb-cfa4981cdba6' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>species :</span></dt><dd>PM2.5 Aerosol</dd><dt><span>units :</span></dt><dd>µg/m3</dd><dt><span>value :</span></dt><dd>hourly values</dd><dt><span>standard_name :</span></dt><dd>mass_concentration_of_pm2p5_ambient_aerosol_in_air</dd></dl></div><div class='xr-var-data'><pre>array([[[[nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         ...,
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan]]],


       [[[nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         ...,
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan]]],


       [[[nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
...
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan]]],


       [[[nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         ...,
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan]]],


       [[[nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         ...,
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan],
         [nan, nan, nan, ..., nan, nan, nan]]]], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-636d1877-4dde-457d-89b0-ee6a278cd65d' class='xr-section-summary-in' type='checkbox'  checked><label for='section-636d1877-4dde-457d-89b0-ee6a278cd65d' class='xr-section-summary' >Attributes: <span>(7)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>title :</span></dt><dd>PM25 Air Pollutant FORECAST at the Surface</dd><dt><span>institution :</span></dt><dd>Data produced by Meteo France</dd><dt><span>source :</span></dt><dd>Data from ENSEMBLE model</dd><dt><span>history :</span></dt><dd>Model ENSEMBLE FORECAST</dd><dt><span>FORECAST :</span></dt><dd>Europe, 20211222+[0H_96H]</dd><dt><span>summary :</span></dt><dd>ENSEMBLE model hourly FORECAST of PM25 concentration at the Surface from 20211222+[0H_96H] on Europe</dd><dt><span>project :</span></dt><dd>MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</dd></dl></div></li></ul></div></div>



**What happened?**
- You may not see any changes but if you look carefuly to `pm2p5_conc` values, you will see that many `nan`. In fact, we now have `nan` values Whenever the criteria within the `where` statement is not met e.g. when PM2.5 <= 25
- Let's plot one time to better see what happened:


```python
fig = plt.figure(1, figsize=[15,10])

# We're using cartopy to project our data.
# (see documentation on cartopy)
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator())
ax.coastlines(resolution='10m')

# We need to project our data to the new projection and for this we use `transform`.
# we set the original data projection in transform (here PlateCarree)
dset.where(dset['pm2p5_conc'] > 25).isel(time=0)['pm2p5_conc'].plot(ax=ax,
                                                                    transform=ccrs.PlateCarree(),
                                                                    vmin = 0, vmax = 35,
                                                                    cmap=cmc.roma_r)
# One way to customize your title
plt.title("Copernicus Monitoring Service PM2.5, 2 day forecasts\n 24th December 2021 at 12:00 UTC\n only values > 25", fontsize=18)
plt.savefig("CAMS-PM2_5-fc-20211224-25.png")
```


    
![png](xarray-Galaxy_files/xarray-Galaxy_40_0.png)
    



```python
fig = plt.figure(1, figsize=[10,10])

# We're using cartopy to project our data.
# (see documentation on cartopy)
proj_plot = ccrs.Mercator()

# We need to project our data to the new projection and for this we use `transform`.
# we set the original data projection in transform (here PlateCarree)
p = dset.where(dset['pm2p5_conc'] > 25).sel(time=slice(np.timedelta64(1,'h'),np.timedelta64(1,'D'))).sel(latitude=slice(43., 40.),
                                                                          longitude=slice(11.,15.))['pm2p5_conc'].plot(transform=ccrs.PlateCarree(),
                                                                                     vmin = 0, vmax = 35,
                                                                                     subplot_kws={"projection": proj_plot},
                                                                                     col='time', col_wrap=4,
                                                                                     robust=True,
                                                                                     aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
                                                                                     cmap=cmc.roma_r)
# We have to set the map's options on all four axes
for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
    ax.coastlines('10m')
    ax.set_title("PM2.5 > 25 μm.m-3" + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
# Save your figure
plt.savefig("CAMS-PM2_5-fc-multi-Italy-25.png")
```


    <Figure size 720x720 with 0 Axes>



    
![png](xarray-Galaxy_files/xarray-Galaxy_41_1.png)
    


## Reduction operations
- We often want to compute the mean of all our dataset, or along a dimension (for instance time).
- If you do not pass any argument to the operation then it is done over all dimension.


```python
dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).mean()
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     ()
Data variables:
    pm2p5_conc  float32 9.118</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-090f2bd5-c9bf-4113-bcfe-b493b5b00069' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-090f2bd5-c9bf-4113-bcfe-b493b5b00069' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-e0cdd071-9a6b-4af1-bd63-4b6f0827a872' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-e0cdd071-9a6b-4af1-bd63-4b6f0827a872' class='xr-section-summary'  title='Expand/collapse section'>Coordinates: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'></ul></div></li><li class='xr-section-item'><input id='section-7f33b641-3aa1-4d3b-aa23-e8ca1646a0fb' class='xr-section-summary-in' type='checkbox'  checked><label for='section-7f33b641-3aa1-4d3b-aa23-e8ca1646a0fb' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>()</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>9.118</div><input id='attrs-0d991ab8-86b3-4d6a-9ba4-981515eddc4d' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-0d991ab8-86b3-4d6a-9ba4-981515eddc4d' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-da741108-19f6-4ea6-a1f9-949d47379dc7' class='xr-var-data-in' type='checkbox'><label for='data-da741108-19f6-4ea6-a1f9-949d47379dc7' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array(9.118174, dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-21598a10-a9cf-43c3-97c0-7a0101e8ed74' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-21598a10-a9cf-43c3-97c0-7a0101e8ed74' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



## Question 4: What is the maximum forecasted PM2.5 value over the Rome-Naples region?

-  **Answer**: We will select a sub-area: 11. East to 15.0 East and 40. N to 43. N. The maximum PM2.5 value is 59.13694382 μm.m-3


```python
dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).max()
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     ()
Data variables:
    pm2p5_conc  float64 59.14</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-8768bfe3-7c55-4c09-b41f-04e9f9762534' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-8768bfe3-7c55-4c09-b41f-04e9f9762534' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-6bd19da8-c908-47b4-a9e6-e1d66b5a6aee' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-6bd19da8-c908-47b4-a9e6-e1d66b5a6aee' class='xr-section-summary'  title='Expand/collapse section'>Coordinates: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'></ul></div></li><li class='xr-section-item'><input id='section-736ea105-2f8b-49e7-807b-77b6af703e2e' class='xr-section-summary-in' type='checkbox'  checked><label for='section-736ea105-2f8b-49e7-807b-77b6af703e2e' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>()</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>59.14</div><input id='attrs-2e3a751b-9981-4467-93f9-cea4f63409ad' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-2e3a751b-9981-4467-93f9-cea4f63409ad' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-08b6a02a-2008-4f1a-86d0-26aa145cac6e' class='xr-var-data-in' type='checkbox'><label for='data-08b6a02a-2008-4f1a-86d0-26aa145cac6e' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array(59.13694382)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-8121eb7d-c074-4bd8-8316-756eb6db1af0' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-8121eb7d-c074-4bd8-8316-756eb6db1af0' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



## Question 5: When is forecasted the maximum PM2.5 value?

-  **Answer**: We will select a sub-area: 11. East to 15.0 East and 40. N to 43. N and average over the entire selected area and search where the maximum PM2.5 value of 59.13694382 μm.m-3 is found. The maximum PM2.5 value occurs on 2021-12-22 at 20:00 UTC.


```python
dset_tmean = dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).max(dim=('latitude', 'longitude'))
dset_tmean
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (level: 1, time: 97)
Coordinates:
  * level       (level) float32 0.0
  * time        (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
Data variables:
    pm2p5_conc  (time, level) float32 42.3 38.29 30.49 ... 15.63 14.94 15.39</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-713b2517-e991-430f-918e-5067e04f3fe0' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-713b2517-e991-430f-918e-5067e04f3fe0' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>level</span>: 1</li><li><span class='xr-has-index'>time</span>: 97</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-1a9742d4-2c2e-4f05-a9e2-b92b07ed2991' class='xr-section-summary-in' type='checkbox'  checked><label for='section-1a9742d4-2c2e-4f05-a9e2-b92b07ed2991' class='xr-section-summary' >Coordinates: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-69005bf0-34b1-433d-9a6f-8443e84da5e8' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-69005bf0-34b1-433d-9a6f-8443e84da5e8' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-dd4ab07b-4613-45c4-b982-967080e43b77' class='xr-var-data-in' type='checkbox'><label for='data-dd4ab07b-4613-45c4-b982-967080e43b77' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>00:00:00 ... 4 days 00:00:00</div><input id='attrs-49c611ef-062a-4ea5-b6d8-426ea84f72e8' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-49c611ef-062a-4ea5-b6d8-426ea84f72e8' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-09d761aa-8a35-4e07-bc7d-b43794e2552e' class='xr-var-data-in' type='checkbox'><label for='data-09d761aa-8a35-4e07-bc7d-b43794e2552e' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array([              0,   3600000000000,   7200000000000,  10800000000000,
        14400000000000,  18000000000000,  21600000000000,  25200000000000,
        28800000000000,  32400000000000,  36000000000000,  39600000000000,
        43200000000000,  46800000000000,  50400000000000,  54000000000000,
        57600000000000,  61200000000000,  64800000000000,  68400000000000,
        72000000000000,  75600000000000,  79200000000000,  82800000000000,
        86400000000000,  90000000000000,  93600000000000,  97200000000000,
       100800000000000, 104400000000000, 108000000000000, 111600000000000,
       115200000000000, 118800000000000, 122400000000000, 126000000000000,
       129600000000000, 133200000000000, 136800000000000, 140400000000000,
       144000000000000, 147600000000000, 151200000000000, 154800000000000,
       158400000000000, 162000000000000, 165600000000000, 169200000000000,
       172800000000000, 176400000000000, 180000000000000, 183600000000000,
       187200000000000, 190800000000000, 194400000000000, 198000000000000,
       201600000000000, 205200000000000, 208800000000000, 212400000000000,
       216000000000000, 219600000000000, 223200000000000, 226800000000000,
       230400000000000, 234000000000000, 237600000000000, 241200000000000,
       244800000000000, 248400000000000, 252000000000000, 255600000000000,
       259200000000000, 262800000000000, 266400000000000, 270000000000000,
       273600000000000, 277200000000000, 280800000000000, 284400000000000,
       288000000000000, 291600000000000, 295200000000000, 298800000000000,
       302400000000000, 306000000000000, 309600000000000, 313200000000000,
       316800000000000, 320400000000000, 324000000000000, 327600000000000,
       331200000000000, 334800000000000, 338400000000000, 342000000000000,
       345600000000000], dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-34892255-fb21-4a0a-974d-e66311e2b9e3' class='xr-section-summary-in' type='checkbox'  checked><label for='section-34892255-fb21-4a0a-974d-e66311e2b9e3' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(time, level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>42.3 38.29 30.49 ... 14.94 15.39</div><input id='attrs-22a592ec-40ed-499c-851f-c0e81324ac4a' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-22a592ec-40ed-499c-851f-c0e81324ac4a' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-b5d951aa-fd79-46c1-9f83-8378e4ce1425' class='xr-var-data-in' type='checkbox'><label for='data-b5d951aa-fd79-46c1-9f83-8378e4ce1425' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[42.296505 ],
       [38.28903  ],
       [30.48666  ],
       [28.831423 ],
       [28.261066 ],
       [26.068716 ],
       [34.163006 ],
       [34.15143  ],
       [37.76875  ],
       [38.484173 ],
       [34.02689  ],
       [30.031181 ],
       [29.249487 ],
       [28.92403  ],
       [30.3602   ],
       [36.182774 ],
       [39.95514  ],
       [42.54687  ],
       [47.64081  ],
       [52.986134 ],
...
       [16.213467 ],
       [16.752007 ],
       [17.096655 ],
       [17.479673 ],
       [17.755537 ],
       [17.728249 ],
       [17.551815 ],
       [17.511532 ],
       [16.69771  ],
       [16.567703 ],
       [16.201101 ],
       [16.107779 ],
       [16.097507 ],
       [16.159868 ],
       [17.916952 ],
       [17.410202 ],
       [17.213158 ],
       [15.629351 ],
       [14.938957 ],
       [15.392316 ]], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-58f4c80c-4233-400e-8ff4-3d58c73622e5' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-58f4c80c-4233-400e-8ff4-3d58c73622e5' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>




```python
np.datetime64('2021-12-22') + dset_tmean.time
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.DataArray &#x27;time&#x27; (time: 97)&gt;
array([&#x27;2021-12-22T00:00:00.000000000&#x27;, &#x27;2021-12-22T01:00:00.000000000&#x27;,
       &#x27;2021-12-22T02:00:00.000000000&#x27;, &#x27;2021-12-22T03:00:00.000000000&#x27;,
       &#x27;2021-12-22T04:00:00.000000000&#x27;, &#x27;2021-12-22T05:00:00.000000000&#x27;,
       &#x27;2021-12-22T06:00:00.000000000&#x27;, &#x27;2021-12-22T07:00:00.000000000&#x27;,
       &#x27;2021-12-22T08:00:00.000000000&#x27;, &#x27;2021-12-22T09:00:00.000000000&#x27;,
       &#x27;2021-12-22T10:00:00.000000000&#x27;, &#x27;2021-12-22T11:00:00.000000000&#x27;,
       &#x27;2021-12-22T12:00:00.000000000&#x27;, &#x27;2021-12-22T13:00:00.000000000&#x27;,
       &#x27;2021-12-22T14:00:00.000000000&#x27;, &#x27;2021-12-22T15:00:00.000000000&#x27;,
       &#x27;2021-12-22T16:00:00.000000000&#x27;, &#x27;2021-12-22T17:00:00.000000000&#x27;,
       &#x27;2021-12-22T18:00:00.000000000&#x27;, &#x27;2021-12-22T19:00:00.000000000&#x27;,
       &#x27;2021-12-22T20:00:00.000000000&#x27;, &#x27;2021-12-22T21:00:00.000000000&#x27;,
       &#x27;2021-12-22T22:00:00.000000000&#x27;, &#x27;2021-12-22T23:00:00.000000000&#x27;,
       &#x27;2021-12-23T00:00:00.000000000&#x27;, &#x27;2021-12-23T01:00:00.000000000&#x27;,
       &#x27;2021-12-23T02:00:00.000000000&#x27;, &#x27;2021-12-23T03:00:00.000000000&#x27;,
       &#x27;2021-12-23T04:00:00.000000000&#x27;, &#x27;2021-12-23T05:00:00.000000000&#x27;,
       &#x27;2021-12-23T06:00:00.000000000&#x27;, &#x27;2021-12-23T07:00:00.000000000&#x27;,
       &#x27;2021-12-23T08:00:00.000000000&#x27;, &#x27;2021-12-23T09:00:00.000000000&#x27;,
       &#x27;2021-12-23T10:00:00.000000000&#x27;, &#x27;2021-12-23T11:00:00.000000000&#x27;,
       &#x27;2021-12-23T12:00:00.000000000&#x27;, &#x27;2021-12-23T13:00:00.000000000&#x27;,
       &#x27;2021-12-23T14:00:00.000000000&#x27;, &#x27;2021-12-23T15:00:00.000000000&#x27;,
...
       &#x27;2021-12-24T10:00:00.000000000&#x27;, &#x27;2021-12-24T11:00:00.000000000&#x27;,
       &#x27;2021-12-24T12:00:00.000000000&#x27;, &#x27;2021-12-24T13:00:00.000000000&#x27;,
       &#x27;2021-12-24T14:00:00.000000000&#x27;, &#x27;2021-12-24T15:00:00.000000000&#x27;,
       &#x27;2021-12-24T16:00:00.000000000&#x27;, &#x27;2021-12-24T17:00:00.000000000&#x27;,
       &#x27;2021-12-24T18:00:00.000000000&#x27;, &#x27;2021-12-24T19:00:00.000000000&#x27;,
       &#x27;2021-12-24T20:00:00.000000000&#x27;, &#x27;2021-12-24T21:00:00.000000000&#x27;,
       &#x27;2021-12-24T22:00:00.000000000&#x27;, &#x27;2021-12-24T23:00:00.000000000&#x27;,
       &#x27;2021-12-25T00:00:00.000000000&#x27;, &#x27;2021-12-25T01:00:00.000000000&#x27;,
       &#x27;2021-12-25T02:00:00.000000000&#x27;, &#x27;2021-12-25T03:00:00.000000000&#x27;,
       &#x27;2021-12-25T04:00:00.000000000&#x27;, &#x27;2021-12-25T05:00:00.000000000&#x27;,
       &#x27;2021-12-25T06:00:00.000000000&#x27;, &#x27;2021-12-25T07:00:00.000000000&#x27;,
       &#x27;2021-12-25T08:00:00.000000000&#x27;, &#x27;2021-12-25T09:00:00.000000000&#x27;,
       &#x27;2021-12-25T10:00:00.000000000&#x27;, &#x27;2021-12-25T11:00:00.000000000&#x27;,
       &#x27;2021-12-25T12:00:00.000000000&#x27;, &#x27;2021-12-25T13:00:00.000000000&#x27;,
       &#x27;2021-12-25T14:00:00.000000000&#x27;, &#x27;2021-12-25T15:00:00.000000000&#x27;,
       &#x27;2021-12-25T16:00:00.000000000&#x27;, &#x27;2021-12-25T17:00:00.000000000&#x27;,
       &#x27;2021-12-25T18:00:00.000000000&#x27;, &#x27;2021-12-25T19:00:00.000000000&#x27;,
       &#x27;2021-12-25T20:00:00.000000000&#x27;, &#x27;2021-12-25T21:00:00.000000000&#x27;,
       &#x27;2021-12-25T22:00:00.000000000&#x27;, &#x27;2021-12-25T23:00:00.000000000&#x27;,
       &#x27;2021-12-26T00:00:00.000000000&#x27;], dtype=&#x27;datetime64[ns]&#x27;)
Coordinates:
  * time     (time) timedelta64[ns] 00:00:00 01:00:00 ... 4 days 00:00:00
Attributes:
    long_name:  FORECAST time from 20211222</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.DataArray</div><div class='xr-array-name'>'time'</div><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 97</li></ul></div><ul class='xr-sections'><li class='xr-section-item'><div class='xr-array-wrap'><input id='section-dc03091d-66ac-4449-a5eb-b32ae2e9c230' class='xr-array-in' type='checkbox' checked><label for='section-dc03091d-66ac-4449-a5eb-b32ae2e9c230' title='Show/hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-array-preview xr-preview'><span>2021-12-22 2021-12-22T01:00:00 ... 2021-12-25T23:00:00 2021-12-26</span></div><div class='xr-array-data'><pre>array([&#x27;2021-12-22T00:00:00.000000000&#x27;, &#x27;2021-12-22T01:00:00.000000000&#x27;,
       &#x27;2021-12-22T02:00:00.000000000&#x27;, &#x27;2021-12-22T03:00:00.000000000&#x27;,
       &#x27;2021-12-22T04:00:00.000000000&#x27;, &#x27;2021-12-22T05:00:00.000000000&#x27;,
       &#x27;2021-12-22T06:00:00.000000000&#x27;, &#x27;2021-12-22T07:00:00.000000000&#x27;,
       &#x27;2021-12-22T08:00:00.000000000&#x27;, &#x27;2021-12-22T09:00:00.000000000&#x27;,
       &#x27;2021-12-22T10:00:00.000000000&#x27;, &#x27;2021-12-22T11:00:00.000000000&#x27;,
       &#x27;2021-12-22T12:00:00.000000000&#x27;, &#x27;2021-12-22T13:00:00.000000000&#x27;,
       &#x27;2021-12-22T14:00:00.000000000&#x27;, &#x27;2021-12-22T15:00:00.000000000&#x27;,
       &#x27;2021-12-22T16:00:00.000000000&#x27;, &#x27;2021-12-22T17:00:00.000000000&#x27;,
       &#x27;2021-12-22T18:00:00.000000000&#x27;, &#x27;2021-12-22T19:00:00.000000000&#x27;,
       &#x27;2021-12-22T20:00:00.000000000&#x27;, &#x27;2021-12-22T21:00:00.000000000&#x27;,
       &#x27;2021-12-22T22:00:00.000000000&#x27;, &#x27;2021-12-22T23:00:00.000000000&#x27;,
       &#x27;2021-12-23T00:00:00.000000000&#x27;, &#x27;2021-12-23T01:00:00.000000000&#x27;,
       &#x27;2021-12-23T02:00:00.000000000&#x27;, &#x27;2021-12-23T03:00:00.000000000&#x27;,
       &#x27;2021-12-23T04:00:00.000000000&#x27;, &#x27;2021-12-23T05:00:00.000000000&#x27;,
       &#x27;2021-12-23T06:00:00.000000000&#x27;, &#x27;2021-12-23T07:00:00.000000000&#x27;,
       &#x27;2021-12-23T08:00:00.000000000&#x27;, &#x27;2021-12-23T09:00:00.000000000&#x27;,
       &#x27;2021-12-23T10:00:00.000000000&#x27;, &#x27;2021-12-23T11:00:00.000000000&#x27;,
       &#x27;2021-12-23T12:00:00.000000000&#x27;, &#x27;2021-12-23T13:00:00.000000000&#x27;,
       &#x27;2021-12-23T14:00:00.000000000&#x27;, &#x27;2021-12-23T15:00:00.000000000&#x27;,
...
       &#x27;2021-12-24T10:00:00.000000000&#x27;, &#x27;2021-12-24T11:00:00.000000000&#x27;,
       &#x27;2021-12-24T12:00:00.000000000&#x27;, &#x27;2021-12-24T13:00:00.000000000&#x27;,
       &#x27;2021-12-24T14:00:00.000000000&#x27;, &#x27;2021-12-24T15:00:00.000000000&#x27;,
       &#x27;2021-12-24T16:00:00.000000000&#x27;, &#x27;2021-12-24T17:00:00.000000000&#x27;,
       &#x27;2021-12-24T18:00:00.000000000&#x27;, &#x27;2021-12-24T19:00:00.000000000&#x27;,
       &#x27;2021-12-24T20:00:00.000000000&#x27;, &#x27;2021-12-24T21:00:00.000000000&#x27;,
       &#x27;2021-12-24T22:00:00.000000000&#x27;, &#x27;2021-12-24T23:00:00.000000000&#x27;,
       &#x27;2021-12-25T00:00:00.000000000&#x27;, &#x27;2021-12-25T01:00:00.000000000&#x27;,
       &#x27;2021-12-25T02:00:00.000000000&#x27;, &#x27;2021-12-25T03:00:00.000000000&#x27;,
       &#x27;2021-12-25T04:00:00.000000000&#x27;, &#x27;2021-12-25T05:00:00.000000000&#x27;,
       &#x27;2021-12-25T06:00:00.000000000&#x27;, &#x27;2021-12-25T07:00:00.000000000&#x27;,
       &#x27;2021-12-25T08:00:00.000000000&#x27;, &#x27;2021-12-25T09:00:00.000000000&#x27;,
       &#x27;2021-12-25T10:00:00.000000000&#x27;, &#x27;2021-12-25T11:00:00.000000000&#x27;,
       &#x27;2021-12-25T12:00:00.000000000&#x27;, &#x27;2021-12-25T13:00:00.000000000&#x27;,
       &#x27;2021-12-25T14:00:00.000000000&#x27;, &#x27;2021-12-25T15:00:00.000000000&#x27;,
       &#x27;2021-12-25T16:00:00.000000000&#x27;, &#x27;2021-12-25T17:00:00.000000000&#x27;,
       &#x27;2021-12-25T18:00:00.000000000&#x27;, &#x27;2021-12-25T19:00:00.000000000&#x27;,
       &#x27;2021-12-25T20:00:00.000000000&#x27;, &#x27;2021-12-25T21:00:00.000000000&#x27;,
       &#x27;2021-12-25T22:00:00.000000000&#x27;, &#x27;2021-12-25T23:00:00.000000000&#x27;,
       &#x27;2021-12-26T00:00:00.000000000&#x27;], dtype=&#x27;datetime64[ns]&#x27;)</pre></div></div></li><li class='xr-section-item'><input id='section-85e9a235-4094-4f44-956c-dc18065f61f9' class='xr-section-summary-in' type='checkbox'  checked><label for='section-85e9a235-4094-4f44-956c-dc18065f61f9' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>00:00:00 ... 4 days 00:00:00</div><input id='attrs-248899d2-1ea1-407d-948f-241565392372' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-248899d2-1ea1-407d-948f-241565392372' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-68f4b6b6-71c9-4472-841e-750ee4e8b33d' class='xr-var-data-in' type='checkbox'><label for='data-68f4b6b6-71c9-4472-841e-750ee4e8b33d' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array([              0,   3600000000000,   7200000000000,  10800000000000,
        14400000000000,  18000000000000,  21600000000000,  25200000000000,
        28800000000000,  32400000000000,  36000000000000,  39600000000000,
        43200000000000,  46800000000000,  50400000000000,  54000000000000,
        57600000000000,  61200000000000,  64800000000000,  68400000000000,
        72000000000000,  75600000000000,  79200000000000,  82800000000000,
        86400000000000,  90000000000000,  93600000000000,  97200000000000,
       100800000000000, 104400000000000, 108000000000000, 111600000000000,
       115200000000000, 118800000000000, 122400000000000, 126000000000000,
       129600000000000, 133200000000000, 136800000000000, 140400000000000,
       144000000000000, 147600000000000, 151200000000000, 154800000000000,
       158400000000000, 162000000000000, 165600000000000, 169200000000000,
       172800000000000, 176400000000000, 180000000000000, 183600000000000,
       187200000000000, 190800000000000, 194400000000000, 198000000000000,
       201600000000000, 205200000000000, 208800000000000, 212400000000000,
       216000000000000, 219600000000000, 223200000000000, 226800000000000,
       230400000000000, 234000000000000, 237600000000000, 241200000000000,
       244800000000000, 248400000000000, 252000000000000, 255600000000000,
       259200000000000, 262800000000000, 266400000000000, 270000000000000,
       273600000000000, 277200000000000, 280800000000000, 284400000000000,
       288000000000000, 291600000000000, 295200000000000, 298800000000000,
       302400000000000, 306000000000000, 309600000000000, 313200000000000,
       316800000000000, 320400000000000, 324000000000000, 327600000000000,
       331200000000000, 334800000000000, 338400000000000, 342000000000000,
       345600000000000], dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-23786b62-946d-47a8-936e-397612ddb2c4' class='xr-section-summary-in' type='checkbox'  checked><label for='section-23786b62-946d-47a8-936e-397612ddb2c4' class='xr-section-summary' >Attributes: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div></li></ul></div></div>




```python
dset_tmean_max = dset_tmean.where(dset_tmean['pm2p5_conc'] == 59.13694382, drop=True)
dset_tmean_max
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (time: 1, level: 1)
Coordinates:
  * level       (level) float32 0.0
  * time        (time) timedelta64[ns] 20:00:00
Data variables:
    pm2p5_conc  (time, level) float32 59.14</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-ead6b4b8-a54b-47b8-88b9-ef117643f60d' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-ead6b4b8-a54b-47b8-88b9-ef117643f60d' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 1</li><li><span class='xr-has-index'>level</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-9bedc97b-b8c2-4c00-9775-3d38e432fc6f' class='xr-section-summary-in' type='checkbox'  checked><label for='section-9bedc97b-b8c2-4c00-9775-3d38e432fc6f' class='xr-section-summary' >Coordinates: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-d83ea6a0-7117-4b65-98b3-53ac27302006' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-d83ea6a0-7117-4b65-98b3-53ac27302006' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-23103088-6fc7-4c8a-befa-67e12b5ff948' class='xr-var-data-in' type='checkbox'><label for='data-23103088-6fc7-4c8a-befa-67e12b5ff948' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>20:00:00</div><input id='attrs-3910a1e1-edd2-4805-a73c-63aad1f73230' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-3910a1e1-edd2-4805-a73c-63aad1f73230' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-f4be306c-6989-4605-b70f-80f69b96bf92' class='xr-var-data-in' type='checkbox'><label for='data-f4be306c-6989-4605-b70f-80f69b96bf92' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array([72000000000000], dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-f2735ac6-6e6f-4759-bc71-5fbf0daf855f' class='xr-section-summary-in' type='checkbox'  checked><label for='section-f2735ac6-6e6f-4759-bc71-5fbf0daf855f' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(time, level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>59.14</div><input id='attrs-7320f179-eab3-4ad4-880b-8c0cbed2eae3' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-7320f179-eab3-4ad4-880b-8c0cbed2eae3' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-87bea97e-b845-4463-bbd6-3550d6e87a32' class='xr-var-data-in' type='checkbox'><label for='data-87bea97e-b845-4463-bbd6-3550d6e87a32' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[59.136944]], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-8f73f8a4-f74c-4bac-8cb7-d825b64faffd' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-8f73f8a4-f74c-4bac-8cb7-d825b64faffd' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



**Remark**: We average over a relatively small so we do not make a weighted average. Use weighted averages when averaging over the entire globe or over a large area where the pixel sizes may vary (depending on the latitude).

## Resample
- We often want to resample (on time dimension) our data:
    - if your resampling frequency is lower than your original data, you would need to apply an operation on the data you group together, such as mean, min, max;
    - if your resampling frequency is higher than your original data, you would need to indicate how to fill the gaps, for instance interpolate and indicate which interpolation method to apply or select nearest values, etc.


```python
dset.resample(time='1D').mean()
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (time: 5, longitude: 700, latitude: 400, level: 1)
Coordinates:
  * time        (time) timedelta64[ns] 0 days 1 days 2 days 3 days 4 days
  * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
  * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
  * level       (level) float32 0.0
Data variables:
    pm2p5_conc  (time, level, latitude, longitude) float32 0.4298 ... 7.501</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-e3e7e2d0-4208-41fc-99ca-f22cf3f209a4' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-e3e7e2d0-4208-41fc-99ca-f22cf3f209a4' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 5</li><li><span class='xr-has-index'>longitude</span>: 700</li><li><span class='xr-has-index'>latitude</span>: 400</li><li><span class='xr-has-index'>level</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-41d74fb3-3778-4dbd-a35b-0b282612c883' class='xr-section-summary-in' type='checkbox'  checked><label for='section-41d74fb3-3778-4dbd-a35b-0b282612c883' class='xr-section-summary' >Coordinates: <span>(4)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>0 days 1 days 2 days 3 days 4 days</div><input id='attrs-4f7a8571-1469-4e0c-9eff-9dcd0f7948d5' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-4f7a8571-1469-4e0c-9eff-9dcd0f7948d5' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-996f2d95-7c00-429e-a559-e830b64fb011' class='xr-var-data-in' type='checkbox'><label for='data-996f2d95-7c00-429e-a559-e830b64fb011' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([              0,  86400000000000, 172800000000000, 259200000000000,
       345600000000000], dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>-24.95 -24.85 ... 44.85 44.95</div><input id='attrs-25a692cf-2fa1-421c-adcc-c069f408b5f8' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-25a692cf-2fa1-421c-adcc-c069f408b5f8' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-20b4c384-47d9-4c3c-8fc3-7f8e5c43c8b4' class='xr-var-data-in' type='checkbox'><label for='data-20b4c384-47d9-4c3c-8fc3-7f8e5c43c8b4' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
        44.949997], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>69.95 69.85 69.75 ... 30.15 30.05</div><input id='attrs-2a2e74ea-21f5-4eb5-9bba-2669bcf3ffcf' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-2a2e74ea-21f5-4eb5-9bba-2669bcf3ffcf' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-d7602d6b-419e-4408-9113-84e853ff2bed' class='xr-var-data-in' type='checkbox'><label for='data-d7602d6b-419e-4408-9113-84e853ff2bed' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd></dl></div><div class='xr-var-data'><pre>array([69.95, 69.85, 69.75, ..., 30.25, 30.15, 30.05], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-57f3d1a1-01ce-4410-938d-d11b7cb95826' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-57f3d1a1-01ce-4410-938d-d11b7cb95826' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-5c61f7d0-6feb-4cb2-a6d6-127c726c74bc' class='xr-var-data-in' type='checkbox'><label for='data-5c61f7d0-6feb-4cb2-a6d6-127c726c74bc' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-12e290c3-2092-4d98-b5db-08d9cfd1bdef' class='xr-section-summary-in' type='checkbox'  checked><label for='section-12e290c3-2092-4d98-b5db-08d9cfd1bdef' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(time, level, latitude, longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.4298 0.4305 ... 7.329 7.501</div><input id='attrs-b7b5a427-c908-4a6b-834d-960e70595e11' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-b7b5a427-c908-4a6b-834d-960e70595e11' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-730a7e5f-d5f3-4b6e-a3c3-8592cb8a6aed' class='xr-var-data-in' type='checkbox'><label for='data-730a7e5f-d5f3-4b6e-a3c3-8592cb8a6aed' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([[[[ 0.42978784,  0.43048358,  0.42945918, ...,  1.3107316 ,
           1.2517792 ,  1.1982225 ],
         [ 0.4370046 ,  0.43546566,  0.4325211 , ...,  1.246305  ,
           1.1844547 ,  1.1439868 ],
         [ 0.43792644,  0.4366493 ,  0.43614247, ...,  1.1642256 ,
           1.1276717 ,  1.1353835 ],
         ...,
         [ 4.587089  ,  4.5268383 ,  4.511112  , ...,  8.950371  ,
           9.012288  ,  9.123664  ],
         [ 4.589078  ,  4.523675  ,  4.533488  , ...,  8.9588995 ,
           9.020508  ,  9.137817  ],
         [ 4.949938  ,  4.8722715 ,  4.8262234 , ...,  8.991719  ,
           9.740005  , 10.128865  ]]],


       [[[ 0.4713411 ,  0.47189176,  0.4671572 , ...,  1.6874713 ,
           1.6840715 ,  1.6898432 ],
         [ 0.48606887,  0.4831296 ,  0.4789001 , ...,  1.5788873 ,
           1.5715913 ,  1.6546102 ],
         [ 0.4892663 ,  0.48747218,  0.48571062, ...,  1.4654318 ,
...
          13.56842   , 13.506455  ],
         [ 4.36068   ,  4.3648267 ,  4.3712196 , ..., 13.8911085 ,
          13.727001  , 13.636907  ],
         [ 4.3599424 ,  4.3731194 ,  4.360737  , ..., 14.034665  ,
          13.876995  , 13.998788  ]]],


       [[[ 0.500001  ,  0.500001  ,  0.500001  , ...,  2.7360647 ,
           2.7907553 ,  2.8610635 ],
         [ 0.51731694,  0.51731694,  0.51731694, ...,  2.696999  ,
           2.7751305 ,  2.8610635 ],
         [ 0.500001  ,  0.500001  ,  0.500001  , ...,  2.5329418 ,
           2.6266909 ,  2.728256  ],
         ...,
         [ 3.9130645 ,  3.8914356 ,  3.8592055 , ...,  6.9082866 ,
           6.61145   ,  6.622023  ],
         [ 3.9357026 ,  3.7816854 ,  3.7221985 , ...,  7.2744718 ,
           7.0306134 ,  7.076543  ],
         [ 3.9141445 ,  3.7201238 ,  3.7152212 , ...,  7.353484  ,
           7.3293824 ,  7.501277  ]]]], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-d70f1912-5bec-4ba3-b4b5-7b196ed6d5b1' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-d70f1912-5bec-4ba3-b4b5-7b196ed6d5b1' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>




```python
dset.resample(time='30min').interpolate('linear')
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (longitude: 700, latitude: 400, level: 1, time: 193)
Coordinates:
  * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
  * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
  * level       (level) float32 0.0
  * time        (time) timedelta64[ns] 00:00:00 00:30:00 ... 4 days 00:00:00
Data variables:
    pm2p5_conc  (time, level, latitude, longitude) float64 0.4202 ... 7.501
Attributes:
    title:        PM25 Air Pollutant FORECAST at the Surface
    institution:  Data produced by Meteo France
    source:       Data from ENSEMBLE model
    history:      Model ENSEMBLE FORECAST
    FORECAST:     Europe, 20211222+[0H_96H]
    summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
    project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-625516ff-bf8e-43c9-a110-aef5d09c8396' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-625516ff-bf8e-43c9-a110-aef5d09c8396' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>longitude</span>: 700</li><li><span class='xr-has-index'>latitude</span>: 400</li><li><span class='xr-has-index'>level</span>: 1</li><li><span class='xr-has-index'>time</span>: 193</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-4eb0af58-9023-4f6e-a1e2-2ebe0990e223' class='xr-section-summary-in' type='checkbox'  checked><label for='section-4eb0af58-9023-4f6e-a1e2-2ebe0990e223' class='xr-section-summary' >Coordinates: <span>(4)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>-24.95 -24.85 ... 44.85 44.95</div><input id='attrs-1e80877f-3b37-4b87-827b-c8942651b549' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-1e80877f-3b37-4b87-827b-c8942651b549' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-0733905a-1a91-4e57-9e1f-4223c1e0e4b0' class='xr-var-data-in' type='checkbox'><label for='data-0733905a-1a91-4e57-9e1f-4223c1e0e4b0' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
        44.949997], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>69.95 69.85 69.75 ... 30.15 30.05</div><input id='attrs-920b5722-840b-4825-a9c4-287d22138cf3' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-920b5722-840b-4825-a9c4-287d22138cf3' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-1dd2a440-1bb5-4348-9a2a-2734f9200342' class='xr-var-data-in' type='checkbox'><label for='data-1dd2a440-1bb5-4348-9a2a-2734f9200342' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd></dl></div><div class='xr-var-data'><pre>array([69.95, 69.85, 69.75, ..., 30.25, 30.15, 30.05], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-e66f8108-8143-4e34-b0ff-6003aeb7cbd2' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-e66f8108-8143-4e34-b0ff-6003aeb7cbd2' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-928dcae7-e340-4d13-a150-63a3d02c6a31' class='xr-var-data-in' type='checkbox'><label for='data-928dcae7-e340-4d13-a150-63a3d02c6a31' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>00:00:00 ... 4 days 00:00:00</div><input id='attrs-ddf96300-3370-4d05-b969-5dbd2adbcc9b' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-ddf96300-3370-4d05-b969-5dbd2adbcc9b' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-1ec2de1f-c679-4a08-b8b2-d4190e425d01' class='xr-var-data-in' type='checkbox'><label for='data-1ec2de1f-c679-4a08-b8b2-d4190e425d01' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array([              0,   1800000000000,   3600000000000,   5400000000000,
         7200000000000,   9000000000000,  10800000000000,  12600000000000,
        14400000000000,  16200000000000,  18000000000000,  19800000000000,
        21600000000000,  23400000000000,  25200000000000,  27000000000000,
        28800000000000,  30600000000000,  32400000000000,  34200000000000,
        36000000000000,  37800000000000,  39600000000000,  41400000000000,
        43200000000000,  45000000000000,  46800000000000,  48600000000000,
        50400000000000,  52200000000000,  54000000000000,  55800000000000,
        57600000000000,  59400000000000,  61200000000000,  63000000000000,
        64800000000000,  66600000000000,  68400000000000,  70200000000000,
        72000000000000,  73800000000000,  75600000000000,  77400000000000,
        79200000000000,  81000000000000,  82800000000000,  84600000000000,
        86400000000000,  88200000000000,  90000000000000,  91800000000000,
        93600000000000,  95400000000000,  97200000000000,  99000000000000,
       100800000000000, 102600000000000, 104400000000000, 106200000000000,
       108000000000000, 109800000000000, 111600000000000, 113400000000000,
       115200000000000, 117000000000000, 118800000000000, 120600000000000,
       122400000000000, 124200000000000, 126000000000000, 127800000000000,
       129600000000000, 131400000000000, 133200000000000, 135000000000000,
       136800000000000, 138600000000000, 140400000000000, 142200000000000,
       144000000000000, 145800000000000, 147600000000000, 149400000000000,
       151200000000000, 153000000000000, 154800000000000, 156600000000000,
       158400000000000, 160200000000000, 162000000000000, 163800000000000,
       165600000000000, 167400000000000, 169200000000000, 171000000000000,
       172800000000000, 174600000000000, 176400000000000, 178200000000000,
       180000000000000, 181800000000000, 183600000000000, 185400000000000,
       187200000000000, 189000000000000, 190800000000000, 192600000000000,
       194400000000000, 196200000000000, 198000000000000, 199800000000000,
       201600000000000, 203400000000000, 205200000000000, 207000000000000,
       208800000000000, 210600000000000, 212400000000000, 214200000000000,
       216000000000000, 217800000000000, 219600000000000, 221400000000000,
       223200000000000, 225000000000000, 226800000000000, 228600000000000,
       230400000000000, 232200000000000, 234000000000000, 235800000000000,
       237600000000000, 239400000000000, 241200000000000, 243000000000000,
       244800000000000, 246600000000000, 248400000000000, 250200000000000,
       252000000000000, 253800000000000, 255600000000000, 257400000000000,
       259200000000000, 261000000000000, 262800000000000, 264600000000000,
       266400000000000, 268200000000000, 270000000000000, 271800000000000,
       273600000000000, 275400000000000, 277200000000000, 279000000000000,
       280800000000000, 282600000000000, 284400000000000, 286200000000000,
       288000000000000, 289800000000000, 291600000000000, 293400000000000,
       295200000000000, 297000000000000, 298800000000000, 300600000000000,
       302400000000000, 304200000000000, 306000000000000, 307800000000000,
       309600000000000, 311400000000000, 313200000000000, 315000000000000,
       316800000000000, 318600000000000, 320400000000000, 322200000000000,
       324000000000000, 325800000000000, 327600000000000, 329400000000000,
       331200000000000, 333000000000000, 334800000000000, 336600000000000,
       338400000000000, 340200000000000, 342000000000000, 343800000000000,
       345600000000000], dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-864e69b9-b84c-4e0e-92fb-00b55f7478e4' class='xr-section-summary-in' type='checkbox'  checked><label for='section-864e69b9-b84c-4e0e-92fb-00b55f7478e4' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(time, level, latitude, longitude)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.4202 0.4331 ... 7.329 7.501</div><input id='attrs-746e1cc4-6d06-4fd2-8241-122ad080d74c' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-746e1cc4-6d06-4fd2-8241-122ad080d74c' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-b37c5b0f-128b-41d2-a395-74fa5e44141a' class='xr-var-data-in' type='checkbox'><label for='data-b37c5b0f-128b-41d2-a395-74fa5e44141a' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>species :</span></dt><dd>PM2.5 Aerosol</dd><dt><span>units :</span></dt><dd>µg/m3</dd><dt><span>value :</span></dt><dd>hourly values</dd><dt><span>standard_name :</span></dt><dd>mass_concentration_of_pm2p5_ambient_aerosol_in_air</dd></dl></div><div class='xr-var-data'><pre>array([[[[ 0.42024356,  0.43306175,  0.43306175, ...,  1.12402189,
           1.11009526,  1.17111671],
         [ 0.42376786,  0.42228994,  0.42031461, ...,  1.07814932,
           1.07138491,  1.04034841],
         [ 0.43263543,  0.42045674,  0.42055619, ...,  1.0738008 ,
           1.06749117,  1.04071796],
         ...,
         [ 4.76865149,  4.74206305,  4.68561745, ..., 13.32007599,
          13.32007599, 13.32007599],
         [ 4.76052284,  4.59739637,  4.59644413, ..., 13.66893768,
          13.85706139, 14.00489712],
         [ 4.72977066,  4.64504528,  4.59644413, ..., 13.73429394,
          14.11473274, 14.21715069]]],


       [[[ 0.42118728,  0.43297519,  0.43492918, ...,  1.14479685,
           1.13278872,  1.13722253],
         [ 0.42522317,  0.42370971,  0.4216065 , ...,  1.11134458,
           1.10068637,  1.08161545],
         [ 0.43242806,  0.42122281,  0.42144307, ...,  1.10227805,
...
           7.02822757,  7.09734917],
         [ 3.95157754,  3.86564445,  3.83766317, ...,  7.58870363,
           7.54283071,  7.66191769],
         [ 3.93363619,  3.83662581,  3.83417451, ...,  7.77208042,
           7.92680812,  8.15759277]]],


       [[[ 0.50000101,  0.50000101,  0.50000101, ...,  2.73606467,
           2.79075527,  2.86106348],
         [ 0.51731694,  0.51731694,  0.51731694, ...,  2.69699907,
           2.77513051,  2.86106348],
         [ 0.50000101,  0.50000101,  0.50000101, ...,  2.53294182,
           2.62669086,  2.72825599],
         ...,
         [ 3.91306448,  3.89143562,  3.85920548, ...,  6.90828657,
           6.6114502 ,  6.62202311],
         [ 3.93570256,  3.78168535,  3.72219849, ...,  7.27447176,
           7.03061342,  7.07654285],
         [ 3.91414452,  3.72012377,  3.71522117, ...,  7.35348415,
           7.32938242,  7.50127697]]]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-c63f7856-ab6a-417b-a4d1-3b34ddc4881a' class='xr-section-summary-in' type='checkbox'  checked><label for='section-c63f7856-ab6a-417b-a4d1-3b34ddc4881a' class='xr-section-summary' >Attributes: <span>(7)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>title :</span></dt><dd>PM25 Air Pollutant FORECAST at the Surface</dd><dt><span>institution :</span></dt><dd>Data produced by Meteo France</dd><dt><span>source :</span></dt><dd>Data from ENSEMBLE model</dd><dt><span>history :</span></dt><dd>Model ENSEMBLE FORECAST</dd><dt><span>FORECAST :</span></dt><dd>Europe, 20211222+[0H_96H]</dd><dt><span>summary :</span></dt><dd>ENSEMBLE model hourly FORECAST of PM25 concentration at the Surface from 20211222+[0H_96H] on Europe</dd><dt><span>project :</span></dt><dd>MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</dd></dl></div></li></ul></div></div>




```python
dset.resample(time='30min').nearest()
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:     (time: 193, longitude: 700, latitude: 400, level: 1)
Coordinates:
  * time        (time) timedelta64[ns] 00:00:00 00:30:00 ... 4 days 00:00:00
  * longitude   (longitude) float32 -24.95 -24.85 -24.75 ... 44.75 44.85 44.95
  * latitude    (latitude) float32 69.95 69.85 69.75 69.65 ... 30.25 30.15 30.05
  * level       (level) float32 0.0
Data variables:
    pm2p5_conc  (time, level, latitude, longitude) float32 0.4202 ... 7.501
Attributes:
    title:        PM25 Air Pollutant FORECAST at the Surface
    institution:  Data produced by Meteo France
    source:       Data from ENSEMBLE model
    history:      Model ENSEMBLE FORECAST
    FORECAST:     Europe, 20211222+[0H_96H]
    summary:      ENSEMBLE model hourly FORECAST of PM25 concentration at the...
    project:      MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</pre><div class='xr-wrap' hidden><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-03c36c27-b673-42b5-bd48-09019e31ce5e' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-03c36c27-b673-42b5-bd48-09019e31ce5e' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 193</li><li><span class='xr-has-index'>longitude</span>: 700</li><li><span class='xr-has-index'>latitude</span>: 400</li><li><span class='xr-has-index'>level</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-71becb6f-0ab1-4776-82ea-9586f076a056' class='xr-section-summary-in' type='checkbox'  checked><label for='section-71becb6f-0ab1-4776-82ea-9586f076a056' class='xr-section-summary' >Coordinates: <span>(4)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>timedelta64[ns]</div><div class='xr-var-preview xr-preview'>00:00:00 ... 4 days 00:00:00</div><input id='attrs-9f8fd211-8921-49e0-8f37-2ba2ca67e403' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-9f8fd211-8921-49e0-8f37-2ba2ca67e403' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-fd54d728-e307-4a9a-9a5c-807cc66d8bf5' class='xr-var-data-in' type='checkbox'><label for='data-fd54d728-e307-4a9a-9a5c-807cc66d8bf5' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>FORECAST time from 20211222</dd></dl></div><div class='xr-var-data'><pre>array([              0,   1800000000000,   3600000000000,   5400000000000,
         7200000000000,   9000000000000,  10800000000000,  12600000000000,
        14400000000000,  16200000000000,  18000000000000,  19800000000000,
        21600000000000,  23400000000000,  25200000000000,  27000000000000,
        28800000000000,  30600000000000,  32400000000000,  34200000000000,
        36000000000000,  37800000000000,  39600000000000,  41400000000000,
        43200000000000,  45000000000000,  46800000000000,  48600000000000,
        50400000000000,  52200000000000,  54000000000000,  55800000000000,
        57600000000000,  59400000000000,  61200000000000,  63000000000000,
        64800000000000,  66600000000000,  68400000000000,  70200000000000,
        72000000000000,  73800000000000,  75600000000000,  77400000000000,
        79200000000000,  81000000000000,  82800000000000,  84600000000000,
        86400000000000,  88200000000000,  90000000000000,  91800000000000,
        93600000000000,  95400000000000,  97200000000000,  99000000000000,
       100800000000000, 102600000000000, 104400000000000, 106200000000000,
       108000000000000, 109800000000000, 111600000000000, 113400000000000,
       115200000000000, 117000000000000, 118800000000000, 120600000000000,
       122400000000000, 124200000000000, 126000000000000, 127800000000000,
       129600000000000, 131400000000000, 133200000000000, 135000000000000,
       136800000000000, 138600000000000, 140400000000000, 142200000000000,
       144000000000000, 145800000000000, 147600000000000, 149400000000000,
       151200000000000, 153000000000000, 154800000000000, 156600000000000,
       158400000000000, 160200000000000, 162000000000000, 163800000000000,
       165600000000000, 167400000000000, 169200000000000, 171000000000000,
       172800000000000, 174600000000000, 176400000000000, 178200000000000,
       180000000000000, 181800000000000, 183600000000000, 185400000000000,
       187200000000000, 189000000000000, 190800000000000, 192600000000000,
       194400000000000, 196200000000000, 198000000000000, 199800000000000,
       201600000000000, 203400000000000, 205200000000000, 207000000000000,
       208800000000000, 210600000000000, 212400000000000, 214200000000000,
       216000000000000, 217800000000000, 219600000000000, 221400000000000,
       223200000000000, 225000000000000, 226800000000000, 228600000000000,
       230400000000000, 232200000000000, 234000000000000, 235800000000000,
       237600000000000, 239400000000000, 241200000000000, 243000000000000,
       244800000000000, 246600000000000, 248400000000000, 250200000000000,
       252000000000000, 253800000000000, 255600000000000, 257400000000000,
       259200000000000, 261000000000000, 262800000000000, 264600000000000,
       266400000000000, 268200000000000, 270000000000000, 271800000000000,
       273600000000000, 275400000000000, 277200000000000, 279000000000000,
       280800000000000, 282600000000000, 284400000000000, 286200000000000,
       288000000000000, 289800000000000, 291600000000000, 293400000000000,
       295200000000000, 297000000000000, 298800000000000, 300600000000000,
       302400000000000, 304200000000000, 306000000000000, 307800000000000,
       309600000000000, 311400000000000, 313200000000000, 315000000000000,
       316800000000000, 318600000000000, 320400000000000, 322200000000000,
       324000000000000, 325800000000000, 327600000000000, 329400000000000,
       331200000000000, 333000000000000, 334800000000000, 336600000000000,
       338400000000000, 340200000000000, 342000000000000, 343800000000000,
       345600000000000], dtype=&#x27;timedelta64[ns]&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>-24.95 -24.85 ... 44.85 44.95</div><input id='attrs-df001f2f-9354-4b9b-8af8-5bc5d4d64783' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-df001f2f-9354-4b9b-8af8-5bc5d4d64783' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-1fafa579-93e8-450d-85df-0f8c5f956b7b' class='xr-var-data-in' type='checkbox'><label for='data-1fafa579-93e8-450d-85df-0f8c5f956b7b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([-24.950012, -24.849976, -24.75    , ...,  44.75    ,  44.850006,
        44.949997], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>69.95 69.85 69.75 ... 30.15 30.05</div><input id='attrs-0da16bf6-69a8-48f1-b835-8904826c9dec' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-0da16bf6-69a8-48f1-b835-8904826c9dec' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-984f4a41-2306-4e60-9dcd-d9ad0647c5c4' class='xr-var-data-in' type='checkbox'><label for='data-984f4a41-2306-4e60-9dcd-d9ad0647c5c4' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>latitude</dd><dt><span>units :</span></dt><dd>degrees_north</dd></dl></div><div class='xr-var-data'><pre>array([69.95, 69.85, 69.75, ..., 30.25, 30.15, 30.05], dtype=float32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>level</span></div><div class='xr-var-dims'>(level)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.0</div><input id='attrs-3a8c018b-4c8a-466c-9968-e0bb9f889370' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-3a8c018b-4c8a-466c-9968-e0bb9f889370' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-d9a6d5a2-2d2d-4571-920a-265574d16b6d' class='xr-var-data-in' type='checkbox'><label for='data-d9a6d5a2-2d2d-4571-920a-265574d16b6d' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>long_name :</span></dt><dd>level</dd><dt><span>units :</span></dt><dd>m</dd></dl></div><div class='xr-var-data'><pre>array([0.], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-a7fea4c1-ef75-42ce-a4c4-9f8e289ba9cc' class='xr-section-summary-in' type='checkbox'  checked><label for='section-a7fea4c1-ef75-42ce-a4c4-9f8e289ba9cc' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>pm2p5_conc</span></div><div class='xr-var-dims'>(time, level, latitude, longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>0.4202 0.4331 ... 7.329 7.501</div><input id='attrs-a133801f-42bf-4ed1-9659-a4b0b9294542' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-a133801f-42bf-4ed1-9659-a4b0b9294542' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-2b66b70f-c804-4a86-a357-8ee6eb4b6ff4' class='xr-var-data-in' type='checkbox'><label for='data-2b66b70f-c804-4a86-a357-8ee6eb4b6ff4' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>species :</span></dt><dd>PM2.5 Aerosol</dd><dt><span>units :</span></dt><dd>µg/m3</dd><dt><span>value :</span></dt><dd>hourly values</dd><dt><span>standard_name :</span></dt><dd>mass_concentration_of_pm2p5_ambient_aerosol_in_air</dd></dl></div><div class='xr-var-data'><pre>array([[[[ 0.420244, ...,  1.171117],
         ...,
         [ 4.729771, ..., 14.217151]]],


       ...,


       [[[ 0.500001, ...,  2.861063],
         ...,
         [ 3.914145, ...,  7.501277]]]], dtype=float32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-95005d7b-a2b9-44ef-bc55-97e435aa80b5' class='xr-section-summary-in' type='checkbox'  checked><label for='section-95005d7b-a2b9-44ef-bc55-97e435aa80b5' class='xr-section-summary' >Attributes: <span>(7)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>title :</span></dt><dd>PM25 Air Pollutant FORECAST at the Surface</dd><dt><span>institution :</span></dt><dd>Data produced by Meteo France</dd><dt><span>source :</span></dt><dd>Data from ENSEMBLE model</dd><dt><span>history :</span></dt><dd>Model ENSEMBLE FORECAST</dd><dt><span>FORECAST :</span></dt><dd>Europe, 20211222+[0H_96H]</dd><dt><span>summary :</span></dt><dd>ENSEMBLE model hourly FORECAST of PM25 concentration at the Surface from 20211222+[0H_96H] on Europe</dd><dt><span>project :</span></dt><dd>MACC-RAQ (http://macc-raq.gmes-atmosphere.eu)</dd></dl></div></li></ul></div></div>



**Remark**: increasing the frequency of your data e.g. artificialy creating data may not be scientifically relevant. Please use carefuly!

## Question 5: Using a Multi-plot between Rome and Naples, and making averages per day, can you tell us if forecasted PM2.5 will increase or decrease?

-  **Answer**: We will select a sub-area: 11. East to 15.0 East and 40. N to 43. N. PM2.5 will overall decrease.


```python
fig = plt.figure(1, figsize=[10,10])

# We're using cartopy to project our data.
# (see documentation on cartopy)
proj_plot = ccrs.Mercator()

sub_dset = dset.sel(latitude=slice(43., 40.), longitude=slice(11.,15.)).resample(time='1D').mean()
# We need to project our data to the new projection and for this we use `transform`.
# we set the original data projection in transform (here PlateCarree)
p = sub_dset['pm2p5_conc'].plot(transform=ccrs.PlateCarree(), 
                                vmin = 0, vmax = 35,
                                subplot_kws={"projection": proj_plot},
                                col='time', col_wrap=5,
                                robust=True,
                                aspect=dset.dims["longitude"] / dset.dims["latitude"],  # for a sensible figsize
                                cmap=cmc.roma_r)
# We have to set the map's options on all four axes
for ax,i in zip(p.axes.flat,  (np.datetime64('2021-12-22') + dset.time.sel(time=slice(np.timedelta64(0),np.timedelta64(1,'D')))).values):
    ax.coastlines('10m')
    ax.set_title("CAMS PM2.5 " + pd.to_datetime(i).strftime("%d %b %H:%S UTC"), fontsize=12)
# Save your figure
plt.savefig("CAMS-PM2_5-fc-multi-Italy-mean-per-day.png")
```



    
![png](xarray-Galaxy_files/xarray-Galaxy_57_2.png)
    


**Remark**:  `Groupby`
- Use `groupby` instead of `resample` when you wish to group over a dimension that is not `time`. `groupby` is very similar to resample but can be applied to any coordinates and not only to time. 
