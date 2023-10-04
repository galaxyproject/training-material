---
layout: tutorial_hands_on

title: QGIS Web Feature Services
questions:
- how to use WFS webservice to access and import GIS data through Galaxy
objectives:
- Deal with QGIS execution management
- Access, filter and import GIS data using WFS webservice
time_estimation: 0H30
key_points:
- From WFS webservice, you can now access GIS remote data, build a query to select subpart of the data, and import it to Galaxy through QGIS Galaxy interactive tool.
tags:
  - GIS
  - Geographical Information System
  - WFS
  - Spatial data
  - QGIS
  - Maps
  - OGC
contributors:
- colineroyaux
- Marie59
- yvanlebras

---


# Introduction


Based on this [QGIS official tutorial](https://docs.qgis.org/2.18/en/docs/training_manual/online_resources/wfs.html), you will learn here how to access, filter and import GIS data through WFS webservice using QGIS Galaxy interactive tool:

In the Geographical Information System landscape, there is existing standards to help users deal with remote data. Most common webservices are Web Map Services (WMS) and Web Feature Services (WFS). If WMS allows users only to access and display maps stored remotely, WFS is giving access to the features of data so you can modify it and create your own data and maps.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Managing QGIS Galaxy interactive tool
QGIS is now integrated in Galaxy as an interactive tool. This kind of tools is working differently than classical tools as it allows the user to interact with a dedicated graphical interface. This kind of tools is used to give access to Jupyter notebooks, RStudio or R Shiny apps for example. 

To use QGIS, you need to use the {% tool [dedicated form](https://ecology.usegalaxy.eu/root?tool_id=interactive_tool_qgis) %}, you can specify input datasets from your hisrtory you want to use in QGIS, or not ;), then press the execute button to launch a QGIS instance. When the graphical user interface of QGIS is ready to be used, a URL will be displayed at the top of the Galaxy center panel. If you don't see it, you can see and access it through the "Active InteractiveTools" space of the "User" menu.

Once you finished your work on QGIS, if you want to retrieve data and/or entire project, you need to save files in /output, then quit QGIS properly through the "Project" Menu tab.

> <hands-on-title>Deploy your own QGIS instance</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Deploy a {% tool [QGIS instance](https://ecology.usegalaxy.eu/root?tool_id=interactive_tool_qgis) %}
> 3. Access QGIS
> 4. Save your project in /working/output folder
>
{: .hands_on}

# Web Feature Services
## Loading WFS layer

> <hands-on-title>Loading WFS layer</hands-on-title>
>
> 1. Click the "open data source manager" button at the top left of QGIS
> 2. Select "WFS / OGC API - Features"
> 3. Click the "New" button
> 4. In the dialog that appears, enter the Name as `nsidc` and the URL as `http://nsidc.org/cgi-bin/atlas_south?version=1.1.0`.
>
>    > <tip-title>copy pasting from computer to QGIS</tip-title>
>    >
>    > * You can expand the QGIS left panel (where there is 3 dots, verticaly) to access the "clipboard" menu, and paste the content you want to paste on a QGIS form. Then, click outside of this panel to collapse it, and you can click for example on the `url` field to paste the url from your clipbaord
>    >
>    > ![QGIS clipboard](../../images/QGIS/qgistuto1.PNG)
>    {: .tip}
>
> 5. Click OK, then you can create the connection with the "connect" button to see a list of available layers
> 6. Find and select the `antarctica_country_border` layer
> 7. Find and select, holding the "CTRL" button, the `south_poles_wfs`layer
> 9. Unselect "only request features overlapping the view extent" option then click "add"
>    > ![QGIS clipboard](../../images/QGIS/qgistuto3.PNG)
> 10. You now have Antarctica border displayed with a symbol showing the south pole
>    > ![QGIS clipboard](../../images/QGIS/qgistuto4.PNG)
> 11. You can double click each layer at the bottom left to modify symbology, notably colour for "south_poles_wfs" modifying the symbol selecting for example "effect drop shadow" and add a label `Geographic South Pole`
>    > ![QGIS clipboard](../../images/QGIS/qgistuto5.PNG)
>    > ![QGIS clipboard](../../images/QGIS/qgistuto6.PNG)
> 12. And to "antarctica_country_border".
>    > ![QGIS clipboard](../../images/QGIS/qgistuto7.PNG)
> 
{: .hands_on}

## Querying a WFS layer

Even if you can select, download and display entire WFS layers, it is often more efficient to interrogate a layer before load it to QGIS. This is a major interest of the use of such web services, as you can save internet bandwith selecting only part of layers you want in your own QGIS instance.

> <hands-on-title>Querying a WFS layer</hands-on-title>
>
> 1. In the WFS server already created in the first step of this tutorial, there is layer called `countries (excluding Antarctica)`. If we want to know where is South Africa related to the `south_poles_wfs`, there is several manner to operate. One can load the entire layer `countries (excluding Antarctica)` and then using it locally, or we can save bandwith and only load locally the needed informations concerning South Africa. We will here use the second manner, querying the WFS layer to obtain only information we will use in our QGIS instance.
> 2. Click the "open data source manager" button at the top left of QGIS
> 3. Select "WFS / OGC API - Features" then the already connected `nsidc` server
> 4. Select the `countries (excluding Antarctica)` layer and click "Build query" button
>    > ![QGIS clipboard](../../images/QGIS/qgistuto8.PNG)
> 6. On the new dialog box, you can copy/paste this query: `SELECT * FROM country_borders_excluding_antarctica WHERE "Countryeng" = 'South Africa'` using the clipboard functionnality from the left QGIS panel to make a "bridge" between your local clipboard and the virtualized QGIS one.
>    > ![QGIS clipboard](../../images/QGIS/qgistuto9.PNG)
> 7. Clicking ok, you can now see the SQL query on a dedicated column of the layers table.
>    > ![QGIS clipboard](../../images/QGIS/qgistuto10.PNG)
> 9. Clicking the "Add" button, you now have the South Africa layer added and diplayed on QGIS.
>    > ![QGIS clipboard](../../images/QGIS/qgistuto11.PNG)
{: .hands_on}

> <question-title></question-title>
>
> 1. Are you seeing how GGIS show the fact that a layer is resulting from a query made on a larger layer?
>
> > <solution-title></solution-title>
> >
> > 1. There is a dedicated "filtre" icon next to the name of the layer.
> >
> {: .solution}
>
{: .question}








# Conclusion


You just did a classical GIS operation,accessing remote data using WFS webservice through QGIS, retrieving features of interest you can share broadly through Galaxy!
