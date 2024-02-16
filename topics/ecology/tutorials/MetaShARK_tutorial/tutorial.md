---
layout: tutorial_hands_on

title: Creating metadata using Ecological Metadata Language (EML) standard
  with EML Assembly Line functionalities
zenodo_link: https://zenodo.org/records/10663465
questions:
- How to generate detailled metadata easily from biodiversity datasets ?
- How to use international metadata standard?
- How to update metadata informations ?
objectives:
- Explain the necessity of using such tools when producing ecological metadata
- Learn how to use the interactive tool MetaShARK
- Understand the challenges MetaShARK is trying to respond to
- Learn how to create rich metadata using Ecological Metadata Language (EML) standard
- Learn how to update EML metadata
time_estimation: 30M
key_points:
- This tool aims to improve FAIR quality of metadata focusing on user exeprience and automatic inferences
- Creating metadata as FAIR as possible is a must
- Be carefull of the format and standard of metadata used only EML metadata will work
tags:
  - Metadata
  - EML
  - Ecology
  - Biodiversity
  - FAIR
  - Data Paper
contributions:
    authorship:
        - yvanlebras
    editing:
        - yvanlebras
        - hexylena
    funding:
        - fnso2019
        - pndb

---



This tutorial aims to teach how to use functionalities of the EML Assembly Line R package to produce rich metadata using the Ecological Metadata Language (EML) international metadata standard. Here, we will notably propose a concrete example on how to use Galaxy Ecology tools to create, evaluate and modify EML metadata content using both EML Assemby Line metadata template tabular files, easily readable and editable by humans, and XML file, devoted to machine.


> <comment-title>What does FAIR mean?</comment-title>
> [FAIR](https://www.go-fair.org/fair-principles/) stand for **Findable, Accessible, Interoperable, Reusable**. 
>
>![FAIR Data Principles](./Images/FAIR_data_principles.jpg){:width="500"}
>These principles were to improve the access and usabiliy of data by the machine and to help making data reusable and shareable for users.
>Metadata is the data used to describe and explain all the context behind the production of data. It is necessary to produce a rich and FAIR metadata in order 
>to permit external users to understand and reuse data for their own studies.
{:  .comment}
> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}




# How can EML Assembly Line functionalities help producing rich metadata easily?

A major gap when a researcher is writing metadata documents is the fact that metadata international standards often use formats not really human readable and/or editable as XML or JSON. To answer this issue, Environmental Data Initiaitve (EDI) through the EML Assembly Line R package propose to generate intermediate metadata template files using classical tabular text format.

Another major issue regarding metadata fill in, is the fact that one need to take a lot of time to write, and often rewrite, metadata elements who can be already filled using automatic inferences or use of webservices. Here again, Environmental Data Initiaitve (EDI) through the EML Assembly Line R package propose to generate automatically information related to data attributes, geographic coverage, taxonomic coverage, using the content of provided datafiles.

Finally, through the MetaShARK R Shiny app created by the french biodiversity data hub research infrastructure (Pôle national de données de Biodiversité (PNDB)), user can use a graphical user interface to apply the EML Assembly Line workflow and benefit from some additionnal functionnalities as:
- The capacity to associate terminological resources terms coming from Bioportal ontologies to data attributes as keywords using CEDAR API,
- Automatic fill in of personal information using ORCID API,
- Automatic production of a data paper draft.


> <comment-title>What is a Data Paper?</comment-title>
> According to the [GBIF](https://www.gbif.org/data-papers) (Global Biodiversity Information Facility), 
>A data paper is a peer reviewed document describing a dataset, published in a peer reviewed journal. It takes effort to prepare, curate and describe data. 
>Data papers provide recognition for this effort by means of a scholarly article.
{:  .comment}

# Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial. You can name it "EML assembly Line tutorial" for example
> 2. Download all files on your local computer from Zenodo: https://zenodo.org/api/records/10663465/files-archive. It is neccessary as for now, MetaShARK is deployed from Galaxy but without having a possibility to directly populate MetaShARK app with datafiles from Galaxy. You will thus then have to upload manually data files from your local computer to MetaShARK.
> 3. Unzip the donwloaded archive so you can access each file independently
> 4. Import tsv, netcdf and geotiff data files directly from [Zenodo]({{ page.zenodo_link }}) so it can be used on some steps of the tutorial.
>     -> Training Data for "Creating metadata using Ecological Metadata Language (EML) standard with EML Assembly Line functionalities" Galaxy tutorial:
>    ```
>    https://zenodo.org/records/10663465/files/datafile_1.tsv
>    https://zenodo.org/records/10663465/files/datafile_2.tsv
>    https://zenodo.org/records/10663465/files/LakeGeneva_phytoplankton_1974-2004.nc
>    https://zenodo.org/records/10663465/files/Present.Surface.pH.tif
>    ```
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    ![Upload datafiles from URLs](./Images/0_Upload_datafiles.png){:width="500"}
> 6. Import shapefile related files into Galaxy using the Galaxy upload tool, then on the "composite" tab, specifying "shp" composite type, then uploading .dbf, .shp and .shx files on the dedicated spaces.
>    ![Upload datafiles from URLs](./Images/1_upload_shapefile.png){:width="500"}
{: .hands_on}


# MetaShARK

To deploy a MetaShARK app, you can go to the Galaxy tool [MetaShARK](https://ecology.usegalaxy.eu/root?tool_id=interactive_tool_metashark) and click "Execute". Then, you have to wait the launch of the app, and when ready to be used, you will see the message "There is an InteractiveTool result view available," with an hyperlink on the "Open" statement allowing you to reach the app clicking on it.

WARNING! Note that it is a beta version of MetaShARK R shiny app. You can thus encounter issues using it!

When oppening MetaShARK, you will have an interface looking like this :

![Interface of MetaShARK](./Images/2_MetaShARK.png){:width="500"}

To start creating metadata, you need to reach the "Fill in EML" module, then specify or complete the automatic "Data package name" and mention the "Dataset title". Here a title can be "Manage heterogeneous data files through EML". Finally, you can choose an open licence between CC-BY-4.0, default, or CC0 then click on "Create". 

![Create a data package](./Images/3_MetaShARK_create.png){:width="500"}

Then you can upload datafiles. Here, you can import these files from the downloaded Zenodo archive:
- datafile_1.tsv
- datafile_2.tsv
- LakeGeneva_phytoplankton_1974-2004.nc
- Present.Surface.pH.tif
- 02_Ref.shp
- 02_Ref.shx
- 02_Ref.dbf

![Upload data files](./Images/4_MetaShARK_upload.png){:width="500"}

MetaShARK will normally guess that the three `02_Ref` files are representing a uniq shapefile. MetaShARK will normally guess each data type and infer list of attributes for each file but the geotiff `Present.Surface.pH.tif` one. So now you need to select this datafile and upload the `attributes_Present.Surface.pH.txt` metadata template file so MetaShARK can fill attributes of this file (here the attribute is named "Present.Surface.pH").

![Upload a metadata template for geotiff attributes](./Images/5_MetaShARK_geotiffattributes.png){:width="500"}

Then you can provide a description for this attribute, for example "Present surface pH", then look at each attribute information of each data file so you can click on the "Next" button and go to the next step, to give informations on categorical variables!

![Look at categorical variables](./Images/6_MetaShARK_catvars.png){:width="500"}

Clicking "Next" button will then allows you to fill spatial informations about all GIS recognized datafiles, here the `Present.Surface.pH.tif` geotiff raster file and the `02_Ref` shapefile vector file. Geotiff is in pixel, accuracy unknown and shepfile is in Point, both are in `GCS_WGS_1984`spatial reference.

![Specify spatial informations for GIS datafiles](./Images/7_MetaShARK_spatialinfo.png){:width="500"}

Next step is devoted to specifying geographic coverage. You can use a method between "columns" or "custom". "Custom" allows you to create one to several geographical sites using a map widget where you can draw limits of each site or enter directly latitude and longitude coordinates. "Columns" method, used here, allows you to specify an attribute conataining site names then associated latitude and longitudes attributes.

![Upload a metadata template for geotiff attributes](./Images/8_MetaShARK_geocov.png){:width="500"}


![Upload a metadata template for geotiff attributes](./Images/9_MetaShARK_taxcov.png){:width="500"}


![Upload a metadata template for geotiff attributes](./Images/10_MetaShARK_personal.png){:width="500"}


![Upload a metadata template for geotiff attributes](./Images/11_MetaShARK_abstract.png){:width="500"}

![Upload a metadata template for geotiff attributes](./Images/12_MetaShARK_methods.png){:width="500"}

![Upload a metadata template for geotiff attributes](./Images/13_MetaShARK_keywords.png){:width="500"}

![Upload a metadata template for geotiff attributes](./Images/14_MetaShARK_makeeml.png){:width="500"}



# Conclusion

Here is the end of this short tutorial aiming in explaining the purpose of MetaShARK and how to use it.
Don't hesitate to contact us if you have any questions :)
