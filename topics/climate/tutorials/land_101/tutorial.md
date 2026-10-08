---
layout: tutorial_hands_on
title: Where to Start with Land Data in Galaxy
zenodo_link: ''  # Add Zenodo DOI if available
questions:
  - What land data is available in Galaxy for climate and ecology studies?
  - How can I monitor land surface changes and degradation?
  - What tools are available for analyzing land cover, vegetation, and soil data?
  - How is land data connected to climate change?
objectives:
  - Understand the role of land surfaces in the Earth's climate system
  - Discover key land datasets (Sentinel-2, QGIS, etc.)
  - Learn how to use Galaxy tools for land data analysis
  - Explore the connections between land degradation and climate change
  - Access tutorials for visualizing and analyzing land data
time_estimation: 1H
key_points:
  - Land as a critical component of the Earth's climate system
  - Monitoring land surface changes with Sentinel-2
  - Assessing land degradation (erosion, salinization, desertification)
  - Using NDVI and other indices for vegetation and climate studies
  - GIS tools for land data analysis and visualization
tags:
  - earth-system
  - land
  - climate
contributions:
  authorship:
    - Marie59
  funding:
    - gaiadata
    - gallantries
    - fairease
    - eurosciencegateway
---

# 🏜️ Where to Start with Land Data in Galaxy

> <comment-title>Climate 101 Introduction</comment-title>
> Before diving into land data, it's important to understand the basics of climate science and its connection to land surfaces.
> Check out our **[Climate 101 Tutorial]( {% link topics/climate/tutorials/climate-101/tutorial.md %} )** to learn about:
> - The difference between weather and climate
> - Essential Climate Variables (ECVs) related to land
> - How land surfaces interact with climate systems
{: .comment}

The **land** component encompasses the solid part of Earth, including rocks, soils, and surface processes. It plays a **critical role in the Earth's climate system** through:
- **Carbon storage** in soils and vegetation
- **Albedo effects** (reflectivity of land surfaces)
- **Water cycle** interactions
- **Land-atmosphere exchanges** of energy, water, and greenhouse gases

This component aims at assessing **land and soil degradation**, such as:
- Erosion
- Loss of organic matter and soil biodiversity
- Compaction
- Salinization
- Landslides
- Contamination
- Sealing
- Desertification
- **Climate change impacts**


## 📌 Key Data Sources and Tools

### Copernicus Data Space Ecosystem
- **Data**: The **Copernicus Sentinel-2 mission** consists of two polar-orbiting satellites in sun-synchronous orbit (180° phase difference).
  - **Capabilities**: Wide swath width (290 km) and high revisit time for monitoring land surface changes.
  - **Applications**: Supports **land monitoring, agriculture, emergency management, risk mapping, forestry, climate change studies, disaster control, and humanitarian relief operations**.
  - **Strong Climate Connection**: Sentinel-2 data is essential for studying land cover changes that affect climate patterns.

- **Scene Classification Algorithm**:
  - Uses reflective properties to detect clouds and snow.
  - Based on threshold tests using TOA reflectance, band ratios, and indices like:
    - **Normalised Difference Vegetation Index (NDVI)**: Measures vegetation health and density.
    - **Normalised Difference Snow and Ice Index (NDSI)**: Detects snow and ice cover.
  - Produces probabilistic cloud and snow mask quality indicators.
  - **Climate Relevance**: NDVI is crucial for studying vegetation changes related to climate variability.

- **Tool**: [Copernicus Data Space Ecosystem](https://usegalaxy.eu/root?tool_id=interactive_tool_copernicus_notebook)
- **Tutorial**: [From NDVI data with OpenEO to time series visualisation with Holoviews]( {% link topics/ecology/tutorials/ndvi_openeo/tutorial.md %})


### QGIS (Geographical Information System)
- **Functionality**: Access, filter, and import GIS data through **WFS (Web Feature Services)** and **WMS (Web Map Services)**.
  - **WMS**: Allows users to access and display maps stored remotely.
  - **WFS**: Provides access to data features for modification and custom map creation.
  - **Strong Climate Connection**: GIS tools help analyze spatial patterns of land degradation and their climate impacts.

- **Tool**: [QGIS Galaxy Interactive Tool](https://usegalaxy.eu/root?tool_id=interactive_tool_qgis)
- **Tutorial**: [QGIS Web Feature Services]( {% link topics/ecology/tutorials/QGIS_Web_Feature_Services/tutorial.md %} )


### FATES (Functionally Assembled Terrestrial Ecosystem Simulator)
- **Purpose**: Advanced modeling of terrestrial ecosystems and their interactions with climate.
- **Tutorials**:
  - [FATES JupyterLab Tutorial]( {% link topics/climate/tutorials/fates-jupyterlab/tutorial.md %} )
  - [FATES Tutorial]( {% link topics/climate/tutorials/fates/tutorial.md %} )
- **Strong Climate Connection**: FATES helps study how land ecosystems respond to and influence climate change.


### 3D Trees Analysis
A suite of tools for processing and visualizing **LiDAR point cloud data** to study forest structure, biomass, and ecosystem dynamics.

- **{% tool [3DTrees: LAS/LAZ Standardization](toolshed.g2.bx.psu.edu/repos/ecology/las_laz_standardization/las_laz_standardization/1.0.0) %}**: Standardize LAS/LAZ files or validate collections for consistency.
- **{% tool [3DTrees: Potree Converter](toolshed.g2.bx.psu.edu/repos/ecology/potree_converter/potree_converter/1.0.0) %}**: Convert LAS/LAZ point clouds to Potree octree format for web visualization.

🔜 **Upcoming**: Public workflow and training materials for 3D tree analysis.

## 🔗 Cross-Domain Connections

The following tools and resources are **also relevant to other Earth System domains**:

### Biodiversity Tools (Terrestrial Focus)
- **{% tool [Get Species Occurrences Data](toolshed.g2.bx.psu.edu/repos/ecology/gbif_data/gbif_data/0.9.0) %}**
  *Also listed in: [Biodiversity]( {% link topics/climate/tutorials/earth_system/tutorial.md%} )*
  - Access GBIF's global biodiversity records for terrestrial species.
