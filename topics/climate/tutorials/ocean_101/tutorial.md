---
layout: tutorial_hands_on
title: Where to Start with Ocean Data in Galaxy
zenodo_link: ''  # Add Zenodo DOI if available
questions:
  - What ocean data is available in Galaxy?
  - How can I access and analyze ocean datasets?
  - What tools are available for oceanography, marine chemistry, and climate interactions?
objectives:
  - Discover key ocean datasets (Argo, Copernicus, EMODnet, etc.)
  - Learn how to use Galaxy tools for ocean data analysis
  - Explore tutorials and workflows for ocean data visualization
time_estimation: 1H
key_points:
  - Ocean as a key component of the Earth’s climate system
  - Real-time monitoring with Argo, Sentinel-3, and EMODnet
  - Tools for ocean data access, analysis, and visualization
tags:
  - earth-system
  - ocean
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

# 🌊 Where to Start with Ocean Data in Galaxy

The ocean is a key component of the Earth’s climate system, requiring **continuous real-time monitoring** to understand its dynamics and predict its evolution. Galaxy provides access to a variety of ocean datasets, tools, and tutorials to help you explore oceanography, marine chemistry, and climate interactions.


## 📌 Key Data Sources and Tools

### Argo Program
- **Data**: Latitude, longitude, time, water temperature, pressure, salinity, chlorophyll, oxygen, and nitrate.
- **Tool**: {% tool [Argo Data Access (Galaxy version 0.1.15+galaxy0)](toolshed.g2.bx.psu.edu/repos/ecology/argo_getdata/argo_getdata/0.1.15+galaxy0) %}
- **Tutorial**: [Analyse Argo Data]( {% link topics/climate/tutorials/argo_pangeo/tutorial.md %} )
- **Workflow**: [Full Analyse Argo Data](https://earth-system.usegalaxy.eu/u/marie.josse/w/full-analyse-argo-data)

### Nitrate DMQC for Argo Floats
- **Description**: Focuses on **Quality Control (QC) and calibration** of nitrate sensors mounted on Argo floats. This includes **Data Management Quality Control (DMQC)** and **Qualified Control Validation (QCV)** analysis to ensure accurate nitrate measurements.
- **Tutorial**: [Nitrate DMQC for Autonomous Platforms such as Argo Floats]( {% link topics/climate/tutorials/ocean_qcv_analysis/tutorial.md %} )
- **Workflow**: [Argo-Glider Nitrate QCV](https://usegalaxy.eu/published/workflow?id=44827462c065bae3)
- **Tools Used**:
  - **QCV Harmonizer**: Harmonizes vocabulary and format for input files.
  - **ODV Collection Manager**: Creates ODV spreadsheet collections for visualization and qualification.
  - **Ocean Data View (ODV)**: Interactive tool for visualizing and qualifying data.
  - **ODV History Manager**: Reports QC flags and data edits.
  - **Biogeochemical Calibration**: Calibrates nitrate using reference datasets and methods.


### Copernicus Data Space Ecosystem
- **Data**: Sea surface topography, sea and land surface temperature, ocean and land surface color (Sentinel-3 mission).
- **Tool**: [Copernicus Data Space Ecosystem (JupyterLab)](https://usegalaxy.eu/root?tool_id=interactive_tool_copernicus_notebook)
- **Tutorial**: [Sentinel-5P Data Visualisation]( {% link topics/climate/tutorials/sentinel5_data/tutorial.md %} )
- **Example**: In JupyterLab, navigate to `samples > openeo > Sentinel_3.ipynb` for a Sentinel-3 data example.


### EMODnet Chemistry
- **Data**: Standardized, harmonized, and validated marine chemical data for European sea basins (eutrophication, contaminants).
- **Tool**: [Ocean Data View (ODV)](https://usegalaxy.eu/root?tool_id=interactive_tool_odv)
- **Tutorials & Workflows**:
  - [Ocean’s Variables Study]{% link topics/climate/tutorials/ocean-variables/tutorial.md %})
  - [Workflow: Ocean’s Variables Study](https://usegalaxy.eu/u/marie.josse/w/oceans-var-v2)

### Copernicus Marine Data Store
- **Data**: Ocean product catalog with ~15 variables (hindcast, current, and forecast data).
- **Tool**: {% tool [Copernicus Marine Data Store retrieve marine data ](toolshed.g2.bx.psu.edu/repos/ecology/copernicusmarine/copernicusmarine/2.3.0+galaxy0) %}

## 🔗 Cross-Domain Connections

The following tools and resources are **also relevant to other Earth System domains**:

### Biodiversity Tools (Marine Focus)
- **{% tool [OBIS Occurrences](toolshed.g2.bx.psu.edu/repos/ecology/obis_data/obis_data/0.0.2) %}**
  *Also listed in: [Biodiversity]( {% link topics/climate/tutorials/earth_system/tutorial.md %} )*
  - Access global marine biodiversity data for science and conservation.

- **[Marine Omics Visualization Workflow](https://usegalaxy.eu/u/marie.josse/w/marine-omics-visualisation)**
  *Also listed in: [Biodiversity]( {% link topics/climate/tutorials/earth_system/tutorial.md %} )*
  - Analyze marine microbial metagenomics datasets.

