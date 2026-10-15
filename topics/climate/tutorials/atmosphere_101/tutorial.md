---
layout: tutorial_hands_on
title: Where to Start with Atmosphere Data in Galaxy
zenodo_link: ''  # Add Zenodo DOI if available
questions:
  - What atmospheric data is available in Galaxy?
  - How can I access and analyze atmospheric datasets for climate and air quality?
  - What tools are available for monitoring atmospheric composition ?
  - How are atmospheric data and climate connected?
objectives:
  - Understand the role of the atmosphere in the Earth's climate system
  - Discover key atmospheric datasets (Sentinel-5P, ERA5, etc.)
  - Learn how to use Galaxy tools for atmospheric data analysis
  - Explore the strong connections between atmospheric data and climate studies
  - Access tutorials for visualizing and analyzing atmospheric data
time_estimation: 1H
key_points:
  - Atmosphere as a layer of gas around the Earth, strongly connected to climate
  - Monitoring atmospheric composition with Sentinel-5P
  - Accessing ERA5 reanalysis data for climate studies
  - Volcanic activity monitoring (SO₂, Aerosol Index) and its climate impacts
  - Using Pangeo for advanced atmospheric data analysis
tags:
  - earth-system
  - atmosphere
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

# ☁️ Where to Start with Atmosphere Data in Galaxy

> <comment-title>Climate 101 Introduction</comment-title>
> Before diving into atmospheric data, it's important to understand the basics of climate science.
> Check out our **[Climate 101 Tutorial]( {% link topics/climate/tutorials/climate-101/tutorial.md %} )** to learn about:
> - The difference between weather and climate
> - Essential Climate Variables (ECVs)
> - Observations, reanalysis, predictions, and projections
{: .comment}

The **atmosphere** is the layer of gas around the Earth, playing a **critical role in regulating climate, weather, and air quality**. Atmospheric data is **strongly connected to climate studies**, as it helps us understand:
- How greenhouse gases affect global temperatures
- The impact of aerosols on climate patterns
- The relationship between atmospheric composition and extreme weather events

Galaxy provides access to **atmospheric datasets, tools, and tutorials** to help you explore these connections.


## 📌 Key Data Sources and Tools

### Copernicus Data Space Ecosystem
- **Data**: Atmospheric composition data from the **Sentinel-5 Precursor (Sentinel-5P)** mission, the first Copernicus mission focused on monitoring the atmosphere.
  - **Partnership**: Developed by ESA, the European Commission, the Netherlands Space Office, industry, data users, and scientists.
  - **Strong Climate Connection**: Sentinel-5P data is essential for studying atmospheric components that drive climate change.
- **Tool**: [Copernicus Data Space Ecosystem](https://usegalaxy.eu/root?tool_id=interactive_tool_copernicus_notebook)
- **Tutorials**:
  - [Sentinel-5P Data Visualisation]( {% link topics/climate/tutorials/sentinel5_data/tutorial.md %} )
  - [Pangeo Tutorial]( {% link topics/climate/tutorials/pangeo/tutorial.md %} )
  - [Pangeo Notebook Tutorial]( {% link topics/climate/tutorials/pangeo-notebook/tutorial.md %} )

#### Volcano 🌋 Use Case
- **Data**: Subset of Sentinel-5P L2 data (e.g., April 1 to May 30, 2021) for the **Antilles region**, focusing on **La Soufrière Saint Vincent** (eruption on April 9, 2021).
  - **Variables**: Sulfur dioxide (SO₂) and Aerosol Index (AI).
  - **Purpose**: High-temporal-resolution data for studying volcanic impacts on **air quality, air traffic, and climate**.
  - **Climate Connection**: Volcanic eruptions can significantly affect climate by injecting aerosols into the stratosphere, which can reflect sunlight and cool the planet.
  - **Relevance**: Useful for scientists and institutions monitoring **air quality, aviation safety, and climate impacts**.


### Climate Data Store
- **Data**: **ERA5 reanalysis data** providing hourly snapshots of the atmosphere, land surface, and ocean waves from **1959 onwards** (preliminary data available for 1950–1958).
  - **Features**: Includes uncertainty estimates to highlight the evolution of observing systems.
  - **Strong Climate Connection**: ERA5 is one of the most comprehensive datasets for studying past and present climate.
- **Tool**: {% tool [Climate Data Store (Galaxy version 0.2.0)](toolshed.g2.bx.psu.edu/repos/climate/c3s/c3s/0.2.0) %}
- **Tutorial**: [Panoply Tutorial]( {% link topics/climate/tutorials/panoply/tutorial.md %} ) (for visualizing atmospheric and climate data).

## 🔗 Cross-Domain Connections

The following atmospheric tools and resources **support climate studies** that may overlap with other domains:

### Climate Tools (Relevant to All Domains)
- **{% tool [Climate Data Store (ERA5)](toolshed.g2.bx.psu.edu/repos/climate/c3s/c3s/0.2.0) %}**
  *Relevant for: Ocean, Land, Biodiversity*
  - Hourly atmospheric, land, and ocean wave data for climate studies.

- **[Sentinel-5P Data](https://usegalaxy.eu/root?tool_id=interactive_tool_copernicus_notebook)**
  *Relevant for: Ocean (marine atmosphere), Land (air quality over land)*
  - Atmospheric composition data (SO₂, NO₂, O₃, aerosol index).
