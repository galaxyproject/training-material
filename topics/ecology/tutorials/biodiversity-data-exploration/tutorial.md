---
layout: tutorial_hands_on

title: Biodiversity data exploration
zenodo_link: https://zenodo.org/record/6107457/files/IMOS-National_Reef_Monitoring_Network_Sub-Facility-Global_reef_fish_abundance_and_biomass.csv?download=1
questions:
- How to explore biodiversity data?
- How to look at Homoscedasticity, normality or collinearity of presences-absence or abundance data?
- How to compare beta diversity taking into account space, time and species components?
objectives:
- Explore Biodiversity data with taxonomic, temporal and geographical informations
- Have an idea about quality content of the data regarding statistical tests like normality or homoscedasticity and coverage like temporal or geographical coverage
time_estimation: 1H
key_points:
- Explore your data before diving into deep analysis
tags:
  - taxonomic data
  - data quality
contributors:
- onorvez
- Marie59
- colineroyaux
- yvanlebras

---


# Introduction


This tutorial will guide you on the exploration of biodiversity data having taxonomic, spatial and temporal informations.

We'll be using Reef Life Survey (RLS) data extracted from the Australian Ocean Data Network (AODN) portal. We'll use a subset done directly on the AODN data portal (https://portal.aodn.org.au/) on this dataset "IMOS - National Reef Monitoring Network Sub-Facility - Global reef fish abundance and biomass". We decided to use data only on the Mollusca phylum from the east coast of Australia between 2008 and 2021. We'll explore this dataset in the view of making statistical analyses so we will check the homoscedasticity and normality of the variables, see if some variables are correlated, how the data is distributed through space and time, etc ... And finally, we'll explore Beta diversity through the computation of the SCBD and LCBD (Species and Local Contribution to Beta Diversity).

> <details-title>Definitions of SCBD and LCBD</details-title>
>
> Species Contribution to Beta Diversity: degree of variation for individual species across the study area.
>
> Local Contribution to Beta Diversity: comparative indicators of the ecological uniqueness of the sites.
>
{: .details}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data preparation

First step is to download biodiversity data on your Galaxy history. Here we will use a "classical" (containing taxonomic, spatial and temporal informations) biodiversity dataset from the well known ["Reef life survey" initiative](https://reeflifesurvey.com/).

## Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial and give it a name (example: "RLS for biodiversity data exploration tutorial")
>    for you to find it again later if needed.
> 2. Import the files from [Zenodo](https://zenodo.org/record/6107457/files/IMOS-National_Reef_Monitoring_Network_Sub-Facility-Global_reef_fish_abundance_and_biomass.csv?download=1)
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets "reef_life_molluscs" for example and preview your dataset
>
>    You can see that the dataset hasn't been detected to be a CSV dataframe, it is because RLS data directly puts the
>    metadata of the dataframe in the first lines before the dataframe so you'll have to remove these lines using the {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>        - {% icon param-text %} *"Remove first"*: `72`
>        - {% icon param-file %} *"from"*: reef_life_molluscs data file
>
>    Then, verify if your new file hasn't got hashtags in the first lines and then ask Galaxy to autodetect datatype (click on the pencil, then "Datatypes" then click on "Auto-detect" button). Galaxy will normally detect it as csv.
>
> 4. Convert datatype CSV to tabular
>
>    {% snippet faqs/galaxy/datasets_convert_datatype.md conversion="Convert CSV to tabular" %}
>
{: .hands_on}

## Customize your dataset

In order to clean unnecessary informations from the table we will now cut a few columns and change the format of time information.

> <hands-on-title>Clean your data</hands-on-title>
>
> 1. Use {% tool [Advanced cut columns from a table (cut)](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cut_tool/1.1.0) %}
>    with following parameters :
>    - {% icon param-files %} *"File to cut"*: Convert CSV to tabular data file
>    - {% icon param-select %} *"Operation"*: `Keep`
>    - {% icon param-select %} *"Delimited by"*: `Tab`
>    - {% icon param-select %} *"Cut by"*: `fields`
>      - {% icon param-select %} *"List of Fields"*: `Column: 8` `Column: 10` `Column: 11` `Column: 12` `Column: 25` `Column: 28`
>
> 2. Use {% tool [Column Regex Find And Replace](toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regexColumn1/1.0.0) %} with following parameters:
>    - {% icon param-files %} *"Select cells from"*: Advanced Cut data file
>    - {% icon param-select %} *"using column"*: `Column: 4`
>    - {% icon param-repeat %} Click *"+ Insert Check"*:
>         - {% icon param-text %} *"Find Regex"*: `([0-9]{4})-[0-9]{2}-[0-9]{2}`
>         - {% icon param-text %} *"Replacement"*: `\1`
>
>
{: .hands_on}

# Data checking

## Homoscedasticity and normality analysis

> <hands-on-title>Here we will check homogeneity of variances (Levene test) for every species and represent it through multiple boxplots and the normal distribution (Kolmogorov-Smirnov test) represented by a distribution histogram and a Q-Q plot.</hands-on-title>
>
> 1. {% tool [Homoscedasticity and normality](toolshed.g2.bx.psu.edu/repos/ecology/ecology_homogeneity_normality/ecology_homogeneity_normality/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: Column Regex Find and Replace data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Select column containing temporal data (year, date, ...)"*: `c4`
>    - {% icon param-select %} *"Select column containing species"*: `c5`
>    - {% icon param-select %} *"Select column containing numerical values (like abundances)"*: `c6`
>
>    You have to get three outputs: the Levene Test for homoscedasticity dataset, the Kolmogrov-Smirnov test for normality  and 9 PNG files in a data collection representing the homogeneity of variances for each species at each time point of the study.
>    If the levene test is significant (P-value in column Pr < 0.5 and at least one * at the end of the 4th line), variances aren't homogeneous, the hypothesis of homoscedasticity is rejected.
>    If the K-S test is significant (p-value < 0.5), your numerical variable isn't normally distributed, the hypothesis of normality is rejected.
>    The two tests have to be significant so variances aren't homogenous and data isn't normally distributed.
>
{: .hands_on}

![Homoscedasticity and normality_example on Sepioteuthis australis](../../images/BiodivExplo/Homoscedasticity_and_normality_of_%20Sepioteuthis%20australis.png)

## Autocorrelation in your data

> <hands-on-title>Autocorrelation</hands-on-title>
>
> 1. {% tool [Variables exploration](toolshed.g2.bx.psu.edu/repos/ecology/ecology_link_between_var/ecology_link_between_var/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: Column Regex Find and Replace data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Variables links exploration"*: `Autocorrelation of one selected numerical variable`
>        - {% icon param-select %} *"Select columns containing numerical values"*: `c6`
>
>     You have to get two outputs, one text file containing the Autocorrelation function values and one PNG file in the data collection showing the autocorrelation for a variable.
>     If the bars of the histogram are strictly confined between the dashed lines (representing 95% confidence interval without white noise), there is auto-correlation.
>
>    Here, we don't see there is autocorrelation.
>
{: .hands_on}

![Variable_exploration_autocorrelation example](../../images/BiodivExplo/Variables_exploration_autocorrelation.png)

## Check collinearity in your data

> <hands-on-title>Collinearity between numerical variables</hands-on-title>
>
> 1. {% tool [Variables exploration](toolshed.g2.bx.psu.edu/repos/ecology/ecology_link_between_var/ecology_link_between_var/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: formatted biodiversity data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Variables links exploration"*: `Collinearity between selected numerical variables for each species`
>        - {% icon param-select %} *"Select column containing species"*: `c5`
>        - {% icon param-select %} *"Select columns containing numerical values"*: `c['4', '6']`
>
>    You have to get two outputs, one describing species we couldn't evaluate and one PNG file with one plot containing multiple correlation plots and the correlation values between each variables.
>
{: .hands_on}

![Variable_exploration_collinarity example](../../images/BiodivExplo/Variables_exploration_collinarity_of_Sepioteuthis%20australis.png)

# Data exploration
## Visualize abundance repartition through space

> <hands-on-title>Abundance map in the environment</hands-on-title>
>
> 1. {% tool [Presence-absence and abundance](toolshed.g2.bx.psu.edu/repos/ecology/ecology_presence_abs_abund/ecology_presence_abs_abund/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: formatted biodiversity data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Variables presence, absence and abundance"*: `Abundance map in the environment `
>        - {% icon param-select %} *"Select column containing latitudes "*: `c2`
>        - {% icon param-select %} *"Select column containing longitudes"*: `c3`
>        - {% icon param-text %} *"What do you study in this analysis ?"*: `Molluscs of Australian east coast`
>        - {% icon param-select %} *"Select column containing taxon "*: `c5`
>    - {% icon param-select %} *"Select column containing abundances "*: `c6`
>
>    You have to get two outputs, one with the map of the abundance through space with the coordinates and one text file to inform you about the geographical extent of your map.
>
{: .hands_on}

![Presence-absence-example](../../images/BiodivExplo/Presence-absence_and_abundance_mappy.png)

## Visualize the number of locations where each taxons are present

> <hands-on-title>Presence count of taxons (barplot)</hands-on-title>
>
> 1. {% tool [Presence-absence and abundance](toolshed.g2.bx.psu.edu/repos/ecology/ecology_presence_abs_abund/ecology_presence_abs_abund/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: formatted biodiversity data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Variables presence, absence and abundance"*: `Presence count of taxons (barplot)`
>        - {% icon param-select %} *"Select column containing your separation variable"*: `c1`
>        - {% icon param-select %} *"Select column containing taxon"*: `c5`
>    - {% icon param-select %} *"Select column containing abundances "*: `c6`
>
>    You have to get two outputs, one with 120 PNG files (one for each site) representing the number of locations where each taxons are present and one text file to inform you about the used locations.
>
{: .hands_on}


## Visualize the rarefaction curves of your species

> <hands-on-title>Rarefaction curve of species</hands-on-title>
>
> 1. {% tool [Presence-absence and abundance](toolshed.g2.bx.psu.edu/repos/ecology/ecology_presence_abs_abund/ecology_presence_abs_abund/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: formatted biodiversity data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Variables presence, absence and abundance"*: `Rarefaction curve of species`
>        - {% icon param-text %} *"Size of subsamples"*: `200`
>        - {% icon param-select %} *"Select column containing species"*: `c5`
>    - {% icon param-select %} *"Select column containing abundances "*: `c6`
>
>    You have to get two outputs, one data collection with one PNG files representing the rarefaction curves of each species in one graph and one tabular file with log informations.
>
{: .hands_on}


![Variable_exploration_rarefaction curves](../../images/BiodivExplo/Presence-absence_and_abundance_rarefaction.png)


## Visualize the dispersion of a numeric variable and the correlation between species absence

> <hands-on-title>Boxplot of dispersion and correlation of absence</hands-on-title>
>
> 1. {% tool [Statistics on presence-absence](toolshed.g2.bx.psu.edu/repos/ecology/ecology_stat_presence_abs/ecology_stat_presence_abs/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: formatted biodiversity data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Select a column containing numerical values (such as the abundance) "*: `c6`
>    - {% icon param-select %} *"Select the column of the x-axis : most commonly species"*: `c5`
>    - {% icon param-select %} *"Select column containing locations "*: `c1`
>    - {% icon param-select %} *"Select column containing temporal data (year, date, ...) "*: `c4`
>
>    You have to get two outputs, one PNG file with the boxplot and dispersion plot of the abundance and one plot representing wether the absence of several species is correlated. In the second plot if you see there is a cross on the round representation the two related species haven't got their absences correlated, the other are correlated and seem to be co-absent.
>
{: .hands_on}


![Boxplot dispersion_example](../../images/BiodivExplo/Boxplot_and_dispersion_plot.png)
![Absence correlation_example](../../images/BiodivExplo/Absence-correlation_plot.png)

# Beta diversity
## Local and Species Contribution to Beta Diversity (SCBD and LCBD)

> <hands-on-title>Task description</hands-on-title>
>
> 1. {% tool [Local Contributions to Beta Diversity (LCBD)](toolshed.g2.bx.psu.edu/repos/ecology/ecology_beta_diversity/ecology_beta_diversity/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: formatted biodiversity data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Select column with abundances"*: `c6`
>    - {% icon param-select %} *"Select column with locations"*: `c1`
>    - {% icon param-select %} *"Select column containing taxon"*: `c5`
>    - {% icon param-select %} *"Select column containing dates"*: `c4`
>    - {% icon param-select %} *"Other LCBD : spatialized representation or xy plot."*: `Spatialized representation`
>        - {% icon param-select %} *"Select column containing latitudes in decimal degrees"*: `c2`
>        - {% icon param-select %} *"Select column containing longitudes in decimal degrees"*: `c3`
>
>    You have to get three outputs. Two text file containing a table with information on the beta diversity and one text file with the list of species that has a SCBD larger than the mean SCBD. One data collection with PNG files showing multiple plots according to one type of variable in order to vizualize the beta diversity.
>
>
{: .hands_on}


![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_LCBD_Beta_diversity_through_space.png)
![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_LCBD_LCBD_sites_time.png)
![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_LCBD_Mean_LCBD_through_time.png)
![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_LCBD_SCBD_Species_Radar_plot.png)


# Conclusion


Here, you just explored your biodiversity dataframe properly and you know a lot more about your data. You can now peacefully make your statiscal analyses as most of the red flags you can get have been revealed by this toolsuite ! Enjoy !

# Bonus: Want to spatially anoymize your data?

A step of this tutorial could be to show you how you can simply apply spatial coordinates anonymization if you want to share data without spatial context, particularly of interest if you want to share endangered species oriented data.

## Bonus! Spatial coordinates anonymization

> <hands-on-title>Task description</hands-on-title>
>
> 1. {% tool [Spatial coordinates anonymization](toolshed.g2.bx.psu.edu/repos/ecology/tool_anonymization/tool_anonymization/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: Column Regex Find and Replace data file
>    - {% icon param-select %} *"First line is a header line"*: `Yes`
>    - {% icon param-select %} *"Select column containing latitudes in decimal degrees"*: `c2`
>    - {% icon param-select %} *"Select column containing longitudes in decimal degrees"*: `c3`
>
{: .hands_on}

