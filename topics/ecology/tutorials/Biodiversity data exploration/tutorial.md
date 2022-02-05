---
layout: tutorial_hands_on

title: Biodiversity data exploration
zenodo_link: https://zenodo.org/record/5930763/files/Reel_life_survey_fish_sample.tabular?download=1
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
contributors:
- onorvez
- Marie59
- colineroyaux
- yvanlebras

---


# Introduction
{:.no_toc}

This tutorial will guide you on the exploration of biodiversity data having taxonomic, spatial and temporal informations.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Downloading biodiversity data

First step is to download biodiversity data on your Galaxy history. Here we will use a "classical" (containing taxonomic, spatial and temporal informations) biodiversity dataset from the well known ["Reef life survey" initiative](https://reeflifesurvey.com/).

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/5930763/files/Reel_life_survey_fish_sample.tabular?download=1)
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets "reef_life_fish" for example
> 4. Check that the datatype is tabular
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Want to spatially anoymize your data 

A first step of this tutorial will show you how you can simply apply spatial coordinates anonymization if you want to share data without spatial context.


## Sub-step with **Spatial coordinates anonymization**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Spatial coordinates anonymization](toolshed.g2.bx.psu.edu/repos/ecology/tool_anonymization/tool_anonymization/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select column containing latitudes in decimal degrees"*: `c9`
>    - *"Select column containing longitudes in decimal degrees"*: `c10`
>
>
>
{: .hands_on}

# Descriptive statistical testing

## Sub-step with **Homoscedasticity and normality**

> ### {% icon hands_on %} Hands-on: Here we will check homogeneity of variances (Levene test) for every species and represent it through multiple boxplots and the normal distribution (Kolmogorov-Smirnov test) represented by a distribution histogram and a Q-Q plot.

If the levene test is significant (P-value in column Pr < 0.5 and at least one * at the end of the 4th line), variances aren't homogeneous, the hypothesis of homoscedasticity is rejected.

If the K-S test is significant (p-value < 0.5), your numerical variable isn't normally distributed, the hypothesis of normality is rejected.
>
> 1. {% tool [Homoscedasticity and normality](toolshed.g2.bx.psu.edu/repos/ecology/ecology_homogeneity_normality/ecology_homogeneity_normality/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select column containing temporal data (year, date, ...)"*: `c11`
>    - *"Select column containing species"*: `c16`
>    - *"Select column containing numerical values (like abundances)"*: `c18`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}



![Boxplot dispersion_example](../../images/BiodivExplo/Boxplot_and_dispersion_plot.png)
![Homoscedasticity and normality_example](../../images/BiodivExplo/Homoscedasticity_and_normality_Homogeneity_of_%20Amblygobius%20phalaena.png)


## Sub-step with **Variables exploration**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Variables exploration](toolshed.g2.bx.psu.edu/repos/ecology/ecology_link_between_var/ecology_link_between_var/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Variables links exploration"*: `Collinearity between selected numerical variables for each species`
>        - *"Select column containing species"*: `c16`
>        - *"Select columns containing numerical values"*: `c['12', '17', '18']`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

![Variable_exploration_example](../../images/BiodivExplo/Variables_exploration_collinarity_of_Amblygobius%20phalaena.png)  

## Sub-step with **Presence-absence and abundance**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Presence-absence and abundance](toolshed.g2.bx.psu.edu/repos/ecology/ecology_presence_abs_abund/ecology_presence_abs_abund/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Variables presence, absence and abundance"*: `Abundance map in the environment `
>        - *"Select column containing latitudes "*: `c9`
>        - *"Select column containing longitudes"*: `c10`
>        - *"What do you study in this analysis ?"*: `fishes`
>        - *"Select column containing taxon "*: `c16`
>    - *"Select column containing abundances "*: `c18`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

![Presence-absence-example](../../images/BiodivExplo/Presence-absence_and_abundance_mappy.png)

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Statistics on presence-absence**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Statistics on presence-absence](toolshed.g2.bx.psu.edu/repos/ecology/ecology_stat_presence_abs/ecology_stat_presence_abs/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select a column containing numerical values (such as the abundance) "*: `c18`
>    - *"Select the column of the x-axis : most commonly species"*: `c16`
>    - *"Select column containing locations "*: `c8`
>    - *"Select column containing temporal data (year, date, ...) "*: `c11`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Local Contributions to Beta Diversity (LCBD)**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Local Contributions to Beta Diversity (LCBD)](toolshed.g2.bx.psu.edu/repos/ecology/ecology_beta_diversity/ecology_beta_diversity/0.0.0) %} with the following parameters:
>    - {% icon param-file %} *"Input table"*: `output` (Input dataset)
>    - *"Select column with abundances"*: `c18`
>    - *"Select column with locations"*: `c8`
>    - *"Select column containing taxon"*: `c16`
>    - *"Select column containing dates"*: `c11`
>    - *"Other LCBD : spatialized representation or xy plot."*: `Spatialized representation`
>        - *"Select column containing latitudes in decimal degrees"*: `c9`
>        - *"Select column containing longitudes in decimal degrees"*: `c10`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_(LCBD)_Beta_diversity_through_space.png) 
![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_(LCBD)_LCBD_sites_time.png) 
![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_(LCBD)_Mean_LCBD_through_time.png) 
![Variable_exploration_example](../../images/BiodivExplo/Local_Contributions_to_Beta_Diversity_(LCBD)_SCBD_Species_Radar_plot.png) 

Final absence correlation plot:
![Absence correlation_example](../../images/BiodivExplo/Absence-correlation_plot.png)

## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
