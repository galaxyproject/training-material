---
layout: tutorial_hands_on
title: Analyzing CMIP6 data with Galaxy Climate JupyterLab
zenodo_link: ''
requirements:
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
        - jupyterlab
  -
    type: "external"
    title: "Programming with Python"
    link: "https://swcarpentry.github.io/python-novice-inflammation/"
  -
    type: "external"
    title: "The Unix Shell"
    link: "http://swcarpentry.github.io/shell-novice/"

questions:
- Why using Climate data from the Coupled Model Intercomparison Project (CMIP) Phase 6? 
- How to use CMIP6 data with Galaxy Climate JupyterLab?
- What is pangeo software ecosystem and how to use it to analyze climate data? 
objectives:
- Learn about Climate data from the Coupled Model Intercomparison Project
- Learn to use Climate historical and Climate projections data
- Learn about Pangeo Software Ecosystem
- Learn to analyze CMIP6 with Galaxy Climate JupyterLab
time_estimation: 3H
key_points:
- The Coupled Model Intercomparison Project
- Climate historical data and Climate projections
- Pangeo ecosystem
contributors:
- annefou

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->
## What is climate prediction?

According to the [Glossary of Meteorology](https://glossary.ametsoc.org/wiki/Climate_prediction), _climate prediction_ is defined as 

~~~  
The prediction of various aspects of the climate of a region during some future period of time.

Climate predictions are generally in the form of probabilities of anomalies of climate variables (e.g., temperature, precipitation), with lead times up to several seasons (
see climate anomaly). The term "climate projection" rather than "climate prediction" is now commonly used for longer-range predictions with a higher degree of uncertainty and a lesser degree of specificity. For example, this term is often used for "predictions" of climate change that depend on uncertain consequences of anthropogenic influences such as land use and the burning of fossil fuels.
~~~

## Motivation
* **Cloud feedbacks** are the major contributor to **high climate sensitivity** in global climate models (GCMs) [(Flato et al., 2013)](https://www.cambridge.org/core/books/climate-change-2013-the-physical-science-basis/evaluation-of-climate-models/94BC2268C864F2C6A18436DB22BD1E5A).
* In a warming climate, the **cloud phase changes** (reduction of ice phase), which has significant **implications for radiation, glacier, and ice sheet mass balances**.
* GCMs **underestimate** the ice water path and **overestimate** the liquid water path compared to satellite measurements [(Komurcu et al., 2014)](https://doi.org/10.1002/2013JD021119).
* The over-/underestimation is likely **related to the parameterization** of ice nucleation and growth processes in GCMs.
* Which will **impact** the precipitation microphysics and **the amount of solid and liquid precipitation** reaching the ground. 
* The **majority** of the total precipitation reaching the ground **originates from mixed-phase processes** [(Mülmenstädt et al. 2015)](https://doi.org/10.1002/2015GL064604)



This tutorial will use Coupled Model Intercomparison Project 6 (CMIP6) output of surface snowfall rate and cloud phase. We are interested in a correlation between the already known cloud phase bias and the surface snowfall rate in GCMS.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# What is CMIP6?
CMIP6 is the Coupled Model Intercomparison Project in its 6th phase. [CMIP6](https://esgf-node.llnl.gov/projects/cmip6/) coordinates independent model intercomparisons and their experiments. Further, they organize and distribute outputs from these standardized experiments. The results by CMIP6 have been widely used to understand the climate and project future climate scenarios, such as in the [IPCC Assessment Reports](https://www.ipcc.ch/).

# Outcome
* We will learn how to retrieve CMIP6 data through the [Pangeo CMIP6](https://pangeo-data.github.io/pangeo-cmip6-cloud/)
* How to make a simple plot of a variable
* Regridd the CMIP6 variables to the exact horizontal resolution with [`xesmf`](https://xesmf.readthedocs.io/en/latest/)
* Calculate an ensemble mean of all used models
* Calculate and plot the seasonal mean 
  
## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="netcdf" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Plot the snowfall for January between 1985 - 2014 from NorESM2-MM
For this task, the CMIP6 experiment_id = 'historical'
> ### {% icon hands_on %} Hands-on: Open CMIP6 online catalog with Pangeo CMIP6
>
>   
>    > ### {% icon solution %} Solution
>    > ```
>    > cat_url = "https://storage.googleapis.com/cmip6/pangeo-cmip6.json"
>    > col = intake.open_esm_datastore(cat_url)
>    > col
>    > ```
> ### {% icon hands_on %} Hands-on: Search CMIP6 high resolution models (~100 km) 
> The variable we later want to plot is snowfall. 
> | shortname     |             Long name                   |      Units    |  levels |
> | ------------- |:---------------------------------------:| -------------:|--------:|
> |  prsn         |    Snowfall Flux                        | [kg m-2 s-1]  | surface |
> 
>    > ### {% icon solution %} Solution
>    > ```
>    > list_models = ['NorESM2-MM', 'TaiESM1', 'EC-Earth3-AerChem',
>    >                'GFDL-ESM4', 'SAM0-UNICON', 'CAS-ESM2-0',
>    >                'MPI-ESM1-2-HR', 'BCC-CSM2-MR', 'E3SM-1-0',
>    >                'E3SM-1-1', 'CMCC-CM2-SR5', 'CESM2-WACCM-FV2',
>    >                'CESM2', 'E3SM-1-1-ECA', 'GFDL-CM4', 'MRI-ESM2-0']
>    > variable_id = ['prsn]
>    > cat = col.search(source_id=list_models, experiment_id=['historical'], variable_id=variable_id[0], member_id=['r1i1p1f1'])
>    > cat.df
>    > ```
> ### {% icon hands_on %} Hands-on: Create dictonary from the list of datasets we found
>    > ### {% icon solution %} Solution
>    > ```
>    > dset_dict = cat.to_dataset_dict(zarr_kwargs={'use_cftime':True})
>    > ```
> ### {% icon hands_on %} Hands-on: Make a plot of mean surface snowfall for January 1985 - 2014
>    > ### {% icon solution %} Solution
>    > ```
>    > # open dataset from dictonary
>    > ds = dset_dict['CMIP.NCC.NorESM2-MM.historical.Amon.gn']
>    > ds
>    > ```
>    > Create 30-year mean for surface snowfall in January and convert the snowfall from [kg m-2 s-1] to [mm day-1]
>    > ```
>    > prsn_month = ds.prsn.groupby('time.month').mean('time')
>    > ```
>    > Plot 30-year for January for the Northern Hemisphere above 30$^o$N
>    > ```
>    > fig, ax = plt.subplots(1,1, figsize=[10,10], subplot_kw={'projection':ccrs.Orthographic(30, 90)})
>    > fig.suptitle('CMIP6 - high resolution (1985 - 2014)', fontsize=16, fontweight="bold")
>    > # Plot cosmetics 
>    > ax.coastlines()
>    > gl = ax.gridlines()
>    > ax.add_feature(cy.feature.BORDERS);
>    > gl.top_labels = False
>    > # Plot variable
>    > im = prsn_month.sel(month = 1).plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar = False,levels = np.arange(0.00, 7.25, 0.25), extend = 'max')
>    > 
>    > #add colorbar
>    > cb = fig.colorbar(im, orientation="vertical", fraction=0.046, pad=0.04)
>    > cb.set_label(label='Snowfall (mm$\,$day$^{-1}$)', weight='bold') 
>    > plt.tight_layout()
>    > fig.subplots_adjust(top=1)
>    > ```



# Regrid CMIP6 model output to NorESM2-MM grid using `xesmf`

> ### {% icon hands_on %} Hands-on: 
>
>   
>    > ### {% icon solution %} Solution
>    > ```



It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:


The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **My Tool**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **My Tool** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: File
>    - *"Parameter"*: `a value`
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


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
