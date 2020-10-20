---
layout: tutorial_hands_on
title: Functionally Assembled Terrestrial Ecosystem Simulator (FATES) with Galaxy Climate JupyterLab
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

questions:
- How to start Galaxy Climate JupyterLab in Galaxy?
- How to upload input data for running CLM-FATES for Nordic sites?
- How to create CLM-FATES case in Galaxy Climate JupyterLab?
- How to customize your runs?
- How to analyze your model outputs?
- How to upload observations to compare with your model outputs?
- How to save your model results into a Galaxy history?
- How to share your results?
objectives:
- Learn to use Galaxy Climate JupyterLab for setting up CLM-FATES
- Running CLM-FATES in Galaxy for single-point locations where in-situ measurements are available.
- Analyzing CLM-FATES output and comparing it with observations.
- Sharing results from the simulations along with corresponding in-situ data.
- Composing, executing and publishing the corresponding Galaxy workflow. 
time_estimation: 12H
key_points:
- Galaxy Climate JupyterLab
- CLM-FATES
- Model analysis
- Comparison with observations.
contributors:
- annefou

---


# Introduction
{:.no_toc}


The practical aims at familiarzing you with running CLM-FATES within Galaxy Climate JupyterLab. 

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

> ### {% icon comment %} Background
>
> Background about ESMs, CTSM and FATES. We can have slides too (TO DO separately).
>
> ##### Earth System Modelling (ESM)
>
> ##### The Community Terrestrial Systems Model
>
> ##### The Community Land Model
>
> ##### Functionally Assembled Terrestrial Ecosystem Simulator (FATES)
>
{:  .comment}

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *fates*.
>    {% include snippets/create_new_history.md %}
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>     TODO: input data for running FATES (so it can be run anywhere if not in data library).
>    https://zenodo.org/record/
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Check the datatype is **tar**
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 4. Rename Datasets
>
>    As "https://zenodo.org/record/?????/files/fates_emerald_inputdata.tar" is not a beautiful name and can give errors for some tools, it is a good practice to change the dataset name by something more meaningfull. For example by removing `https://zenodo.org/record/?????/files/` to obtain `fates_emerald_inputdata.tar`, respectively.
>
>    {% include snippets/rename_dataset.md %}
>
> 5. Add a tag to the dataset corresponding to `fates`
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

## Opening up Climate JupyterLab

> ### {% icon hands_on %} Hands-on: Launch JupyterLab for Ocean / Atmosphere / Land / Climate Python ecosystem in Galaxy
>
> Currently JupyterLab for Ocean / Atmosphere / Land / Climate Python ecosystem in Galaxy is available on [Live.useGalaxy.eu](https://live.usegalaxy.eu) only. JupyterLab for Ocean / Atmosphere / Land / Climate Python ecosystem and not the default JupyterLab in Galaxy contains all the python packages and additional software we need for running Earth System Model, including Functionally Assembled Terrestrial Ecosystem Simulator (FATES). The default JupyterLab in Galaxy would not be sufficient for executing all the tasks in th
is tutorial.
>
> 1. Open the JupyterLab tool {% icon tool %} by clicking [here](https://live.usegalaxy.eu/?tool_id=interactive_tool_climate_notebook){:target="_blank"}​ with the following par
ameters:
> 2. Click Execute
> 3. The tool will start running and will stay running permanently
> 4. Click on the "User" menu at the top and go to "Active Interactive Tools" and locate the JupyterLab instance you started.
> 5. Click on your JupyterLab instance (please not that it may take a few minutes before you can click on the link to your jupyterLab instance).
>
{: .hands_on}


You should now be looking at a page with the JupyterLab interface:

![Jupyterlab climate session](../../images/jupyterlab_climate_session.png)


# Import input data to JupyterLab

We will be using a JupyterLab Terminal for most of this tutorial.

> ### {% icon hands_on %} Hands-on: Open a JupyterLab Terminal
> To open a new terminal, select theterminal in the new Launcher tab. 
> More information can be found on the [JupyterLab documentation on Terminals](https://jupyterlab.readthedocs.io/en/stable/user/terminal.html).
>
> Import the FATES input dataset from your history:
> ```
> get -i fates_emerald_inputdata.tar -t name
> ```
>
> Then untar this file:
> ```
> tar xf /import/fates_emerald_inputdata.tar --directory $HOME/
> ```
{: .hands_on}

# CLM-FATES single point simulations

## Get CLM-FATES EMERALD code

> ### {% icon hands_on %} Hands-on: Clone CLM-FATES for Nordic sites
> 
> ```
> git clone -b release-emerald2.0.0 https://github.com/NordicESMhub/ctsm.git
> cd ctsm
> ./manage_externals/checkout_externals
> ```
>
{: .hands_on}

## Create CLM-FATES new case

> ### {% icon hands_on %} Hands-on: Create CLM-FATES new case for ALP1 site
>
> ```
> source activate cesm
>
> cd cime/scripts
> ./create_newcase --case ../../../ctsm_cases/fates_alp1 --compset 2000_DATM%1PTGSWP3_CLM50%FATES_SICE_SOCN_MOSART_SGLC_SWAV --res 1x1_ALP1 --machine espresso --run-unsupported
> ```
>
{: .hands_on}

## Setup, build and submit your first simulation

> ### {% icon hands_on %} Hands-on: Setup, build and submit
>
> ```
> cd ctsm_cases/fates_alp1
> ./case.setup
> ./case.build
> ./xmlchange DOUT_S=FALSE
> ./case.submit > case_submit.out 2>&1
> ```
{: .hands_on}

## Customize your run

> ### {% icon hands_on %} Hands-on: Run 5 years
>
> ```
> ./xmlchange --file env_run.xml --id RUN_STARTDATE --val 0001-01-01      # set up the starting date of your simulation 
> ./xmlchange --file env_run.xml --id STOP_OPTION --val nyears            # set the simulation periods to "years"
> ./xmlchange --file env_run.xml --id STOP_N --val 5                      # set the length of simulation, i.e, how many years
> ./xmlchange --file env_run.xml --id CONTINUE_RUN --val TRUE             # if you want to continue your simulation from restart file, set it to TRUE
> ./xmlchange --file env_run.xml --id RESUBMIT --val 10                   # set up how many times you want to resubmit your simulation.
>                                                                         # e.g, STOP_N=5, RESUBMIT=10, you will have simulation for 5+5*10=55 
> ./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_START --val 1901    # set up the start year of the atmospheric forcing 
> ./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_END --val 1950      # set up the end year of the atmospheric forcing
> ./xmlchange DOUT_S=TRUE
> ./case.submit > case_submit_sontinue_run.out 2>&1
> ```
>  
{: .hands_on}

# Analysis

## Analyzing FATES-CLM model outputs

> ### {% icon hands_on %} Hands-on: Open a new Python notebook
> Create a notebook by clicking the `+` button in the file browser and then selecting a kernel in the new Launcher tab:
> Get more information online at [JupyterLab notebooks](https://jupyterlab.readthedocs.io/en/stable/user/notebook.html).
{: .hands_on}

### Use `xarray` to read and plot

In this section, we give an example on how to visualize your results using `xarray`:

```
import xarray as xr

xr.set_options(display_style="html")
%matplotlib inline

dset = xr.open_dataset("x.nc")
## Comparisons with observations

> ### {% icon hands_on %} Hands-on: Import observations into your JupyterLab session
>
>
{: .hands_on}

# Save your results to your Galaxy history

Open a JupyterLab Terminal as the following commands will be executed from the command line.

> ### {% icon hands_on %} Hands-on: Put your data to your Galaxy history
>
> ```
> cd /home/jovyan/
> tar cvf archive_emerald_fates_test.tar archive
> ```
> Then you are now ready to put your dataset into Galaxy. As it can be large, we recommend to use FTP:
> ```
> curl -T {"archive_emerald_fates_test.tar"} ftp://ftp.usegalaxy.eu --user USER:PASSWORD --ssl
> ```
> Where you replace `USER` by your galaxy username (what you used to login to Galaxy and `PASSWORD` by your Galaxy login password.
>
{: .hands_on}

# Share your work

One of the most important features of Galaxy comes at the end of an analysis. When you have published striking findings, it is important that other researchers are able to reproduce your in-silico experiment. Galaxy enables users to easily share their workflows and histories with others.

To share a history, click on the {% icon galaxy-gear %} icon in the history panel and select `Share or Publish`. On this page you can do 3 things:

1. **Make History Accessible via Link**. This generates a link that you can give out to others. Anybody with this link will be able to view your history.
2. **Make History Accessible and Publish**. This will not only create a link, but will also publish your history. This means your history will be listed under `Shared Data → Histories` in the top menu.
3. **Share with a user**. This will share the history only with specific users on the Galaxy instance.

> ### {% icon comment %} Permissions
> Different servers have different default permission settings. Some servers create all of your datasets completely private to you, while others make them accessible if you know the secret ID.
>
> Be sure to select **Also make all objects within the History accessible** whenever you make a history accessible via link, otherwise whomever you send your link to might not be able to see your history.
{: .comment}

> ### {% icon hands_on %} Hands-on: Share history
>
> 1. Share your history with your neighbour.
> 2. Find the history shared by your neighbour. Histories shared with specific users can be accessed by those users under their top masthead "User" menu under `Histories shared with me`.
{: .hands_on}



# Conclusion

{:.no_toc}

We have learnt to run single-point simulations with FATES-CLM through the Galaxy Climate JupyterLab.
