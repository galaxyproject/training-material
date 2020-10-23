---
layout: tutorial_hands_on
title: Functionally Assembled Terrestrial Ecosystem Simulator (FATES)
zenodo_link: 'https://doi.org/10.5281/zenodo.4108341'
requirements:
  -
    type: "internal"
    topic_name: introduction
    tutorials:
        - galaxy-intro-101-everyone

questions:
- Who to run CLM-FATES with CLM-FATES Galaxy tool?
- How to upload input data for running CLM-FATES?
- How to customize your runs?
- How to analyze your model outputs?
- How to create a workflow for running multi-site configurations?
- How to share your workflow?
objectives:
- Setting up CLM-FATES case.
- Customizing your runs.
- Interactive visualization with Panoply
- Automating your analyzes and visualisations of your CLM-FATES case.
- Creating multi-case scenarios.
- Composing, executing and publishing CML-FATES workflow..
time_estimation: 12H
key_points:
- CLM-FATES
- Quick visualization of your results with Panoply
- Create multi-site simulations with a Galaxy workflow
contributors:
- annefou

---


# Introduction
{:.no_toc}


The practical aims at familiarzing you with running CLM-FATES in Galaxy and analyzing the model results. 
It will also teach you on how to create Galaxy workflow for your CLM-FATES simulations to make your research fully reproducible.

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
> FATES is the “Functionally Assembled Terrestrial Ecosystem Simulator”.
> FATES needs what we call a "Host Land Model" (HLM) to run and in this tutorial
> we will be using the [Community Land Model](http://www.cesm.ucar.edu/models/clm/)
> of the [Community Terrestrial Systems Model](https://github.com/ESCOMP/CTSM) (CLM-CTSM).
> FATES was derived from the CLM Ecosystem Demography model (CLM(ED)), which was documented in 
> {% cite Fisher2015 %}.
> and this technical note was first published as an appendix to [that paper](https://pdfs.semanticscholar.org/396c/b9f172cb681421ed78325a2237bfb428eece.pdf).
> The [FATES documentation](https://fates-docs.readthedocs.io/en/latest/index.html) will provide some more insight on FATES too.
>
{:  .comment}

## Step 1: Get CLM-FATES input data

Preparing CLM-FATES input data is out of scope for this tutorial. We assume the input data tarball contains the following folders:

```
atm   cpl   lnd   share
```

Each sub-folder will then contain all the necessary inputs for running your CLM-FATES case.
For the purpose of this tutorial, input data for a single point location ALP1 (61.0243N, 8.12343E) has been prepared and is ready to use.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial. If you are not inspired, you can name it *fates*.
>    {% include snippets/create_new_history.md %}
> 2. Import the files from [Zenodo](https://doi.org/10.5281/zenodo.4108341) or from the shared data library
>
>    ```
>    https://zenodo.org/record/4108341/files/inputdata_version2.0.0_ALP1.tar
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
>    As "https://zenodo.org/record/4108341/files/inputdata_version2.0.0_ALP1.tar" is not a beautiful name and can give errors for some tools,
>    it is a good practice to change the dataset name by something more meaningfull. For example by removing `https://zenodo.org/record/4108341/files/` to obtain `inputdata_version2.0.0_ALP1.tar`, respectively.
>
>    {% include snippets/rename_dataset.md %}
>
> 5. Add a tag to the dataset corresponding to `fates`
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Step 2: Setting up a CLM-FATES simulation

We will be using the CTSM/FATES-EMERALD Galaxy tool to evaluate ???

> ## {% icon comment %} Tip: Finding your tool
>
> Different Galaxy servers may have tools available under different sections, therefore it is often useful to use the **search bar** at the top of the tool panel to find your tool.
>
> Additionally different servers may have multiple, similarly named tools which accomplish similar functions. When following tutorials, you should use precisely the tools that they describe. For real analyses, however, you will need to search among the various options to find the one that works for you.
>
{: .comment}

> ### {% icon hands_on %} Hands-on: Creating a new CTSM/FATES-EMERALD case
>
> 1. {% tool [CTSM/FATES-EMERALD](ctsm_fates) %} with the following parameters:
>    - {% icon param-file %} *"inputdata for running FATES EMERALD"*: select the **inputdata_version2.0.0_ALP1.tar** file from your history
>    - *Name of your case*: ALP1_exp
>    - *Expand 'Customize the model run period' to change the default values:
>        - **Determines the model run initialization type.**: select **hybrid**
>        - **Reference case for hybrid or branch runs**: ALP1_refcase
>        - **Reference date for hybrid or branch runs (yyyy-mm-dd)**: 2100-01-01
>        - **Run start date (yyyy-mm-dd). Only used for startup or hybrid runs.**: 2100-01-01
>        - **restart for running FATES EMERALD**: ALP1_refcase_2100-01-01.tar
>        - **Provides a numerical count for STOP_OPTION.**: 5
>        - **Sets the run length along with STOP_N and STOP_DATE**: nyears
>    - Click **Execute**
>
>    > ### {% icon comment %} Tip: search for the tool
>    >
>    > Use the **tools search box** at the top of the tool panel to find **Remove beginning** {% icon tool %}.
>    {: .tip}
>
>    > ## {% icon comment %} startup versus hybrid
>    >
>    >  TODO: Explain here why we want to start from an hybrid and not startup
>    >
>    {: .comment}
>
> 2. Check that the datatype is **netcdf**
>
>    Files you uploaded are in netcdf format. In Galaxy, Datatypes are, by default, automatically guessed. Here, as necdf is a derivative of the h5 format, Galaxy automatically affect the h5 datatype to netcdf files. To cope with that, one can change the datatype manually, once datasets uploaded (as shown below) OR you can directly specify datatype on the upload tool form so Galaxy will not try to automatically guess it.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 3. **Rename** {% icon galaxy-pencil %} the output dataset to `ALP1.nc`
>
>    {% include snippets/rename_dataset.md %}
>
> 4. Click on the new history item to expand it
>
>    > ### {% icon question %} Questions
>    >
>    > 1. Which tags are present on this resulting dataset? (You may have to refresh the history panel to see the tags)
>    > 2. 
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Step 3: Quick visualization with Panoply

## Opening up Panoply

> ### {% icon hands_on %} Hands-on: Launch Panoply
>
>  Panoply is available as a Galaxy interactive environment and may not be available on all Galaxy servers.
>
> > ### {% icon tip %} Tip: Launch Panoply in Galaxy
> > Currently Panoply in Galaxy is available on useGalaxy.eu instance, on the "Interactive tools" tool panel section or, as all interactive tools, from the dedicated usGalaxy.eu subdomain: [Live.useGalaxy.eu](https://live.usegalaxy.eu)
> >
> > 1. Open the Panoply tool {% icon tool %} by clicking [here](https://live.usegalaxy.eu/?tool_id=interactive_tool_panoply){:target="_blank"}
> > 2. Check **ALP1.nc** dataset selected in the netcdf input field
> > 3. Click Execute
> > 4. The tool will start running and will stay running permanently
> > 5. Click on the "User" menu at the top and go to "Active Interactive Tools" and locate the Panoply instance you started.
> > 6. Click on your Panoply instance
> >    ![Panoply dataset selection](../../images/select_dataset.png "Select dataset")
> > 7. Click on **ALP1.nc** dataset
> {: .tip}
{: .hands_on}

## Inspect metadata

> ### {% icon hands_on %} Hands-on: Inspect dataset
>
> 1. Inspect dataset content
> 
>    Here you can look at the dataset (ALP1.nc) and related variables (FSDS, FSA, AREA_TREE, BIOMASS_CANOPY, etc) 
>
>    > ### {% icon question %} Question
>    >
>    > ?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > {: .solution}
>    {: .question}
>
>
> 2. Inspect the area occupied by woody plants (AREA_TREE) variable
>
>    > ### {% icon question %} Question
>    >
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

# Step 4: Using Galaxy tools for analysing your CLM-FATES simulation

# Step 5: Generate a Galaxy workflow

# Step 6: Rerun your workflow with a different location 

Use any of the available input datasets (ALP2, etc.).

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

We have learnt to run single-point simulations with FATES-CLM and generate workflows for multi-site scenarios.
