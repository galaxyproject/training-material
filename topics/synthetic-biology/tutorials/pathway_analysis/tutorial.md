---
layout: tutorial_hands_on

title: Evaluating and ranking a set of pathways based on multiple metrics
zenodo_link: https://zenodo.org/record/6628296
questions:
- How to evaluate a set of heterologous pathways ?
objectives:
- Calculate the production flux of the lycopene target molecule using _Flux Balance Analysis_ tool.
- Compute thermodynamics values to optimize the yield of the reaction producing the lycopene using _Thermo_ tool.
- Compute the global score for the previous annotated pathways.
- Rank the computed heterologous pathways depending on their score.
time_estimation: 20M
contributors:
- kenza12
- tduigou
- breakthewall
- guillaume-gricourt
- ioanagry
- jfaulon

---


# Introduction


Progress in synthetic biology is enabled by powerful bioinformatics tools such as those aimed to design metabolic pathways for the production of chemicals. These tools are available in SynBioCAD portal which is the first Galaxy set of tools for synthetic biology and metabolic engineering ({% cite Hrisson2022 %}).

In this tutorial, we will use a set of tools from the **Pathway Analysis workflow** which will enable you to evaluate a set of heterelogous pathways previously produced by the RetroSynthesis workflow in a chassis organism (_E. coli_ iML1515). These workflows are available in the [Galaxy SynbioCAD platform](https://galaxy-synbiocad.org). The goal is to inform the user of the theoretically best performing pathways by ranking them based on the four following criteria: target product flux, thermodynamic feasibility, pathway length and enzyme availability.

We recommend that you follow the Retrosynthesis tutorial before starting the current tutorial which will enable you to find pathways to synthesize heterologous compounds producing Lycopene in the _E. coli_ chassis organism.

Four main steps will be run using the following workflow:

To rank the computed heterologous pathways, we need to calculate some metrics. This is why an in-house Flux Balance Analysis (FBA) was developed to calculate the production flux of a given target (e.g. lycopene). The method forces a fraction of its maximal flux through the biomass reaction while optimizing for the target molecule. This is achieved by the _Flux Balance Analysis_ tool.

Secondly, we will use the _Thermo_ tool to estimate thermodynamics values (based on Gibbs free energies) for each pathway to know whether a producing pathway is feasible in physiological conditions. The contribution of individual reactions to the final pathway thermodynamic is balanced solving a linear equation system.

After that, the _Score Pathway_ tool is used to calculate a global score combining target flux, pathway thermodynamics, pathway length and enzyme availability.

Finally, the pathway are ranked based on the global score using the _Rank Pathways_ tool.

![This image shows the pathway analysis workflow viewed from the Galaxy editor interface. The pathway analysis workflow will evaluate and rank pathways based on multiple metrics. First, Flux Balance Analysis (FBA) node will take as inputs a collection of SBMLs representing the heterologous pathways to be evaluated, another SBML representing the chassis, and 2 additional inputs are expected for indicating the biomass reaction ID and the cell compartment ID that will be used during the FBA analysis. The FBA method will optimize the production of the target under the constraint that the model should still produce some level biomass. Then, the Thermo tool is used to estimate thermodynamics values for each pathway to know whether it is favorable towards target production in physiological conditions. After that the Score Pathway tool will calculate a global score combining target flux, pathway thermodynamics, pathway length and enzyme availability. Finally, the pathways are ranked based on the global score using the Rank Pathways tool. All exchange files are SBMLs. The output is a CSV file which contains the pathway IDs and their corresponding global score.](../../images/pathway_analysis_workflow.jpg)

Note that we will run the steps of this workflow individually so as not to neglect the understanding of the intermediate steps. Then, we will run the workflow automatically so that it itself retrieves the outputs from the previous step and gives them as input to the next tool.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Data Preparation

First we need to upload and prepare the following inputs to analyze:

- A set of pathways provided in the SBML format (Systems Biology Markup Language) to be ranked, modeling heterologous pathways such as those outputted by the **RetroSynthesis workflow** (available in [Galaxy SynbioCAD platform](https://galaxy-synbiocad.org)).

- The GEM (Genome-scale Metabolic Models) which is a formalized representation of the metabolism of the host organism (the model is _E. coli_ iML1515), provided in the SBML format.

## Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial named *Pathway Analysis*.
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) :
>
>    ```
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_001_0001.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_001_0006.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_001_0011.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_002_0001.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_002_0006.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_002_0011.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_003_0001.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_003_0116.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/rp_003_0231.xml
>    https://zenodo.org/api/files/5db78fa1-b8cb-4046-b57c-8a9d00806f42/SBML_Model_iML1515.xml
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Create a list or collection named `Heterologous pathways` and composed of the 9 rpSBML pathways.
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
{: .hands_on}

# Compute the target product flux

Notice that the starting compounds (in other words, _the precursors_) of the predicted pathways (also referred as the _heterologous pathways_) are compounds that have been initially extracted from the genome-scale metabolic model (GEM) of the organism we are interested in (also referred as _chassis_). While this step is out of the scope of the present Pathway Analysis tutorial, this means that the precursors of predicted pathways are also present in the chassis model. Hence, predicted pathways and the chassis organism model can be merged to construct "augmented" whole-cell models, enabling flux analysis of these metabolic systems. This is what we'll do here to predict the production flux of a compound of interest. 

Within the frame of this tutorial, we'll use the _E. coli_ iML1515 GEM (downloaded from the [BiGG database](http://bigg.ucsd.edu/)) to model the chassis metabolism of _E. coli_ and the target compound is the lycopene. The provided _E. coli_ model is in the SBML. The extraction of precursor compounds and the pathway prediction have already been performed during the RetroSynthesis workflow (available in [Galaxy SynbioCAD platform](https://galaxy-synbiocad.org)). 

The FBA (Flux Balance Analysis) method used to calculate the flux is a mathematical approach (as decribed in section Methods in {% cite Hrisson2022 %}) which uses the COBRApy package ({% cite Ebrahim2013 %}) and proposes 3 different analysis methods (standard FBA, parsimonious FBA, fraction of reaction). The first two methods are specific to the COBRApy package and the last one `Fraction of Reaction` is an in-house analysis method (as decribed in section Methods in {% cite Hrisson2022 %}) to consider the cell needs for its own maintenance while producing the target compound.

Within the workflow, the purpose of the _Flux Balance Analysis_ tool is to predict the production flux of the targeted compound, while considering the cellular needs. Under such simulation conditions, the analysis that returns a low production flux may be due to some precursor compounds having a limiting production flux, nor cofactor fluxes involved not being sufficiently balanced by the chassis native metabolism. Pathways with high flux would be caused by both the precursor compounds and the cofactors being in abundance. In either case, bottlenecks that limit the flux of the pathway may be investigated (this is outside of the scope of the workflow) and pathways that do not theoretically generate high yields can be filtered out.

We first perform an FBA (with COBRApy) optimizing the biomass reaction and record its maximal theoretical flux. The upper and lower bounds of the biomass reaction are then set to a same amount, equals to a fraction of its previously recorded optimum (default is 75% of its optimum). The method then performs a second FBA where biomass flux is enforced to this fraction of its optimum while optimizing the target production flux. Simulated fluxes are recorded directly into the SBML file and all changed flux bounds are reset to their original values before saving the output file.

![This picture describes the process to obtain annotated SBML pathways with calculated fluxes. First, FBA tool takes as input one SBML representing the heterologous pathway and another SBML representing the chassis. The two SBMLs are merged to render an augmented model containing both the reactions of the heterologous pathway and the chassis. Then FBA tool uses the COBRApy package to optimize the producing flux of the target reaction, under the constraint that the flux of the biomass reaction should be equals to 75% of its maximal theoretical value. At the end, the calculated fluxes are recorded as annotations into the SBML of the heterologous pathway.](../../images/fba_calculation.png)

> <details-title>Comment</details-title>
>
> Blocking compounds that cannot provide any flux are temporarily removed from heterologous reactions for the FBA evaluation. Such cases can happen due to side substrates or products of predicted reactions that do not match any chassis compound, representing dead-end paths.
>
{: .details}

> <hands-on-title>Calculating the flux of a target using Flux Balance Analysis (FBA)</hands-on-title>
>
> 1. Run {% tool [Flux balance analysis](toolshed.g2.bx.psu.edu/repos/iuc/rpfba/rpfba/5.12.1) %} with the following parameters:
>    - {% icon param-collection %} *"Pathway (rpSBML)"*: Select `Heterologous pathways` (Input dataset collection) from your current history.
>    - {% icon param-file %} *"Model (SBML)"*: Select `SBML_Model_iML1515.xml` (Input dataset) from your current history.
>    - *"SBML compartment ID"*: Leave the default value `c`.
>
>    > <comment-title>Choose a compartment corresponding to your model</comment-title>
>    >
>    > You can specify the compartment from which the chemical species were extracted.
>    > The default is `c`, the BiGG code for the cytoplasm.
>    {: .comment}
>
>    - *"biomass reaction ID"*: Specify the biomass reaction ID that will be used for the "fraction" simulation, type `R_BIOMASS_Ec_iML1515_core_75p37M`.
>
>    > <comment-title>How to select the biomass reaction ID ?</comment-title>
>    >
>    > The biomass reaction ID objective is extracted from the current model *E.Coli iML1515*. You can search the term `biomass` in your XML model and pick the ID where the term `core` appears.
>    {: .comment}
>
>    - *"Constraint based simulation type"*: `Fraction of Reaction`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the FBA score for `rp_003_0001` pathway ?
>
> > <solution-title></solution-title>
> >
> > 1. View the SBML rp_003_0001 file and look for `fba_fraction` value in `<groups:listOfGroups>` section: value= `0.23693089430893874`.
> >
> {: .solution}
>
{: .question}

## Compute thermodynamics values

The goal of the thermodynamic analysis is to estimate the feasibility of the predicted pathways toward target production, in physiological conditions. The eQuilibrator libraries ({% cite Flamholz_2011 %}) are used to calculate the formation energy of compounds by either using public database IDs (when referenced within the tools internal database) or by decomposing the chemical structure and calculating its energy of formation using the component contribution method.

The reaction Gibbs energy is estimated by combining the energy of formation of the compounds involved in the reaction (with consideration for the stoichiometric coefficients).

Finally, the thermodynamic of a pathway is estimated by combining the Gibbs energy of reactions involved in it. The contribution of individual reactions to the final pathway thermodynamic is balanced using a linear equation system, according to the relative uses of intermediate compounds across the pathway (See Thermodynamics in Methods section for further details: {% cite Hrisson2022 %}). A pathway Gibbs energy below zero indicates that the thermodynamic is favorable toward the production of the target.

![This image describe the process of computing thermodynamics values of a pathway, to assess whether a pathway is favorable towards target production in physiological conditions. rpThermo tool takes as input the output of rpFBA, which corresponds to the pathway to be analyzed (as an SBML file). rpThermo relies on eQuilibrator libraries to calculate the formation energy of compounds. The reaction Gibbs energy is estimated by combining the energy of formation of the compounds involved in the reaction (with consideration for the stoichiometric coefficients). From this, the thermodynamic of a pathway is estimated by combining the Gibbs energy of all reactions taking part in it. The individual contribution of reactions to the final pathway thermodynamic is weighted using a linear equation system, to take into consideration of the relative needs of intermediate compounds across the pathway. For the record, a pathway Gibbs energy below zero indicates that the thermodynamic is favorable toward the production of the target.](../../images/rpthermo.png)

Secondly, we will use the _Thermo_ tool to estimate thermodynamics values (based on Gibbs free energies) for each pathway to know whether a producing pathway is feasible in physiological conditions

> <hands-on-title>Compute thermodynamics values for each pathway using rpThermo tool</hands-on-title>
>
> 1. {% tool [Thermo](toolshed.g2.bx.psu.edu/repos/tduigou/rpthermo/rpthermo/5.12.1) %} with the following parameters:
>    - {% icon param-file %} *"Input File"*: `pathway_with_fba` (output of **Flux balance analysis** {% icon tool %})
>
>    > <comment-title></comment-title>
>    >
>    > The tool takes as input pathways in SBML format and returns annotated pathways (with thermodynamics information for each reaction) in SBML format too.
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the thermodynamic score attributed to the reaction with the following *EC (Enzyme Commission) number* *2.5.1.29* for `rp_001_0001` pathway ?
>
> > <solution-title></solution-title>
> >
> > 1. View the SBML rp_001_0001 file and search the reaction ID `2.5.1.29` contained in `<listOfReactions>`. The corresponding value is indicated in `thermo_dGm_prime` field : `-242.348`.
> >
> {: .solution}
>
{: .question}

## Compute the global score of pathways

The _Pathway Score_ tool provides a global score for a given pathway previously annotated by the _Flux Balance Analysis_ and _Thermo_ tools. This score is computed by a machine learning (ML) model (cf. Machine Learning Global Scoring in {% cite Hrisson2022 %}). The model takes as input features describing the pathway (thermodynamic feasibility, target flux with fixed biomass, length) and the reactions within the pathway (reaction SMARTS, Gibbs free energy, enzyme availability score) and prints out the probability for the pathway to be a valid pathway. The ML model has been trained on literature data (cf. section Benchmarking with literature data in {% cite Hrisson2022 %}) and by a validation trial (cf. section Benchmarking by expert validation trial in {% cite Hrisson2022 %}).

> <hands-on-title>Compute the global score using the _Pathway Score_ tool</hands-on-title>
>
> 1. {% tool [Score Pathway](toolshed.g2.bx.psu.edu/repos/tduigou/rpscore/rpscore/5.12.1) %} with the following parameters:
>    - {% icon param-file %} *"Pathway (rpSBML)"*: `pathway_with_thermo` (output of **Thermo** {% icon tool %})
>
>    > <comment-title></comment-title>
>    >
>    > This tool will output a new annotated SBML file representing the pathway, containing the `global_score` annotation.
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the computed global score for the `rp_001_0001` pathway ?
>
> > <solution-title></solution-title>
> >
> > 1. View the SBML file `rp_001_0001` and search for `global_score` : value=`0.975147980451584`.
> >
> {: .solution}
>
{: .question}

## Rank annotated pathways

Finally, _Rank Pathways_ ranks the previous set of heterologous pathways, based on their global score, to reveal what are the most likely pathways to produce the target molecule (here it is _lycopene_) in a given organism of interest (_E. coli_ in this tutorial).

> <hands-on-title>Rank annotated pathways using rpRanker tool</hands-on-title>
>
> 1. {% tool [Rank Pathways](toolshed.g2.bx.psu.edu/repos/tduigou/rpranker/rpranker/5.12.1) %} with the following parameters:
>    - {% icon param-file %} *"Pathways"*: `scored_pathway` (output of **Score Pathway** {% icon tool %})
>
>    > <comment-title></comment-title>
>    >
>    > This tool will output a CSV file which contains the pathway IDs and their corresponding global score.
>    {: .comment}
>
{: .hands_on}


> <question-title></question-title>
>
> 1. What are the 3 top-ranked pathways ?
>
> > <solution-title></solution-title>
> >
> > 1. `002_0011`, `001_0011` ,`002_0006`.
> >
> {: .solution}
>
{: .question}

# Run the **Pathway Analysis Workflow**

In this section, you can run the Pathway Analysis Workflow more easily and fastly following these instructions:

> <hands-on-title>Execute the entire workflow in one go.</hands-on-title>
>
> 1. Import your **Pathway Analysis Workflow** by uploading the [**workflow file**](https://training.galaxyproject.org/training-material/topics/synthetic-biology/tutorials/pathway_analysis/workflows/main_workflow.ga).
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
> 2. Click on *Workflow* on the top menu bar of Galaxy. You will see **Pathway Analysis Workflow**
> 3. Click on the {% icon workflow-run %} (*Run workflow*) button next to your workflow
> 4. Provide the workflow with the following parameters:
>    - {% icon param-file %} *"Heterologous pathways":* Select `Heterologous pathways` (Input dataset collection) from your current history.
>    - {% icon param-file %} *"Chassis where to produce target from"*: Select `SBML_Model_iML1515.xml` (Input dataset) from your current history.
>    - *"Cell compartment ID"*: Enter value `c`.
>    - *"biomass reaction ID"*: Specify the biomass reaction ID that will be restricted in the "fraction" simulation type `R_BIOMASS_Ec_iML1515_core_75p37M`.
>
>    > <comment-title></comment-title>
>    >
>    > All the outputs will be automatically generated and identical to the previous ones. 
>    {: .comment}
{: .hands_on}

# Conclusion


To select the best pathways for producing the lycopene in *E. coli*, some metrics have to be estimated, namely production flux of the target and pathway thermodynamics. A global score is then computed by combining these criteria with others (pathway length, enzyme availability score, reaction SMARTS) using a machine learning model. These steps achieved using the tools of the presented Pathway Analysis workflow.

![This scheme represents the pathway analysis workflow enabling the identification of the best pathways for producing a molecule of interest. To do that the workflow takes as input the collection of pathways to be scored and the metabolic model of the chassis (all files are SBMLs). Iteratively, each pathway is merged with the chassis model, and several metrics are evaluated such as the target production flux (using the rpFBA tool) and the thermodynamics of the pathway (using the rpThermo tool). A global score of each pathway combining these metrics with others is computed using a machine learning method (with the rpScore tool). Finally, the pathways are ranked from best to worst according to their global score (using the rpRanker tool). During the workflow, all metrics are stored as SBML annotations. The final output is a CSV file which contains the pathway IDs and their corresponding global score.](../../images/pathway_analysis_scheme.png)