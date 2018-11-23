---
layout: tutorial_hands_on

title: Running Molecular Dynamics simulations
zenodo_link: ''
questions:
- Which molecular dynamics engine can I use to run simulations?
- How to start molecular dynamics simulations?
- Which bioinformatics techniques is important to know for this type of data?
objectives:
- Learn about available MD engines provided in BRIDGE
requirements:
  -
    title: "Setting up molecular systems"
    type: "internal"
    link: "/computational-chemistry/tutorials/setting-up-molecular-systems/tutorial.html"
time_estimation: 3H
key_points:
- Several MD engines are available in BRIDGE
- Workflows are available for common simulations tasks such as equilibration and production dynamics for various ensembles (NVT, NPT)
contributors:
- chrisbarnettster
- tsenapathi

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

***TODO***: *Improve the introduction*

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Simulation with NAMD (imported from CHARMM-GUI)

This workflow uses NAMD as a molecular dynamics engine.An NVT simulation is followed by an NPT simulation.



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
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

If you've imported from CHARMM-GUI read through the theory then skip to [namd_nvt](#namd_nvt) otherwise carry on and start with [setup](#setup). The workflows as shown below for [CHARMM and NAMD](#figure-1) or [only NAMD](#figure-2) will be used.



![Snapshot of CHARMM and NAMD analysis workflow](images/NAMD_CHARMMGUI_workflow.png"A simple simulation workflow starting with CHARMM for setup and NAMD to continue the production simulation")
![Snapshot of NAMD analysis workflow](images/NAMD_workflow.png "A simple NAMD simulation workflow")




***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} NVT theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

> ### {% icon details %} NPT theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

## **setup**
{: #setup}

Prepare a protein ligand system in CHARMM. 
This tool will:
- solvate the protein-ligand comples in TIP3P waters
- neutralise the system by using 0.05M NaCl
- conduct a short minimisation 

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **setup** {% icon tool %} with the following parameters:
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
> 1. Why use 0.05M NaCl?
> 2. Why use TIP3P water?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## **minimizer**

Reduces the energy of the system by running a minimisation algorithm. 
This tools will:
- Minimise energy using steepest descent followed by Adopted Newton Raphson (using the define number of steps_
- Setup periodic boundaries and generate Particle Mesh Ewald (PME)
- generate reference structures for restraints in NAMD (if selected)

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **minimizer** {% icon tool %} with the following parameters:
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

## **namd_nvt**
{: #namd_nvt}
Classical NVT dynamics.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **namd_nvt** {% icon tool %} with the following parameters:
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

## **namd_npt**
Classical NPT dynamics.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **namd_npt** {% icon tool %} with the following parameters:
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
  1. Is this converged?
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

# Conclusion
{:.no_toc}

