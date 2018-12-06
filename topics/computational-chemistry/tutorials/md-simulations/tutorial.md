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
- You've run an equilibrium and production MD simulation using NAMD
contributors:
- chrisbarnettster
- tsenapathi

---


# Introduction
{:.no_toc}


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Simulation with NAMD (imported from CHARMM-GUI)

This process can be accomplished by selected each tools from the tools menu or by importing the workflow. The workflow method is most efficient and the individual tools used in the workflow are discussed below. 


If you've imported from CHARMM-GUI read through the theory then skip to [NAMD_MD](#namd_md) otherwise carry on and start with [setup](#setup). The workflows as shown below for [CHARMM and NAMD](#figure-1) or [only NAMD](#figure-2) will be used.



![Snapshot of CHARMM and NAMD analysis workflow](images/NAMD_CHARMMGUI_workflow.png"A simple simulation workflow starting with CHARMM for setup and NAMD to continue the production simulation")
![Snapshot of NAMD analysis workflow](images/NAMD_workflow.png "A simple NAMD simulation workflow")



> ### {% icon details %} NVT, NPT and statistical mechanics theory
>
> [See Statistical Mechanics, McQuarrie for more in depth theory ISBN:9781891389153](https://books.google.co.za/books/about/Statistical_Mechanics.html?id=itcpPnDnJM0C&redir_esc=y)
>
{: .details}

## Additional prep with CHARMM (SKIP)

### **setup**
{: #setup}

Prepare a protein ligand system in CHARMM. 
This tool will:
- solvate the protein-ligand comples in TIP3P waters
- neutralise the system by using 0.05M NaCl
- conduct a short minimisation 

> ### {% icon hands_on %} Hands-on:
>
> 1. **setup** {% icon tool %} with the following parameters:
>
>
{: .hands_on}

### **minimizer**

Reduces the energy of the system by running a minimisation algorithm. 
This tools will:
- Minimise energy using steepest descent followed by Adopted Newton Raphson (using the define number of steps_
- Setup periodic boundaries and generate Particle Mesh Ewald (PME)
- generate reference structures for restraints in NAMD (if selected)

> ### {% icon hands_on %} Hands-on: 
>
> 1. **minimizer** {% icon tool %} with the following parameters:
>
>
{: .hands_on}


## NAMD MD workflow
{: #namd_md}

This workflow uses NAMD as a molecular dynamics engine. An NVT simulation is followed by an NPT simulation.

> ### {% icon hands_on %} Hands-on: Access the workflow
> Access the published workflows
> ![List of published workflows](images/published_workflows_update.png "Published workflows")
> Choose to import the NAMD MD workflow from published workflows
> ![Import workflow](images/import_workflow.png "Import workflow")
> Choose to run a workflow from your available workflows
> ![Your workflows](images/workflows_yours.png "Your workflows")
{: .hands_on}



### **namd_nvt**
{: #namd_nvt}
Classical NVT dynamics.

> ### {% icon hands_on %} Hands-on: NVT dynamics
>
> 1. **namd_nvt** {% icon tool %} with the following parameters:
>
> This tool runs classical molecular dynamics simulations in NAMD using an NVT ensemble. User can run the simulation in small time intervals. The coordinates, velocities and the extended system files can be use to restart the simulations. Harmonic restraints can be used. NAMD collective variable module is used to give RMSD harmonic restraints.
> ![Snapshot of NAMD NVT tool parameters part 2](images/namd_nvt_tool_params.png "NAMD NVT parameters")
{: .hands_on}



### **namd_npt**
Classical NPT dynamics.

> ### {% icon hands_on %} Hands-on: NPT dynamics
>
> 1. **namd_npt** {% icon tool %} with the following parameters:
>
> This tool runs classical molecular dynamics simulations in NAMD using an NPT ensemble. User can run the simulation in small time intervals. The coordinates, velocities and the extended system files can be use to restart the simulations. Harmonic restraints can be used. NAMD collective variable module is used to give RMSD harmonic restraints.
> ![Snapshot of NAMD NPT tool parameters part 1](images/namd_npt_part1.png "NAMD NPT parameters 1")
> ![Snapshot of NAMD NPT tool parameters part 2](images/namd_npt_part2.png "NAMD NPT parameters 2")
{: .hands_on}


# Conclusion
{:.no_toc}

