---
layout: tutorial_hands_on
enable: false

title: Running molecular dynamics simulations using NAMD
zenodo_link: 'https://zenodo.org/record/3234841'
level: Intermediate
questions:
- How do I use the NAMD engine in Galaxy?
- What is the correct procedure for performing a simple molecular dynamics simulation of a protein?
objectives:
- Learn about the NAMD engine provided in BRIDGE
requirements:
  -
    type: "internal"
    topic_name: computational-chemistry
    tutorials:
      - setting-up-molecular-systems
time_estimation: 3H
key_points:
- Several MD engines are available in BRIDGE
- Workflows are available for common simulations tasks such as equilibration and production dynamics for various ensembles (NVT, NPT)
- You've run an equilibrium and production MD simulation using NAMD
follow_up_training:
  -
    type: "internal"
    topic_name: computational-chemistry
    tutorials:
      - analysis-md-simulations
contributors:
- chrisbarnettster
- tsenapathi
- simonbray

---


# Introduction



> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


In this tutorial we will perform a simulation with the popular [NAMD](http://www.ks.uiuc.edu/Research/namd/) molecular dynamics software. Please note NAMD tools are not currently available on a public Galaxy server due to licensing issues. If you are interested in following this tutorial, you will need to download the [BRIDGE docker container](https://github.com/scientificomputing/BRIDGE) and [download NAMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD) yourself.

This tutorial is made up of two parts. In the first section, we will look at preparation of a system (solvation, charge neutralisation, energy minimisation) using CHARMM. In the second section, we will perform an equilibration and production simulation, using NAMD. If you already completed the [Setting up molecular systems]({% link topics/computational-chemistry/tutorials/setting-up-molecular-systems/tutorial.md %}) tutorial, which covers the use of the CHARMM graphical user interface (GUI), you have already prepared your system, so go straight to the [second section](#md-simulations-with-namd), using the files you prepared earlier.

The process can be accomplished by selecting each tool from the tools menu, or by importing the workflow. The workflow method is most efficient and the individual tools used in the workflows are discussed below. The entire workflow (preparation + simulation) is shown below for [CHARMM and NAMD](#workflows).



> <details-title>NVT, NPT and statistical mechanics theory</details-title>
>
> [See Statistical Mechanics, McQuarrie for more in depth theory ISBN:9781891389153](https://books.google.co.za/books/about/Statistical_Mechanics.html?id=itcpPnDnJM0C&redir_esc=y)
>
{: .details}

# System preparation with CHARMM

If you already prepared your system using the CHARMM-GUI, and saved the output files, you can skip this section.

## Setup

Initially, we need to prepare a protein-ligand system in CHARMM. 

This tool will:
- solvate the protein-ligand complex, using the TIP3P water model
- neutralise the system, using 0.05M NaCl
- conduct a short energy minimisation 

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from the Zenodo link provided.
>    ```
>    https://zenodo.org/record/3234841/files/cbh1test.crd
>    https://zenodo.org/record/3234841/files/cbh1test.psf
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Check that the datatype is correct. The crd file should have the CRD datatype and the psf file the PSF datatype.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

> <hands-on-title>initial processing</hands-on-title>
>
> Run **System Setup** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"psf input"*: Protein structure file of the protein-ligand system. (Input dataset)
>    - {% icon param-file %} *"crd input"*: Coordinate file of the protein-ligand system. (Input dataset)
>    - *"Buffer"*: Edge to edge distance between the protein and the edge of the water box in angstroms
>    - *"Custom topology and parameter files for ligands"* (optional)
>
{: .hands_on}

### Energy minimisation

The setup provides us with the CRD and PSF files needed to perform a simulation. In addition, we have a PRM file which defines the parameters of the unit cell. Now, we need to perform energy minimisation. This relaxes the system, removing any steric clashes or unusual geometry which could artificially raise the energy.

This tools will:
- Minimise energy using a steepest descent algorithm followed by Adopted Newton Raphson (using the defined number of steps)
- Set up periodic boundaries and generate Particle Mesh Ewald (PME)
- generate reference structures for restraints in NAMD (if selected)

> <hands-on-title>energy minimization</hands-on-title>
>
> Run **Energy Minimizer** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"system_setup_crd"*: Coordinate file generated by the setup tool. (Input dataset)
>    - {% icon param-file %} *"system_setup_psf"*: Protein structure file generated by the setup tool. (Input dataset)
>    - {% icon param-file %} *"waterbox parameters input"*: Water box parameter file generated by the setup tool. (Input dataset)
>    - *"Minimization steps"*: `1000`
>    - *"Create reference structures for RMSD restraints for NAMD?"*: `No`
>    - *"Custom topology and parameter files for ligands"*: `No`
>
{: .hands_on}


# MD simulations with NAMD

At this point we are ready to run the simulation, which uses NAMD as a molecular dynamics engine. An NVT simulation is followed by an NPT simulation.


### NVT

Classical NVT dynamics, maintaining constant **n**umber of particles, **v**olume and **t**emperature.

This tool runs classical molecular dynamics simulations in NAMD using an NVT ensemble. User can run the simulation in small time intervals. The coordinates, velocities and the extended system files can be use to restart the simulations. If required, harmonic restraints can be used to maintain the protein shape. These restraints, in particular RMSD harmonic restraints can be added with the NAMD collective variable module.

> <hands-on-title>NVT dynamics</hands-on-title>
>
> 1. **NAMD MD Simulator (NVT)** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"xplor psf input"*: PSF file from CHARMM preparation (Input dataset)
>    - {% icon param-file %} *"pdb input"*: PDB file from CHARMM preparation (Input dataset)
>    - {% icon param-file %} *"PME grid specs"*: Generated by the setup tool (Input dataset)
>    - {% icon param-file %} *"waterbox prm input"*: Water box parameter file generated by the setup tool. (Input dataset)
>    - *"Temperature (K)"*: `300`
>    - *"Are you restarting a simulation?"*: `No`
>    - *"Use Harmonic restraints simulation?"*: `No`
>    - *"Custom parameter file for ligands"*: `No`
>    - *"DCD Frequency (ps)"*: `1` (Frequency to record frames in the DCD trajectory)
>    - *"Simulation Time (ps)"*: `10`
>    - *"Number of processors"*: `4`
{: .hands_on}


### NPT
Classical NPT dynamics, maintaining constant **n**umber of particles, **p**ressure and **t**emperature.

This tool runs classical molecular dynamics simulations in NAMD using an NPT ensemble. User can run the simulation in small time intervals. The coordinates, velocities and the extended system files can be use to restart the simulations. Harmonic restraints can be used. NAMD collective variable module is used to give RMSD harmonic restraints.

> <hands-on-title>NPT dynamics</hands-on-title>
>
> 1. **NAMD MD Simulator (NPT)** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"xplor psf input"*: PSF file from CHARMM preparation (Input dataset)
>    - {% icon param-file %} *"pdb input"*: PDB file from CHARMM preparation (Input dataset)
>    - {% icon param-file %} *"PME grid specs"*: Generated by the setup tool (Input dataset)
>    - {% icon param-file %} *"waterbox prm input"*: Water box parameter file generated by the setup tool. (Input dataset)
>    - *"Temperature (K)"*: `300`
>    - *"Pressure (bar)"*: `1.01325`
>    - *"Are you restarting a simulation?"*: `Yes`
>    - *"Coordinates from the previous simulation"*:  from the NVT simulation
>    - *"Velocities from the previous simulation"*:  from the NVT simulation
>    - *"Extended system of the previous simulation"*: from the NVT simulation
>    - *"Use Harmonic restraints simulation?"*: `No`
>    - *"Custom parameter file for ligands"*: `No`
>    - *"DCD Frequency (ps)"*: `1` (Frequency to record frames in the DCD trajectory)
>    - *"Simulation Time (ps)"*: `15`
>    - *"Number of processors"*: `4`
{: .hands_on}

# Workflows

Both the CHARMM preparatory workflow and the NAMD simulation workflow are available as an alternative to executing individual tools.

![Snapshot of CHARMM and NAMD analysis workflow]({% link topics/computational-chemistry/images/NAMD_CHARMMGUI_workflow.png %} "A simple simulation workflow starting with CHARMM for setup and NAMD to continue the production simulation")

![Snapshot of NAMD analysis workflow]({% link topics/computational-chemistry/images/NAMD_workflow.png %} "A simple NAMD simulation workflow")

# Conclusion

After completing the steps, or running the workflow, we have successfully produced a trajectory (the xtc file) which describes the atomic motion of the system. This can be viewed using molecular visualization software or analysed further; please visit the visualization and [analysis]({% link topics/computational-chemistry/tutorials/analysis-md-simulations/tutorial.md %}) tutorials for more information.


