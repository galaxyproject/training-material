---
layout: tutorial_hands_on

title: Running molecular dynamics simulations using GROMACS
zenodo_link: 'https://zenodo.org/record/2598415#.XJDyvoUo9hE'
questions:
- How do I use the GROMACS engine in Galaxy?
- What is the correct procedure for performing a simple molecular dynamics simulation of a protein?
objectives:
- Learn about the GROMACS engine provided via Galaxy.
- Understand the structure of the molecular dynamics workflow.
- Prepare a short (1 ns) trajectory for a simulation of a protein.
time_estimation: 2H
key_points:
- Molecular dynamics produces a trajectory describing the atomic motion of a system.
- Preparation of the system is required; setup, solvation, minimization, equilibration.
contributors:
- simonbray

---


# Introduction
{:.no_toc}

Molecular dynamics (MD) is a method to simulate molecular motion by iterative application of Newton's laws of motion. It is often applied to large biomolecules such as proteins or nucleic acids.

Multiple packages exist for performing MD simulations. One of the most popular is the open-source GROMACS, which is the subject of this tutorial. Other MD packages which are also provided in Galaxy are NAMD and CHARMM (link).

This is a introductory guide to using GROMACS in Galaxy to prepare and perform molecular dynamics on a small protein (lysozyme). It is based on the GROMACS tutorial provided by Justin Lemkul [here](http://www.mdtutorials.com/gmx/lysozyme/index.html).

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Simulation with GROMACS

A GROMACS workflow is provided for this tutorial; we will discuss the tools that make up each of the steps.

## **preliminaries**
To perform simulation, an initial PDB file is required. This should be 'cleaned' of solvent and any other non-protein atoms. Solvent will be re-added in a subsequent step.

A prepared file is available via Zenodo. Alternatively, you can prepare the file yourself. Download a PDB structure file from the [Protein Data Bank](https://www.rcsb.org/) and remove the unwanted atoms using [grep](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool).


> ### {% icon hands_on %} Hands-on (optional):
>
> 1. Go the the PDB website (LINK) and search for the code 1AKI. Download the structure and upload to Galaxy.
> 2. Use [grep](https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool) to remove all lines that refer to non-protein atoms. Select 'Don't Match', enter 'HETATM' under 'Regular Expression'.
>
{: .hands_on}

## **setup**
{: #setup}

Now we prepare a topology for a protein structure. The topology file contains all the information required to describe the molecule for the purposes of simulation. A force field and water model must be selected for topology calculation. Multiple choices are available for each.

In addition, a 'box' is created, centered on the structure, prior to the next step (solvation). Options include box dimensions and shape; here, rhombic dodecahedron is the best option, as it can contain the protein using in the smallest volume, thus reducing the simulation resources devoted to the solvent.

Finally, a 'position restraint file' is created which will be used for NVT/NPT equilibration.

This tool will:
- create a 'topology' file
- convert a PDB protein structure into a GRO file
- center the structure in a simulation box
- create a position restraint file

> ### {% icon hands_on %} Hands-on:
>
> **setup** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"PDB input file"*: `1AKI_clean.pdb` (Input dataset)
>    - *"Water model"*: `SPC/E`
>    - *"Force field"*: `OPLS/AA`
>    - *"Ignore hydrogens"*: `no`
>    - *"Box type"*: `Rectangular box with all sides equal`
>    - *"Generate detailed log"*: `yes`
>    - Press **Execute**
>
{: .hands_on}

## **solvate**
{: #solvate}

The next stage is protein solvation. 

At this stage ions are also automatically added to neutralize the charge of the system. In our case, as lysozyme has a charge of +8, 8 chloride anions are added.

This tool will:
- add water molecules to fill the box defined in the setup


> ### {% icon hands_on %} Hands-on:
>
> **solvate** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file produced by setup
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology produced by setup
>    - *"Water model for solvation"*: `SPC`
>    - *"Generate detailed log"*: `yes`
>    - Press **Execute** 
>
{: .hands_on}

## **EM**
{: #EM}

To obviate any steric clashes or unusual geometry which would artificially raise the energy of the system, we must relax the structure by running an energy minimization (EM) algorithm.

Here, and in the later steps, two options are presented under 'Parameter input'. Firstly, the default setting, which we will use for this tutorial, requires options to be selected through the Galaxy interface. Alternatively, you can choose to upload an MDP (molecular dynamics parameters) file to define the simulation parameters. Using your own MDP file will allow greater customization, as not all parameters are implemented in Galaxy (yet); however, it requires a more advanced knowledge of GROMACS. Description of all parameters can be found [here](http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html).


This tool will:
- Run an energy minimization algorithm on the system.

> ### {% icon hands_on %} Hands-on: 
>
> **EM** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - *"Generate detailed log"*: `yes`
>    - *"Parameter input"*: `Use default (partially customisable) setting`
>    - *"Choice of integrator"*: `Steepest descent algorithm` (most common choice for EM)
>    - *"Neighbor searching"*: `Generate a pair list with buffering` (the '[Verlet scheme](https://en.wikipedia.org/wiki/Verlet_list)')
>    - *"Electrostatics"*: `Fast smooth Particle-Mesh Ewald (SPME) electrostatics`
>    - *"Distance for the Coulomb cut-off"*: `1.0`
>    - *"Cut-off distance for the short-range neighbor list"*: `1.0` (but irrelevant as we are using the Verlet scheme)
>    - *"Short range van der Waals cutoff"*: `1.0`
>    - *"Number of steps for the MD simulation"*: `50000`
>    - *"EM tolerance"*: `1000.0`
>    - *"Maximum step size"*: `0.01`
>    - Press **Execute** 
>
{: .hands_on}

## **nvt**
{: #nvt}

At this point equilibration of the solvent around the solute (i.e. the protein) is necessary. This is performed in two stages: equilibration under an NVT ensemble, followed by an NPT ensemble. Use of the NVT ensemble entails maintaining constant number of particles, volume and temperature.

For equilibration, the protein must be held in place while the solvent is allowed to move freely around it. This is done using the position restraint file created in system setup.

This tool will:
 - Perform equilibration using classical NVT dynamics.

> ### {% icon hands_on %} Hands-on: NVT dynamics
>
> **nvt** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - {% icon param-file %} *"Position restraint file"*: Position restraint file produced by setup
>    - *"Parameter input"*: `Use default (partially customisable) setting`
>    - *"Choice of integrator"*: `A leap-frog algorithm for integrating Newton’s equations of motion` (A basic leap-frog integrator)
>    - *"Bond constraints"*: `Bonds with H-atoms` (bonds involving H are constrained)
>    - *"Neighbor searching"*: `Generate a pair list with buffering` (the '[Verlet scheme](https://en.wikipedia.org/wiki/Verlet_list)')
>    - *"Electrostatics"*: `Fast smooth Particle-Mesh Ewald (SPME) electrostatics`
>    - *"Temperature"*: `300`
>    - *"Step length in ps"*: `0.002`
>    - *"Number of steps that elapse between saving data points (velocities, forces, energies)"*: `5000`
>    - *"Distance for the Coulomb cut-off"*: `1.0`
>    - *"Cut-off distance for the short-range neighbor list"*: `1.0` (but irrelevant as we are using the Verlet scheme)
>    - *"Short range van der Waals cutoff"*: `1.0`
>    - *"Number of steps for the NVT simulation"*: `50000`
>    - *"Trajectory output"*: `Return no trajectory output` (we are not interested in how the system evolves to the equilibrated state, merely the final structure)
>    - *"Structure output"*: `Return the .gro structure`
>    - *"Generate detailed log"*: `yes`
>    - Press **Execute** 

> {: .hands_on}

## **npt**
{: #npt}
Having stabilized the temperature of the system with NVT equilibration, we also need to stabilize the pressure of the system. We therefore equilibrate again using the NPT (constant number of particles, pressure, temperature) ensemble.

> ### {% icon hands_on %} Hands-on: NPT dynamics
>
> **npt** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - {% icon param-file %} *"Checkpoint (TOP) file"*: Checkpoint file produced by NVT equilibration
>    - {% icon param-file %} *"Position restraint file"*: Position restraint file produced by setup
>    - *"Parameter input"*: `Use default (partially customisable) setting`
>    - *"Choice of integrator"*: `A leap-frog algorithm for integrating Newton’s equations of motion` (A basic leap-frog integrator)
>    - *"Bond constraints"*: `Bonds with H-atoms` (bonds involving H are constrained)
>    - *"Neighbor searching"*: `Generate a pair list with buffering` (the '[Verlet scheme](https://en.wikipedia.org/wiki/Verlet_list)')
>    - *"Electrostatics"*: `Fast smooth Particle-Mesh Ewald (SPME) electrostatics`
>    - *"Temperature"*: `300`
>    - *"Step length in ps"*: `0.002`
>    - *"Number of steps that elapse between saving data points (velocities, forces, energies)"*: `5000`
>    - *"Distance for the Coulomb cut-off"*: `1.0`
>    - *"Cut-off distance for the short-range neighbor list"*: `1.0` (but irrelevant as we are using the Verlet scheme)
>    - *"Short range van der Waals cutoff"*: `1.0`
>    - *"Number of steps for the NPT simulation"*: `50000`
>    - *"Trajectory output"*: `Return no trajectory output` (we are not interested in how the system evolves to the equilibrated state, merely the final structure)
>    - *"Structure output"*: `Return the .gro structure`
>    - *"Generate detailed log"*: `yes`
>    - Press **Execute** 
> {: .hands_on}

## **mdrun**
{: #mdrun}
Now that equilibration is complete, we can release the position restraints and perform a production MD simulation.

> ### {% icon hands_on %} Hands-on: Production simulation
>
> 1. **mdrun** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - {% icon param-file %} *"Checkpoint (TOP) file"*: Checkpoint file produced by NPT equilibration
>    - *"Parameter input"*: `Use default (partially customisable) setting`
>    - *"Choice of integrator"*: `A leap-frog algorithm for integrating Newton’s equations of motion` (A basic leap-frog integrator)
>    - *"Bond constraints"*: `Bonds with H-atoms` (bonds involving H are constrained)
>    - *"Neighbor searching"*: `Generate a pair list with buffering` (the '[Verlet scheme](https://en.wikipedia.org/wiki/Verlet_list)')
>    - *"Electrostatics"*: `Fast smooth Particle-Mesh Ewald (SPME) electrostatics`
>    - *"Temperature"*: `300`
>    - *"Step length in ps"*: `0.002`
>    - *"Number of steps that elapse between saving data points (velocities, forces, energies)"*: `5000`
>    - *"Distance for the Coulomb cut-off"*: `1.0`
>    - *"Cut-off distance for the short-range neighbor list"*: `1.0` (but irrelevant as we are using the Verlet scheme)
>    - *"Short range van der Waals cutoff"*: `1.0`
>    - *"Number of steps for the simulation"*: `500000`
>    - *"Trajectory output"*: `Return .xtc file (reduced precision)` (this time, we save the trajectory)
>    - *"Structure output"*: `Return the .gro structure`
>    - *"Generate detailed log"*: `yes`
>    - Press **Execute** 
> {: .hands_on}

# Conclusion

After completing the steps, or running the workflow, we have successfully produced a trajectory whcih describes the atomic motion of the system. This can be viewed using molecular visualization software or analysed further; please visit the visualization and analysis tutorials for more information.


{:.no_toc}

