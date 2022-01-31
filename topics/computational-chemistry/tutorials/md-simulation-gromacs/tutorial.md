---
layout: tutorial_hands_on

title: Running molecular dynamics simulations using GROMACS
zenodo_link: 'https://zenodo.org/record/2598415'
level: Intermediate
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
follow_up_training:
  -
    type: "internal"
    topic_name: computational-chemistry
    tutorials:
      - analysis-md-simulations
      - htmd-analysis
contributors:
  - simonbray

---


# Introduction
{:.no_toc}

Molecular dynamics (MD) is a method to simulate molecular motion by iterative application of Newton's laws of motion. It is often applied to large biomolecules such as proteins or nucleic acids.

Multiple packages exist for performing MD simulations. One of the most popular is the open-source GROMACS, which is the subject of this tutorial. Other MD packages which are also wrapped in Galaxy are [NAMD]({% link topics/computational-chemistry/tutorials/md-simulation-namd/tutorial.md %}) and CHARMM (available in the [docker container](https://github.com/scientificomputing/BRIDGE)).

This is a introductory guide to using GROMACS ({% cite abraham15 %}) in Galaxy to prepare and perform molecular dynamics on a small protein. For the tutorial, we will perform our simulations on hen egg white lysozyme.

> ### {% icon comment %} More information
> This guide is based on the GROMACS tutorial provided by Justin Lemkul [here](http://www.mdtutorials.com/gmx/lysozyme/index.html) - please consult it if you are interested in a more detailed, technical guide to GROMACS.
{: .comment}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Process

Prior to performing simulations, a number of preparatory steps need to be executed.

The process can be divided into multiple stages:
 1. Setup (loading data, solvation i.e. addition of water and ions)
 2. Energy minimization of the protein
 3. Equilibration of the solvent around the protein (with two ensembles, NVT and NPT)
 4. Production simulation, which produces our trajectory.

The trajectory is a binary file that records the atomic coordinates at multiple time steps, and therefore shows the dynamic motion of the molecule. Using visualization software, we can display this trajectory as a film displaying the molecular motion of the protein. We will discuss each step making up this workflow in more detail.


# Getting data
To perform simulation, an initial PDB file is required. This should be 'cleaned' of solvent and any other non-protein atoms. Solvent will be re-added in a subsequent step.

A prepared file is available via Zenodo. Alternatively, you can prepare the file yourself. Download a PDB structure file from the [Protein Data Bank](https://www.rcsb.org/) and remove the unwanted atoms using the grep text processing tool. This simply removes the lines in the PDB file that refer to the unwanted atoms.


> ### {% icon hands_on %} Hands-on: Upload an initial structure
> 1. First of all, create a new history and give it a name.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Use the {% tool [Get PDB](toolshed.g2.bx.psu.edu/repos/bgruening/get_pdb/get_pdb/0.1.0) %} tool to download a PDB file for simulation:
>    - *"PDB accession code"*: `1AKI`
>
> 3. Use the {% tool [Search in textfiles (grep)](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} text processing tool to remove all lines that refer to non-protein atoms.
>    - *"Select lines from"*: uploaded PDB file
>    - *"that"*: `Don't Match`
>    - *"Regular Expression"*: `HETATM`
>
{: .hands_on}

> ### {% icon comment %} Alternative upload
> As an alternative option, if you prefer to upload the cleaned file directly from Zenodo, you can do so with the following link:
> ```
> https://zenodo.org/record/2598415/files/1AKI_clean.pdb
> ```
>
> {% snippet faqs/galaxy/datasets_import_via_link.md %}
{: .comment}

> ### {% icon details %} Background: What is the PDB (Protein Data Bank) and format?
>
> The Protein Data Bank (PDB) format contains atomic coordinates of biomolecules and provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies. Each structure is stored under a four-letter accession code. For example, the PDB file we will use is assigned the code [1AKI](https://www.rcsb.org/pdb/explore/explore.do?structureId=1AKI).
>
> More resources:
>
>  -  Multiple structures are stored and can be queried at [https://www.rcsb.org/](https://www.rcsb.org/)
>  - Documentation describing the PDB file format is available from the wwPDB at [http://www.wwpdb.org/documentation/file-format.php](http://www.wwpdb.org/documentation/file-format.php).
{: .details}

## Lysozyme
The protein we will look at in this tutorial is hen egg white [lysozyme](https://en.wikipedia.org/wiki/Lysozyme), a widely studied enzyme which is capable of breaking down the polysaccharides of many bacterial cell walls. It is a small (129 residues), highly stable globular protein, which makes it ideal for our purposes.

![Structure of lysozyme openly available from https://commons.wikimedia.org/wiki/File:Lysozyme.png]({% link topics/computational-chemistry/images/Lysozyme.png %} "Structure of lysozyme")

# Setup

The **GROMACS initial setup** {% icon tool %} tool uses the PDB input to create three files which will be required for MD simulation.

Firstly, a topology for the protein structure is prepared. The topology file contains all the information required to describe the molecule for the purposes of simulation - atom masses, bond lengths and angles, charges. Note that this automatic construction of a topology is only possible if the building blocks of the molecules (i.e. amino acids in the case of a protein) are precalculated for the given force field. A force field and water model must be selected for topology calculation. Multiple choices are available for each; we will use the OPLS/AA force field and SPC/E water model. Secondly, a GRO structure file is created, storing the structure of the protein. Finally, a 'position restraint file' is created which will be used for NVT/NPT equilibration. We will return to this later.

In summary, the initial setup tool will:
- create a 'topology' file
- convert a PDB protein structure into a GRO file, with the structure centered in a simulation box (unit cell)
- create a position restraint file

After these files have been generated, a further step is required to define a simulation box (unit cell) in which the simulation can take place. This can be done with the **GROMACS structure configuration** {% icon tool %} tool. It also defines the unit cell 'box', centered on the structure. Options include box dimensions and shape; here, while a cuboidal box may be most intuitive, rhombic dodecahedron is the most efficient option, as it can contain the protein using the smallest volume, thus reducing the simulation resources devoted to the solvent.

> ### {% icon hands_on %} Hands-on: perform initial processing
>
> 1. Run {% tool [GROMACS initial setup](toolshed.g2.bx.psu.edu/repos/chemteam/gmx_setup/gmx_setup/2020.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"PDB input file"*: `1AKI_clean.pdb` (Input dataset)
>    - *"Water model"*: `SPC/E`
>    - *"Force field"*: `OPLS/AA`
>    - *"Ignore hydrogens"*: `No`
>    - *"Generate detailed log"*: `Yes`
>
> > ### {% icon question %} Question
> >
> > Why is it necessary to provide an input structure containing no non-protein molecules?
> >
> > > ### {% icon solution %} Solution
> > > Automatic topology construction only succeeds if the components of the structure are recognized. For example, providing a structure of a protein in complex with a non-protein ligand or cofactor will result in an error.
> > {: .solution}
> {: .question}
>
> 2. Run {% tool [GROMACS structure configuration](toolshed.g2.bx.psu.edu/repos/chemteam/gmx_editconf/gmx_editconf/2020.4+galaxy0) %} with the following parameters:
>    - *"Input structure"*: GRO output from initial setup tool
>    - *"Output format"*: `GRO file`
>    - *"Configure box?"*: `Yes`
>    - *"Box dimensions in nanometers"*: `1.0`
>    - *"Box type"*: `Rectangular box with all sides equal`
>    - *"Generate detailed log"*: `Yes`
{: .hands_on}


# Solvation

The next stage is protein solvation, performed using **GROMACS solvation and adding ions** {% icon tool %}. Water molecules are added to the structure and topology files to fill the unit cell. At this stage sodium or chloride ions are also automatically added to neutralize the charge of the system. In our case, as lysozyme has a charge of +8, 8 chloride anions are added.

![Solvated protein]({% link topics/computational-chemistry/images/solvated_protein.png %} "Solvated protein in a cubic unit cell")

> ### {% icon hands_on %} Hands-on: solvation
>
> {% tool [GROMACS solvation and adding ions](toolshed.g2.bx.psu.edu/repos/chemteam/gmx_solvate/gmx_solvate/2020.4+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file produced by the structure configuration tool
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology produced by setup
>    - *"Water model for solvation"*: `SPC`
>    - *"Add ions to neutralise system?"*: `Yes, add ions`
>    - *"Specify salt concentration (sodium chloride) to add, in mol/liter"*: `0`
>    - *"Generate detailed log"*: `Yes`
>
{: .hands_on}

# Energy minimization

To remove any steric clashes or unusual geometry which would artificially raise the energy of the system, we must relax the structure by running an energy minimization (EM) algorithm.

Here, and in the later steps, two options are presented under 'Parameter input'. Firstly, the default setting, which we will use for this tutorial, requires options to be selected through the Galaxy interface. Alternatively, you can choose to upload an MDP (molecular dynamics parameters) file to define the simulation parameters. Using your own MDP file will allow greater customization, as not all parameters are implemented in Galaxy (yet); however, it requires a more advanced knowledge of GROMACS. Description of all parameters can be found [here](http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html).

> ### {% icon hands_on %} Hands-on: energy minimization
>
> {% tool [GROMACS energy minimization](toolshed.g2.bx.psu.edu/repos/chemteam/gmx_em/gmx_em/2020.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file produced by solvation tool
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - *"Parameter input"*: `Use default (partially customisable) setting`
>    - *"Choice of integrator"*: `Steepest descent algorithm` (most common choice for EM)
>    - *"Neighbor searching"*: `Generate a pair list with buffering` (the '[Verlet scheme](https://en.wikipedia.org/wiki/Verlet_list)')
>    - *"Electrostatics"*: `Fast smooth Particle-Mesh Ewald (SPME) electrostatics`
>    - *"Distance for the Coulomb cut-off"*: `1.0`
>    - *"Cut-off distance for the short-range neighbor list"*: `1.0` (but irrelevant as we are using the Verlet scheme)
>    - *"Short range van der Waals cutoff"*: `1.0`
>    - *"Number of steps for the MD simulation"*: `50000`
>    - *"EM tolerance"*: `1000`
>    - *"Maximum step size"*: `0.01`
>    - *"Generate detailed log"*: `Yes`
>
{: .hands_on}

# Equilibration

At this point equilibration of the solvent around the solute (i.e. the protein) is necessary. This is performed in two stages: equilibration under an NVT ensemble, followed by an NPT ensemble. Use of the NVT ensemble entails maintaining constant **n**umber of particles, **v**olume and **t**emperature, while the NPT ensemble maintains constant **n**umber of particles, **p**ressure and **t**emperature. (The NVT ensemble is also known as the isothermal-isochoric ensemble, while the NPT ensemble is also known as the isothermal-isobaric ensemble).

During the first equilibration step (NVT), the protein must be held in place while the solvent is allowed to move freely around it. This is achieved using the position restraint file we created in system setup. When we specify this restraint, protein movement is not totally forbidden, but is energetically punished. During the second NPT step, we remove the restraints.

## NVT equilibration
Firstly, we perform equilibration using classical NVT dynamics.

> ### {% icon hands_on %} Hands-on: NVT dynamics
>
> {% tool [GROMACS simulation](toolshed.g2.bx.psu.edu/repos/chemteam/gmx_sim/gmx_sim/2020.4+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - *"Use a checkpoint (CPT) file"*: `No CPT input`
>    - *"Produce a checkpoint (CPT) file"*: `Produce CPT output`
>    - *"Apply position restraints"*: `Apply position restraints`
>    - {% icon param-file %} *"Position restraint file"*: Position restraint file produced by 'Setup' tool.
>    - *"Ensemble"*: `Isothermal-isochoric ensemble (NVT).`
>    - *"Trajectory output"*: `Return no trajectory output` (we are not interested in how the system evolves to the equilibrated state, merely the final structure)
>    - *"Structure output"*: `Return .gro file`
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
>    - *"Generate detailed log"*: `Yes`
{: .hands_on}

## NPT equilibration
Having stabilized the temperature of the system with NVT equilibration, we also need to stabilize the pressure of the system. We therefore equilibrate again using the NPT (constant number of particles, pressure, temperature) ensemble.

Note that we can continue where the last simulation left off (with new parameters) by using the checkpoint (CPT) file saved at the end of the NVT simulation.

> ### {% icon hands_on %} Hands-on: NPT dynamics
>
> {% tool [GROMACS simulation](toolshed.g2.bx.psu.edu/repos/chemteam/gmx_sim/gmx_sim/2020.4+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - *"Use a checkpoint (CPT) file"*: `Continue simulation from a CPT file.`
>    - {% icon param-file %} *"Checkpoint (CPT) file"*: Checkpoint file produced by NVT equilibration
>    - *"Produce a checkpoint (CPT) file"*: `Produce CPT output`
>    - *"Apply position restraints"*: `No position restraints`
>    - {% icon param-file %} *"Position restraint file"*: None
>    - *"Ensemble"*: `Isothermal-isobaric ensemble (NPT).`
>    - *"Trajectory output"*: `Return no trajectory output`
>    - *"Structure output"*: `Return .gro file`
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
>    - *"Generate detailed log"*: `Yes`
{: .hands_on}

> ### {% icon question %} Question
>
> Why is the position of the protein restrained during equilibration?
>
> > ### {% icon solution %} Solution
> > The purpose of equilibration is to stabilize the temperature and pressure of the system; these are overwhelmingly dependent on the solvent. Structural changes in the protein are an additional complicating variable, which can more simply be removed by restraining the protein.
> {: .solution}
{: .question}


# Production simulation
Now that equilibration is complete, we can release the position restraints. We are now finally ready to perform a production MD simulation.

> ### {% icon hands_on %} Hands-on: Production simulation
>
> 1. {% tool [GROMACS simulation](toolshed.g2.bx.psu.edu/repos/chemteam/gmx_sim/gmx_sim/2020.4+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO structure file
>    - {% icon param-file %} *"Topology (TOP) file"*: Topology
>    - *"Use a checkpoint (CPT) file"*: `Continue simulation from a CPT file.`
>    - {% icon param-file %} *"Checkpoint (CPT) file"*: Checkpoint file produced by NPT equilibration
>    - *"Produce a checkpoint (CPT) file"*: `No CPT output`
>    - *"Apply position restraints"*: `No position restraints`
>    - *"Ensemble"*: `Isothermal-isobaric ensemble (NPT).`
>    - *"Trajectory output"*: `Return .xtc file (reduced precision)` (this time, we save the trajectory)
>    - *"Structure output"*: `Return .gro file`
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
>    - *"Generate detailed log"*: `Yes`
{: .hands_on}

# Workflow

A GROMACS workflow is provided for this tutorial. Overall, the workflow takes a PDB (Protein Data Bank) structure file as input and returns a MD trajectory.

![GROMACS workflow]({% link topics/computational-chemistry/images/workflow_gromacs.png %} "The basic GROMACS workflow")

# Conclusion
{:.no_toc}

After completing the steps, or running the workflow, we have successfully produced a trajectory (the xtc file) which describes the atomic motion of the system. This can be viewed using molecular visualization software or analysed further; please visit the visualization and [analysis]({% link topics/computational-chemistry/tutorials/analysis-md-simulations/tutorial.md %}) tutorials for more information.

![Trajectory]({% link topics/computational-chemistry/images/traj.gif %} "Trajectory produced using the GROMACS workflow, visualized with the NGL viewer")




