---
layout: tutorial_hands_on

title: High Throughput Molecular Dynamics and Analysis
level: Advanced
zenodo_link: ''
questions:
- How are protein-ligand systems parameterized for molecular dynamics simulation?
- What kind of analysis can be carried out on molecular trajectories?
- How can high-throughput MD be used to study multiple ligands?
objectives:
- Learn about force-fields and MD parameterization
- Learn how to conduct MD simulation and analysis for a protein-ligand system
- Understand how different molecular interactions contribute to the binding affinity of various ligands for the Hsp90 protein.
time_estimation: 3H
key_points:
- Simulating protein-ligand systems is more complex than simply simulating protein-only systems
- ....
contributors:
- simonbray
- tsenapathi
- chrisbarnettster
- bgruening

---


# Introduction
{:.no_toc}

This tutorial provides an introduction to using high-throughput molecular dynamics to study protein-ligand interaction, as applied to N-terminus of Hsp90 (heat shock protein 90).


<!-- This is a comment. -->

<!-- General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial. -->


<!-- **Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)** -->

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background

## What is high-throughput molecular dynamics?
Molecular dynamics (MD) is a method to simulate molecular motion by iterative application of Newton’s laws of motion. It is often applied to large biomolecules such as proteins or nucleic acids. A common application is to assess the interaction between these macromolecules and a number of small molecules (e.g. potential drug candidates). This tutorial provides a guide to setting up and running a high-throughput workflow for screening multiple small molecules, using the open-source GROMACS tools provided through the Galaxy platform.


## Why is Hsp90 interesting to study?
The 90 kDa heat shock protein (Hsp90) is a chaperone protein responsible for catalyzing the conversion of a wide variety of proteins to a functional form; examples of the Hsp90 clientele, which totals several hundred proteins, include nuclear steroid hormone receptors and protein kinases. The mechanism by which Hsp90 acts varies between clients, as does the client binding site; the process is dependent on post-translational modifications of Hsp90 and the identity of co-chaperones which bind and regulate the conformational cycle.

Due to its vital biochemical role as a chaperone protein involved in facilitating the folding of many client proteins, Hsp90 is an attractive pharmaceutical target. In particular, as protein folding is a potential bottleneck to slow cellular reproduction and growth, blocking Hsp90 function using inhibitors which bind tightly to the ATP binding site could assist in treating cancer; for example, the antibiotic geldanamycin and its analogs are under investigation as possible anti-tumor agents.

![Hsp90 structure, with a ligand bound]({% link topics/computational-chemistry/images/hsp90lig.png %} "Structure of Hsp90, with a ligand bound. <a href="https://usegalaxy.eu/u/sbray/v/hsp90-lig">Click to view</a> in NGL. ({% cite ngl %})")


## Get data

First of all, download the required data.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Search Galaxy for the 'Get PDB' tool. Request the accession code ```6hhr```.
> 3. Rename the dataset to 'Hsp90 structure'
> 4. Check that the datatype is correct (PDB file).
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
{: .hands_on}


# Simulation

## Topology generation

Now we have downloaded a PDB structure of the protein we wish to study, we will start parameterizing it for MD simulation.

GROMACS distinguishes between constant and dynamic attributes of the atoms in the system. The constant attributes (e.g. atom charges, bonds connecting atoms) are listed in the topology (TOP file), while dynamic attributes (attributes that can change during a simulation, e.g. atom position, velocities and forces) are stored in structure (PDB or GRO) and trajectory (XTC and TRR) files.

The PDB file we start from only explicitly states atom type and position. Therefore, before beginning simulation, we need to calculate the rest of the information contained within the topology file. There are a range of force fields which perform these calculations in slightly different ways.

Parameterization needs to be done separately for the ligand and protein. Therefore, the first step is to separate the PDB file into two sets of coordinates - one for the ligand and one for the protein.

> ### {% icon question %} Question
>
> 1. Why do protein and ligand need to be parameterized separately?
>
> > ### {% icon solution %} Solution
> >
> > 1. Protein and small molecules are constructed differently. A protein is made up of 20 different building blocks (amino acids) - therefore, to construct a protein topology, amino acid topologies simply need to be combined appropriately. By contrast, the structure of small molecules is far more flexible and needs to be calculated for each different structure.
> >
> {: .solution}
>
{: .question}

### Extract protein and ligand coordinates

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Search in textfiles** {% icon tool %} with the following parameters:
>    - *"Select lines from"*: 'Hsp90 structure'   
>    - *"that"*: `Don't Match`
>    - *"Regular Expression"*: `HETATM`
> 2. Rename output to 'Protein (PDB)'
> 3. **Search in textfiles** {% icon tool %} with the following parameters:
>    - *"Select lines from"*: 'Hsp90 structure'   
>    - *"that"*: `Match`
>    - *"Regular Expression"*: `AG5E`
> 4. Rename output to 'Ligand (PDB)'
>
{: .hands_on}


### Set up protein topology

Firstly, we need to calculate the topology for the protein file. We will use the **GROMACS initial setup** {% icon tool %} tool.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GROMACS initial setup** {% icon tool %} with the following parameters:
>    - *"PDB input file"*: 'Protein (PDB)' file
>    - *"Force field"*: `gaff`
>    - *"Water model"*: `TIP3P`
>    - *"Generate detailed log"*: `Yes`
>
>    > ### {% icon comment %} Comment
>    > A force field is essentially a function to calculate the potential energy of a system, based on various empirical parameters (for the atoms, bonds, charges, dihedral angles and so on). There are a number of families of forcefields; some of the most commonly used include CHARMM, AMBER, GROMOS and OPLS. Here, we use GAFF (general AMBER force field), which is a generalized AMBER force field which can be applied to almost any small organic molecule, not just macromolecules such as proteins.
>    >
>    >
>    > A wide range of models exist for modeling water. Here we are using the common TIP3P model, which is an example of a 'three-site model' - so-called because the molecule is modeled using three points, corresponding to the three atoms of water. (Four- and five-site models include additional 'dummy atoms' representing the negative charges of the lone pairs of the oxygen atom).
>    {: .comment}
>
{: .hands_on}

The tool produces four outputs: a GRO file (containing the coordinates of the protein), a TOP file (containing other information, including on charges, masses, bonds and angles), an ITP file (which will be used to restrain the protein position in the equilibration step later on), and a log for the tool.

Please note all GROMACS tools output a log. Generally, you only need to look at this when a job fails. It provides useful information for debugging if we encounter any problems.


### Generate a topology for the ligand

To generate a topology for the ligand, we will use the **acpype** {% icon tool %} tool. This provides a convenient interface to the AmberTools suite and allows us to easily create the ligand topology in the format required by GROMACS.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Generate MD topologies for small molecules** {% icon tool %} with the following parameters:
>    - *"Input file"*: 'Ligand (PDB)'
>    - *"Charge of the molecule"*: `0`
>    - *"Multiplicity"*: `1`
>    - *"Force field to use for parameterization"*: `AMBER14SB`
>    - *"Save GRO file?"*: `Yes`
>
{: .hands_on}



## Solvation and energy minimization

Having generated topologies, we now need to combine them, define the box which contains the system, add solvent and ions, and perform an energy minimization step.

### Combine topology and GRO files

> ### {% icon hands_on %} Hands-on: Combine GRO files
>
> 1. On the `Structure file (GRO format)` created by the **acpype** tool,click on the `Visualize this data` icon. Select `Editor` to open the file using the text editor integrated into Galaxy. Select all the lines starting with `1 GSE` and copy your selection.
> 2. Open the Protein GRO file by clicking on the `Visualize this data` button on the dataset.
> 3. Paste the lines from the ligand GRO file just before the last line.
> 4. If you scroll back to the top, you will see that the total number of atoms in the system is given in the second line (`3280`). You have just added 21 new atoms, so increase the value by 21 to `3301`.
> 5. Click `Export` to save your changes as a new dataset. Make sure the datatype of the new file is still `GRO`. Rename to `System GRO file`.
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Combine topology files
>
> 1. On the ligand `Topology` created by the **acpype** tool, right-click on the `Visualize this data` icon and open the link in a new tab. Select the first section in the file, starting with `[ atomtypes ]`, and copy the selection.
> 2. Returning to the first tab, open the protein TOP file using the text editor integrated into Galaxy by clicking on the `Visualize this data` button on the dataset.
> 3. Paste the lines from the ligand ITP file near to the top of the file, just after the line `#include "amber99sb.ff/forcefield.itp"`.
> 4. Go back to the ligand ITP file and select the rest of the file (from `[ moleculetypes ]`) onwards. Copy the selection.
> 5. In the protein TOP file, paste the selection near to the bottom of the file, before the line `; Include water topology` (and just after the position restraint file). Notice that the `[ moleculetype ]` section you just copied starts with `base` - this is the name acpype has given to the ligand. Feel free to change this to whatever you prefer - `ligand`, or `GSE`.
> 6. Finally, we need to state in the topology that we have included a new kind of molecule. Go to the final section (`[ molecules ]`) and add a new line `base` (or whatever name you gave the ligand in step 5), with a 1 in the `#mols` column.
> 5. Click `Export` to save your changes as a new dataset. Make sure the datatype of the new file is still `TOP`. Rename to `System topology`.
{: .hands_on}

If this procedure was too complicated, you can download the combined files here: LINK. However, you will find it useful to understand the information contained within topology files and learn how to make changes to it.

### Create the simulation box with **GROMACS structure configuration**

The next step, once combined coordinate (GRO) and topology (TOP) files have been created, is to create a simulation box in which the system is situated.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GROMACS structure configuration** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input structure"*: `System GRO file` (Input dataset)
>    - *"Configure box?"*: `Yes`
>        - *"Box dimensions in nanometers"*: `1.0`
>        - *"Box type"*: `Triclinic`
>    - *"Generate detailed log"*: `Yes`
>
>
>    > ### {% icon comment %} Comment
>    >
>    > This tool simply adds a box of the required dimensions to the GRO file. A distance of at least 1.0 nm is recommended to avoid interactions between the protein and its mirror image. On the other hand, increasing the box size too far will increase the simulation time, due to the greater number of solvent molecules which need to be treated.
>    {: .comment}
>
{: .hands_on}

### Solvation

The next step is solvation of the newly created simulation box. Note that the system is charged (depending on the pH) - the solvation tool also adds sodium or chloride ions as required to neutralise this.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GROMACS solvation and adding ions** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: `output` (output of **GROMACS structure configuration** {% icon tool %})
>    - {% icon param-file %} *"System topology"*: `output`
>    - *"Generate detailed log"*: `Yes`
>
>
{: .hands_on}


### Energy minimization

The next step is energy minimization, which can be carried out using the **GROMACS energy minimization** {% icon tool %} tool.

> ### {% icon question %} Question
>
> 1. What is the purpose of energy minimization?
>
> > ### {% icon solution %} Solution
> >
> > 1. Running an energy minimization (EM) algorithm relaxes the structure by removing any steric clashes or unusual geometry which would artificially raise the energy of the system.
> >
> {: .solution}
>
{: .question}


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GROMACS energy minimization** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file."*: `output1` (output of **GROMACS solvation and adding ions** {% icon tool %})
>    - {% icon param-file %} *"Topology (TOP) file."*: `output2` (output of **GROMACS solvation and adding ions** {% icon tool %})
>    - *"Parameter input"*: `Use default (partially customisable) setting`
>        - *"Number of steps for the MD simulation"*: `50000`
>        - *"EM tolerance"*: `1000.0`
>    - *"Generate detailed log"*: `Yes`
>    - Rename output to `Minimized GRO file`
>
{: .hands_on}


## Equilibration

We now carry out equilibration in two stages: NVT and NPT. This is discussed at greater length in the basic GROMACS tutorial. Equilibration requires restraining the protein structure - we use the ITP file produced by the initial setup tool for this.

> ### {% icon comment %} More detail about equilibration
>
> At this point equilibration of the solvent around the solute (i.e. the protein) is necessary. This is performed in two stages: equilibration under an NVT ensemble, followed by an NPT ensemble. Use of the NVT ensemble entails maintaining constant **n**umber of particles, **v**olume and **t**emperature, while the NPT ensemble maintains constant **n**umber of particles, **p**ressure and **t**emperature. (The NVT ensemble is also known as the isothermal-isochoric ensemble, while the NPT ensemble is also known as the isothermal-isobaric ensemble).
>
> For equilibration, the protein must be held in place while the solvent is allowed to move freely around it. This is achieved using the position restraint file we created in system setup. When we specify this restraint, protein movement is not totally forbidden, but is energetically punished.
>
{: .comment}


> ### {% icon hands_on %} Hands-on: NVT equilibration
>
> 1. **GROMACS simulation** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: `Minimized GRO file` (from energy minimization step)
>    - {% icon param-file %} *"Topology (TOP) file"*: TOP file produced by solvation step.
>    - In *"Inputs"*:
>        - {% icon param-file %} *"Position restraint (ITP) file"*: ITP file produced by initial setup step.
>    - In *"Outputs"*:
>        - *"Trajectory output"*: `Return .xtc file (reduced precision)`
>        - *"Structure output"*: `Return .gro file`
>        - *"Produce a checkpoint (CPT) file"*: `Produce CPT output`
>    - In *"Settings"*:
>        - *"Parameter input"*: `Use default (partially customisable) setting`
>            - *"Bond constraints (constraints)"*: `All bonds (all-bonds).`
>            - *"Temperature /K"*: `300`
>            - *"Step length in ps"*: `0.0002`
>            - *"Number of steps that elapse between saving data points (velocities, forces, energies)"*: `1000`
>            - *"Number of steps for the simulation"*: `50000`
>    - *"Generate detailed log"*: `Yes`
>
>
{: .hands_on}

Having stabilized the temperature of the system with NVT equilibration, we also need to stabilize the pressure of the system. We therefore equilibrate again using the NPT (constant number of particles, pressure, temperature) ensemble.

Note that we can continue where the last simulation left off (with new parameters) by using the checkpoint (CPT) file saved at the end of the NVT simulation.

> ### {% icon hands_on %} Hands-on: NPT equilibration
>
> 1. **GROMACS simulation** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: GRO output of **GROMACS simulation** {% icon tool %} (NVT equilibration)
>    - {% icon param-file %} *"Topology (TOP) file"*: TOP file produced by solvation step.
>    - In *"Inputs"*:
>        - {% icon param-file %} *"Checkpoint (CPT) file"*: Output of **GROMACS simulation** {% icon tool %} (NVT equilibration))
>        - {% icon param-file %} *"Position restraint (ITP) file"*: ITP file produced by initial setup step.
>    - In *"Outputs"*:
>        - *"Trajectory output"*: `Return .xtc file (reduced precision)`
>        - *"Structure output"*: `Return .gro file`
>        - *"Produce a checkpoint (CPT) file"*: `Produce CPT output`
>    - In *"Settings"*:
>        - *"Ensemble"*: `Isothermal-isobaric ensemble (NPT)`
>        - *"Parameter input"*: `Use default (partially customisable) setting`
>            - *"Bond constraints (constraints)"*: `All bonds (all-bonds).`
>            - *"Temperature /K"*: `300`
>            - *"Step length in ps"*: `0.002`
>            - *"Number of steps that elapse between saving data points (velocities, forces, energies)"*: `1000`
>            - *"Number of steps for the simulation"*: `50000`
>    - *"Generate detailed log"*: `Yes`
>
>
{: .hands_on}


## Main simulation

We can now remove the restraints and continue with the simulation. The simulation will run for 1 million steps, with a step size of 1 fs, so wil have a total length of 1 ns.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GROMACS simulation** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"GRO structure file"*: Output of **GROMACS simulation** {% icon tool %} (NPT equilibration)
>    - {% icon param-file %} *"Topology (TOP) file"*: Output of the solvation step
>    - In *"Inputs"*:
>        - {% icon param-file %} *"Checkpoint (CPT) file"*: Output of **GROMACS simulation** {% icon tool %} (NPT simulation))
>    - In *"Outputs"*:
>        - *"Trajectory output"*: `Return .xtc file (reduced precision)`
>        - *"Structure output"*: `Return .gro file`
>        - *"Produce a checkpoint (CPT) file"*: `Produce CPT output`
>    - In *"Settings"*:
>        - *"Ensemble"*: `Isothermal-isobaric ensemble (NPT)`
>        - *"Parameter input"*: `Use default (partially customisable) setting`
>            - *"Temperature /K"*: `300`
>            - *"Step length in ps"*: `0.001`
>            - *"Number of steps that elapse between saving data points (velocities, forces, energies)"*: `1000`
>            - *"Number of steps for the simulation"*: `1000000`
>    - *"Generate detailed log"*: `Yes`
>
>
{: .hands_on}


# Analysis

An analysis of the GROMACS simulation outputs (structure and trajectory file) will be carried out using Galaxy tools developed for computational chemistry~\cite{senapathi_biomolecular_2019} based on popular analysis software, such as MDAnalysis~\cite{michaudagrawal_mdanalysis_2011}, MDTraj~\cite{mcgibbon_mdtraj_2015}, and  Bio3D~\cite{skjaerven_integrating_2014}. These tools output both tabular files as well as a variety of attractive plots.

![Analysis workflow for Hsp90 with ligand](../../images/workflow_htmd_analysis.png "Analysis workflow for Hsp90 with a ligand")

## Create PDB file needed by most analysis tools

Before beginning a detailed analysis, the structure and trajectory files generated previously need to be converted into different formats. First, convert the structural coordinates of the system in GRO format into PDB format. This PDB file will be used by most analysis tools as a starting structure.  Next, convert the trajectory from XTC to DCD format, as a number of tools (particularly those based on Bio3D) only accept trajectories in DCD format.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **GROMACS structure configuration** {% icon tool %} with the following parameters:
>    - *"Output format"*: `PDB file`
>    - *"Configure box?"*: `No`
>
>    > ### {% icon comment %} Comment
>    >
>    > This tool can also be used to carry out initial setup (as discussed in the simulation methods section) for GROMACS simulations and convert from PDB to GRO format.
>    {: .comment}
>
{: .hands_on}


## Convert trajectory to DCD format

Convert from XTC to DCD format. A number of the analysis tools being used have been built to analyse trajectories in CHARMM's DCD format.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MDTraj file converter** {% icon tool %} with the following parameters:
>    - *"Output format"*: `DCD file`
>    >
>    > ### {% icon comment %} Comment
>    >
>    > This tool can also be used to interconvert between several trajectory formats.
>    {: .comment}
>
{: .hands_on}


## RMSD analysis - protein

The Root Mean Square Deviation (RMSD) and Root Mean Square Fluctuation (RMSF) are calculated to check the stability and conformation of the protein and ligand through the course of the simulation.
RMSD is a standard measure of structural distance between coordinate
sets that measures the average distance between a group of atoms. The
RMSD of the C$\alpha$ atoms of the protein backbone is calculated here and
is a measure of how much the protein conformation has changed between different time points in the trajectory.


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RMSD Analysis** {% icon tool %} with the following parameters:
>    - *"Select domains"*: `C-alpha`
>
>    > ### {% icon comment %} Comment
>    >
>    > Note that for more complex systems, you may need to consider a more focused selection. For example, if you have a ligand that is a protein consider modifying this selection.
>    {: .comment}
>
{: .hands_on}

![RMSD timeseries Hsp90](../../images/htmd_analysis_rmsd1_series.png "RMSD timeseries for the Hsp90 Cα atoms")

![RMSD histogram Hsp90](../../images/htmd_analysis_rmsd1_histo.png "RMSD histogram for the Hsp90 Cα atoms")

The RMSD time series for the protein shows a thermally stable and equilibrated structure that plateaus at 1.0{\AA} with an average RMSD between 0.8{\AA} and 1.0{\AA}. There are no large conformational changes during the simulation. The RMSD histogram confirms this, see Figure \ref{fig:rmsdprotein}. Note these graphs are automatically created by Galaxy as part of the tool's outputs.


## RMSD analysis - ligand

Calculating the RMSD of the ligand is necessary to check if it is stable in the active site and to identify possible binding modes. If the ligand is not stable, there will be large fluctuations in the RMSD.

For the RMSD analysis of the ligand, the `Select domains' parameter of the tool can for convenience be set to `Ligand'; however, this automatic selection sometimes fails. The other alternative, which we apply here, is to specify the `Residue ID' in the textbox provided. In this example the ligand's Residue ID is `G5E'. The output is the requested RMSD data as a time series, the RMSD plotted as a time series and as a histogram (for example, see Figure \ref{fig:rmsdprotein} in the results section).


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RMSD Analysis** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"DCD trajectory input"*: `output` (output of **MDTraj file converter** {% icon tool %})
>    - {% icon param-file %} *"PDB input"*: `output` (output of **GROMACS structure configuration** {% icon tool %})
>    - *"Select domains"*: `Residue ID`
>        - *"Residue ID"*: `G5E`
>
>
{: .hands_on}

In our case the ligand is stable with a single binding mode. The RMSD fluctuates around 0.3{\AA}, with a slight fluctuation near the end of the simulation. This is more clearly seen in the histogram, see Figure \ref{fig:rmsdligand}. The conformation seen during simulation is very similar to that in the crystal structure and the ligand is stable in the active site.

![RMSD timeseries Hsp90 ligand](../../images/htmd_analysis_rmsd2_series.png "RMSD timeseries for the Hsp90 Residue ID G5E (ligand)")

![RMSD histogram Hsp90 ligand](../../images/htmd_analysis_rmsd2_histo.png "RMSD histogram for the Hsp90 Residue ID G5E (ligand)")

![Hsp90 ligand binding poses](../../images/htmd_analysis_bindingposes.png "Two binding poses seen for this ligand")


## RMSF analysis

The Root Mean Square Fluctuation (RMSF) is valuable to consider, as it represents the deviation at a reference position over time. The fluctuation in space of particular amino acids in the protein are considered. The C$\alpha$ of the protein, designated by \texttt{C-alpha}, is a good selection to understand the change in protein structure. Depending on the system these fluctuations can be correlated to experimental techniques including Nuclear Magnetic Resonance (NMR) and M\"{o}ssbauer spectroscopy~\cite{berjanskii_nmr_2006,kuzmanic_determination_2010}. The output from the tools is the requested RMSF data and the RMSF plotted as a time series (for example, see Figure \ref{fig:rmsf} in the results section).

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **RMSF Analysis** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"DCD trajectory input"*: `output` (output of **MDTraj file converter** {% icon tool %})
>    - {% icon param-file %} *"PDB input"*: `output` (output of **GROMACS structure configuration** {% icon tool %})
>    - *"Select domains"*: `C-alpha`
>
>
{: .hands_on}


When considering the RMSF (Figure \ref{fig:rmsf}), fluctuations greater than 1.0{\AA} are of interest; for example see the fluctuations near residue positions 50, 110 and 160.  Inspecting the structure with molecular visualization software such as VMD, these can be seen to correspond to flexible loop regions on the protein surface. In addition, very large fluctuations are seen for the C-terminus; this is common and no investigation is needed.

Note that the first few residues of this protein are missing in the PDB, and therefore residue position 0 in the RMSF corresponds to position 17 in the Hsp90 FASTA primary sequence. This is a fairly common problem that can occur with molecular modeling of proteins, where there may be missing residues at the beginning or within the sequence.

![RMSF Hsp90](../../images/htmd_analysis_rmsf.png "RMSF for Hsp90")


## PCA analysis

Principal component analysis (PCA) converts a set of correlated
observations (movement of selected atoms in protein) to a set of principal
components (PCs) which are linearly independent (or uncorrelated). Here several related tools are used.
The PCA tool calculates the PCA in order to determine the relationship between statistically meaningful conformations (major global motions) sampled during the trajectory. The C$\alpha$ carbons of the protein backbone are again a good selection for this purpose.  Outputs include the PCA raw data and figures of the relevant principal components (PCs) as well as an eigenvalue rank plot (see Figure \ref{fig:pca}) which is used to visualize the proportion of variance due to each principal component (remembering that the PCs are ranked eigenvectors based on the variance).
Having discovered the principal components usually these are visualized. The PCA visualization tool will create trajectories of specific principal components which can be viewed in a molecular viewer such as VMD~\cite{hump_vmd_1996} or NGL viewer~\cite{Rose2018ngl}. We also consider the PCA cosine content which when close to 1 indicates that the simulation is not converged and a longer simulation is needed. For values below 0.7, no statement can be made about convergence or lack thereof.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PCA** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"DCD trajectory input"*: `output` (output of **MDTraj file converter** {% icon tool %})
>    - {% icon param-file %} *"PDB input"*: `output` (output of **GROMACS structure configuration** {% icon tool %})
>    - *"Select domains"*: `C-alpha`
>
>
{: .hands_on}

The first three principal components are responsible for 32.8\% of the total variance, as seen in the eigenvalue rank plot (Figure \ref{fig:pca}). The first principal component (PC1) accounts for 15.4\% of the variance (see PC1 vs PC2 and eigenvalue rank plots in Figure \ref{fig:pca}). Visualization of PC1 using VMD shows a rocking motion and wagging of the C-terminus.

![PCA Hsp90](../../images/htmd_analysis_pca.png "PCA for Hsp90")


![PC1 RMSF Hsp90](../../images/htmd_analysis_pc1_rmsf.png "PC1 vs RMSF for Hsp90")


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Cosine Content** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"DCD/XTC trajectory input"*: `output` (output of **MDTraj file converter** {% icon tool %})
>    - {% icon param-file %} *"PDB/GRO input"*: `output` (output of **GROMACS structure configuration** {% icon tool %})
>
>
{: .hands_on}


The PCA cosine content of the dominant motion related to PC1 is 0.93, indicating that the simulation is not fully converged. This is expected due to the short simulation length. For production level simulations, it is the norm to extend simulations to hundreds of nanoseconds in length, if not microseconds. As this tutorial is designed to be carried out on public webservers, we limit simulations to 1 ns, as we cannot provide a large amount of computational resources for training purposes.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **PCA visualization** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"DCD trajectory input"*: `output` (output of **MDTraj file converter** {% icon tool %})
>    - {% icon param-file %} *"PDB input"*: `output` (output of **GROMACS structure configuration** {% icon tool %})
>    - *"Select domains"*: `C-alpha`
>
>
{: .hands_on}

![PC1 Hsp90 gif](../../images/htmd_analysis_pc1_hsp90.gif "PC1 motion for Hsp90")

## Hydrogen bond analysis

Hydrogen bonding interactions contribute to binding and are worth investigating, in particular persistent hydrogen bonds. All possible hydrogen bonding interactions between the two selected regions, here the protein and the ligand, are investigated over time using the VMD hydrogen bond analysis tool included in Galaxy. Hydrogen bonds are identified and in the output the total number of hydrogen bonds and  occupancy over time is returned.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Hydrogen Bond Analysis using VMD** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"DCD/XTC trajectory input"*: `output` (output of **MDTraj file converter** {% icon tool %})
>    - {% icon param-file %} *"PDB/GRO input"*: `output` (output of **GROMACS structure configuration** {% icon tool %})
>    - *"Selection 1"*: `protein`
>    - *"Selection 2"*: `resname G5E`
>
>
{: .hands_on}

The active site of this protein is quite hydrophobic, yet multiple hydrogen bonds were identified. The hydrogen bond between aspartate-93 and the ligand (as identified in the crystal structure) was found to be persistent, meeting the hydrogen bond criteria for 89.22\% of the simulation. A hydrogen bond between the ligand and the carbonyl group of glycine-97 was found to have a 15.27\% occupancy. Hydrogen bonding interactions with threonine-184, asparagine-51 and lysine-58 were also observed but these are not persistent and only present for a minority of the simulation. These values can be accessed from the 'Percentage occupancy of the H-bond' output of the hydrogen bond analysis tool.

## Optional: Estimate the binding free energy

We can estimate the binding free energy between a ligand and a receptor using Molecular Mechanics Poisson-Boltzman Surface Area (MMPBSA). For the simplest form of calculation, a Single Trajectory Estimate, a simulation of the complex in water is run, this was done with GROMACS previously. The trajectory of this complex is used to estimate the MMPBSA or MMGBSA depending on the options chosen. A General Born (GB) calculation is recommended here as it completes in reasonable time.

### Converting parameters

The binding free energy between a ligand and a receptor can be estimated using Molecular Mechanics Poisson-Boltzman or General Born Surface Area (MMPBSA or MMGBSA). For
the simplest form of calculation, a Single Trajectory Estimate, a
simulation of the complex in water is run. The trajectory of this complex is used to estimate the
MMPBSA or MMGBSA depending on the options chosen. A General Born (GB)
calculation is recommended here as it completes in reasonable time.

The MMPBSA_MMGBSA tool supports the AMBER topology format only. To convert from GROMACS
'.top' and '.gro' to AMBER topologies, use the `Convert Parameters' tool which uses
ParmEd\cite{Swails2016} to carry out the conversion. The result will be four AMBER `.prmtop'
files which are required for the binding free energy estimate
calculation. Selections are provided to the tool in order to remove unneccesary molecules. For example, to create the ligand `.prmtop' output, the ligand selection \texttt{!:G5E} selects and removes all molecules that are not the ligand.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Convert Parameters** {% icon tool %} with the following parameters:
>    - *"Force Field format"*: `GROMACS`
>        - {% icon param-file %} *"Input topology (top) file"*: `output` (Input dataset)
>        - {% icon param-file %} *"Input structure (gro) file"*: `output` (Input dataset)
>    - *"Ligand selection"*: `!:G5E`
>    - *"Receptor selection"*: `:NA,CL,SOL,G5E`
>    - *"Complex selection"*: `:NA,CL,SOL`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    > The example selection makes several assumptions. The ligand is named G5E, the solvent SOL and that the ions are NA and CL. These selection will depend on the system that was built.
>    >
>    {: .comment}
>

## Binding free energy estimate


> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **mmpbsa mmgbsa** {% icon tool %} with the following parameters:
>    - In *"Input"*:
>        - *"Single or Multiple Trajectories"*: `Single Trajectory Protocol (STP)`
>            - {% icon param-file %} *"AMBER prmtop input for Ligand"*: `output` (Input dataset)
>    - In *"General parameters"*:
>        - *"Final frame to analyse"*: `2000`
>        - *"interval between frames analysed"*: `5`
>        - *"quasi-harmonic entropy calculation"*: `Yes`
>        - *"Strip mask"*: `:NA:SOL:CL`
>    - In *"Details of calculation and parameters"*:
>        - *"General Born calculation"*: `yes`
>        - *"Poisson Boltzman calculation"*: `no`
>        - *"Decomposition Analysis"*: `yes`
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

TODO
An MMGBSA calculation using the single trajectory protocol yields a binding estimate of XXXX kcal/mol (YYYY kJ/mol) for this ligand. The negative value suggests the complex is thermodynamically stable, as we expect. Note that this value does not include the entropy contribution, but this is likely to be relatively small.

# Optional: Automating high throughput calculations
Up until this step, Galaxy tools have been applied sequentially to datasets. This is useful to gain an understanding of the steps involved, but becomes tedious if the workflow needs to be run on multiple protein-ligand systems. Fortunately, Galaxy allows entire workflows to be executed with a single mouse-click, enabling straightforward high-throughput analyses.

We will demonstrate the high-throughput capabilities of Galaxy by running the workflow detailed so far on a further three ligands.


<!-- \begin{tcolorbox}
\textbf{High-throughput MD}
%\textbf{Hands-on: Data upload}

\begin{enumerate}
\def\labelenumi{\arabic{enumi}.}
\tightlist
\item
  Create a new history for running the high-throughput workflow and name it `Hsp90 HTMD simulation'

\item
  Upload the SD-file containing the new ligand structures from Zenodo (LINK...) and rename it `Ligands (SDF)'
\item
  Import the simulation workflow from the European~\cite{eu_htmd_simulation_workflow} or the South African Galaxy server~\cite{za_htmd_simulation_workflow}.
\item
  Run the imported workflow with the following parameters:

  \begin{itemize}
  \tightlist
  \item
    \emph{``SDF file with (docked) ligands''}: `Ligands (SDF)' file.
  \end{itemize}
\item
  Import the analysis workflow from the European~\cite{eu_htmd_analysis_workflow} or the South African Galaxy server~\cite{za_htmd_analysis_workflow} (also available through Zenodo).
\item
  Run the imported workflow with the following parameters:
  \begin{itemize}
  \tightlist
  \item
    \emph{``Send results to a new history''}: `Yes'
  \item
    \emph{``History name''}: `Hsp90 HTMD analysis'
  \item
    \emph{``GRO input''}: Collection of GRO files produced by simulation workflow
  \item
    \emph{``XTC input''}: Collection of XTC files produced by simulation workflow
  \end{itemize}
\end{enumerate}
\end{tcolorbox} -->

This process runs the entire simulation and analysis procedure described so far on the new set of ligands. It uses Galaxy's collection~\cite{gtn_collections} feature to organize the data; each item in the history is a collection (essentially a directory containing multiple individual datasets) containing one file corresponding to each of the input ligands.

Note that the SD-file needs to contain ligands with the correct 3D coordinates for MD simulation. The easiest way to obtain these is using a molecular docking tool such as Autodock Vina~\cite{Trott2009} or rDock~\cite{Ruiz2014}; tutorials and workflows are available for both of these from the Galaxy Training Network. As an example, the history in which the SD-file used in the HTMD workflow is generated (using AutoDock Vina) is provided~\cite{eu_6hhr}.


Apart from manual setups or collections, there are several other alternatives which are helpful in scaling up workflows. Galaxy supports and provides training material for converting histories to workflows~\cite{gtn_toworkflow}, using multiple histories~\cite{gtn_multiple}, and the Galaxy Application Programming Interface (API)~\cite{gtn_api}. For beginners and users who prefer a visual interface, automation can be done using multiple histories and collections with the standard Galaxy user interface.

If you are able to write small scripts, you can automate everything you have learned here with the Galaxy API. This allows you to interact with the server to automate repetitive tasks and create more complex workflows (which may have repetition or branching). The simplest way to access the API is through the Python library BioBlend~\cite{sloggett_bioblend}. An example Python script, which uses BioBlend to run the GROMACS simulation workflow for each of a list of ligands, is given in the hands-on box below.

<!-- \begin{tcolorbox}
\textbf{BioBlend script}
\begin{minted}[linenos=false, breaklines, breakafter=d, fontsize=\small]{python}
from bioblend import galaxy

# Server and account details
API_KEY = 'YOUR USEGALAXY.EU API KEY'
gi = galaxy.GalaxyInstance(key=API_KEY,
    url='https://usegalaxy.eu/')

# ID for GROMACS workflow
workflow_id = 'adc6d049e9283789'

# Dataset IDs for ligands to dock
ligands = {
# ligand_name: dataset ID,
'lig1': '11ac94870d0bb33a79c5fa18b0fd3b4c',
# ...
}

# Loop over ligands, invoking workflow
for name, _id in ligands.items():
    inv = gi.workflows.invoke_workflow(
        workflow_id,
        inputs={
            '1': {'src': 'hda', 'id': _id}
        },
        history_name=f'HTMD run on {name}'
    )
\end{minted}
\end{tcolorbox} -->


# Conclusion
{:.no_toc}

This tutorial provides a guide on how to study protein-ligand interaction using molecular dynamics in Galaxy. Performing such analyses in Galaxy makes it straightforward to set up, schedule and run workflows, removing much of the difficulty from MD simulation. Thus, the technical barrier to performing high-throughput studies is greatly reduced. Results are structured in the form of Galaxy histories or collections, and include ready-plotted diagrams, which ensure data can be easily understood and reproduced if necessary. Apart from streamlining the process for existing MD users, this tutorial should also prove useful as a pedagogical guide for educating students or newcomers to the field.

After completing the tutorial, the user will be familiar at a basic level with a range of MD analysis techniques, and understand the steps required for a typical MD simulation. Thus, they will be equipped to apply these tools to their own problems.
