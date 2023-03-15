---
layout: tutorial_hands_on

title: Virtual screening of the SARS-CoV-2 main protease with rxDock and pose scoring
level: Intermediate
zenodo_link: 'https://zenodo.org/record/3730474'
questions:
- How can candidate ligands be generated and docked to a protein in Galaxy?
- How can the poses of the docked ligands be evaluated?
- How can a workflow for drug virtual screening be constructed in Galaxy?
objectives:
- Understand how Galaxy was used to perform docking and pose scoring on the SARS-CoV-2 main protease (MPro).
- Replicate the study on a (very) small scale
- Gain familiarity with the docking and scoring techniques involved.
time_estimation: 2H  # Just 1 week (if you have 5000 CPUs) ;)
key_points:
- Galaxy can support large, rapid studies in computational chemistry
- Protein-ligand docking contributes to the discovery of new drugs
requirements:
  -
    type: "internal"
    topic_name: computational-chemistry
    tutorials:
        - cheminformatics
tags:
- covid19
contributors:
- simonbray

---

# Introduction


This tutorial provides a companion to the work performed in March 2020 by InformaticsMatters, the Diamond Light Source, and the European Galaxy Team to perform virtual screening on candidate ligands for the SARS-CoV-2 main protease (MPro). This work is described [here](https://covid19.galaxyproject.org/cheminformatics).

In this tutorial, you will perform protein-ligand docking to MPro using rxDock ({% cite rdock %}), a version of the popular rDock software, and score the results using two different methods. The same tools will be used as in the original study, but with a smaller dataset.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background

Early in March 2020, the Diamond Light Source completed a successful fragment screen on MPro, which provided 55 fragment hits (see their [press release](https://www.diamond.ac.uk/covid-19/for-scientists/Main-protease-structure-and-XChem.html) here). In an effort to identify candidate molecules for binding, InformaticsMatters, the XChem group and the European Galaxy team joined forces to construct and execute a Galaxy workflow for performing and evaluating molecular docking on a massive scale.

An initial list of 41,000 candidate molecules was assembled by using the Fragalysis fragment network to elaborate from the initial fragment hits, as described [here](https://diamondlightsource.atlassian.net/wiki/spaces/FRAG/pages/8323192/The+Astex+Fragment+network). These were used as inputs for the docking and scoring workflow. The workflow consists of the following steps, each of which was carried out using tools installed on the European Galaxy server:
1. Charge enumeration of the 41,000 candidate molecules selected based on the fragment hits.
2. Generation of 3D conformations based on SMILES strings of the candidate molecules.
3. Docking of molecules into each of the MPro structures using rxDock.
4. Evaluation of the docking poses using a TransFS, a deep learning approach ({% cite transfs %}) developed by the XChem group and collaborators, and SuCOS scoring ({% cite sucos %}), which compares the poses with the structures of the original fragment hits.

The original study required almost 20 years of CPU time, not counting GPU resources consumed. This is obviously not reproducible as a tutorial. Therefore, we will repeat the workflow with a small library of just 100 molecules, on a single MPro fragment structure. Links will be provided to original Galaxy histories, with notes to explain where and why things were done differently to the tutorial.

![MPro structure, with a fragment bound]({% link topics/computational-chemistry/images/mpro.png %} "Structure of MPro, with a fragment bound. <a href="https://usegalaxy.eu/u/sbray/v/mpro-x0072">Click to view</a> in NGL. ({% cite ngl %})")

# Get data

We require three datasets for the simulation and analysis:
1. A list of 100 ligand candidates. These are the molecules which will be docking into the protein binding site.
2. A PDB file of the receptor MPro protein (without ligand or solvent).
3. A list of fragment hits (17 in total) in SDF format.

> <details-title>Differences with the original study</details-title>
>
> Of the initial 55 fragment hits, 17 were chosen for further study. From these, 41,587 compounds were generated using the Fragalysis fragment network for further study. The 100 compounds used in the tutorial are taken from this list.
>
>
> Starting data is available from this Galaxy history: [https://usegalaxy.eu/u/sbray/h/mpro-raw-data](https://usegalaxy.eu/u/sbray/h/mpro-raw-data).
> 
> 
> This history contains 103 files. One of these (`Initial candidates for docking`) contains the 41k candidate compounds in SMILES format. The remaining 102 files (all with names beginning with `Mpro-x...`) provide structural information on the fragment hits - 6 files per hit (hence 17 x 6 = 102). 
>
>
> The identity of the files is as follows:
>
> - the `*_0.mol` files contain the fragment structure in mol format.
> - the `*_0.pdb` files contain the fragment structure in pdb format.
> - the `*_0_apo.pdb` files contain the protein with solvent, but without ligand
> - the `*_0_apo-desolv.pdb` files contain the protein without either solvent or ligand
> - the `*_0_apo-solv.pdb` files contain only solvent
> - the `*_0_bound.pdb` file contain everything (protein, ligand and solvent)
>
> The PDB file of the receptor that we are using is `Mpro-x0195_0_apo-desolv.pdb`. In other words, the structure is derived from just one fragment hit. In the original study, however, all compounds were docked against all of the fragment hit structures.
{: .details}

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo](https://zenodo.org/record/3730474):
>
>    ```
>    https://zenodo.org/record/3730474/files/candidates.smi
>    https://zenodo.org/record/3730474/files/Mpro-x0195_0_apo-desolv.pdb
>    https://zenodo.org/record/3730474/files/hits.sdf
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets `Candidates SMILES`, `Receptor PDB` and `Hits SDF` respectively.
> 4. Check that the datatypes (`smi`, `pdb`, and `sdf` respectively) are correct. In particularly, check the `Candidates SMILES` file, as the SMILES datatype is not detected automatically by Galaxy.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
{: .hands_on}

# Preparation for docking

Before docking, the candidate ligands need to be prepared for docking with the following steps: 1) charge enumeration, 2) generation of three-dimensional structures, and 3) splitting the SD-file into a collection.

> <details-title>Differences with the original study</details-title>
>
> This stage is carried out as described here, except of course with the full set of 42,000 compounds. See [here](https://covid19.galaxyproject.org/cheminformatics/1-DockingPrep/) for more details.
{: .details}

## Charge enumeration

Many of the compounds may contain functional groups which can exist in multiple charge states, and this will affect the quality of binding to the receptor dramatically. Therefore, we perform 'charge enumeration', which means that we generate all charge forms of the compounds within a certain pH range.

> <hands-on-title>Charge enumeration</hands-on-title>
>
> 1. {% tool [Enumerate charges](toolshed.g2.bx.psu.edu/repos/bgruening/enumerate_charges/enumerate_charges/2020.03.4+galaxy0) %} with the following parameters:
>    - *"Input molecule data"*: `Candidate SMILES`
>    - *"Minimum pH"*: `4.4`
>    - *"Maximum pH"*: `10.4`
> 2. Rename the output file `Enumerated candidates SMILES`.
>
>    > <comment-title></comment-title>
>    >
>    > The **Enumerate charges** {% icon tool %} tool is based on the Dimorphite-DL program. ({% cite Ropp2019 %})
>    {: .comment}
>
{: .hands_on}

The output is another SMILES file, with several hundred entries.

## Generate three-dimensional conformations

So far our list of enumerated candidate compounds is still in SMILES format; we need to produce three-dimensional structures in SDF format for docking. This can be done with the **Compound conversion** {% icon tool %} tool.

If you are not familiar with SMILES and SDF formats, consult the introductory [protein-ligand docking tutorial](../cheminformatics/tutorial.html) for more details.

> <hands-on-title>Convert to SDF format</hands-on-title>
>
> 1. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - *"Molecular input file"*: `Enumerated candidates` dataset.
>    - *"Output format"*: `MDL MOL format (sdf, mol)`
>        - *"Generate 3D coordinates"*: `Yes`
> 2. Rename the output file `Enumerated candidates SDF`.
>
>    > <comment-title></comment-title>
>    >
>    > The **Compound conversion** {% icon tool %} tool is based on the OpenBabel toolkit. ({% cite OBoyle2011 %})
>    {: .comment}
>
{: .hands_on}

## Splitting the SD-file into a collection

The next stage is to split the SD-file with the candidate ligands into a set of smaller SD-files.


> <question-title></question-title>
>
> Why is splitting the file necessary?
>
> > <solution-title></solution-title>
> >
> > The rxDock tool performs one docking at a time (more technically: the task is not parallelized, as it uses only a single CPU). Therefore, splitting the large SD-file into many small files allows the work to be carried out by multiple Galaxy jobs in parallel, so it completes faster.
> >
> > In the original study, this kind of parallelization was absolutely essential because of the enormous dataset; at some points, there were 5,000 docking jobs running concurrently on the European Galaxy server. Even on the much smaller scale of this tutorial, we can speed things up considerably using this trick.
> >
> {: .solution}
{: .question}

> <hands-on-title>Split the SD-files</hands-on-title>
>
> 1. {% tool [Split file to dataset collection](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_to_collection/split_file_to_collection/0.5.0) %} with the following parameters:
>    - *"Select the file type to split"*: `SD-files`
>    - *"SD-file to split"*: `Enumerated candidates SDF`
>        - *"Specify number of output files or number of records per file?"*: `Number of output files`
>            - *"Number of new files"*: `10`
>        - *"Method to allocate records to new files"*: `Alternate output files`
>
>
{: .hands_on}

# Active site preparation

The active site also needs to be prepared for docking, using the following steps: 1) conversion to MOL2 format, and 2) generation of the active site using the **rbcavity** {% icon tool %} tool.

> <details-title>Differences with the original study</details-title>
>
> This stage was carried out as described here. However, it was repeated for each of the fragment hit structures, not just the `Mpro-x0195_0_apo-desolv.pdb` file used here. See [here](https://covid19.galaxyproject.org/cheminformatics/2-ActiveSitePrep/) for more details.
{: .details}

## Convert protein structure to MOL2 format

The receptor file we are using is in PDB format, but the rxDock tool we use for docking requires an input in MOL2 format. Therefore, we first convert the file.

> <hands-on-title>Conversion to MOL2 format</hands-on-title>
>
> 1. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - *"Molecular input file"*: `Receptor PDB` dataset.
>    - *"Output format"*: `Sybyl Mol2 format (mol2)`
> 2. Rename the output file `Receptor MOL2`.
>
>
{: .hands_on}

## Generate Frankenstein ligand

For docking with rxDock, a file needs to be created defining the active site. This requires two input files - one for the protein and one for the ligand. We want an active site generation that takes into account the features of all 17 fragments, and therefore need to generate a 'Frankenstein ligand' which possesses the properties of all the fragments. A very simple Galaxy tool is available for this.

> <question-title></question-title>
>
> What is a 'Frankenstein ligand' and why do we need it?
>
> > <solution-title></solution-title>
> >
> > The Frankenstein ligand combines the atoms of all the fragments and therefore occupies the entire space of the binding site. Therefore, it is the best choice for cavity definition. See the [information provided by InformaticsMatters](https://www.informaticsmatters.com/blog/2018/11/23/cavities-and-frankenstein-molecules.html) for more details.
> >
> {: .solution}
>
{: .question}

> <hands-on-title>Generate Frankenstein ligand</hands-on-title>
>
> 1. {% tool [Create Frankenstein ligand](toolshed.g2.bx.psu.edu/repos/bgruening/ctb_frankenstein_ligand/ctb_frankenstein_ligand/2013.1-0+galaxy0) %} with the following parameters:
>    - *"Input file"*: `Hits SDF`
> 2. Rename the file to `Frankstein SDF`.
>
{: .hands_on}

## Generate active site definition

The active site can now be generated using the **rbcavity** {% icon tool %} tool, which requires the receptor in MOL2 format as input as well as a single reference ligand in Mol/SDF format. We use the Frankenstein ligand as the reference.

> <hands-on-title>Active site preparation</hands-on-title>
>
> 1. {% tool [rxDock cavity definition](toolshed.g2.bx.psu.edu/repos/bgruening/rxdock_rbcavity/rxdock_rbcavity/0.1.5) %} with the following parameters:
>    - *"Receptor"*: `Receptor MOL2`
>    - *"Reference ligand"*: `Frankenstein SDF`
>    - *"Mapper sphere radius"*: `3.0`
>    - *"Mapper small sphere radius"*: `1.0`
>    - *"Mapper minimum volume"*: `100`
>    - *"Mapper volume increment"*: `0`
>    - *"Mapper grid step"*: `0.5`
>    - *"Cavity weight"*: `1.0`
>
> 2. Rename the output file `Active site`.
>    > <comment-title></comment-title>
>    >
>    > The meanings of these parameters are too complex to go into in this tutorial. If you are interested, see the [rDock documentation](http://rdock.sourceforge.net/wp-content/uploads/2015/08/rDock_User_Guide.pdf) for more details.
>    {: .comment}
>
{: .hands_on}

# Docking and scoring

Docking and scoring are now performed, using the following steps: 1) docking using rxDock, 2) recombining the results into a single SDF file, 3) TransFS scoring, and 4) SuCOS scoring.

> <details-title>Differences with the original study</details-title>
>
> This section in the original study differed from this tutorial in the following ways:
> 1. Docking was performed on over 100,000 enumerated candidates, rather than the 300 used here.
> 2. 25 different poses were generated per candidate, rather than 5, as in this tutorial.
> 3. Because of the large number of poses to score (more than a million), the scoring steps were parallelized by splitting into collections. This is skipped in the tutorial.
> 4. The entire process was repeated 17 times, using a different fragment hit as the receptor structure each time.
> See [here](https://covid19.galaxyproject.org/cheminformatics/3-Docking/) and [here](https://covid19.galaxyproject.org/cheminformatics/3-Docking/) for more details. A full list of Galaxy histories generated is listed [here](https://covid19.galaxyproject.org/cheminformatics/Histories/).
{: .details}

## Docking with rxDock

> <hands-on-title>Docking</hands-on-title>
>
> 1. {% tool [rxDock docking](toolshed.g2.bx.psu.edu/repos/bgruening/rxdock_rbdock/rxdock_rbdock/0.1.5) %} with the following parameters:
>    - *"Receptor"*: `Receptor MOL2`
>    - *"Active site"*: `Active site`
>    - *"Ligands"*: `Split file` collection
>    - *"Number of dockings"*: `5`
>    - *"Number of best poses"*: `5`
>
>    > <comment-title></comment-title>
>    >
>    > For more information about docking, check out the [introductory tutorial](../cheminformatics/tutorial.html). It uses a different tool, AutoDock Vina, rather than rxDock, but the general principles are the same.
>    {: .comment}
>
{: .hands_on}

## Collapse collection to a single file 

Having created a collection to parallelize the docking procedure, we can now recombine the results to a single file.

> <hands-on-title>Collapse collection</hands-on-title>
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.2) %} with the following parameters:
>    - *"Collection of files to collapse into single dataset"*: Output of docking step
> 2. Rename to `Docked poses SDF`.
>
{: .hands_on}

The output file should contain around 1,900 docked poses in SDF format.

## Pose scoring with TransFS

In this step, we carry out scoring of the poses using TransFS. This is a deep learning approach developed at the University of Oxford, employing augmentation of training data with incorrectly docked ligands to prompt the model to learn from protein-ligand interactions. ({% cite transfs %})

The TransFS scoring returns a value (saved as `<TransFSScore>` in the SDF file) between 0 (poor) and 1 (good).

> <hands-on-title>TransFS scoring</hands-on-title>
>
> 1. {% tool [XChem TransFS pose scoring](toolshed.g2.bx.psu.edu/repos/bgruening/xchem_transfs_scoring/xchem_transfs_scoring/0.4.0) %} with the following parameters:
>    - *"Receptor"*: `Receptor PDB`
>    - *"Ligands"*: `Docked poses SDF`
>    - *"Distance to waters"*: `2`
>
{: .hands_on}

## Pose scoring with SuCOS

This step involves scoring of the poses from each molecule against the original fragment screening hit ligands using the SuCOS MAX shape and feature overlay algorithm. ({% cite sucos %}) The conformation and position of the poses are compared to known structures (i.e. the fragment hits) to determine a score.

SuCOS scoring returns a value (saved as `<Max_SuCOS_Score>` in the SDF file) between 0 (poor) and 1 (good).

> <hands-on-title>SuCOS scoring</hands-on-title>
>
> 1. {% tool [Max SuCOS score](toolshed.g2.bx.psu.edu/repos/bgruening/sucos_max_score/sucos_max_score/2020.03.4+galaxy0) %} with the following parameters:
>    - *"Ligands to be scored"*: Output of the TransFS step
>    - *"Set of clusters to score against"*: `Hits SDF`
> 2. Rename the output file to `Scored poses`.
{: .hands_on}

# Compound selection

The aim of the original study was to select 500 candidate molecules for synthesis and experimental study. In order to do this, the data for all fragment hits had to be combined (i.e. so that each compound was assigned the lowest score from all the fragments). The resulting table was then compared with a list of compounds available from [Enamine](https://enamine.net/) and [Chemspace](https://chem-space.com/) and the 500 highest scoring matching compounds were selected for purchase.

This step is skipped in the tutorial, as only 100 compounds were tested, using only a single fragment hit structure. If you want, though, check out the [history](https://usegalaxy.eu/u/timdudgeon/h/top-500-enamine--chemspace-bb) and [workflow](https://usegalaxy.eu/u/timdudgeon/w/filter-results) used.

# Conclusion


This tutorial guided you through docking and scoring of a small set of compounds to the MPro protein. Hopefully, you have a better understanding of how docking can be practically used, as well as how the original study was performed.
