---
layout: tutorial_hands_on

title: Protein docking
zenodo_link: ''
questions:
- What is docking?
- How can I perform a simple docking workflow in Galaxy?
objectives:
- Dock a variety of ligands to the active site of the Hsp90 protein
time_estimation: 3H
key_points:
- Docking allows 'virtual screening' of drug candidates
contributors:
- simonbray

---


# Introduction
{:.no_toc}

The aim of protein-ligand docking is to find the optimal binding between a small molecule and a protein. Generally, the goal is to search for a potential drug candidate. Firstly, a target protein is identified which is involved in a disease. Secondly, a 'library' of ligands which may be able to bind to this protein and interfere with its function is assembled. Each of the compounds is then 'docked' into the protein to find the optimal binding position and energy.

Docking is a form of molecular modelling, but several simplifications are made in comparison to methods such as molecular dynamics. Most significantly, the receptor is generally considered to be rigid, with covalent bond lengths and angles held constant. Charges and protonation states are also not permitted to change. While these approximations lower accuracy to some extent, they increase computational speed, which is necessary to screen a large compound library in a realistic amount of time.

In this tutorial, you will perform docking of ligands into the N-terminus of the Hsp90 protein. The tools used for docking are based on the open-source software [Autodock Vina](http://vina.scripps.edu/).

> ### {% icon details %} Biological background
>
> The 90 kDa heat shock protein (Hsp90) is a chaperone protein responsible for catalyzing the conversion of a wide variety of proteins to a functional form; examples of the Hsp90 clientele, which totals several hundred proteins, include nuclear steroid hormone receptors and protein kinases. The mechanism by which Hsp90 acts varies between clients, as does the client binding site; the process is dependent on post-translational modifications of Hsp90 and the identity of co-chaperones which bind and regulate the conformational cycle. 
>
>
> Due to its vital biochemical role as a chaperone protein involved in facilitating the folding of many client proteins, Hsp90 is an attractive pharmaceutical target. In particular, as protein folding is a potential bottleneck to slow cellular reproduction and growth, blocking Hsp90 function using inhibitors which bind tightly to the ATP binding site could assist in treating cancer; for example, the antibiotic geldanamycin and its analogs are under investigation as possible anti-tumor agents. 
>
{: .details}



> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Download data

For this exercise, we need two datasets: a protein structure and a library of compounds. We will download the former directly from the Protein Data Bank; the latter is available from Zenodo.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Search Galaxy for the 'Get PDB' tool. Request the accession code ```2crb```.
> 2. Now, import the compound library from [Zenodo]() or from the shared data library
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

# Separating protein and ligand structures

You can view the contents of the downloaded PDB file by pressing the 'View data' icon in the history pane. Scrolling past the first few slides, the lines containing the atomic coordinates of the protein begin with ```ATOM```, while those describing the atoms of the ligand and the solvent water molecules are labelled ```HETATM```. We can use the grep tool to separate these molecules into separate files.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Find the grep tool in Galaxy
> 2. **grep** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: Downloaded PDB file
>    - {% icon param-file %} *"that"*: `Don't match`
>    - {% icon param-file %} *"Regular Expression"*: `HETATM`
>
>   All other parameters can be left as their defaults. The result is a file with all non-protein atoms removed. Rename the dataset 'Protein'
> 2. Now repeat the previous step, but this time enter `ATOM` under *"Regular Expression"*. This produces a file with only non-protein atoms. Rename the dataset 'Ligand'. (Note: the file does also contain water molecules, but these will be removed in the next step).
> 3. **openbabel_compound_convert** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: Ligand PDB file created in step 2.
>    - {% icon param-file %} *"Output format"*: `MDL MOL format (sdf, mol)`
>    - {% icon param-file %} *"Generate 3D coordinates"*: `Yes`
>
>   All other parameters can be left as their defaults. Applying this tool will generate a representation of the structure of the ligand in MOL format.
> 4. **autodock_vina_prepare_box** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input ligand"*: Ligand MOL file created in step 3.
>    - {% icon param-file %} *"x-axis buffer"*: `5`
>    - {% icon param-file %} *"y-axis buffer"*: `5`
>    - {% icon param-file %} *"z-axis buffer"*: `5`
>    - {% icon param-file %} *"Exhaustiveness"*: `1`
> 5. **autodock_vina_prepare_receptor** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select a PDB file"*: Protein PDB file created in step 1.
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

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.



***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:

# Creating and processing the compound library

The compound library is provided as a SMILES file. For docking, we would like to convert to SDF format.

> ### {% icon details %} What are SMILES and SDF formats?
>
> SMILES and SD-files both represent chemical structures. A SMILES file represents the 2D structure of a molecule as a chemical graph. By contrast, the SDF format stores the three-dimensional coordinates of a structure, similarly to a PDB file.
>
> For docking, we need the three-dimensional coordinates of the ligand; thus we want to convert from SMILES to SDF format.
{: .details}

