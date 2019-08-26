---
layout: tutorial_hands_on

title: Protein-ligand docking
zenodo_link: ''
questions:
- What is cheminformatics?
- What is protein-ligand docking?
- How can I perform a simple docking workflow in Galaxy?
objectives:
- Create a miniature compound library using the ChEMBL database
- Dock a variety of ligands to the active site of the Hsp90 protein
time_estimation: 3H
key_points:
- Docking allows 'virtual screening' of drug candidates
- Molecular fingerprints encode features into a bitstring
- The ChemicalToolbox contains many tools for cheminformatics analysis
contributors:
- simonbray

---


# Introduction
{:.no_toc}

The aim of protein-ligand docking is to find the optimal binding between a small molecule and a protein. Generally, the goal is to search for a potential drug candidate. Firstly, a target protein is identified which is involved in a disease. Secondly, a 'library' of ligands which may be able to bind to this protein and interfere with its function is assembled. Each of the compounds is then 'docked' into the protein to find the optimal binding position and energy.

Docking is a form of molecular modelling, but several simplifications are made in comparison to methods such as molecular dynamics. Most significantly, the receptor is generally considered to be rigid, with covalent bond lengths and angles held constant. Charges and protonation states are also not permitted to change. While these approximations lower accuracy to some extent, they increase computational speed, which is necessary to screen a large compound library in a realistic amount of time.

In this tutorial, you will perform docking of ligands into the N-terminus of the Hsp90 protein. The tools used for docking are based on the open-source software [AutoDock Vina](http://vina.scripps.edu/).

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

For this exercise, we need two datasets: a protein structure and a library of compounds. We will download the former directly from the Protein Data Bank; the latter will be created by searching the ChEMBL database.

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Search Galaxy for the 'Get PDB' tool. Request the accession code ```2brc```.
> 3. Rename the dataset to 'Hsp90 structure'
> 4. Check that the datatype is correct (PDB file).
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

![Hsp90 N-terminus structure]({{ site.baseurl }}{% link topics/computational-chemistry/images/2brc.png %} "Structure of Hsp90 N-terminus, as recorded on the PDB. Visualization produced using VMD.")

# Separating protein and ligand structures

You can view the contents of the downloaded PDB file by pressing the 'View data' icon in the history pane. After the header section (about 500 lines), tha atoms of the protein and their coordinates are listed. The lines begin with ```ATOM```. At the end of the file, the atomic coordinates of the ligand and the solvent water molecules are also listed, labelled ```HETATM```. We can use the grep tool to separate these molecules into separate files.

> ### {% icon details %} What is grep?
>
> Grep is a command line tool for searching textfiles for lines which match a search query.
>
> There is a Galaxy tool available based on grep, which we will apply to the downloaded PDB tool.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. **Search in textfiles (grep)** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: Downloaded PDB file
>    - {% icon param-file %} *"that"*: `Don't match`
>    - {% icon param-file %} *"Regular Expression"*: `HETATM`
>    - All other parameters can be left as their defaults. The result is a file with all non-protein (`HETATM`) atoms removed. Rename the dataset 'Protein (PDB)'.
> 2. **Search in textfiles (grep)** {% icon tool %} with the following parameters. Here, we use grep again to produce a file with only non-protein atoms.
>    - {% icon param-file %} *"Select lines from"*: Downloaded PDB file
>    - {% icon param-file %} *"that"*: `Don't match`
>    - {% icon param-file %} *"Regular Expression"*: `ATOM`
>    - All other parameters can be left as their defaults. This produces a file with only non-protein atoms. Rename the dataset 'Ligand (PDB)'. (Note: the file does also contain water molecules, but these will be removed in the next step).
> 3. **Compound conversion** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: Ligand PDB file created in step 2.
>    - {% icon param-file %} *"Output format"*: `MDL MOL format (sdf, mol)`
>    - {% icon param-file %} *"Generate 3D coordinates"*: `Yes`
>    - All other parameters can be left as their defaults. Applying this tool will generate a representation of the structure of the ligand in MOL format. Rename the dataset 'Ligand (MOL)'.
>
{: .hands_on}

At this stage separate protein and ligand files have been created. Next, we want to generate a compound library we can use for docking.


# Creating and processing the compound library

In this step we will create a compound library, using data from the ChEMBL database.

> ### {% icon details %} What chemical databases are available?
>
> Multiple databases are available online which provide access to chemical information, e.g. chemical structures, reactions, or literature. In this tutorial, we use a tool which searches the [ChEMBL](https://www.ebi.ac.uk/chembl/) database. There are also Galaxy tools available for downloading data from [PubChem](https://pubchem.ncbi.nlm.nih.gov/). 
>
> Other notable databases include [Drugbank](https://www.drugbank.ca/) and [ChemSpider](http://www.chemspider.com/).
>
{: .details}

We will generate our compound library by searching ChEMBL for compounds which have a similar structure to the ligand in the PDB file we downloaded in the first step. There is a Galaxy tool for accessing ChEMBL which requires a data input in SMILES format; thus, the first step is to convert the 'Ligand' PDB file to a SMILES file. Then the search is performed, returning a SMILES file. For docking, we would like to convert to SDF format, which we can do once again using the 'Compound conversion' tool.

> ### {% icon hands_on %} Hands-on: Generate compound library
>
> 1. **Compound conversion** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: 'Ligand' PDB file
>    - {% icon param-file %} *"Output format"*: `SMILES format (SMI)`
>    - Leave all other options as default.
> 2. Rename the output of the 'compound conversion' step to 'Ligand SMILES'.
> 3. **Search ChEMBL database** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"SMILES input type"*: File
>    - {% icon param-file %} *"Input file"*: 'Ligand SMILES' file
>    - {% icon param-file %} *"Search type"*: `Similarity`
>    - {% icon param-file %} *"Tanimoto cutoff score"*: `80`
>    - {% icon param-file %} *"Filter for Lipinski's Rule of Five"*: `Yes`
>    - All other parameters can be left as their defaults.
> > ### {% icon question %} Question
> >
> > Why are compounds filtered for Lipinski's Rule of Five?
> >
> > > ### {% icon solution %} Solution
> > > Lipinski's rule of five is an empirical rule which can be used to determine the 'druglikeness' of a molecule. The rule consists of four criteria which relate to the pharmokinetics of the molecule.
> > {: .solution}
> {: .question}
>
> 5. Optional: experiment with different combinations of options - adding different filters, adjusting the Tanimoto coefficient.
> 6. Check that the datatype is correct (smi). This step is essential, as Galaxy does not automatically detect the datatype for SMILES files.
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
> 7. Rename dataset 'Compound library'.
>
{: .hands_on}

> ### {% icon details %} What are SMILES and SDF formats?
>
> SMILES and SD-files both represent chemical structures. A SMILES file represents the 2D structure of a molecule as a chemical graph. In other words, it states only the atoms and the connectivity between them. An example of a SMILES string (taken from the ligand in the PDB file) is `c1c2OCCOc2ccc1c1c(C)[nH]nc1c1cc(CC)c(O)cc1O`. For more information on how the notation works, please consult the [OpenSMILES specification](http://opensmiles.org/opensmiles.html) or the description provided by [Wikipedia](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system). A more comprehensive alternative to the SMILES system is the International Chemical Identifier (InChI).
>
> Neither SMILES nor InChIs display the three-dimensional structure of a molecule. By contrast, the SDF (structure data file) format encodes three-dimensional atomic coordinates of a structure, similarly to a PDB file.
>
> In a previous step, we also generated a MOL file - this format is closely related to the SDF format. The difference is that MOL files can store only a single molecule, whereas SD-files can encode single or multiple molecules. Multiple molecules are separated by lines containing four dollar signs (`$$$$`).
>
> For docking, we need the three-dimensional coordinates of the ligand; thus, we want to convert from SMILES to SDF format.
{: .details}

# Aside: other options for the compound library

<!-- Add hydrogen atoms at a certain pH value
Calculate molecular descriptors with Mordred
Change title to metadata value.
Compound search - an advanced molecular search program using SMARTS
Compute physico-chemical properties for a set of molecules
Conformer calculation for molecules (confab)
Descriptors calculated with RDKit
Filter a set of molecules from a file
Fragmenter splits a molecule in fragments
Merging fragmented molecules
Molecules to Fingerprints with different fingerprint types
Multi Compound Search an advanced molecular grep program using SMARTS
Multi Compound Search an advanced molecular grep program using SMARTS
Natural Product likeness calculator
NxN Clustering of molecular fingerprints
Open Molecule Generator Exhaustive generation of chemical structures
Remove counterions and fragments
Remove protonation state of every atom
Remove small molecules from a library of compounds
SDF to Fingerprint extract fingerprints from sdf files metadata
Similarity Search of fingerprint data sets
Substructure Search of fingerprint data sets
Substructure Search of fingerprint data sets
Taylor-Butina Clustering of molecular fingerprints
Visualisation of compounds -->

The ChemicalToolBox contains a large number of cheminformatics tools. This section will demonstrate some of the useful functionalities available. If you are just interested in docking, feel free to skip this section - or, just try out the tools which look particularly interesting.

### Visualization

It can be useful to visualize the compounds generated. There is a tool available for this in Galaxy based on OpenBabel, an open-source library for analyzing chemical data.

> ### {% icon hands_on %} Hands-on: Visualization of chemical structures
>
> 1. **Visualisation** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: Compound library
>    - {% icon param-file %} *"Embed molecule as CML"*: `No`
>    - {% icon param-file %} *"Draw all carbon atoms"*: `No`
>    - {% icon param-file %} *"Use thicker lines"*: `No`
>    - {% icon param-file %} *"Property to display under the molecule"*: `Molecule title`
>    - {% icon param-file %} *"Sort the displayed molecules by"*: `Molecular weight`
>    - {% icon param-file %} Format of the resultant picture"*: `SVG`
>
{: .hands_on}

This produces an SVG image of all the structures generated ordered by molecular weight.

![Hsp90 N-terminus structure]({{ site.baseurl }}{% link topics/computational-chemistry/images/compound_library.png %} "Structures of the compounds from ChEMBL.")

### Calculation of fingerprints and clustering

In this step we will assess the similarity of the compounds in our small library to each other and cluster them into groups accordingly. A key tool in cheminformatics for measuring molecular similarity is fingerprinting, which entails extracting chemical properties of molecules and storing them as a bitstring. These bitstrings can easily be compared computationally.

Initially, we will add an second column to the SMILES compound library containing a label for each molecule, and concatenate (join together) the library file with the original SMILES file for the ligand from the PDB file.

> ### {% icon hands_on %} Hands-on: Calculate molecular fingerprints
> 1. **Add column** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Add this value"*: `mol`.
>    - {% icon param-file %} *"to Dataset"*: 'Compound library' file.
>    - {% icon param-file %} *"Iterate"*: `Yes`
> 2. **Concatenate datasets** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: Output of previous step, 'Ligand SMILES'.
>    - {% icon param-file %} *"to Dataset"*: 'Compound library' file.
>    - {% icon param-file %} *"Iterate"*: `Yes`
>    - Rename dataset 'Labelled compound library'
> 3. **Molecules to Fingerprints** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"molecule file"*: 'Labelled compound library' file.
>    - {% icon param-file %} *"Type of fingerprint"*: `Open Babel FP2 fingerprints`
>    - Rename to 'Fingerprints'.
> 2. **Taylor-Butina clustering** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fingerprint dataset"*: 'Fingerprints' file.
>    - {% icon param-file %} *"threshold"*: `0.8`
> 3. **NxN Clustering** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Fingerprint dataset"*: 'Ligand (MOL)' file.
>    - {% icon param-file %} *"threshold"*: `0.0`
>    - {% icon param-file %} *"Format of the resulting picture"*: `SVG`
>    - {% icon param-file %} *"Output options"*: `Image`
{: .hands_on}

Taylor-Butina clustering provides a classification of the compounds into different groups or clusters, based on their structural similarity. We can therefore see how similar the compounds are to the original ligand, and after docking, we can compare the results to the proposed clusters to observe if there is any correlation.

![Fingerprinting]({{ site.baseurl }}{% link topics/computational-chemistry/images/fingerprints.png %} "A simple fingerprinting system. Each 1 or 0 in the bitstring corresponds to the presence or absence of a particular feature in the molecule. In this case the presence of phenyl, amine and carboxylic acid groups are encoded.")

The image produced by the NxN clustering shows the clustering in the form of a dendrogram.

![NxN clustering]({{ site.baseurl }}{% link topics/computational-chemistry/images/nxn.png %} "Dendrogram produced by NxN clustering. The library used to produce this image is generated with a Tanimoto cutoff of 70, hence 51 search results are shown, plus the original ligands contained in the PDB file.")

> ### {% icon details %} Further investigation (optional)
>
> * Try generating fingerprints using some of the other nine different protocols available and monitor how this affects the clustering.
>
> * For both the Taylor-Butina and NxN clustering, a threshold has to be set. Try varying this value and observe how the clustering results are affected.
>
{: .details}

# Prepare files for docking

A processing step now needs to be applied to the protein structure and the docking candidates - each of the structures needs to be converted to PDBQT format prior to application of the AutoDock Vina docking tool.

In addition, docking requires the coordinates of a binding site to be defined. Effectively, this defines a cuboidal volume in which the docking software attempts to define an optimal binding site. In this case, we already know the location of the binding site, since the downloaded PDB structure contained a bound ligand. There is a tool in Galaxy which can be used to automatically create a configuration file for docking when ligand coordinates are already known.

> ### {% icon details %} How to find the binding site of an apoprotein?
>
> If the structure contains no ligand in complex with the protein (i.e. apoprotein), identifying the binding site is not so trivial as in this example. Fortunately, software is available for automatic detection of pockets which may be promising candidates for ligand binding sites. If you want, check out the fpocket tool and try running it on the Hsp90 structure.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Generate PDBQT and config files for docking
>
> 1. **Compound conversion** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: 'Protein' file.
>    - {% icon param-file %} *"Output format"*: `PDBQT input format`.
>    - Leave all other options unchanged.
>    - Rename to 'Protein PDBQT'.
> 2. **Compound conversion** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: 'Compound library' file.
>    - {% icon param-file %} *"Output format"*: `PDBQT`
>    - {% icon param-file %} *"Split multi-molecule files into a collection"*: `Yes`
>    - Leave all other options unchanged.
> 3. **Calculate the box parameters for an AutoDock Vina job** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input ligand"*: 'Ligand (MOL) file.
>    - {% icon param-file %} *"x-axis buffer"*: `5`
>    - {% icon param-file %} *"y-axis buffer"*: `5`
>    - {% icon param-file %} *"z-axis buffer"*: `5`
>    - {% icon param-file %} *"Exhaustiveness"*: `1`
>    - The 'Random seed' parameter should be left empty.
>    - Rename to 'Docking config file'
{: .hands_on}

# Docking

Finally, the docking itself can be performed.

> ### {% icon hands_on %} Hands-on: Perform docking
>
> 1. **Docking** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Receptor"*: 'Protein PDBQT' file.
>    - {% icon param-file %} *"Ligand"*: 'Prepared ligands' collection.
>    - {% icon param-file %} *"Specify parameters"*: 'Docking config file'
>    - {% icon param-file %} *"Exhaustiveness"*: leave blank (it was specified in the previous step)
>    - {% icon param-file %} *"Output format"*: `PDBQT (and separate file with binding scores)`
{: .hands_on}

The output consists of two collections, containing respectively structural files (PDBQT format) and scoring files for each of the ligands.

View one of the scoring files (click on the eye icon in the history pane). You will see something that looks like this:

```
#################################################################
# If you used AutoDock Vina in your work, please cite:          #
#                                                               #
# O. Trott, A. J. Olson,                                        #
# AutoDock Vina: improving the speed and accuracy of docking    #
# with a new scoring function, efficient optimization and       #
# multithreading, Journal of Computational Chemistry 31 (2010)  #
# 455-461                                                       #
#                                                               #
# DOI 10.1002/jcc.21334                                         #
#                                                               #
# Please see http://vina.scripps.edu for more information.      #
#################################################################

Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 1
Performing search ... done.
Refining results ... done.

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1         -9.2      0.000      0.000
   2         -8.9      2.290      5.660
   3         -7.9      2.782      5.767
   4         -7.6      3.125      7.494
   5         -7.5      2.628      7.444
   6         -7.2      3.453      7.822
   7         -6.9      2.625      7.542
   8         -6.8      3.596      6.902
   9         -6.7      3.293      8.765
  10         -6.0      3.263      7.244
  11         -5.4      3.513      7.115
  12         -4.7      3.316      5.864
  13         -4.7      3.125      6.982
  14         -3.1      3.911      8.546
  15         37.0      4.403      6.901
Writing output ... done.

```

In the table, fifteen different binding modes are listed, from most to least energetically favorable. The second column shows the binding affinity in kcal/mol. The first row shows the most favorable binding mode. Thus we can see that the optimal binding energy of the ligand is -9.2 kcal/mol (or in SI units, -38.5 kJ/mol).

Compare the optimal binding energies for several ligands, selecting examples from different NxN clusters.