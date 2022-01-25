---
layout: tutorial_hands_on

title: Protein-ligand docking
level: Intermediate
zenodo_link: ''
questions:
- What is cheminformatics?
- What is protein-ligand docking?
- How can I perform a simple docking workflow in Galaxy?
objectives:
- Create a small compound library using the ChEMBL database
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

Cheminformatics is the use of computational techniques and information about molecules to solve problems in chemistry. This involves a number of steps: retrieving data on chemical compounds, sorting data for properties which are of interest, and extracting new information. This tutorial will provide a brief overview of all of these, centered around protein-ligand docking, a molecular modelling technique. The purpose of protein-ligand docking is to find the optimal binding between a small molecule (ligand) and a protein. It is generally applied to the drug discovery and development process with the aim of finding a potential drug candidate. First, a target protein is identified. This protein is usually linked to a disease and is known to bind small molecules. Second, a 'library' of possible ligands is assembled. Ligands are small molecules that bind to a protein and may interfere with protein function. Each of the compounds in the library is then 'docked' into the protein to find the optimal binding position and energy.

Docking is a form of molecular modelling, but several simplifications are made in comparison to methods such as molecular dynamics. Most significantly, the receptor is generally considered to be rigid, with covalent bond lengths and angles held constant. Charges and protonation states are also not permitted to change. While these approximations reduce accuracy to some extent, they increase computational speed, which is necessary to screen a large compound library in a realistic amount of time.

In this tutorial, you will perform docking of ligands into the N-terminus of Hsp90 (heat shock protein 90). The tools used for docking are based on the open-source software [AutoDock Vina](http://vina.scripps.edu/) ({% cite Trott2009 %}).

> ### {% icon details %} Biological background
>
> The 90 kDa heat shock protein (Hsp90) is a chaperone protein responsible for catalyzing the conversion of a wide variety of proteins to a functional form; examples of the Hsp90 clientele, which totals several hundred proteins, include nuclear steroid hormone receptors and protein kinases. The mechanism by which Hsp90 acts varies between clients, as does the client binding site; the process is dependent on post-translational modifications of Hsp90 and the identity of co-chaperones which bind and regulate the conformational cycle. 
>
>
> Due to its vital biochemical role as a chaperone protein involved in facilitating the folding of many client proteins, Hsp90 is an attractive pharmaceutical target. In particular, as protein folding is a potential bottleneck to slow cellular reproduction and growth, blocking Hsp90 function using inhibitors which bind tightly to the ATP binding site could assist in treating cancer; for example, the antibiotic geldanamycin and its analogs are under investigation as possible anti-tumor agents. 
>
> See {% cite Pearl2006 %} for a review on the structure and mechanism of the Hsp90 protein. Alternatively, read more at [PDB-101](https://pdb101.rcsb.org/motm/108) and [Wikipedia](https://en.wikipedia.org/wiki/Hsp90).
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

For this exercise, we need two datasets: a protein structure and a library of compounds. We will download the former directly from the Protein Data Bank; the latter will be created by searching the ChEMBL database ({% cite Gaulton2016 %}).

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Search Galaxy for the {% tool [Get PDB](toolshed.g2.bx.psu.edu/repos/bgruening/get_pdb/get_pdb/0.1.0) %} tool. Request the accession code ```2brc```.
> 3. Rename the dataset to 'Hsp90 structure'
> 4. Check that the datatype is correct (PDB file).
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

![Hsp90 N-terminus structure]({% link topics/computational-chemistry/images/2brc.png %} "Structure of Hsp90 N-terminus, as recorded on the PDB. Visualization produced using VMD  ({% cite Humphrey1996 %}).")

# Separating protein and ligand structures

You can view the contents of the downloaded PDB file by pressing the 'View data' icon in the history pane. After the header section (about 500 lines), the atoms of the protein and their coordinates are listed. The lines begin with ```ATOM```. At the end of the file, the atomic coordinates of the ligand and the solvent water molecules are also listed, labelled ```HETATM```. We will use the grep tool to separate these molecules into separate files, and then convert the ligand file into SDF/MOL format using the 'Compound conversion' tool, which is based on OpenBabel, an open-source library for analyzing chemical data ({% cite OBoyle2011 %}).

> ### {% icon details %} What is grep?
>
> Grep is a command-line tool for searching text files for lines which match a search query.
>
> There is a Galaxy tool available based on grep, which we will apply to the downloaded PDB tool.
>
{: .details}

> ### {% icon hands_on %} Hands-on: Separate protein and ligand
>
> 1. {% tool [Search in textfiles (grep)](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: Downloaded PDB file 'Hsp90 structure'
>    - {% icon param-file %} *"that"*: `Don't match`
>    - {% icon param-file %} *"Regular Expression"*: `HETATM`
>    - All other parameters can be left as their defaults. 
>    - Rename the dataset **'Protein (PDB)'**.
>
>    The result is a file with all non-protein (`HETATM`) atoms removed.
> 2. {% tool [Search in textfiles (grep)](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters. Here, we use grep again to produce a file with only non-protein atoms.
>    - {% icon param-file %} *"Select lines from"*: Downloaded PDB file 'Hsp90 structure' 
>    - {% icon param-file %} *"that"*: `Match`
>    - {% icon param-file %} *"Regular Expression"*: `CT5` (the name of the ligand in the PDB file)
>    - All other parameters can be left as their defaults. 
>    - Rename the dataset **'Ligand (PDB)'**.
>
>    This produces a file which only contains ligand atoms.
> 3. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: Ligand PDB file created in step 2.
>    - {% icon param-file %} *"Output format"*: `MDL MOL format (sdf, mol)`
>    - {% icon param-file %} *"Add hydrogens appropriate for pH"*: `7.4`
>    - All other parameters can be left as their defaults. 
>    - Change the datatype to 'mol' and rename the dataset **'Ligand (MOL)'**.
>
>    Applying this tool will generate a representation of the structure of the ligand in MOL format.
>
{: .hands_on}

At this stage, separate protein and ligand files have been created. Next, we want to generate a compound library we can use for docking.


## Creating and processing the compound library

In this step we will create a compound library, using data from the ChEMBL database.

> ### {% icon details %} What chemical databases are available?
>
> Multiple databases are available online which provide access to chemical information, e.g. chemical structures, reactions, or literature. In this tutorial, we use a tool which searches the [ChEMBL](https://www.ebi.ac.uk/chembl/) database. There are also Galaxy tools available for downloading data from [PubChem](https://pubchem.ncbi.nlm.nih.gov/). 
>
> Other notable databases include [Drugbank](https://www.drugbank.ca/) and [ChemSpider](http://www.chemspider.com/).
>
{: .details}

We will generate our compound library by searching ChEMBL for compounds which have a similar structure to the ligand in the PDB file we downloaded in the first step. There is a Galaxy tool for accessing ChEMBL which requires data input in SMILES format; thus, the first step is to convert the 'Ligand' PDB file to a SMILES file. Then the search is performed, returning a SMILES file. For docking, we would like to convert to SDF format, which we can do once again using the 'Compound conversion' tool.

> ### {% icon hands_on %} Hands-on: Generate compound library
>
> 1. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: 'Ligand' PDB file
>    - {% icon param-file %} *"Output format"*: `SMILES format (SMI)`
>    - Leave all other options as default.
> 2. Rename the output of the 'compound conversion' step to **'Ligand SMILES'**.
> 3. {% tool [Search ChEMBL database](toolshed.g2.bx.psu.edu/repos/bgruening/chembl/chembl/0.10.1+galaxy4) %} with the following parameters:
>    - {% icon param-file %} *"SMILES input type"*: File
>    - {% icon param-file %} *"Input file"*: 'Ligand SMILES' file
>    - {% icon param-file %} *"Search type"*: `Similarity`
>    - {% icon param-file %} *"Tanimoto cutoff score"*: `40`
>    - {% icon param-file %} *"Filter for Lipinski's Rule of Five"*: `Yes`
>    - All other parameters can be left as their defaults.
>
>    > ### {% icon question %} Question
>    >
>    > Why are compounds filtered for Lipinski's Rule of Five?
>    >
>    > > ### {% icon solution %} Solution
>    > > Lipinski's rule of five is an empirical rule which can be used to determine the 'druglikeness' of a molecule. The rule consists of four criteria which relate to the pharmacokinetics of the molecule. The rule is discussed in ({% cite Lipinski2004 %}), or you can also read more on [Wikipedia](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five).
>    > {: .solution}
>    {: .question}
>
>    Optional: experiment with different combinations of options - adding different filters, adjusting the Tanimoto coefficient.
> 4. Check that the datatype is correct (smi). This step is essential, as Galaxy does not automatically detect the datatype for SMILES files.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
> 5. Rename dataset **'Compound library'**.
>
{: .hands_on}

> ### {% icon tip %} Tip: Problems using the ChEMBL tool?
> A number of users encounter issues with the ChEMBL tool - sometimes the tool fails, or the output is returned successfully but is empty. If this happens to you, try the following:
> > Rerun the tool - if a transient error on the ChEMBL server was at fault, this might be enough to fix it.
> > Try modifying some of the parameters. For example, reducing the Tanimoto coefficient should increase the number of compounds returned.
> > If all else fails, you can use the following list of SMILES:
```
Cc1n[nH]c(c2ccc(O)c(Cl)c2)c1c3ccc4OCCOc4c3	CHEMBL187670
COc1ccc(cc1)c2c(C)[nH]nc2c3ccc(O)cc3O	CHEMBL192894
COc1ccc(c(O)c1)c2onc(C)c2c3ccc4OCCOc4c3	CHEMBL1541585
CCOc1ccc(c(O)c1)c2[nH]nc(C)c2c3ccccc3OC	CHEMBL1504505
CN(CCc1c(C)n[nH]c1C)Cc2cn(C)nc2c3ccc4OCCOc4c3	CHEMBL1560480
COc1ccc(c(O)c1)c2[nH]nc(C)c2c3ccc4OCCCOc4c3	CHEMBL362893
CCCc1c(OCCCOc2cc(O)c(cc2CC)c3cc[nH]n3)ccc4CCC(Oc14)C(=O)O	CHEMBL81401
Cc1cccc(n1)c2[nH]nc(C)c2c3ccnc4ccccc34	CHEMBL129153
COc1ccc(cc1OC)c2c(C)n[nH]c2c3ccc(O)cc3O	CHEMBL1595327
COc1ccc(c(O)c1)c2noc(C)c2c3ccc4ccccc4n3	CHEMBL1486235
CCc1cc(c(O)cc1O)c2[nH]nc(C)c2c3ccc4OCCOc4c3	CHEMBL399530
Cc1[nH]nc(c2cc(Cl)c(O)cc2O)c1c3ccc4OCCOc4c3	CHEMBL191074
COc1ccc(c(O)c1)c2[nH]ncc2c3ccc4OCCOc4c3	CHEMBL1415374
Cc1n[nH]c(c2cc(Cl)ccc2O)c1c3ccc4OCCOc4c3	CHEMBL187678
Cc1noc(c2ccc(O)cc2O)c1c3ccc4OCCOc4c3	CHEMBL582320
CCOC(=O)c1oc(cc1)c2c(C)[nH]nc2c3cc(CC)c(O)cc3O	CHEMBL3932805
NC(=O)c1ccc2[nH]nc(c3ccc4OCCOc4c3)c2c1	CHEMBL3900406
COc1ccc(cc1OC)c2cc([nH]n2)c3c(O)c(OC)c4occc4c3OC	CHEMBL1351838
CCCc1cc(c(O)cc1OC)c2[nH]ncc2c3ccc4OCCCOc4c3	CHEMBL1443258
Oc1cc(O)c(cc1Cl)c2[nH]ncc2c3ccc4OCCOc4c3	CHEMBL191228
CCc1cc(c(O)cc1O)c2n[nH]cc2c3ccc4OCCOc4c3	CHEMBL3187010
Oc1ccccc1c2cc([nH]n2)c3ccc4OCCOc4c3	CHEMBL1567097
CCCc1cc(c(O)cc1O)c2[nH]ncc2c3ccc4OCCCOc4c3	CHEMBL1578064
COc1ccc(c(O)c1)c2[nH]nc(C)c2c3ccc(OC)c(OC)c3	CHEMBL1335688
COc1ccc(cc1)c2c(N)n[nH]c2c3cc(OC)c4OCCOc4c3	CHEMBL2408971
CCC(C)c1cc(c(O)cc1O)c2[nH]ncc2c3ccc4OCOc4c3	CHEMBL187674
Cc1noc(c2ccc(O)cc2O)c1c3ccc4OCCCOc4c3	CHEMBL587334
OCc1cn(nc1c2ccc3OCCOc3c2)c4ccccc4	CHEMBL1549407
Oc1ccc(c(O)c1)c2[nH]ncc2c3cccc4cccnc34	CHEMBL1305951
Cc1n[nH]c(c2ccc(O)cc2O)c1c3ccc4OCCOc4c3	CHEMBL188965
Cc1[nH]nc(c2ccc3OCC(=O)Nc3c2)c1c4ccc(F)cc4	CHEMBL3337723
CCOc1ccc(c(O)c1)c2n[nH]c(C)c2c3ccc(OC)cc3	CHEMBL1698243
COc1ccc(cc1)c2cc([nH]n2)c3c(O)c(OC)c4occc4c3OC	CHEMBL1402615
CCc1cc(c(O)cc1O)c2[nH]ncc2c3ccc4OCOc4c3	CHEMBL1412538
CCOc1cc(O)c(cc1CC)c2nc(N)ncc2c3ccc4OCCOc4c3	CHEMBL547662
COc1ccc(cc1)c2c(N)onc2c3cc(OC)c4OCCOc4c3	CHEMBL3113121
CCCc1cc(c(O)cc1O)c2n[nH]cc2c3ccc4OCCOc4c3	CHEMBL3956397
CCC(C)c1cc(c(O)cc1O)c2[nH]nc(C)c2c3ccc4OCCOc4c3	CHEMBL190919
CCC(C)c1cc(c(O)cc1OC)c2[nH]ncc2c3ccc4OCCOc4c3	CHEMBL435501
Cn1cc(CNCc2c[nH]nc2c3ccc(F)cc3)c(n1)c4ccc5OCCOc5c4	CHEMBL1537178
Oc1ccc(F)cc1c2cc([nH]n2)c3ccc4OCCOc4c3	CHEMBL1451528
CCc1cc(c(O)cc1OCC(=O)O)c2n[nH]c(C)c2c3ccc4OCCCOc4c3	CHEMBL3952001
Cc1[nH]nc(c2ccc(O)c(O)c2O)c1c3ccc(Cl)cc3	CHEMBL1092945
COc1cc(cc(OC)c1OC)c2n[nH]nc2c3ccc4OCCOc4c3	CHEMBL3740841
Oc1ccc(c(O)c1)c2n[nH]cc2c3ccc4OCOc4c3	CHEMBL577176
```
> Don't worry if you can't get it to work - successfully generating this list is a very minor part of the tutorial!
{: .tip}


There are some other tools available, which will not be used in this tutorial, which help to develop a more focused compound library. For example, the 'Natural product likeness calculator' and 'Drug-likeness' tools assign a score to compounds based on how similar they are to typical natural products and drugs respectively, which could then be used to filter the library. If you are interested, you can try testing them out on the library just generated.

> ### {% icon tip %} Tip: Generating a compound library
> If you try using this tutorial using your own data, you might encounter some issues. Important things to remember:
> * If you encounter an error, check the SMILES file only has a single column. Additional columns can be removed using the 'Cut' tool.
> * If the output file is empty, it may be that the ChEMBL database doesn't have any compounds similar to the input. Consider lowering the Tanimoto coefficient to 70 if this is the case and removing filters (including the Lipinski RO5 filter). If this doesn't help, you will have to use another source of chemical data (e.g. PubChem).
> * Finally, please remember this step is totally optional if you already have a list of compounds for docking (in SMILES or another format). In this case you can upload them to Galaxy and continue with the next step.
{: .tip}

> ### {% icon details %} What are SMILES and SDF formats?
>
> SMILES and SD-files both represent chemical structures. A SMILES file represents the 2D structure of a molecule as a chemical graph. In other words, it states only the atoms and the connectivity between them. An example of a SMILES string (taken from the ligand in the PDB file) is `c1c2OCCOc2ccc1c1c(C)[nH]nc1c1cc(CC)c(O)cc1O`. For more information on how the notation works, please consult the [OpenSMILES specification](http://opensmiles.org/opensmiles.html) or the description provided by [Wikipedia](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system). A more comprehensive alternative to the SMILES system is the International Chemical Identifier (InChI).
>
> Neither SMILES nor InChI format contain the three-dimensional structure of a molecule. By contrast, the SDF (structure data file) format encodes three-dimensional atomic coordinates of a structure, similar to a PDB file.
>
> In a previous step, we also generated a MOL file - this format is closely related to the SDF format. The difference is that MOL files can store only a single molecule, whereas SD-files can encode single or multiple molecules. Multiple molecules are separated by lines containing four dollar signs (`$$$$`).
>
> For docking, we need the three-dimensional coordinates of the ligand; thus, we want to convert from SMILES to SDF format.
{: .details}

# Prepare files for docking

A processing step now needs to be applied to the protein structure and the docking candidates - each of the structures needs to be converted to PDBQT format before using the AutoDock Vina docking tool.

Further, docking requires the coordinates of a binding site to be defined. Effectively, this defines a 'box' in which the docking software attempts to define an optimal binding site. In this case, we already know the location of the binding site, since the downloaded PDB structure contained a bound ligand. There is a tool in Galaxy which can be used to automatically create a configuration file for docking when ligand coordinates are already known.

> ### {% icon hands_on %} Hands-on: Generate PDBQT and config files for docking
>
> 1. {% tool [Prepare receptor](toolshed.g2.bx.psu.edu/repos/bgruening/autodock_vina_prepare_receptor/prepare_receptor/1.5.7+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Select a PDB file"*: 'Protein' PDB file.
> 2. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: 'Compound library' file.
>    - {% icon param-file %} *"Output format"*: `SDF`
>    - {% icon param-file %} *"Generate 3D coordinates"*: `Yes`
>    - {% icon param-file %} *"Add hydrogens appropriate for pH"*: `7.4`
>    - Leave all other options unchanged.
>    - Rename to 'Prepared ligands'
> 3. {% tool [Calculate the box parameters for an AutoDock Vina job](toolshed.g2.bx.psu.edu/repos/bgruening/autodock_vina_prepare_box/prepare_box/2021.03.4+galaxy0) %}  with the following parameters:
>    - {% icon param-file %} *"Input ligand or pocket"*: `Ligand (MOL)` file.
>    - {% icon param-file %} *"x-axis buffer"*: `5`
>    - {% icon param-file %} *"y-axis buffer"*: `5`
>    - {% icon param-file %} *"z-axis buffer"*: `5`
>    - {% icon param-file %} *"Random seed"*: `1`
>    - Rename to 'Docking config file'.
{: .hands_on}

Perhaps you are interested in a system which does not have a ligand within the binding site (an apoprotein). In this case you need to run the fpocket tool to identify potential binding sites in the protein structure - take a look at the following (optional) section.


> ### {% icon details %} How to find the binding site of an apoprotein?
>
> If the structure contains no ligand in complex with the protein (i.e. apoprotein), there is an additional step of identifying the binding site. Software is available for automatic detection of pockets which may be promising candidates for ligand binding sites. For example, let's try out the fpocket tool ({% cite LeGuilloux2009 %}) on the Hsp90 structure, imagining we don't have access to the file containing the ligand coordinates.
> > ### {% icon hands_on %} Hands-on: Finding binding sites using fpocket
> >
> > 1. {% tool [fpocket](toolshed.g2.bx.psu.edu/repos/bgruening/fpocket/fpocket/3.1.4.2+galaxy0) %} with the following parameters:
> >   - {% icon param-file %} *"Input file"*: `Protein (PDB)` file.
> >   - {% icon param-file %} *"Type of pocket to detect"*: `Small molecule binding sites`
> >   - {% icon param-file %} *"Output files"*: select `PDB files containing the atoms in contact with each pocket`, `Log file containing pocket properties`.
> >
> > 2. {% tool [Calculate the box parameters for an AutoDock Vina job](toolshed.g2.bx.psu.edu/repos/bgruening/autodock_vina_prepare_box/prepare_box/2021.03.4+galaxy0) %} with the following parameters:
> >    - {% icon param-file %} *"Input ligand or pocket"*: `pocket2` PDB file from the `Atoms in contact with each pocket` collection.
> >    - {% icon param-file %} *"x-axis buffer"*: `5`
> >    - {% icon param-file %} *"y-axis buffer"*: `5`
> >    - {% icon param-file %} *"z-axis buffer"*: `5`
> >    - {% icon param-file %} *"Exhaustiveness (optional)"*: `1`
> >    - {% icon param-file %} *"Random seed"*: `1`
> >    - Rename to 'Docking config file derived from pocket'.
> >
> >  The fpocket tool generates two different outputs: a `Pocket properties` log file containing details of all the pockets which fpocket found in the protein. The second output is a collection (a list) containing one PDB file for each of the pockets. Each of the PDB files contains only the atoms in contact with that particular pocket. Note that fpocket assigns a score to each pocket, but you should not assume that the top scoring one is the only one where compounds can bind! For example, the pocket where the ligand in the `2brc` PDB file binds is ranked as the second-best according to fpocket.
> {: .hands_on}
>
> You can compare the config file generated with the one generated from the ligand directly - if you picked the right pocket (`pocket2`) the coordinates should be pretty similar.
>
{: .details}


# Docking

Now that the protein and the ligand library have been correctly prepared and formatted, docking can be performed.

> ### {% icon hands_on %} Hands-on: Perform docking
>
> 1. {% tool [Docking](toolshed.g2.bx.psu.edu/repos/bgruening/autodock_vina/docking/1.1.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Receptor"*: 'Protein PDBQT' file.
>    - {% icon param-file %} *"Ligands"*: 'Prepared ligands' file. 
>    - {% icon param-file %} *"Specify pH value for ligand protonation"*: `7.4` 
>    - {% icon param-file %} *"Specify parameters"*: 'Upload a config file to specify parameters' 
>    - {% icon param-file %} *"Box configuration"*: 'Docking config file'
>    - {% icon param-file %} *"Exhaustiveness"*: leave blank (it was specified in the previous step)
{: .hands_on}

The output consists of a collection, which contains an SDF output file for each ligand, containing multiple docking poses and scoring files for each of the ligands. We will now perform some processing on these files which extracts scores from the SD-files and selects the best score for each. 

# Optional: cheminformatics tools applied to the compound library

The ChemicalToolbox contains a large number of cheminformatics tools. This section will demonstrate some of the useful functionalities available. If you are just interested in docking, feel free to skip this section - or, just try out the tools which look particularly interesting.

(This section can also be completed while waiting for the docking, which can take some time to complete.)

### Visualization

It can be useful to visualize the compounds generated. There is a tool available for this in Galaxy based on OpenBabel.

> ### {% icon hands_on %} Hands-on: Visualization of chemical structures
>
> 1. {% tool [Visualisation](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_svg_depiction/openbabel_svg_depiction/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: Compound library
>    - {% icon param-file %} *"Embed molecule as CML"*: `No`
>    - {% icon param-file %} *"Draw all carbon atoms"*: `No`
>    - {% icon param-file %} *"Use thicker lines"*: `No`
>    - {% icon param-file %} *"Property to display under the molecule"*: `Molecule title`
>    - {% icon param-file %} *"Sort the displayed molecules by"*: `Molecular weight`
>    - {% icon param-file %} *"Format of the resultant picture"*: `SVG`
>
{: .hands_on}

This produces an SVG image of all the structures generated ordered by molecular weight.

![Image showing structures of compounds from ChEMBL]({% link topics/computational-chemistry/images/compound_library.png %} "Structures of the compounds from ChEMBL.")

### Calculation of fingerprints and clustering

In this step, we will group similar molecules together. A key tool in cheminformatics for measuring molecular similarity is fingerprinting, which entails extracting chemical properties of molecules and storing them as a bitstring. These bitstrings can be easily compared computationally, for example with a clustering method. The fingerprinting tools in Galaxy are based on the Chemfp tools ({% cite Dalke2013 %}).

Before clustering, let's label each compound. To do so add a second column to the SMILES compound library containing a label for each molecule. The ```Ligand SMILES``` file is also labelled something like ```/data/dnb02/galaxy_db/files/010/406/dataset_10406067.dat``` (the exact name will vary) and we would like to give it a more useful name. When labelling is complete, we can concatenate (join together) the library file with the original SMILES file for the ligand from the PDB file.

> ### {% icon hands_on %} Hands-on: Calculate molecular fingerprints
> 1. {% tool [Replace](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_find_and_replace/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Ligand SMILES`.
>    - {% icon param-file %} *"Find pattern"*: add the current label of the SMILES here. You can find it by clicking the 'view' button next to the `Ligand SMILES` dataset - it will look something like `/data/dnb02/galaxy_db/files/010/406/dataset_10406067.dat`.
>    - {% icon param-file %} *"Replace with"*: `ligand`
> 2. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: Output of the previous step.
>    - Click on **Insert Dataset** and in the new selection box which appears, select 'Compound library'.
>    - Run the step and rename the output dataset 'Labelled compound library'. 
> 3. {% tool [Molecule to fingerprint](toolshed.g2.bx.psu.edu/repos/bgruening/chemfp/ctb_chemfp_mol2fps/1.5) %} with the following parameters:
>    - {% icon param-file %} *"Molecule file"*: 'Labelled compound library' file.
>    - {% icon param-file %} *"Type of fingerprint"*: `Open Babel FP2 fingerprints`
>    - Rename to 'Fingerprints'.
{: .hands_on}

Taylor-Butina clustering  ({% cite Butina1999 %}) provides a classification of the compounds into different groups or clusters, based on their structural similarity. This methods shows us how similar the compounds are to the original ligand, and after docking, we can compare the results to the proposed clusters to observe if there is any correlation.

![Image showing a Fingerprinting System]({% link topics/computational-chemistry/images/fingerprints.png %} "A simple fingerprinting system. Each 1 or 0 in the bitstring corresponds to the presence or absence of a particular feature in the molecule. In this case, the presence of phenyl, amine and carboxylic acid groups are encoded.")

> ### {% icon hands_on %} Hands-on: Cluster molecules using molecular fingerprints
> 1. {% tool [Taylor-Butina clustering](toolshed.g2.bx.psu.edu/repos/bgruening/chemfp/ctb_chemfp_butina_clustering/1.5) %} with the following parameters:
>    - {% icon param-file %} *"Fingerprint dataset"*: 'Fingerprints' file.
>    - {% icon param-file %} *"threshold"*: `0.8`
> 2. {% tool [NxN clustering](toolshed.g2.bx.psu.edu/repos/bgruening/chemfp/ctb_chemfp_nxn_clustering/1.5.1) %} with the following parameters:
>    - {% icon param-file %} *"Fingerprint dataset"*: 'Fingerprints' file.
>    - {% icon param-file %} *"threshold"*: `0.0`
>    - {% icon param-file %} *"Format of the resulting picture"*: `SVG`
>    - {% icon param-file %} *"Output options"*: `Image`
{: .hands_on}

The image produced by the NxN clustering shows the clustering in the form of a dendrogram, where individual molecules are represented as vertical lines and merged into clusters. Merges are represented by horizontal lines. The y-axis represents the similarity of data points to each other; thus, the lower a cluster is merged, the more similar the data points are which it contains. Clusters in the dendogram are colored differently. For example, all molecules connected in red are similar enough to be grouped into the same cluster. 

![NxN clustering]({% link topics/computational-chemistry/images/nxn.png %} "Dendrogram produced by NxN clustering. The library used to produce this image is generated with a Tanimoto cutoff of 80; here 15 search results are shown, plus the original ligand contained in the PDB file.")

> ### {% icon details %} Further investigation (optional)
>
> * Try generating fingerprints using some of the other nine different protocols available and monitor how this affects the clustering.
>
> * For both the Taylor-Butina and NxN clustering, a threshold has to be set. Try varying this value and observe how the clustering results are affected.
>
{: .details}

# Post-processing and plotting

From our collection of SD-files, we first extract all stored values into tabular format and then combine the files together to create a single tabular file.

> ### {% icon hands_on %} Hands-on: Process SD-files
>
> 1. {% tool [Extract values from an SD-file](toolshed.g2.bx.psu.edu/repos/bgruening/sdf_to_tab/sdf_to_tab/2020.03.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input SD-file"*: Collection of SD-files generated by the docking step. (Remember to select the 'collection' icon!)
>    - {% icon param-file %} *"Include the property name as header"*: `Yes` 
>    - {% icon param-file %} *"Include SMILES as column in output"*: `Yes`
>    - {% icon param-file %} *"Include molecule name as column in output"*: `Yes` 
>    - Leave all other paramters unchanged.
> 2. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/4.2) %} with the following parameters:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: Collection of tabular files generated by the previous step.
>    - {% icon param-file %} *"Keep one header line"*: `Yes` 
>    - {% icon param-file %} *"Append File name"*: `No` 
>
>    {% snippet faqs/galaxy/tools_select_collection.md datatype="datatypes" %}
>
> 3. {% tool [Compound conversion](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_compound_convert/openbabel_compound_convert/3.1.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: choose one of the SD-files from the collection generated by the docking step.
>    - {% icon param-file %} *"Output format"*: `Protein Data Bank format (pdb)` 
>    - {% icon param-file %} *"Split multi-molecule files into a collection"*: `Yes`
>    - Leave all other parameters unchanged. 
{: .hands_on}

We now have a tabular file available which contains all poses calculated for all ligands docked, together with scores and RMSD values for the deviation of each pose from the optimum. We also have PDB files for some of the docking poses which can be inspected using the NGLViewer visualization embedded in Galaxy.

> ### {% icon details %} Further investigation (optional)
>
> There are visualizations available in Galaxy for producing various types of plots - for example, a scatter plot of RMSD (compared to optimal docking pose) against docking score can be calculated very easily:
> 
> ![Scatter plot of RMSD against docking score]({% link topics/computational-chemistry/images/scatterplot.png %} "Scatter plot of RMSD against docking score.")
> 
> The plot shows a correlation between the two variables, as expected, though only a slight one, given the narrow range of docking scores in the dataset and the structural similarity of the ligands tested.
> 
> A more advanced exercise would be to generate molecular descriptors for the molecule (there are three tools available for this in Galaxy: RDKit, Mordred and PaDEL) and to identify those which are particularly predictive of docking score, either by inspecting plots or using the statistics tools built into Galaxy.
{: .details}

Use the NGLviewer to inspect the protein and various ligand poses generated by docking. This can be done using either the visualization of NGLViewer in Galaxy, or via the [NGL website](http://nglviewer.org/ngl/).

![Two docking poses for a ligand bound to Hsp90]({% link topics/computational-chemistry/images/activesite.gif %} "Two docking poses for a ligand bound to the active site of Hsp90. One (docking score -8.4) can be seen to be bound deeper in the active site than the other (docking score -5.7), which is reflected in the difference between the docking scores.")