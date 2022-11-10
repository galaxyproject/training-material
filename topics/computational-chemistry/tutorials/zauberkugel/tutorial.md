---
layout: tutorial_hands_on

title: Protein target prediction of a bioactive ligand with Align-it and ePharmaLib
level: Intermediate
zenodo_link: 'https://zenodo.org/record/6055897/files/ePharmaLib_PHARAO_plasmodium.phar'
questions:
- What is a pharmacophore model?
- How can I perform protein target prediction with a multi-step workflow or the one-step Zauberkugel workflow?
objectives:
- Create an SMILES file of a bioactive ligand.
- Screen the query ligand against a pharmacophore library.
- Analyze the results of the protein target prediction.
time_estimation: 2H
key_points:
- A pharmacophore is an abstract description of the molecular features of a bioactive ligand.
- Pharmacophore-based target prediction is an efficient and cost-effective method.
contributors:
- aurelienmoumbock
- simonbray

---


# Introduction


Historically, the pharmacophore concept was formulated in 1909 by the German physician and Nobel prize laureate Paul Ehrlich ({% cite Ehrlich1909 %}). According to the [International Union of Pure and Applied Chemistry (IUPAC)](https://iupac.org/), a pharmacophore is defined as “an ensemble of steric and electronic features that is necessary to ensure the optimal supramolecular interactions with a specific biological target and to trigger (or block) its biological response” ({% cite Wermuth1998 %}). Starting from the cocrystal structure of a non-covalent protein–ligand complex (e.g. Figure 1), pharmacophore perception involves the extraction of the key molecular features of the bioactive ligand at the protein–ligand contact interface into a single model ({% cite Moumbock2019 %}). These pharmacophoric features mainly include: H-bond acceptor (HACC or A), H-bond donor (HDON or D), lipophilic group (LIPO or H), negative center (NEGC or N), positive center (POSC or P), and aromatic ring (AROM or R) moieties. Moreover, receptor-based excluded spheres (EXCL) can be added in order to mimic spatial constraints of the binding pocket (Figure 2). Once a pharmacophore model has been generated, a query can be performed either in a forward manner, using several ligands to search for novel putative hits of a given target, or in a reverse manner, by screening a single ligand against multiple pharmacophore models in search of putative protein targets ({% cite Steindl2006 %}).

![PDB ID: 4MVF]({% link topics/computational-chemistry/images/4MVF-STU.png %} "Crystal Structure of *Plasmodium falciparum* calcium-dependent protein kinase 2 (CDPK2) complexed with staurosporine (STU) with PDB ID: [4MVF](https://www.rcsb.org/structure/4mvf). Image generated using Maestro (Schrödinger LLC, NY).")

Bioactive compounds often bind to several target proteins, thereby exhibiting polypharmacology. However, experimentally determining these interactions is laborious, and structure-based virtual screening of bioactive compounds could expedite drug discovery by prioritizing hits for experimental validation. The recently reported ePharmaLib ({% cite Moumbock2021 %}) dataset is a library of 15,148 e-pharmacophores modeled from solved structures of pharmaceutically relevant protein–ligand complexes of the screening Protein Data Bank (sc-PDB, {% cite Desaphy2014 %}). ePharmaLib can be used for target fishing of phenotypic hits, side effect predictions, drug repurposing, and scaffold hopping.

![STU]({% link topics/computational-chemistry/images/STU.png %} "Depiction of the 2D structure of staurosporine (left) and 3D structure (right) with key pharmacophoric features extracted from the STU–CDPK2 complex (PDB ID: [4MVF](https://www.rcsb.org/structure/4mvf)). Image generated using Maestro (Schrödinger LLC, NY).")

In this tutorial, you will perform pharmacophore-based target prediction of a bioactive ligand known as staurosporine (Figure 2) with the ePharmaLib subset representing *Plasmodium falciparum* protein targets (138 pharmacophore models) and the open-source pharmacophore alignment program [Align-it](https://anaconda.org/bioconda/align_it), formerly known as PHARAO ({% cite Taminau2008 %}).

> <details-title>Pharmacology of staurosporine</details-title>
>
>  Staurosporine (PDB hetID: [STU](https://www.rcsb.org/ligand/STU)) is an indolocarbazole secondary metabolite isolated from several bacteria of the genus Streptomyces. It displays diverse biological activities such as anticancer and antiparasitic activities ({% cite Nakano2009 %}).
{: .details}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Create a history

As a first step, we create a new history for the analysis.

> <hands-on-title>Hands-on 1: Create history</hands-on-title>
>
> 1. Create a new history.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename it to `Staurosporine target prediction`.
>
>    {% snippet faqs/galaxy/histories_rename.md name = "Staurosporine target prediction" %}
>
{: .hands_on}


# Get data

For this exercise, we need two datasets: the ePharmaLib pharmacophore library (PHAR format) and a query ligand structure file (SMI format). 

## Fetching the ePharmaLib dataset

Firstly, we will retrieve the concatenated ePharmaLib subset representing *P. falciparum* protein targets. 

> <hands-on-title>Hands-on 2: Upload ePharmaLib</hands-on-title>
>
> 1. Upload the dataset from the [Zenodo](https://zenodo.org) link provided to your Galaxy history.
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md link="hhttps://zenodo.org/record/6055897/files/ePharmaLib_PHARAO_plasmodium.phar" %}
>
>    > <comment-title>ePharmaLib versions</comment-title>
>    >
>    > Two versions of the ePharmaLib (PHAR & PHYPO formats) have been created for use with the pharmacophore alignment programs [Align-it](https://anaconda.org/bioconda/align_it) and [Phase](https://www.schrodinger.com/products/phase), respectively. Both versions can be broken down into small datasets. e.g. for human targets. They are freely available at [Zenodo](https://zenodo.org) under the link:
>    >	```
>    >	https://zenodo.org/record/6055897
>    >  ```
>    {: .comment}
>
> 2. Change the datatype from `tabular`to `phar`. This step is essential, as Galaxy does not automatically detect the datatype for PHAR files.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="phar" %}
>
> 3. You can view the contents of the downloaded PHAR file by pressing the **eye icon** (View data) for this dataset. 
>
>    > <details-title>What is a PHAR file?</details-title>
>    >
>    > A PHAR file is essentially a series of lines containing the three-dimensional coordinates of pharmacophoric features and excluded spheres. The first column specifies a feature type (e.g. HACC is a hydrogen bond acceptor). Subsequent columns specify the position of the feature center in a three-dimensional space. Individual pharmacophores are separated by lines containing four dollar signs (`$$$$`). The pharmacophores of the ePharmaLib dataset were labeled according to the following three-component code *PDBID-hetID-UniprotEntryName*.
>    {: .details}
{: .hands_on}

## Creating a query ligand structure file

In this step, we will manually create an SMI file containing the SMILES of staurosporine.

> <details-title>What are SMILES and the SMI file format?</details-title>
>
> The simplified molecular-input line-entry system (SMILES) is a string notation for describing the 2D chemical structure of a compound. It only states the atoms present in a compound and the connectivity between them. As an example, the SMILES string of acetone is `CC(=O)C`. SMILES strings can be imported by most molecule editors and converted into either two-dimensional structural drawings or three-dimensional models of the compounds, and vice versa. For more information on how the notation works, please consult the [OpenSMILES specification](http://opensmiles.org/opensmiles.html) or the description provided by [Wikipedia](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system). 
>
{: .details}

> <hands-on-title>Hands-on 3: Create an SMI file</hands-on-title>
>
>
> 1. Create a new file using the Galaxy upload manager, with the following contents. Make sure to select the datatype (with **Type**) as `smi`. This step is essential, as Galaxy does not automatically detect the datatype for SMI files.
>
>    ```
>    C[C@@]12[C@@H]([C@@H](C[C@@H](O1)N3C4=CC=CC=C4C5=C6C(=C7C8=CC=CC=C8N2C7=C53)CNC6=O)NC)OC	staurosporine
>    ```
>
>    {% snippet faqs/galaxy/datasets_create_new_file.md format="smi" %}
>
> > <tip-title>SMILES generation</tip-title>
> >
> > A SMILES string can automatically be generated from a ligand name or 2D structure with a desktop molecule editor such [ChemDraw®](https://perkinelmerinformatics.com/products/research/chemdraw/) and [Marvin®](https://chemaxon.com/products/marvin), or with web-based molecule editors such as [PubChem Sketcher](https://pubchem.ncbi.nlm.nih.gov//edit3/index.html) and [ChemDraw® JS](https://chemdrawdirect.perkinelmer.cloud/js/sample/index.html). Moreover, the pre-computed SMILES strings of a large number of bioactive compounds can be retrieved from chemical databases such as [PubChem](https://pubchem.ncbi.nlm.nih.gov/). e.g.
> >
> >    ```
> >    https://pubchem.ncbi.nlm.nih.gov/compound/44259#section=Isomeric-SMILES&fullscreen=true
> >    ```
> >
> {: .tip}
>
{: .hands_on}

> <question-title></question-title>
>
> Why do we specifically use a so-called isomeric SMILES string?
>
> > <solution-title></solution-title>
> >
> > Staurosporine is a chiral molecule possessing four chiral centers. The SMILES notation allows the specification of configuration at tetrahedral centers and double bond geometry, by marking atoms with `@` or `@@`. These are structural features that cannot be specified by connectivity alone, and therefore SMILES which encode this information are termed isomeric SMILES. A notable feature of these rules is that they allow rigorous partial specification of chirality.
> >
> {: .solution}
>
{: .question}


# Pre-processing

Prior to pharmacophore alignment, the predominant ionization state(s) of the query ligand as well as its 3D conformers should be generated. Also, the pharmacophore dataset will be split into a collection of individual pharmacophore files.
 
## Ligand hydration

More often than not, the bioactive form of a compound is its predominant form at physiological pH (7.4). In this step, we predict the most probable ionization state(s) of the query ligand at pH 7.4 with the cheminformatics toolkit OpenBabel ({% cite OBoyle2011 %}).

> <hands-on-title>Hands-on 4: Add hydrogen atoms</hands-on-title>
>
> 1. {% tool [Add hydrogen atoms](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_addh/openbabel_addh/3.1.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: `staurosporine.smi` (from Hands-on 3)
>    - *"Add hydrogens to polar atoms only (i.e. not to carbon atoms)"*: `Yes`
>
> 2. Rename the output to `staurosporine_hydrated`.
>
>    {% snippet faqs/galaxy/datasets_rename.md name = "staurosporine_hydrated" %}
>
{: .hands_on}

> <question-title></question-title>
>
> Nitrogen-containing functional groups are known to be basic. Which of them present in staurosporine (Figure 2) do you expect to be protonated at pH 7.4, and which not? And why?
>
> > <solution-title></solution-title>
> >
> > Only the secondary N-methylamino group will be protonated because indoles, much like aromatic amides, are typically not basic.
> >
> {: .solution}
>
{: .question}

## Splitting ePharmaLib into individual pharmacophores

The ePharmaLib subset representing *P. falciparum* protein targets (*ePharmaLib_PHARAO_plasmodium.phar*) is a concatenated file containing 148 individual pharmacophore files. To speed up our analysis, it is preferable to split the dataset into individual files in order to perform several pharmacophore alignments in parallel, using Galaxy's collection functionality.

> <hands-on-title>Hands-on 5: Splitting ePharmaLib</hands-on-title>
>
> 1. {% tool [Split file](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_to_collection/split_file_to_collection/0.5.0) %} with the following parameters:
>    - *"Select the file type to split"*: `Generic`
>        - {% icon param-file %} *"File to split"*: `ePharmaLib_PHARAO_plasmodium.phar` (from Hands-on 2)
>        - *"Method to split files"*: `Specify record separator as regular expression`
>            - *"Regex to match record separator"*: `\$\$\$\$`
>            - *"Split records before or after the separator?"*: `After`
>        - *"Specify number of output files or number of records per file?"*: `Number of records per file ('chunk mode')`
>        - *"Base name for new files in collection"*: `epharmalib`
>        - *"Method to allocate records to new files"*: `Maintain record order`
>
> 2. Rename the output to `ePharmaLib_PLAF_split`.
>
>    {% snippet faqs/galaxy/datasets_rename.md name = "ePharmaLib_PLAF_split" %}
>
{: .hands_on}

## Ligand conformational flexibility

To reduce the calculation time, the Align-it ({% cite Taminau2008 %}) tool performs rigid alignment rather than flexible alignment. Conformational flexibility of the ligand is accounted for by introducing a preliminary step, in which a set of energy-minimized conformers for the query ligand are generated with the RDConf ({% cite rdconf %}) tool (using the RDKit ({% cite landrum2013rdkit %}) toolkit). 

> <hands-on-title>Hands-on 6: Low-energy ligand conformer search</hands-on-title>
>
> 1. {% tool [RDConf: Low-energy ligand conformer search](toolshed.g2.bx.psu.edu/repos/bgruening/rdconf/rdconf/2020.03.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `staurosporine_hydrated` (from Hands-on 4)
>    - *"Maximum number of conformers to generate per molecule"*: `100`
>
> 2. Rename the output to `staurosporine_3D_conformers`.
>
>    {% snippet faqs/galaxy/datasets_rename.md name = "staurosporine_3D_conformers" %}
>
>    > <comment-title>RDConf</comment-title>
>    >
>    > It is recommended to use the default settings, except for the number of conformers which should be changed to 100. As a rule of thumb, a threshold of 100 conformers appropriately represents the conformational flexibility of a compound with less than 10 rotatable bonds. The output SDF (structure data file) format encodes three-dimensional atomic coordinates of each conformer, separated by lines containing four dollar signs (`$$$$`).
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> Have a look at the contents of the created collection `staurosporine_3D_conformers`. Why were less than 100 conformers were generated for staurosporine?
>
> > <solution-title></solution-title>
> >
> > Staurosporine is a fused 8-ring system with only two rotatable bonds, due to its planar aromatic 5-ring indolocarbozole scaffold which confers a high structural rigidity upon the compound, i.e. it exists in relatively few energetically distinct 3D conformations.
> >
> {: .solution}
>
{: .question}


# Pharmacophore alignment

In this step, the ligand conformer dataset (SDF format) is converted on-the-fly to a pharmacophore dataset (PHAR format) and simultaneously aligned to the individual pharmacophores of the ePharmaLib dataset in a batch mode with Align-it ({% cite Taminau2008 %}). The pharmacophoric alignments and thus the predicted targets are ranked in terms of a scoring metric: `Tversky index` = [0,1]. The higher the Tversky index, the higher the likelihood of the predicted protein–ligand interaction.

> <hands-on-title>Hands-on 7: Pharmacophore alignment</hands-on-title>
>
> 1. {% tool [Pharmacophore alignment](toolshed.g2.bx.psu.edu/repos/bgruening/align_it/ctb_alignit/1.0.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Defines the database of molecules that will be used to screen"*: `staurosporine_3D_conformers` (from Hands-on 7)
>    - {% icon param-file %} *"Reference molecule"*: `ePharmaLib_PLAF_split` (from Hands-on 5)
>    - *"No normal information is included during the alignment"*: `Yes`
>    - *"Disable the use of hybrid pharmacophore points"*: `Yes`
>    - *"Only structures with a score larger than this cutoff will be written to the files"*: `0.0`
>    - *"Maximum number of best scoring structures to write to the files"*: `1`
>    - *"This option defines the used scoring scheme"*: `TVERSKY_REF`
>
{: .hands_on}

# Post-processing

The above pharmacophore alignment produces three types of outputs: the aligned pharmacophores (PHAR format), aligned structures (SMI format), and alignment scores (tabular format). Of these results, only the alignment scores are of interest and will be post-processed prior to analysis.

## Concatenating the pharmacophore alignment scores

The alignment score of the best ranked ligand conformer aligned against each ePharmaLib pharmacophore is stored in an individual file. In total, this job generates a collection of 138 output files which should be concatenated in a single file, for a better overview of the predictions.

> <hands-on-title>Hands-on 8: Concatenating the scores</hands-on-title>
>
> 1. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: `scores` (from Hands-on 7)
>
> 2. Rename the output to `concatenated_scores`.
>
>    {% snippet faqs/galaxy/datasets_rename.md name = "concatenated_scores" %}
>
{: .hands_on}

## Ranking the predicted protein targets

The resulting `concatenated_scores` needs to be re-sorted according to the alignment metric, the Tversky index, i.e. the 10th column. The pharmacophores of the ePharmaLib dataset were labeled according to the following three-component code *PDBID-hetID-UniprotEntryName*. The contents of the `concatenated_scores` are as follows:

	------    ---------------------------------------------------------------------
	column    Content
	------    ---------------------------------------------------------------------
	     1    Id of the reference structure
	     2    Maximum volume of the reference structure
	     3    Id of the database structure
	     4    Maximum volume of the database structure
	     5    Maximum volume overlap of the two structures
	     6    Overlap between pharmacophore and exclusion spheres in the reference
	     7    Corrected volume overlap between database pharmacophore and reference
	     8    Number of pharmacophore points in the processed pharmacophore
	     9    TANIMOTO score
	    10    TVERSKY_REF score
	    11    TVERSKY_DB score
	------    --------------------------------------------------------------------- 

> <hands-on-title>Hands-on 9: Sort Dataset</hands-on-title>
>
> 1. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `concatenated_scores` (from Hands-on 8)
>    - *"on column"*: `c10`
>
> 2. Rename the output to `final_target_prediction_scores`.
>
>    {% snippet faqs/galaxy/datasets_rename.md name = "final_target_prediction_scores" %}
>
> 3. You can view the contents of the collection `final_target_prediction_scores` by pressing the **eye icon** (View data).
>
> The top-ranked protein of our target prediction experiment is *4mvf-STU-CDPK2_PLAFK* (Figures 1 & 2) with a Tversky index = 0.73. The general observation that can be made from this ranking of protein hits is the high self-retrieval rate of known targets, which demonstrates the high prediction accuracy of the method. The higher the Tversky index, the higher the likelihood of the predicted protein–ligand interaction; with a value of 0.5 corresponding to a 50% likelihood.
>
{: .hands_on}

> <question-title></question-title>
> 
> Why was a perfect pharmacophore alignment (Tversky index =  1) not achieved for the top-ranked protein target for which the cocrystallized ligand is staurosporine (*STU*)?
>
> > <solution-title></solution-title>
> >
> > A perfect pharmacophore alignment because a computational conformer generator (here RDConf in Hands-on 6) is unlikely to be able to reproduce a crystallographic (native) ligand pose with 100% accuracy.
> >
> {: .solution}
>
{: .question}

# One-step Zauberkugel workflow vs. multi-step workflow

For pharmacophore-based protein target prediction, you can choose to use Galaxy tools separately and in succession as described above, or alternatively use the one-step Zauberkugel workflow as described below (Figure 3).

> <hands-on-title>Upload the Zauberkugel workflow</hands-on-title>
>
> Upload the Zauberkugel workflow from the following URL:
>
> ```
> https://github.com/galaxyproject/training-material/blob/main/topics/computational-chemistry/tutorials/zauberkugel/workflows/main_workflow.ga
> ```
>
> {% snippet faqs/galaxy/workflows_import.md %}
>
> The Zauberkugel workflow requires only two inputs; the ligand structure file (SMI format) and the ePharmaLib dataset (PHAR format). The output of the prediction of human targets of staurosporine performed with the ePharmaLib human target subset (<https://zenodo.org/record/6055897>) and this workflow is available as a [Galaxy history](https://usegalaxy.eu/u/aurelien_moumbock/h/zauberkugel).
{: .hands_on}

![Snapshot of Zauberkugel workflow]({% link topics/computational-chemistry/images/zauberkugel.png %} "Zauberkugel — protein target prediction of a bioactive ligand with Align-it and ePharmaLib")

# Further analysis

To obtain a docking pose of a protein–ligand interaction predicted from pharmacophore-based protein target prediction, follow the [Protein–ligand docking](https://training.galaxyproject.org/training-material/topics/computational-chemistry/tutorials/cheminformatics/tutorial.html) Galaxy training.


# Conclusion

