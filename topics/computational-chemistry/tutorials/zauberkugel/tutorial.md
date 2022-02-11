---
layout: tutorial_hands_on

title: Protein target prediction of a bioactive ligand with Align-it and ePharmaLib
level: Intermediate
zenodo_link: 'https://zenodo.org/record/5898113/files/ePharmaLib_PHARAO_human.phar'
questions:
- What is a pharmacophore model?
- How can I perform protein target prediction with individual tools (multiple steps) or the Zauberkugel workflow (single step)?
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

---


# Introduction
{:.no_toc}

Historically, the pharmacophore concept was formulated in 1909 by the German physician and Nobel prize laureate Paul Ehrlich ({% cite Ehrlich1909 %}). According to the [International Union of Pure and Applied Chemistry (IUPAC)](https://iupac.org/), a pharmacophore is defined as “an ensemble of steric and electronic features that is necessary to ensure the optimal supramolecular interactions with a specific biological target and to trigger (or block) its biological response” ({% cite Wermuth1998 %}). Starting from the cocrystal structure of a non-covalent protein--ligand complex (e.g. Figure 1), pharmacophore perception involves the extraction of the key molecular features of the bioactive ligand at the protein--ligand contact interface into a single model ({% cite Moumbock2019 %}). These pharmacophoric features mainly include: H-bond acceptor (HACC or A), H-bond donor (HDON or D), lipophilic group (LIPO or H), negative center (NEGC or N), positive center (POSC or P), and aromatic ring (AROM or R) moieties. Moreover, receptor-based excluded spheres (EXCL) can be added in order to mimic spatial constraints of the binding pocket (Figure 2). Once a pharmacophore model has been generated, a query can be performed either in a forward manner, using several ligands to search for novel putative hits of a given target, or in a reverse manner, by screening a single ligand against multiple pharmacophore models in search of putative protein targets ({% cite Steindl2006 %}).

![PDB ID: 4FR4]({% link topics/computational-chemistry/images/4FR4-STU.png %} "Crystal structure of the human serine/threonine-protein kinase in complex with staurosporine (PDB ID: 4FR4). Image generated using Maestro (Schrödinger LLC, NY).")

Bioactive compounds often bind to several target proteins, thereby exhibiting polypharmacology. However, experimentally determining these interactions is laborious, and structure-based virtual screening of bioactive compounds could expedite drug discovery by prioritizing hits for experimental validation. The recently reported ePharmaLib ({% cite Moumbock2021 %}) dataset is a library of 15,148 e-pharmacophores modeled from solved structures of pharmaceutically relevant protein--ligand complexes of the screening Protein Data Bank (sc-PDB, {% cite Desaphy2014 %}). ePharmaLib can be used for target fishing of phenotypic hits, side effect predictions, drug repurposing, and scaffold hopping.

![STU]({% link topics/computational-chemistry/images/STU.png %} "Overlay of staurosporine onto its key pharmacophoric features (PDB ID: 4FR4). Image generated using Maestro (Schrödinger LLC, NY).")

In this tutorial, you will perform pharmacophore-based target prediction of a bioactive ligand known as staurosporine (Figure 2) with the ePharmaLib subset representing human protein targets (5,010 pharmacophore models) and the open-source pharmacophore alignment program [Align-it](https://anaconda.org/bioconda/align_it), formerly known as PHARAO ({% cite Taminau2008 %}).

> ### {% icon details %} Pharmacology of staurosporine
>
>  Staurosporine (PDB hetID: [STU](https://www.rcsb.org/ligand/STU)) is an indolocarbazole secondary metabolite isolated from several bacteria of the genus Streptomyces. It displays highly potent and broad inhibition of human protein kinases, and serves as the precursor of numerous marketed kinase inhibitors, notably imatinib --- the first FDA-approved small-molecule kinase inhibitor ({% cite mura2018 %}).
{: .details}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Create a history

As a first step, we create a new history for the analysis.

> ### {% icon hands_on %} Hands-on 1: Create history
>
> 1. Create a new history.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename it to `Staurosporine target prediction`.
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}


# Get data

For this exercise, we need two datasets: the ePharmaLib pharmacophore library (PHAR format) and a query ligand structure file (SMI format). 

## Fetching the ePharmaLib dataset

Firstly, we will retrieve the concatenated ePharmaLib subset representing human protein targets. 

> ### {% icon hands_on %} Hands-on 2: Upload ePharmaLib
>
> 1. Upload the dataset from the [Zenodo](https://zenodo.org) link provided to your Galaxy history.
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md link="https://zenodo.org/record/5898113/files/ePharmaLib_PHARAO_human.phar" %}
>
>    > ### {% icon comment %} ePharmaLib versions
>    >
>    > Two versions of the ePharmaLib (PHAR & PHYPO formats) have been created for use with the pharmacophore alignment programs [Align-it](https://anaconda.org/bioconda/align_it) and [Phase](https://www.schrodinger.com/products/phase), respectively. They are freely available at [Zenodo](https://zenodo.org) under the link:
>    >	```
>    >	https://zenodo.org/record/5898113/
>    >  ```
>    {: .comment}
>
> 2. Change the datatype from `tabular`to `phar`. This step is essential, as Galaxy does not automatically detect the datatype for PHAR files.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="phar" %}
>
{: .hands_on}

## Creating a query ligand structure file

In this step, we will manually create an SMI file containing the SMILES of staurosporine.

> ### {% icon details %} What are SMILES and the SMI file format?
>
> The simplified molecular-input line-entry system (SMILES) is a string notation for describing the 2D chemical structure of a compound. It only states the atoms present in a compound and the connectivity between them. As an example, the SMILES string of acetone is `CC(=O)C`. SMILES strings can be imported by most molecule editors and converted into either two-dimensional structural drawings or three-dimensional models of the compounds, and vice versa. For more information on how the notation works, please consult the [OpenSMILES specification](http://opensmiles.org/opensmiles.html) or the description provided by [Wikipedia](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system). 
>
{: .details}

> ### {% icon hands_on %} Hands-on 3: Create an SMI file
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
> > ### {% icon tip %} Tip: SMILES generation
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

> ### {% icon question %} Question
>
> Why do we specifically use a so-called isomeric SMILES string?
>
> > ### {% icon solution %} Solution
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

> ### {% icon hands_on %} Hands-on 4: Add hydrogen atoms
>
> 1. {% tool [Add hydrogen atoms](toolshed.g2.bx.psu.edu/repos/bgruening/openbabel_addh/openbabel_addh/3.1.1+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Molecular input file"*: `output` (Input dataset)
>    - *"Add hydrogens to polar atoms only (i.e. not to carbon atoms)"*: `Yes`
>
{: .hands_on}

> ### {% icon question %} Question
>
> Nitrogen-containing functional groups are known to be basic. Which of them present in staurosporine do you expect to be protonated at pH 7.4, and which not? And why?
>
> > ### {% icon solution %} Solution
> >
> > Only the secondary N-methylamino group will be protonated because indoles much like aromatic amides are typically not basic.
> >
> {: .solution}
>
{: .question}

## Splitting ePharmaLib into individual pharmacophores

The ePharmaLib subset representing human protein targets (*ePharmaLib_PHARAO_human.phar*) is a concatenated file containing 5,010 individual pharmacophore files. In order to speed up our analysis, it is preferable to split the dataset into individual files in order to perform several pharmacophore alignments in parallel, using Galaxy's collection functionality.

> ### {% icon hands_on %} Hands-on 5: Splitting ePharmaLib
>
> 1. {% tool [Split file](toolshed.g2.bx.psu.edu/repos/bgruening/split_file_to_collection/split_file_to_collection/0.5.0) %} with the following parameters:
>    - *"Select the file type to split"*: `Generic`
>        - {% icon param-file %} *"File to split"*: `output` (Input dataset)
>        - *"Method to split files"*: `Specify record separator as regular expression`
>            - *"Regex to match record separator"*: `\$\$\$\$`
>            - *"Split records before or after the separator?"*: `After`
>        - *"Specify number of output files or number of records per file?"*: `Number of records per file ('chunk mode')`
>        - *"Base name for new files in collection"*: `epharmalib`
>        - *"Method to allocate records to new files"*: `Maintain record order`
>
>    > ### {% icon details %} What is a PHAR file?
>    >
>    > A PHAR file is essentially a series of lines containing the three-dimensional coordinates of pharmacophoric features and excluded spheres. The first column specifies a feature type (e.g. HACC is a hydrogen bond acceptor). Subsequent columns specify the position of the feature center in a three-dimensional space. Individual pharmacophores are separated by lines containing four dollar signs (*$$$$*).
>    {: .details}
>
{: .hands_on}

## Ligand conformational flexibility

In order to avoid combinatorial explosion, the Align-it ({% cite Taminau2008 %}) tool performs rigid alignment rather than flexible alignment. Conformational flexibility of the ligand is accounted for by introducing a preliminary step, in which a set of energy-minimized conformers for the query ligand are generated with the RDConf ({% cite rdconf %}) tool (using the RDKit ({% cite landrum2013rdkit %}) toolkit). 

> ### {% icon hands_on %} Hands-on 6: Low-energy ligand conformer search
>
> 1. {% tool [RDConf: Low-energy ligand conformer search](toolshed.g2.bx.psu.edu/repos/bgruening/rdconf/rdconf/2020.03.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `outfile` (output of **Add hydrogen atoms** {% icon tool %})
>    - *"Maximum number of conformers to generate per molecule"*: `100`
>
>    > ### {% icon comment %} RDConf
>    >
>    > It is recommended to use the default settings, except for the number of conformers which should be changed to 100. As a rule of thumb, a threshold of 100 conformers appropriately represents the conformational flexibility of a compound with less than 10 rotatable bonds. The output SDF (structure data file) format encodes three-dimensional atomic coordinates of each conformer, separated by lines containing four dollar signs (*$$$$*).
>    {: .comment}
>
{: .hands_on}

# Pharmacophore alignment

In this step, the ligand conformer dataset (SDF format) is converted on-the-fly to a pharmacophore dataset (PHAR format) and simultaneously aligned to the individual pharmacophores of the ePharmaLib dataset in a batch mode with Align-it ({% cite Taminau2008 %}). The pharmacophoric alignments and thus the predicted targets are ranked in terms of a scoring metric: `Tversky index` = [0,1]. The higher the Tversky index, the higher the likelihood of the predicted protein--ligand interaction. A cut-off of 0.5 is applied to prioritize the most promising targets.

> ### {% icon hands_on %} Hands-on 7: Pharmacophore alignment
>
> 1. {% tool [Pharmacophore alignment](toolshed.g2.bx.psu.edu/repos/bgruening/align_it/ctb_alignit/1.0.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Defines the database of molecules that will be used to screen"*: `outfile` (output of **RDConf: Low-energy ligand conformer search** {% icon tool %})
>    - {% icon param-file %} *"Reference molecule"*: `list_output_generic` (output of **Split file** {% icon tool %})
>    - *"No normal information is included during the alignment"*: `Yes`
>    - *"Disable the use of hybrid pharmacophore points"*: `Yes`
>    - *"Only structures with a score larger than this cutoff will be written to the files"*: `0.5`
>    - *"Maximum number of best scoring structures to write to the files"*: `1`
>    - *"This option defines the used scoring scheme"*: `TVERSKY_REF`
>
{: .hands_on}

This step takes about 45 minutes to an hour to complete.

> ### {% icon tip %} Tip: Problems using the Align-it tool?
> Due to the large number of files contained in the ePharmaLib dataset collection, users occasionally encounter technical issues with this particular Align-it job. Sometimes the job takes too long to successfully complete (more than the estimated 1 hour duration) depending on the available Galaxy computing resources. And in other cases the job completes with errors; if this happens, delete the errored outputs and rerun the job.
{: .tip}

# Post-processing

The above pharmacophore alignment produces three types of outputs: the aligned pharmacophores (PHAR format), aligned structures (SMI format), and alignment scores (tabular format). Of these results, only the alignment scores are of interest and will be post-processed prior to analysis.

## Concatenating the pharmacophore alignment scores

The alignment score of the best ranked ligand conformer aligned against each ePharmaLib pharmacophore is stored in individual files. In total, this job generates 5,010 output files which should be concatenated in a single file, for a better overview of the predictions.

> ### {% icon hands_on %} Hands-on 8: Concatenating the scores
>
> 1. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: `score_result_file` (output of **Pharmacophore alignment** {% icon tool %})
>
{: .hands_on}

## Ranking the predicted protein targets

The resulting single alignment score file needs to be re-sorted according to the alignment metric, the Tversky index, i.e. the 10th column. The pharmacophores of the ePharmaLib dataset were labeled according to the following three-component code *PDBID-hetID-UniprotEntryName*. Given that staurosporine (hetID: *STU*) appears in several PDB cocrystal structures (43 times in the sc-PDB dataset from which ePharmaLib was built), we shall evaluate the predictions results by identifying the predicted proteins whereby the cocrystallized ligand is `STU`.

> ### {% icon hands_on %} Hands-on 9: Sort Dataset
>
> 1. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>    - *"on column"*: `c10`
>
> 2. {% tool [Text reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.2) %} with the following parameters:
>    - *"AWK Program"*: `/STU/{print NR, $0}`
>
>    > ### {% icon comment %} Further investigation (optional)
>    >
>    > The second task performed in the above *Hands-on 9* is optional and should only be performed if the query ligand possess a *PDB hetID*, i.e. is present as a cocrystallized structure in PDB structures.
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
> 
> Output of task 1:
> 1. For predicted targets in the top 10 whose cocrystallized ligands are other than *STU*, are these ligands structurally similar to *STU*?
> Output of task 2:
> 2. Out of the 43 STU-containing targets of the sc-PDB ({% cite Desaphy2014 %}), how many were retrieved by the current target prediction experiment?
> 3. Why a perfect pharmacophore alignment (`Tversky index =  1`) was not achieved for protein targets where the cocrystallized ligand is staurosporine (*STU*)?
>
> > ### {% icon solution %} Solution
> >
> > 1. With the exception of *TIY* (purpurogallin), the other cocrystallized ligands for the top 10 targets, namely: *ITQ* , *UCM*, and *UCN* are *STU*-derived compounds, wherein the indolocarbazole scaffold was preserved.
> > 2. Thirty nine (39).
> > 3. A perfect pharmacophore alignment because a computational conformer generator can hardly reproduce with 100% accuracy a crystallographic (native) ligand pose.
> >
> {: .solution}
>
{: .question}

# Zauberkugel workflow vs. individual tools

For pharmacophore-based protein target prediction, you can choose to use Galaxy tools separately and in succession as described above, or alternatively use the Zauberkugel workflow as described below (Figure 3).

![Snapshot of Zauberkugel workflow]({% link topics/computational-chemistry/images/zauberkugel.png %} "Zauberkugel --- protein target prediction of a bioactive ligand with Align-it and ePharmaLib")

> ### {% icon hands_on %} Upload the Zauberkugel workflow
>
>    {% snippet faqs/galaxy/workflows_import.md %}
>
>    ```
>    https://github.com/galaxyproject/training-material/blob/main/topics/computational-chemistry/tutorials/zauberkugel/workflows/main_workflow.ga
>    ```
>
>    The Zauberkugel workflow requires only two inputs; the ligand structure file (SMI format) and the ePharmaLib dataset (PHAR format).
>
{: .hands_on}


# Further analysis

To obtain a docking pose of a protein--ligand interaction predicted from pharmacophore-based protein target prediction, try out the [Protein-ligand docking](https://github.com/galaxyproject/training-material/blob/main/topics/computational-chemistry/tutorials/cheminformatics/workflows/main_workflow.ga).

# Conclusion
{:.no_toc}
