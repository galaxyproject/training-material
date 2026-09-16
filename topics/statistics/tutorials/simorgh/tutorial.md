---
layout: tutorial_hands_on

title: Comparing ligand-binding site predictions across different protein structure modalities using SIMORGH
questions:
- TO_BE_ADDED
- TO_BE_ADDED
objectives:
- TO_BE_ADDED
- TO_BE_ADDED
- TO_BE_ADDED
time_estimation: 1H
key_points:
- TO_BE_ADDED
- TO_BE_ADDED
contributors:
- Nilchia
- amisteromid
subtopic: TO_BE_ADDED
answer_histories:
    - label: "usegalaxy.eu"
      history: TO_BE_ADDED
      date: TO_BE_ADDED
---

Proteins are not static objects. They constantly move and change shape, and these conformational changes can be especially important when a ligand binds.
A protein structure captured **before ligand binding** is called the **apo** form, while a structure captured **with a ligand bound** is called the **holo** form.

Why does this matter?

In the holo structure, the protein has already adopted a conformation compatible with the ligand. As a result, the binding pocket is often easier to identify: the pocket is already **pre-organized around the ligand**.

But what happens when the ligand is not there?
And what if the structure we use is not an experimentally determined structure at all, but an **Colabfold prediction**?

This tutorial explores exactly that question. We will use SIMORGH ({% cite SIMORGH2026 %}) to address this challenge.

Simorgh is a deep learning framework built on SE(3)-equivariant neural networks for ligand-binding site prediction.
It combines local geometric reasoning with information from diverse protein states.
A message-passing module first encodes the spatial context of individual residues within each structure, followed by an aggregation module that integrates these residue-level representations across the ensemble.
By explicitly modeling structural variability, Simorgh identifies binding sites that are meaningful across multiple protein states, rather than being biased toward a single structure.


> <warning-title>LICENSE</warning-title>
>
> SIMORGH is only available for NON-COMMERCIAL use. Permission is only granted for academic, research, and educational purposes. Before using, be sure to review, agree, and comply with the license.
> For commercial use, please review the SIMORGH license on GitHub and contact the [copyright holders](https://github.com/amisteromid/SIMORGH)
{: .warning}


We will use the **Escherichia coli Signal Recognition Particle Receptor FtsY** as our example.
FtsY binds **Guanosine-5'-diphosphate (GDP)**, and experimentally determined structures are available in both apo and holo conformations.

We will work with two experimental structures; one is 6N5J which is in Apo state ({% cite Ataide2019-en %}), and the other one is 6FQD which is in Holo state ({% cite Mrusek2018-nw %}).

The interesting question is:

If we predict the ligand-binding site from different structural representations of the same protein, how much does the prediction change?

To answer this, we will compare **three structural modalities**:
* Native apo structure — experimentally determined
* Native holo structure — experimentally determined
* Colabfold structure — computationally predicted from sequence


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get the structures

## Experimental structures

> <hands-on-title> Get PDB file </hands-on-title>
>
> 1. {% tool [Get PDB file](toolshed.g2.bx.psu.edu/repos/bgruening/get_pdb/get_pdb/0.1.1) %} with the following parameters:
>    - *"PDB accession code"*: `6N5J`
>    - In *"Aditional options"*:
>       - *"Output Tags"*: `#apo`
>
> 2. Rename the data `6N5J PDB`
>
> 3. {% tool [Get PDB file](toolshed.g2.bx.psu.edu/repos/bgruening/get_pdb/get_pdb/0.1.1) %} with the following parameters:
>    - *"PDB accession code"*: `6FQD`
>    - In *"Aditional options"*:
>       - *"Output Tags"*: `#holo`
>
> 4. Rename the data `6FQD PDB`
>
{: .hands_on}

## Predicted structures

We will use Colabfold ({% cite Mirdita2022-ko %}) to predict the 3D structure of our protein.

> <hands-on-title> Colabfold prediction </hands-on-title>
>
> 1. Create a new **fasta** dataset (`FtsY.fasta`) from the following:
>
>    ```text
>    6N5J_1|Chains A, B|Signal recognition particle receptor FtsY|Escherichia coli (strain K12) (83333)
>    GFARLKRSLLKTKENLGSGFISLFRGKKIDDDLFEELEEQLLIADVGVETTRKIITNLTEGASRKQLRDAEALYGLLKEE
>    MGEILAKVDEPLNVEGKAPFVILMVGVNGVGKTTTIGKLARQFEQQGKSVMLAAGDTFRAAAVEQLQVWGQRNNIPVIAQ
>    HTGADSASVIFDAIQAAKARNIDVLIADTAGRLQNKSHLMEELKKIVRVMKKLDVEAPHEVMLTIDASTGQNAVSQAKLF
>    HEAVGLTGITLTKLDGTAKGGVIFSVADQFGIPIRYIGVGERIEDLRPFKADDFIEALFARED
>    ```
>
>    {% snippet faqs/galaxy/datasets_create_new_file.md name="FtsY.fasta" format="fasta" %}
>
> 2. {% tool [Colabfold MSA](toolshed.g2.bx.psu.edu/repos/iuc/colabfold_msa/colabfold_msa/1.5.5+galaxy1) %} with the following parameters:
>    - *"Data input method"*: `FASTA file`
>    - *"Query sequence fasta"*: `FtsY.fasta`
>
> 3. {% tool [Colabfold Alphafold](toolshed.g2.bx.psu.edu/repos/iuc/colabfold_alphafold/colabfold_alphafold/1.5.5+galaxy1) %} with the following parameters:
>    - *Tar file output from colabfold MSA tool"*: `Output of Colabfold MSA`
>    - In *"Advanced options"*:
>       - *"How many recycles to run?"*: `5`
>       - *"Number of ensembles"*: `1`
>       - *"Set seed"*: `42`
>       - *"Number of models to use for structure prediction"*: `1`
>       - *"Use AMBER"*: `Don't use AMBER`
>
> 3. {% tool [ Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - *Input List"*: `PDB predictions` output of Colabfold Alphafold
>    - *How should a dataset be selected?"*: `The first dataset`
>
> 4. Rename the data `FtsY Colabfold predicted`
>
{: .hands_on}


> <question-title></question-title>
>
> 1. TO_BE_ADDED
>
> > <solution-title></solution-title>
> >
> > 1. TO_BE_ADDED
> >
> {: .solution}
>
{: .question}

# Structure preprocessing

## PDB cleanup

Next we should preprocess the raw exprimental PDB structures to fix structural anomalies. We will use pdbfixer ({% cite pdbfixer %}) to identify and model missing heavy atoms, insert missing flexible loops or internal side chains based on the sequence, remove crystallographic water molecules or unwanted ligand, replace non-standard residue entries with standard amino acids, and add explicit hydrogen atoms at specified pH levels to ensure physical completeness for downstream modeling.

> <hands-on-title>  PDBFixer </hands-on-title>
>
> 1. {% tool [ PDBFixer](toolshed.g2.bx.psu.edu/repos/chemteam/pdbfixer/pdbfixer/1.8.1+galaxy0) %} with the following parameters:
>    - *"PDB input file"*: `6N5J PDB`
>    - *"Missing atoms to be added"*: `Heavy atoms only`
>    - *"Which heterogens to keep"*: `None`
>    - *"Replace nonstandard residues with standard equivalents?"*: `Yes`
>
> 2. Rename the data `6N5J PDB fixed`
>
> 3. {% tool [ PDBFixer](toolshed.g2.bx.psu.edu/repos/chemteam/pdbfixer/pdbfixer/1.8.1+galaxy0) %} with the following parameters:
>    - *"PDB input file"*: `6FQD PDB`
>    - *"Missing atoms to be added"*: `Heavy atoms only`
>    - *"Which heterogens to keep"*: `None`
>    - *"Replace nonstandard residues with standard equivalents?"*: `Yes`
>
> 4. Rename the data `6FQD PDB fixed`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. TO_BE_ADDED
>
> > <solution-title></solution-title>
> >
> > 1. TO_BE_ADDED
> >
> {: .solution}
>
{: .question}

## Retrieve chain A

Since FtsY is a homodimer, we will use a single chain A for simplicity. For this task, we will use *"Text reformating with awk"* tool ({% cite Gruning2018-lh %}) to reformat our pdb files.

The command `rectype = substr($0,1,6)`, gets the characters 1-6 of the current line (such as ATOM, HETATM, or TER).

Then the command `chain = substr($0,22,1)` gets the character 22 of the line which corresponds to either A or B (the chain identifiers).

Then, `if ((rectype ~ /^ATOM/ || rectype ~ /^HETATM/ || rectype ~ /^TER/) && chain == "B") next` skips any line that is `ATOM, HETATM, or TER` and it is chain `B`.

Finally, `print`, prints out any line that is not skipped.

In summary, "Keep everything except ATOM/HETATM/TER records from chain B."

Copy the following script in the ***"AWK Program"*** of the following tool:
```
{
  rectype = substr($0,1,6)
  chain = substr($0,22,1)
  if ((rectype ~ /^ATOM/ || rectype ~ /^HETATM/ || rectype ~ /^TER/) && chain == "B") next
  print
}
```

> <hands-on-title>   Text reformatting with awk </hands-on-title>
>
> 1. {% tool [  Text reformatting with awk](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/9.11+galaxy0) %} with the following parameters:
>    - *"File to process"*: `6N5J PDB fixed`
>    - *"AWK Program"*: `copy the text from above`
>
> 2. Rename the data `6N5J PDB fixed chainA`
>
> 3. {% tool [  Text reformatting with awk](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/9.11+galaxy0) %} with the following parameters:
>    - *"File to process"*: `6FQD PDB fixed`
>    - *"AWK Program"*: `copy the text from above`
>
> 4. Rename the data `6FQD PDB fixed chainA`
>
{: .hands_on}

# Binding Site Prediction with SIMORGH

A single protein structure is only one snapshot of a moving molecule. To capture this motion, we generate a conformational ensemble for each structure using BBFlow ({% cite Wolf2025-dm %}).
BBFlow is a structure-conditioned flow-matching model for generating diverse protein backbone ensembles.
Given a protein backbone as input, it stochastically samples structurally distinct states, capturing the underlying geometric heterogeneity of the protein.
This structural diversity is then leveraged by SIMORGH to predict ligand-binding sites, reducing bias from relying on a single input structure and improving the performance.


> <hands-on-title> SIMORGH </hands-on-title>
>
> 1. {% tool [SIMORGH](toolshed.g2.bx.psu.edu/repos/bgruening/SIMORGH/SIMORGH/1.0.1+galaxy0) %} with the following parameters:
>    - *"I certify that I am NOT using this tool for commercial purposes."*: `Yes`
>    - *"PDB file"*: `6N5J PDB fixed chainA`
>    - *"Run BBflow?"*: `Yes`
>    - *"Number of conformations to sample"*: `8`
>
> 2. {% tool [SIMORGH](toolshed.g2.bx.psu.edu/repos/bgruening/SIMORGH/SIMORGH/1.0.1+galaxy0) %} with the following parameters:
>    - *"I certify that I am NOT using this tool for commercial purposes."*: `Yes`
>    - *"PDB file"*: `6FQD PDB fixed chainA`
>    - *"Run BBflow?"*: `Yes`
>    - *"Number of conformations to sample"*: `8`
>
> 3. {% tool [SIMORGH](toolshed.g2.bx.psu.edu/repos/bgruening/SIMORGH/SIMORGH/1.0.1+galaxy0) %} with the following parameters:
>    - *"I certify that I am NOT using this tool for commercial purposes."*: `Yes`
>    - *"PDB file"*: `FtsY Colabfold predicted`
>    - *"Run BBflow?"*: `Yes`
>    - *"Number of conformations to sample"*: `8`
>
{: .hands_on}


> <question-title></question-title>
>
> 1. TO_BE_ADDED
>
> > <solution-title></solution-title>
> >
> > 1. TO_BE_ADDED
> >
> {: .solution}
>
{: .question}


# What Does SIMORGH Actually Learn?

SIMORGH can incorporate protein dynamics in two stages.

* Stage 1 --> **Conformation Augmentation**
    The first encoder sees different conformations of the same protein.
    Instead of learning from a single static structure, SIMORGH is exposed to the **structural variability of the protein**.
* Stage 2 --> **Learned Aggregation**
    SIMORGH then learns how to combine the embeddings from these different conformations into **one final prediction**.

In other words:
**Multiple conformations** --> **Multiple embeddings** --> **Learned aggregation** --> **Final prediction**

This allows the model to use information from the protein's dynamic ensemble rather than relying on a single structural snapshot.


# Bonus Challenge

Although this would not count as a prediction that incorporates learned protein dynamics, you can run SIMORGH without BBFlow to isolate the effect of aggregation.

This time we run SIMORGH on `6N5J PDB fixed chainA` without using BBflow.

> <hands-on-title> SIMORGH </hands-on-title>
>
> 1. {% tool [SIMORGH](toolshed.g2.bx.psu.edu/repos/bgruening/SIMORGH/SIMORGH/1.0.1+galaxy0) %} with the following parameters:
>    - *"I certify that I am NOT using this tool for commercial purposes."*: `Yes`
>    - *"PDB file"*: `6N5J PDB fixed chainA`
>    - *"Run BBflow?"*: `No`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. TO_BE_ADDED
>
> > <solution-title></solution-title>
> >
> > 1. TO_BE_ADDED
> >
> {: .solution}
>
{: .question}

# Conclusion

TO_BE_ADDED
