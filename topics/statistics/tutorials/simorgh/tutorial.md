---
layout: tutorial_hands_on

title: TO_BE_ADDED
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
And what if the structure we use is not an experimentally determined structure at all, but an **AlphaFold prediction**?

This tutorial explores exactly that question.

We will use the **Escherichia coli Signal Recognition Particle Receptor FtsY** as our example.
FtsY binds **Guanosine-5'-diphosphate (GDP)**, and experimentally determined structures are available in both apo and holo conformations.

The two experimental structures:
* 6N5J | pdb_00006n5j → Apo form
* 6FQD | pdb_00006fqd → Holo form, co-crystallized with GDP


The interesting question is:

If we predict the ligand-binding site from different structural representations of the same protein, how much does the prediction change?

To answer this, we will compare **three structural modalities**:
* Native apo structure — experimentally determined
* Native holo structure — experimentally determined
* AlphaFold structure — computationally predicted from sequence

Since FtsY is a homodimer, we will use a single chain for simplicity.


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
> 1. {% tool [Get PDB file](toolshed.g2.bx.psu.edu/repos/bgruening/get_pdb/get_pdb/0.1.0) %} with the following parameters:
>    - *"PDB accession code"*: `6N5J`
>
> 2. Rename the data `6N5J PDB`
>
> 3. {% tool [Get PDB file](toolshed.g2.bx.psu.edu/repos/bgruening/get_pdb/get_pdb/0.1.0) %} with the following parameters:
>    - *"PDB accession code"*: `6FQD`
>
> 4. Rename the data `6FQD PDB`
>
{: .hands_on}
## Predicted structures

> <hands-on-title> Alphafold 2 prediction </hands-on-title>
>
> 1. {% tool [Alphafold 2](toolshed.g2.bx.psu.edu/repos/galaxy-australia/alphafold2/alphafold/2.3.2+galaxy4) %} with the following parameters:
>    - *"Fasta Input"*: `Paste sequence into textbox`
>    - *"Paste sequence"*: >6N5J_1|Chains A, B|Signal recognition particle receptor FtsY|Escherichia coli (strain K12) (83333)
GFARLKRSLLKTKENLGSGFISLFRGKKIDDDLFEELEEQLLIADVGVETTRKIITNLTEGASRKQLRDAEALYGLLKEEMGEILAKVDEPLNVEGKAPFVILMVGVNGVGKTTTIGKLARQFEQQGKSVMLAAGDTFRAAAVEQLQVWGQRNNIPVIAQHTGADSASVIFDAIQAAKARNIDVLIADTAGRLQNKSHLMEELKKIVRVMKKLDVEAPHEVMLTIDASTGQNAVSQAKLFHEAVGLTGITLTKLDGTAKGGVIFSVADQFGIPIRYIGVGERIEDLRPFKADDFIEALFARED
>    - In *"Advanced Options "*:
>       - *"Limit model outputs"*: `1`
>
> 2. Rename the data `FtsY predicted`
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


# Binding Site Prediction with SIMORGH

TO_BE_ADDED

# What Does SIMORGH Actually Learn?

TO_BE_ADDED

# Bonus Challenge

TO_BE_ADDED

# Conclusion

TO_BE_ADDED
