---
layout: tutorial_hands_on
topic_name: Proteomics
tutorial_name: WIP_proteinQuant_SIL
---

# Introduction
To compare protein amounts in different samples from MS/MS data, two different experiment setups exist. Firstly, unmodified proteins can be measured in separate runs at one sample per MS-run. Secondly, proteins of samples to compare can be labelled with small chemical tags, mixed, and measured side-by-side in a single MS-run. 
There are two types of chemical tags: isobaric tags display the same mass on first hand, but fragment during the generation of the MS/MS spectra to yield reporter ions of different mass. The intensity of those reporter ions can be compared in MS/MS spectra. There are two types of isobaric tags commercially available: tandem mass tags (TMT) and isobaric tags for relative and absolute quantitation (iTRAQ). 
The second type of chemical tags are isotopic. They are chemically identical, but differ in their mass due to incorporated stable isotopes. Examples of different isotopic tags for stable isotope labelling (SIL) are ICAT, SILAC, dimethylation, or heavy oxygen (<sup>18</sup>O).
Incorporation of stable isotopes results in different peptide masses on MS1 level, which give rise to coeluting ion traces in the TIC with a mass difference typical for each different chemical tag.

This tutorial deals with protein quantitation via stable isotope labelling (SIL). We will use tools of the OpenMS suite.  
Because we solely cover *quantitation*, you need to perform peptide and protein ID in beforehand. To learn about protein ID in Galaxy, please consider [this tutorial](./proteinID_SG_PS.md).
The tutorial covers *relative* quantitation only (i.e. comparison of abundances in different samples, no absolute quantitation of peptides / proteins).

If you still are in the planning phase of your quantitative proteomics experiment, you may want to consider our tutorial on different [labelling methods](./labelfree-vs-labelled.md) first.

> In the Hands-on section of this tutorial, we will use a quantitative comparison of **HEK _OR_ E.coli** cell lysate as a test dataset. In this experiment, the very same cell lysate was once labelled with light, once with heavy **dimethyl _OR_ SILAC** and both samples were subsequently mixed in a certain ratio.
> Your objective in this hands-on-tutorial is to find out the correct mixing ratio of the test sample.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 2. [MS1 Feature Detection](#ms1-feature-detection)
> 1. [Importing and Converting Peptide and Protein IDs](#importing-and-converting-peptide-and-protein-ids)
> 3. [Quant to ID matching](#quant-to-id-matching)
> 4. [Evaluation and Optimization of Quantitation Results](#evaluation-and-optimization-of-quantitation-results)
{: .agenda}


# MS1 Feature Detection
Quantitation on MS1 level may in principle be carried out without prior knowledge of peptide / protein IDs. However, some quantitation algorithms take the IDs as an input to make sure that every PSM that was identified is also quantified. This is not the case in our example here. 
In the OpenMS suite, most of the provided tools for MS1 feature detection quantify solely based on mzML files. The advantage of this approach is that quantitations can be made on strict criteria to reduce misquantitations. The drawback is that not all IDs can be matched to a quantitation later on in the workflow. 
The tool settings need to be carefully tested and evaluated manually to obtain optimal results. We will explain this in the section [Evaluation and Optimization of Quantitation Results](#evaluation-and-optimization-of-quantitation-results). 

> ### :pencil2: Hands-on: MS1 Feature Detection
> 
> {: .hands_on}


# Importing and Converting Peptide and Protein IDs

# Quant to ID matching

# Evaluation and Optimization of Quantitation Results