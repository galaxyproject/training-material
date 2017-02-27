---
layout: tutorial_hands_on
topic_name: Peptide and Protein Quantification via Stable Isotope Labelling (SIL)
tutorial_name: proteinQuant_SIL
---

# Introduction
To compare protein amounts in different samples from MS/MS data, two different experiment setups exist. Firstly, unmodified proteins can be measured in separate runs at one sample per MS-run. Secondly, proteins of samples to compare can be labelled with small chemical tags, mixed, and measured side-by-side in a single MS-run. There are two types of chemical tags: isobaric tags display the same mass on first hand, but fragment during the generation of the MS/MS spectra to yield reporter ions of different mass. The intensity of those reporter ions can be compared in MS/MS spectra. There are two types of isobaric tags commercially available: tandem mass tags (TMT) and isobaric tags for relative and absolute quantitation (iTRAQ). The second type of chemical tags are isotopic. They are chemically identical, but differ in their mass due to incorporated stable isotopes. Examples of different isotopic tags for stable isotope labelling (SIL) are ICAT, SILAC, dimethylation, or heavy oxygen (<sup>18</sup>O). Incorporation of stable isotopes results in different peptide masses on MS1 level, which give rise to coeluting ion traces in the TIC with a mass difference typical for each different chemical tag.

This tutorial deals with protein quantitation via stable isotope labelling (SIL). We will use tools of the OpenMS suite. The tutorial covers *relative* quantitation only (i.e. comparison of abundances in different samples, no absolute quantitation of peptides / proteins). Because we solely cover the quantitation, you need to perform peptide and protein ID in beforehand. To learn about protein ID in Galaxy, please consider [this tutorial](./proteinID_SG_PS.md).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 2. [MS1 Feature Detection](#ms1-feature-detection)
> 1. [Importing and Converting Peptide and Protein IDs](#importing-ids)
> 3. [Quant to ID matching](#quant-to-id-matching)
> 4. [Evaluation and Optimization of Quantitation Results](#evaluation-and-optimization-of-quantitation-results)


# MS1 Feature Detection
Quantitation on MS1 level may in principle be carried out without prior knowledge of peptide / protein IDs. However, some quantitation algorithms take the IDs as an input to make sure that every PSM that was identified is also quantified. This is not the case in our example here. In the OpenMS suite, most of the provided tools for MS1 feature detection quantify solely based on mzML files. The advantage of this approach is that quantitations can be made on strict criteria to reduce misquantitations. The drawback is that not all IDs can be matched to a quantitation later on in the workflow. The tool settings need to be carefully tested and evaluated manually to obtain optimal results. We will explain this in the section [Evaluation and Optimization of Quantitation Results](#evaluation-and-optimization-of-quantitation-results). 

> ### :nut_and_bolt: Comment: OpenMS MS1 Quantitation Tools
> The OpenMS suite features quite a lot of different MS1 feature detection tools (called "FeatureFinder"s). Momentarily (as of February 2017), the recommended standard tool is the **FeatureFinderCentroided** :wrench:. It has been thoroughly tested, but is not longer improved and may be replaced by the newer **FeatureFinderMultiplex** :wrench: in future. We will focus on the first tool here, but you may well try the later one if it more suits your need. A recent development is the **FeatureFinderIdentification** :wrench:. Unlike the other FeatureFinders in the OpenMS suite, it takes PSM IDs as an input. However, it is still under development and not recommended as default tool right now (as of February 2017). For more information on OpenMS FeatureFinders, please consult [this discussion](https://github.com/OpenMS/OpenMS/issues/2424#issuecomment-282293381). 

> ### :pencil2: Hands-on: MS1 Feature Detection
> 
> Here, we will use a quantitative comparison of **HEK _OR_ E.coli** cell lysate as a test dataset. In this experiment, the very same cell lysate was once labelled with light, once with heavy **dimethyl _OR_ SILAC** and both samples were subsequently mixed in a certain ratio. Your objective in this hands-on-tutorial is to find out the correct mixing ratio. *And no, the ratio is not mentioned in the description of the dataset on PRIDE. We left out that particular information on purpose. Haha.*
> 
> 

<a name="importing-ids"/></a>
# Importing and Converting Peptide and Protein IDs


# Quant to ID matching
