---
layout: tutorial_hands_on
topic_name: Peptide and Protein Quantification via Stable Isotope Labelling (SIL)
tutorial_name: proteinQuant_SIL
---

# Introduction
To compare protein amounts in different samples from MS/MS data, two different experiment setups exist. Firstly, unmodified proteins can be measured in separate runs at one sample per MS-run. Secondly, proteins of samples to compare can be labelled with small chemical tags, mixed, and measured side-by-side in a single MS-run. There are two types of chemical tags: isobaric tags display the same mass on first hand, but fragment during the generation of the MS/MS spectra to yield reporter ions of different mass. The intensity of those reporter ions can be compared in MS/MS spectra. There are two types of isobaric tags commercially available: tandem mass tags (TMT) and isobaric tags for relative and absolute quantitation (iTRAQ). The second type of chemical tags are isotopic. They are chemically identical, but differ in their mass due to incorporated stable isotopes. Examples of different isotopic tags for stable isotope labelling (SIL) are ICAT, SILAC, dimethylation, or heavy oxygen (<sup>18</sup>O). Incorporation of stable isotopes results in different peptide masses on MS1 level, which give rise to coeluting ion traces in the TIC with a mass difference typical for each different chemical tag.

This tutorial deals with protein quantitation via stable isotope labelling (SIL). We will use tools of the OpenMS suite. The tutorial covers relative quantitation only (i.e. comparison of abundances in different samples, no absolute quantitation of peptides / proteins).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Importing and Converting Peptide and Protein IDs](#importing-ids)
> 2. [MS1 Feature Detection](#ms1-feature-detection)
> 3. [Quant to ID matching](#quant-to-id-matching)

# Importing and Converting Peptide and Protein IDs

# MS1 Feature Detection

# Quant to ID matching
