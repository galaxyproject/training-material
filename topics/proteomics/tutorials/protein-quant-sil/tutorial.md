---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: protein-quant-sil
---

# Introduction
{:.no_toc}

To compare protein amounts in different samples from MS/MS data, two different experiment setups exist. Firstly, unmodified proteins can be measured in separate runs at one sample per MS-run. Secondly, proteins of samples to compare can be labelled with small chemical tags, mixed, and measured side-by-side in a single MS-run.
There are two types of chemical tags: isobaric tags display the same mass on first hand, but fragment during the generation of the MS/MS spectra to yield reporter ions of different mass. The intensity of those reporter ions can be compared in MS/MS spectra. There are two types of isobaric tags commercially available: tandem mass tags (TMT) and isobaric tags for relative and absolute quantitation (iTRAQ).
The second type of chemical tags are isotopic. They are chemically identical, but differ in their mass due to incorporated stable isotopes. Examples of different isotopic tags for stable isotope labelling (SIL) are ICAT, SILAC, dimethylation, or heavy oxygen (<sup>18</sup>O).
Incorporation of stable isotopes results in different peptide masses on MS1 level, which give rise to coeluting ion traces in the TIC with a mass difference typical for each different chemical tag.

This tutorial deals with protein quantitation via stable isotope labelling (SIL). We will use tools of the OpenMS suite.  
Because we solely cover *quantitation*, you need to perform peptide and protein ID in beforehand. To learn about protein ID in Galaxy, please consider [this tutorial](./proteinID_SG_PS.md).

This tutorial covers *relative* quantitation only (i.e. comparison of abundances in different samples, no *absolute* quantitation of peptides / proteins).

If you still are in the planning phase of your quantitative proteomics experiment, you may want to consider our tutorial on different [quantitation methods](./labelfree-vs-labelled.md) first.

> ### :pencil2: Hands-on: Introduction
> In the Hands-on section of this tutorial, we will use a quantitative comparison of **HEK _OR_ E.coli** cell lysate as a test dataset. In this experiment, the very same cell lysate was once labelled with light, once with heavy **dimethyl _OR_ SILAC** and both samples were subsequently mixed in a certain ratio. For a detailed description of the dataset, please refer to the description in the [PRIDE archive]().
> Your objective in this hands-on-tutorial is to find out the correct mixing ratio of the test sample.
> {: .hands_on}

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# MS1 Feature Detection
Quantitation on MS1 level may in principle be carried out without prior knowledge of peptide / protein IDs. However, some quantitation algorithms take the IDs as an input to make sure that every PSM that was identified is also quantified. This is not the case in our example here.
In the OpenMS suite, most of the provided tools for MS1 feature detection quantify solely based on mzML files. The advantage of this approach is that quantitations can be made on strict criteria to reduce misquantitations. The drawback is that not all IDs can be matched to a quantitation later on in the workflow.
The tool settings need to be carefully tested and evaluated manually to obtain optimal results. We will explain this in the section [Evaluation and Optimization of Quantitation Results](#evaluation-and-optimization-of-quantitation-results).

> ### :pencil2: Hands-on: MS1 Feature Detection
>
> 1. Import the test data from [zenodo](). The file type of the data is mzML. The data have not been modified during the conversion from the machine raw file, neither background removal, nor peak picking (centroiding) has been performed.
> 2. Run `FeatureFinderMultiplex`:wrench: on the mzML file. Change **`Labelling`** to `\[ \] \[Arg6,Lys6\]`.
>
>   > ### :bulb: Tip (Expert level): Detecting features of knockouts
>   > In biology, there are rarely cases in which a gene product is completely shut off between two conditions. Rather, most changes are **word missing**. However, in some situations, you will have the situation that a protein is present in only one of the tested conditions and completely lacking in another. A classical example would be comparing a "knockout" mouse with its "wild-type" counterpart.
>   > Due to the feature detection algorithm of `FeatureFinderMultiplex`, those features would normally be disregarded, as they do not look like typical features in labelled samples.
>   >
>   > However, there is a built-in option in `FeatureFinderMultiplex` that enables finding of "knockout features". If you expect one or more proteins to be completely missing in at least one of your conditions, select the advanced option **`knockouts present` **.
>   > Switching on this option is not recommended as a default setting, as it increases the possibility of false positives. When using this option, be advised to check for false positives carefully, as described [below](#expert-level-evaluation-and-optimization-of-quantitation-results).
> {: .hands_on}

# Peptide and Protein Identification and Conversion

> ### :pencil2: Hands-on: Peptide and Protein Identification and Conversion
> 1. Run the workflow "ProteinID_SG_PS" on the test dataset.
> 2. Use **`IDConverter`** to convert the mzid output of `Peptide Shaker` :wrench: to mzidentML.

# Quant to ID matching

> ### Hands-on: Quant to ID matching
>
> 1. Run `ProteinQuantifier` :wrench: on
> 2.
>
>   > ### :question: Questions

>   > {: .questions}
> {: .hands_on}

# Expert level: Evaluation and Optimization of Quantitation Results
`FeatureFinderMultiplex` searches for multiple similar features that elute at the same time, but diverge by a mass shift fitting to the label used. `FeatureFinderMultiplex` uses several parameters that may be used to optimize your search results. Two important parameters are **`Average elution time` and `Averagine similarity`**.
    - `Average elution time`: To improve results, you may look at the mzML file first and find out the average elution time of peaks. **How?**
    - `Averagine similarity`: describes the similarity of two features. Play around with this parameter to optimize the number of features detected. Be careful, reducing it may lead to false positives.

> ### Comment: Benchmarking parameters for opimization - What is a good result?
>
> The quality of the results can be measured by so-called "benchmarking parameters".
>  
> Benchmarking parameters for `FeatureFinderMultiplex`:
>     1. Number of features that can be linked to peptide IDs (and vice-versa).
>     2. Although the first parameter gives a good measure of quality, it does not rule out that false positive features or IDs were matched. To check for false positives, you will have to scan through your data manually. To do so, open both the `FeatureFinderMultiplex` consensusXML output and the mzidentML file into [TOPPView]().
>         - **Caution** Manual evaluation is prone to biases, as you can look solely at small parts of your data. To avoid this, try to look at the *very same* areas / the same features of all different result files.
>     3. If you were using the option **`knockouts present` **, check, if the detected "knockout features" fit to your expectations.

> ### Hands-on: Evaluation and Optimization of Quantitation Results
>
> 1. Run the whole WF again, change ** a single setting (averagine?)** in `FeatureFinderMultiplex`.
> 2. Run `FileInfo` :wrench: on the results -> number of ID-Feature-matches
> 3. Run `dongs` :wrench: on the results -> restrict to a small areas
> 4. Open results in TOPPView.
>
>   > ### :question: Questions
>   > 1. Which setting led to more ID-Feature-matches?
>   > 2. Using the default settings, how many features were not mapped to IDs? How many IDs were not mapped to features?
>   > {: .questions}
> {: .hands_on}
