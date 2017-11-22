---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: protein-quant-sil
---

# Introduction
{:.no_toc}

To compare protein amounts in different samples from MS/MS data, two different experiment setups exist. Firstly, unmodified proteins can be measured in separate runs at one sample per MS-run. Secondly, proteins of samples to compare can be labelled with small chemical tags, mixed, and measured side-by-side in a single MS-run.
There are two types of chemical tags:
  1. Isobaric tags display the same mass on first hand, but fragment during the generation of the MS/MS spectra to yield reporter ions of different mass. The intensity of those reporter ions can be compared in MS/MS spectra. There are two types of isobaric tags commercially available: *tandem mass tags* (TMT) and *isobaric tags for relative and absolute quantitation* (iTRAQ).
  2. Isotopic tags are chemically identical, but differ in their mass due to incorporated stable isotopes. Examples of different isotopic tags for stable isotope labelling (SIL) are ICAT, SILAC, dimethylation, or heavy oxygen (<sup>18</sup>O).

This tutorial deals with protein quantitation via stable isotope labelling (SIL). For isotopic tags, quantitation can be achieved by comparing the intensity of MS1 peptide mass traces. The whole MS1 profile of a peptide, i.e. the intensities of all its isotopic peaks over time, is called a *peptide feature* (Figure 1a). Incorporation of stable isotopes results in different peptide masses on MS1 level, which give rise to coeluting ion traces in the TIC with a mass difference typical for each different chemical tag (Figure 1b). Figure originally published in [Nilse et al, 2015](http://www.ncbi.nlm.nih.gov/pubmed/25931027).

![ms1 feature](../../images/protein-quant-sil_ms1feature.png "MS1 mass traces. A) Two peptide features of co-eluting SIL peptides. B) MS1 spectra at a given RT. C) XIC monoisotopic peak light peptide. D) XIC monoisotopic peak heavy peptide.")

## Prerequisites
{:.no_toc}

If you still are in the planning phase of your quantitative proteomics experiment, you may want to consider our tutorial on different [quantitation methods]({{site.url}}/topics/proteomics/tutorials/labelfree-vs-labelled/tutorial.html) first.

To learn about *protein identification* in Galaxy, please consider [this tutorial]({{site.url}}/topics/proteomics/tutorials/protein-id-oms/tutorial.html).

> ### {% icon hands_on %} Hands-on: Introduction
> In the hands-on section of this tutorial, we will use a quantitative comparison of HEK cell lysate as a test dataset. In this experiment, HEK cells were once labelled with light, once with heavy SILAC. Both cultures were lysed simultaneously and the cell lysates were mixed in a certain ratio. For a detailed description of the dataset, please refer to the description in the [PRIDE archive]().
>
> Your objective in this hands-on-tutorial is to find out the correct mixing ratio of the test sample.
{: .hands_on}

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

> ### {% icon hands_on %} Hands-on: MS1 Feature Detection
>
> 1. Import the test data from [zenodo](https://zenodo.org/record/1051552). The file type of the data is mzML.
[//]: # The data have not been modified during the conversion from the machine raw file, neither background removal, nor peak picking (centroiding) has been performed. **Is this still true? I believe not, check!**
> 2. Run ***FeatureFinderMultiplex*** {% icon tool %} with
>   - the mzML file as **LC-MS dataset in centroid or profile mode**, and
>   - **Labels used for labelling the samples** to `[ ][Arg6,Lys6]`.
{: .hands_on}

# Peptide and Protein Identification

In this tutorial, peptide identification will be performed using the workflow of the previous [Peptide ID Tutorial]({{site.url}}/topics/proteomics/tutorials/protein-id-oms/tutorial.html). 

A common problem in mass spectrometry are misassigned mono-isotopic precursor peaks. Although most search engines allow for some adaptation of the monoisotopic peak, we will instead perform a recalculation of the monoisotopic peaks based on the previously identified features prior to peptide identification.
This step facilitates mapping peptide IDs to identified features [later on](#mapping-features-to-ids). To do so, we will use the OpenMS tool ***HighResPrecursorMassCorrector*** {% icon tool %}.

[//]: # TODO: Read about monoisotopic peak problem, give citation to review!

> ### {% icon hands_on %} Hands-on: Peptide and Protein Identification and Conversion
> 1. Run ***HighResPrecursorMassCorrector*** {% icon tool %} with
>   - the `mzML` file as **Input file**,
>   - and the output of ***FeatureFinderMultiplex*** as **Features used to correct precursor masses**.
> 2. Import the [Protein identification using OpenMS tutorial workflow]({{site.url}}/topics/proteomics/tutorials/protein-id-oms/workflows/workflow.ga) and modify it:
>   - Delete the **PeakPickerHiRes** {% icon tool %} node, as the MS2 data of our test dataset are already centroided.
>   - Connect the `mzML` input directly to the **XTandemAdapter** {% icon tool %} node.
>   - Change the **XTandemAdapter** {% icon tool %} parameters:
>       - Add the variable modifications `Label:13C(6) (K)` and `Label:13C(6) (R)`.
> 3. Run the workflow with
>   - the output of ***HighResPrecursorMassCorrector*** as `1: Input: mzML dataset`
>   - the human FASTA database as `2: Human FASTA database including decoys`
>
>   > ### {% icon tip %} Tip: Using Galaxy Workflows
>   > If you want to learn more about Galaxy workflows, please consult the [Galaxy Introduction]({{site.url}}/topics/introduction/tutorials/galaxy-intro-101/tutorial.html#the-workflow-editor)
>   {: .tip}
{: .hands_on}

# Quant to ID matching

We now have feature quantifications for MS1 elution peaks, peptide identifications for the MS2 spectra and protein identifications. 
The next step is to map the MS2-based peptide identifications to the quantified MS1 precursor peaks ("peptide features"). This will enable the quantification of identified peptides. 

Sometimes several peptide identifications are mapped to a feature. The tool [IDConflictResolver](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_IDConflictResolver.html) filters the mapping so that only the identification with the best score is associated to each feature.

Finally, we will combine the peptide quantifications to protein quantifications.

> ### {% icon hands_on %} Hands-on: Quant to ID matching
>
> 1. Run ***IDMapper*** {% icon tool %} with
>   - the output of ***IDFilter*** as **Protein/peptide identifications file**
>   - the `consensusXML` output of ***FidoAdapter*** as **Feature map/consensus map file**
>   - **Match using RT and m/z of sub-features instead of consensus RT and m/z** set to `Yes`.
> 2. Change the filetype of the ***IDMapper*** output to `consensusXML`.
> 3. Run ***FileFilter*** {% icon tool %} with
>   - **Remove features without annotations** set to `Yes`, and
>   - **Remove unassigned peptide identifications** set to `Yes`.
> 4. Run ***IDConflictResolver*** {% icon tool %}.
> 5. Run ***ProteinQuantifier*** {% icon tool %} with
>   - the output of ***IDConflictResolver*** as **Input file**,
>   - the output of ***IDFilter*** as **Protein inference results [...]**,
>   - **Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all)** set to `0`, 
>   - **Averaging method used to compute protein abundances from peptide abundances** set to `sum`, and
>   - **Add the log2 ratios of the abundance values to the output** set to `Yes`.
>
>   > ### {% icon question %} Questions
>   > 1. How many proteins were successfully quantified?
>   {: .question}
>
>   > ### {% icon comment %} Comment: ProteinQuantifier parameters
>   > Peptide quantitation algorithms are more precise for high abundant peptides. Therefore, it is recommended to base protein quantitations on those peptides. In ProteinQuantifier, you may restrict the calculation of protein abundances to the most abundant peptides by using the option "Calculate protein abundance from this number of proteotypic peptides".
>   > However, we recommend to use the averaging method `sum` instead. By using this option, protein ratios are based on the sum of all peptide abundances. Thus, highly abundant peptides thus have more influence on protein abundance calculation than low abundant peptides. 
>   > A simple Sum-of-Intensities algorithm provided the best estimates of true protein ratios in a comparison of several protein quantitation algorithms ([Carrillo et al., Bioinformatics, 2009](https://www.ncbi.nlm.nih.gov/pubmed/19892804)).
>   {: .comment}
{: .hands_on}

# Evaluation and Optimization of Quantitation Results

Protein quantitation is a multi-step procedure. Many parameters of different steps influence the final results. Therefore, it is recommended to optimize the tool parameters for each dataset and to carefully evaluate quantitation results. While the total number of quantified proteins is a first important parameter for optimization, it is also necessary to visualize the results and check for correct feature finding and ID mapping.

Galaxy does not feature a tool for proteomics visualization, we recommend to use the OpenMS Viewer [TOPPView](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_TOPPView.html). Basic TOPPView tutorials are available as [videos](https://www.openms.de/getting-started/command-line-and-visualisations/) and a more comprehensive tutorial as [PDF](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/TOPP_tutorial.pdf).

For the optimization of tool parameters, it is recommended not to work with a complete LC-MS/MS run. Instead, we will use ***FileFilter*** to extract a small *RT-slice* of our input dataset, i.e. a fraction of the original dataset that was measured during a short period of time. Reducing the test data reduces the time needed for analysis and facilitates visual examination of the data.
Be aware that only very small parts of your data can be checked by visual examination. To minimize biases, try to look at the same areas / features of each result file.

> ### {% icon hands_on %} Hands-on: Data reduction and visual evaluation with TOPPView
>
> 1. Run ***FileFilter*** {% icon tool %} with
>   - **Retention time range to extract** set to `2000:2200`.

## Critical parameters for optimization

[//]: # **Three steps to optimize: 1) FFMult, 2) HighResPrecMassCorr, 3) IdMapper**
[//]: # **Test if HiResMassCorr can be optimized or if it only has litte impact**

Three steps are critical for protein quantitations, the most relevant tool parameters are listed below:
1. Feature finding by ***FeatureFinderMultiplex***:
    - **Typical retention time [s] over which a characteristic peptide elutes**: To improve results, you may look at the mzML file first and find out the average elution time of peaks. **How?**
    - **Lower bound for the intensity of isotopic peaks**:
    - **Lower bound for the intensity of isotopic peaks**: If you already used an intensity cutoff during mzML preprocessing, set this parameter to `0`.
    - **Two peptides in a multiplet are expected to have the same isotopic pattern**: The grade of similarity between isotopic patterns of the light and the heavy peptide.
    - **The isotopic pattern of a peptide should resemble the averagine model at this m/z position**: The grade of similarity between the measured and the theoretical isotopic pattern. The theoretical pattern is estimated by the averagine model based upon peptide mass, as the peptide sequences are unknown.
2. Precursor corrections by ***HighResPrecMassCorr***:
    - **The precursor mass tolerance**
    - **Additional retention time tolerance added to feature boundaries**
3. Mapping of peptide IDs to features by ***IDMapper***:
    - **RT tolerance (in seconds) for the matching of peptide identifications and (consensus) features**:
    - **m/z tolerance (in ppm or Da) for the matching of peptide identifications and (consensus) features**: 

For optimization, it is critical to change **only one parameter at a time**. Also, it is recommended to optimize the tools in the order of their position in the workflow.

> ### {% icon hands_on %} Hands-on: Optimization of Quantitation Results
>
> 1. Extract a workflow out of your history or import the [premade workflow](./workflows/workflow.ga).
> 2. Run the whole WF again, change **a single setting (averagine?)** in ***FeatureFinderMultiplex*** {% icon tool %}.
> 2. Run ***FileInfo*** {% icon tool %} on the results -> number of ID-Feature-matches
> 3. Open results in TOPPView.
>
>   > ### {% icon question %} Questions
>   > 1. Which setting led to more ID-Feature-matches?
>   > 2. Using the default settings, how many features were not mapped to IDs? How many IDs were not mapped to features?
>   {: .question}
{: .hands_on}
