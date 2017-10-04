---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: protein-id-oms-consensus
---

# Introduction
{:.no_toc}

Identifying the proteins contained in a sample is an important step in any proteomic experiment. However, in most settings, proteins are digested to peptides before the LC-MS/MS analysis. In this so-called "bottom-up" procedure, only peptide masses are measured. Therefore, protein identification cannot be performed directly from raw data, but is a multi-step process:

1. Raw data preparations
2. Peptide-to-Spectrum matching
3. Peptide inference
4. Protein inference

It is generally recommended to use more than one peptide search engine and use the combined results for the final peptide inference ([Shteynberg et al., 2013, Mol. Cell. Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23720762)). Again, there are several software solutions for this, e.g. iProphet (TPP) or SearchGUI. In this tutorial we will use the OpenMS tool ***ConsensusID*** {% icon tool %}.
This tutorial covers peptide and protein **identification** only, but you may use the output of this tutorial for the [tutorial on protein quantitation]({{site.url}}/topics/proteomics/tutorials/protein-quant-sil/tutorial.html).

For an alternative ID pipeline using the [Compomics](https://compomics.com/) tools [SearchGUI](https://compomics.github.io/projects/searchgui.html) and [PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html), please consult [this tutorial]({{site.url}}/topics/proteomics/tutorials/protein-id-sg-ps/tutorial.html). However, the latter tutorial does not allow to continue with the tutorial on protein quantitation.

# Input data
{:.no_toc}

As an example dataset, we will use an LC-MS/MS analysis of HeLa cell lysate published
in [Vaudel et al., 2014, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/24678044). Detailed information
about the dataset can be found on [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD000674).
For step 2 we will use a validated human Uniprot FASTA database with appended decoys.
If you already completed the tutorial on [Database Handling]({{site.url}}/topics/proteomics/tutorials/database-handling/tutorial.html)
you can use the constructed database after the **DecoyDatabase** {% icon tool %} step. You can find a prepared database, as well as the input proteomics data in different file formats on [Zenodo](https://zenodo.org/record/796184).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparing raw data

Raw data conversion is the first step of any proteomic data analysis. The most common converter is MSConvert from the [ProteoWizard software suite](http://proteowizard.sourceforge.net/), the format to convert to is mzML.

> ### {% icon hands_on %} Optional Hands-On: Preparing raw data
>
> This part of the Hands-On section is optional, because it cannot be performed on most GalaxyP instances due to licensing reasons. Therefore, we provide the [input data](https://zenodo.org/record/796184) also already converted to `.mzML` file format. If you choose to omit this part of the Hands-On section, please download the file "qExactive01819_vendor-peakPicking.mzml" from [https://zenodo.org/record/546301/files/qExactive01819_vendor-peakPicking.mzml] and directly proceed to [peptide identification](peptide-identification).
>
> 1. Create a new history for this Peptide and Protein ID exercise.
> 2. Load the example dataset into your history from this [link](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/01/PXD000674/qExactive01819.raw).
> 3. Rename the dataset to "Test data".
> 4. Run ***msconvert*** {% icon tool %} on the test data to convert to the mzML format.
>
>   > ### {% icon comment %} Comment: Local Use of MSConvert
>   > The vendor libraries used by MSConvert are only licensed for Windows systems and are therefore rarely implemented in Galaxy instances. If ***msconvert*** {% icon tool %} is not available in your Galaxy instance, please install the software on a Windows computer and run the conversion locally. You can find a detailed description of the necessary steps [here](http://genesis.ugent.be/files/costore/practicals/bioinformatics-for-proteomics/1-Peptide-and-Protein-Identification/1.2-Peak-List-Generation/1.2-Peak-List-Generation.pdf). Afterwards, upload the resulting mzML file to your Galaxy history.
>  {: .comment}
>
> 5. Run ***PeakPickerHiRes*** {% icon tool %} on the resulting mzML file. Use the option `+ Insert param.algorithm_ms_levels` and change the entry to "2". Thus, peak picking will only be performed on MS2 level.
>
>   > ### {% icon comment %} Comment: Peak Picking
>   > Depending on your machine settings, raw data will be generated either in profile mode or centroid mode. For most peptide search engines, the MS2 data have to be converted to centroid mode, a process called "peak picking" or "centroiding". 
>   > Machine vendors offer algorithms to extract peaks from profile raw data. This is implemented in ***msconvert*** {% icon tool %} and can be run in parallel to the mzML conversion. However, the OpenMS tool ***PeakPickerHiRes*** {% icon tool %} is reported to generate better results ([Lange et al., 2006, Pac Symp Biocomput](https://www.ncbi.nlm.nih.gov/pubmed/17094243)) and is therefore recommended for quantitative studies ([Vaudel et al., 2010, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/19953549)).
>   >
>   > While most peptide search engines depend on centroided MS2 peaks, quantitation algorithms used by OpenMS typically use non-centroided MS1 data. Therefore, peak picking is only performed on MS2 level in this tutorial. If you need centroided MS1 and MS2 level, rerun ***PeakPickerHiRes*** {% icon tool %} with standard settings on the output of step 5.
>   >
>   > If your data were generated on a low resolution mass spectrometer, use ***PeakPickerWavelet*** {% icon tool %} instead.
>   {: .comment}
{: .hands_on}

