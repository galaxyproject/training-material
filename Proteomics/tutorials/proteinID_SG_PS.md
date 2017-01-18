---
layout: tutorial_hands_on
topic_name: Peptide and Protein ID
tutorial_name: proteinID_SG_PS
---

# Introduction

Identifying the proteins contained in a sample is an important step in any proteomic experiment. Protein identification cannot be performed directly from raw data, but is a multi-step process: 

1. Raw data preparations
2. MS2 spectrum-to-sequence matching 
3. Peptide inference
4. Protein inference

A plethora of different software solutions exists for each step. In this tutorial, we will show how to use ***msconvert*** :wrench: and ***PeakPickerHiRes*** :wrench: for step 1, ***Search GUI*** :wrench: and ***Peptide Shaker*** :wrench: for the steps 2-4.

As example data, we will use a cell lysate of the HeLa cell line [Vaudel et al., 2014, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/24678044) . Detailed information about the dataset can be found on [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD001907) . For step 2 we will search against the database constructed in the Database Handling tutorial.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Preparing raw data](#Preparing-raw-data)
> 2. [Spectrum-to-Sequence matching](#Peptide-to-Spectrum-matching) 
> 3. [Peptide and Protein Inference](#Peptide-and-Protein-Inference)

# Preparing raw data

Raw data conversion is the first step of any proteomic data analysis. The most common converter is MSConvert, the format to convert to is mzML. SearchGUI takes only mgf format as input, but as we need the mzML format for several other tasks, we will convert to mzML first.

> ### :pencil2: Hands-on: Preparing raw data
>
> 1. Create a new history for this Peptide and Protein ID exercise.
> 2. Load the example dataset into your history from this link: ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/01/PXD000674/qExactive01819.raw
> 3. Rename the dataset to "Test data".
> 4. Run ***msconvert*** :wrench: on the test data to convert to the mzML format
>
>	> ### :nut_and_bolt: Comment: Local Use of MSConvert
>	> The vendor libraries used by MSConvert need a Windows system and is therefore hard to implement in Galaxy. If ***msconvert*** :wrench: is not available in your Galaxy instance, please install the software on a Windows computer and run the conversion locally. You can find a detailed description of the necessary steps (here)[http://genesis.ugent.be/files/costore/practicals/bioinformatics-for-proteomics/1-Peptide-and-Protein-Identification/1.2-Peak-List-Generation/1.2-Peak-List-Generation.pdf] . Afterwards, upload the resulting mzML file to your Galaxy history.
>
> 5. Run ***PeakPickerHiRes*** :wrench: on the resulting mzML file.
> 6. Run ***FileConverter*** :wrench: on the picked mzML to convert to mgf format.

# Spectrum-to-Sequence matching
MS/MS experiments identify peptides by isolating them and subsequently colliding them with a gas for fragmentation. This method generates a spectrum of peptide fragment masses for each isolated peptide - an MS2 spectrum. To find out the sequence of the unfragmented peptide, the MS2 spectrum is compared to a theoretical spectrum generated from a protein database. This step is called spectrum-to-sequence matching. Different peptide search engines have been developed to fulfill the matching procedure. 

It is generally recommended to use more than one peptide search engine and use the combined results for the final peptide inference (Shteynberg et al., 2013, Mol. Cell. Proteomics)[https://www.ncbi.nlm.nih.gov/pubmed/23720762] . Again, there are several software solutions for this step, e.g. iProphet (TPP) or ConsensusID (OpenMS). In this tutorial, we will use the combination of SearchGUI to perform and combine the multiple searches and PeptideShaker for the peptide inference.

> ### :pencil2: Hands-on: Spectrum-to-Sequence matching
>
> 1. Open ***Search GUI*** :wrench: 