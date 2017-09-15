---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: protein-id-oms
---

# Introduction
{:.no_toc}

Identifying the proteins contained in a sample is an important step in any proteomic experiment. However, in most settings, proteins are digested to peptides before the LC-MS/MS analysis. In this so-called "bottom-up" procedure, only peptide masses are measured. Therefore, protein identification cannot be performed directly from raw data, but is a multi-step process:

1. Raw data preparations
2. Peptide-to-Spectrum matching
3. Peptide inference
4. Protein inference

A plethora of different software solutions exists for each step. In this tutorial, we will use ***msconvert*** {% icon tool %}  for raw data conversion and tools from the [OpenMS software suite](https://openms.de) for all other steps. We will use one peptide search engine at first and later on show how to expand the workflow for using multiple search engines. Protein inference will be performed with the Fido algorithm ([Serang et al, JPR, (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20712337)).

For an alternative ID pipeline using the [Compomics](https://compomics.com/) tools [SearchGUI](https://compomics.github.io/projects/searchgui.html) and [PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html), please consult [this tutorial]({{site.url}}/topics/proteomics/tutorials/protein-id-sg-ps/tutorial.html).

# Input data
{:.no_toc}

As an example dataset, we will use an LC-MS/MS analysis of HeLa cell lysate published
in [Vaudel et al., 2014, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/24678044). Detailed information
about the dataset can be found on [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD000674).
For step 2 we will use a validated human Uniprot FASTA database with appended decoys.
If you already completed the tutorial on [Database Handling]({{site.url}}/topics/proteomics/tutorials/database-handling/tutorial.html)
you can use the constructed database after the **DecoyDatabase** {% icon tool %} step. You can find a prepared database, as well as the input proteomics data in different file formats on [Zenodo](https://zenodo.org/record/546301).

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
> This part of the Hands-On section is optional, because it cannot be performed on most GalaxyP instances due to licensing reasons. Therefore, we provide the [input data](https://zenodo.org/record/546301) also already converted to `.mzML` file format. If you choose to omit this part of the Hands-On section, please download the file "qExactive01819.mzml_profile" from [https://zenodo.org/record/546301/files/qExactive01819_profile.mzml] and continue with step 5 below (***PeakPickerHiRes*** {% icon tool %}).
>
> 1. Create a new history for this Peptide and Protein ID exercise.
> 2. Load the example dataset into your history from this [link](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/01/PXD000674/qExactive01819.raw).
> 3. Rename the dataset to "Test data".
> 4. Run ***msconvert*** {% icon tool %} on the test data to convert to the mzML format
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


# Peptide Identification
MS/MS experiments identify peptides by isolating them and subsequently colliding them with a gas for fragmentation. This method generates a spectrum of peptide fragment masses for each isolated peptide - an MS2 spectrum. 
To find out the peptide sequences, the MS2 spectrum is compared to a theoretical spectrum generated from a protein database. This step is called peptide-to-spectrum (also: spectrum-to-sequence) matching. Accordingly, a peptide that is successfully matched to a sequence is termed PSM (Peptide-Spectrum-Match). There can be multiple PSMs per peptide, if the peptide was fragmented several times. Different peptide search engines have been developed to fulfill the matching procedure. 
Here, we will use the search engine [X!Tandem](https://www.ncbi.nlm.nih.gov/pubmed/14976030).

> ### {% icon hands_on %} Hands-On: Peptide Identification
>
> 1. Copy the prepared protein database from the tutorial [Database Handling](../database-handling/tutorial.md) into your current history by using the multiple history view or upload the ready-made database from this [link](https://zenodo.org/record/892005/files/Human_database_including_decoys_%28cRAP_and_Mycoplasma_added%29.fasta).
> 2. Run the tool ***XTandemAdapter*** {% icon tool %} on the centroided mzML input file. Click `+ Insert param_fixed_modifications` and choose `Carbamidomethyl (C)`, then click `+ Insert param_variable_modifications` and choose `Oxidation (M)` and `Execute`.
> 3. FDR filtering
> 5. Run ***FileInfo*** {% icon tool %} to get basic information about the identified peptides.
>
>   > ### {% icon comment %} Comment: Advanced X!Tandem Parameters
>   > The search engine X!Tandem features many more advanced options than the ones reflected in the ***XTandemAdapter*** {% icon tool %}. If you need those advanced options, the ***XTandemAdapter*** {% icon tool %} allows for the optional input of a classic X!Tandem parameter file. Upload your parameter file to the history and use it as input in the field `Default X!Tandem configuration file`. You may set the option `-ignore_adapter_param` to `Yes` to overwrite all options set by the GUI.
>   {: .comment}
>
>   > ### {% icon question %} Questions:
>   > 1. How many peptides were identified?
>   > 2. How many peptides with oxidized methionine were identified?
>   >
>   >  <details>
>   >  <summary>Click to view answers</summary>
>   >   <ol type="1">
>   >     <li> You should have identified 3,325 peptides and 1,170 proteins.</li>
>   >     <li> 328 peptides contain an oxidized methionine (MeO). To get to this number, you can use ***Select*** {% icon tool %} on the Peptide Report and search for either "Oxidation of M" or "M\<ox\>".</li>
>   >   </ol>
>   >  </details>
>   {: .question}
{: .hands_on}

# Using Multiple Peptide Search Engines
It is generally recommended to use more than one peptide search engine and use the combined results for the final peptide inference ([Shteynberg et al., 2013, Mol. Cell. Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23720762)). Again, there are several software solutions for this, e.g. iProphet (TPP) or SearchGUI. In this tutorial we will use the OpenMS tool ***ConsensusID*** {% icon tool %}.

# Protein Inference
In bottom-up proteomics, it is necessary to combine the identified peptides to proteins. This is not a trivial task, as proteins are redundant to some degree. Thus, not every peptide can be assigned to only one protein. 
The OpenMS suite implemented the [Fido](https://www.ncbi.nlm.nih.gov/pubmed/20712337) algorithm for protein inference. Fido uses a Bayesian probabilistic model to group and score proteins based on peptide-spectrum matches.

# FDR Filtering

# Analysis of Contaminants
The FASTA database used for the peptide to spectrum matching contained some entries that were not expected to stem from the HeLa cell lysate, but are common contaminations in LC-MS/MS samples. The main reason to add those is to avoid false assignment of the spectra to other proteins. However, it also enables you to check for contaminations in your samples. **CAVE:** in human samples, many proteins that are common contaminants may also stem from the real sample. The real source can never be determined for sure.

> ### {% icon hands_on %} Hands-On: Analysis of Contaminants
>
> 1. Run ***Select*** {% icon tool %} on the Peptide Shaker Protein Report to select all lines that match the pattern "CONTAMINANT".
> 2. Remove all contaminants from your protein list by running ***Select*** {% icon tool %} on the Peptide Shaker Protein Report. Select only those lines that **do not** match the pattern "CONTAMINANT".
>
>   > ### {% icon question %} Questions:
>   > 1. Which contaminants did you identify? Where do these contaminations come from?
>   > 2. What other sources of contaminants exist?
>   > 2. How many mycoplasma proteins did you identify? Does this mean that the analyzed HeLa cells were infected with mycoplasma?
>   >
>   >  <details>
>   >  <summary>Click to view answers</summary>
>   >   <ol type="1">
>   >       <li> TRY_BOVIN is bovine trypsin. It was used to degrade the proteins to peptides. ALBU_BOVIN is bovine serum albumin. It is added to cell culture medium in high amounts.</li>
>   >       <li> Contaminants often stem from the experimenter, these are typically keratins or other high-abundant human proteins. Basically any protein present in the room of the mass spectrometer might get into the ion source, if it is airborne. As an example, sheep keratins are sometimes found in proteomic samples, stemming from clothing made of sheep wool.</li>
> >   <li> There should be five mycoplasma proteins in your protein list. However, all of them stem from different mycoplasma species. Also, every protein was identified by one peptide only. You can see this in column 17-19 of your output. These observations make it very likely that we are facing false positives here. As we were allowing for a false discovery rate of 1 %, we would expect 12 false positive proteins in our list. False positives are distributed to random peptides in the FASTA database. Our database consists of about 20,000 human proteins and 4,000 mycoplasma proteins. Therefore, we would expect 20 % of all false positives to match to mycoplasma proteins.</li>
>   >   </ol>
>   >  </details>
>   {: .question}
{: .hands_on}

# Evaluation of Peptide and Protein IDs
***Peptide Shaker*** {% icon tool %} provides you with validation results for the identified PSM, peptides and proteins. It classifies all these IDs in the categories "Confident" or "Doubtful". On each level, the meaning of these terms differs to some extent:

- **PSMs** are marked as "Doubtful" when the measured MS2 spectrum did not fit well to the theoretical spectrum.
- **Peptides** have a combined scoring of their PSMs. They are marked as "Doubtful", when the score is below a set threshold. The threshold is defined by the false discovery rate (FDR).
- **Proteins** are marked as "Doubtful", when they were identified by only a single peptide or when they were identified solely by "Doubtful" peptides.

> ### {% icon hands_on %} Hands-On: Evaluation of Peptide and Protein IDs
>
> 1. Remove all "Doubtful" proteins from your protein list by running ***Select*** {% icon tool %} on the Peptide Shaker Protein Report. Select only those lines that **do not** match the pattern "Doubtful".
>
>   > ### {% icon question %} Questions:
>   > 1. How to exclude mycoplasma proteins?
>   > 2. How many "Confident" non-contaminant proteins were identified?
>   >
>   >  <details>
>   >  <summary>Click to view answers</summary>
>   >   <ol type="1">
>   >       <li> Add another ***Select*** {% icon tool %} matching the pattern "HUMAN".</li>
>   >   <li> You should have identified 582 human non-contaminant proteins that were validated to be "Confident".<\li>
>   >   </ol>
>   >  </details>
>   {: .question}
{: .hands_on}

# Premade Workflow

A premade workflow for this tutorial can be found [here](workflows/wf_proteinID_SG_PS.ga)

# Further Reading

- [Search GUI and Peptide Shaker tutorials at Compomics](https://compomics.com/bioinformatics-for-proteomics/)
- [Using Search GUI and Peptide Shaker in Galaxy](https://drive.google.com/file/d/0B6bIeOvjBkbWVnhMLWxXdGVUY3M/view)
