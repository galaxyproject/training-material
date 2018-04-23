---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: protein-id-sg-ps
---

# Introduction
{:.no_toc}

Identifying the proteins contained in a sample is an important step in any proteomic experiment. However, in most settings, proteins are digested to peptides before the LC-MS/MS analysis. In this so-called "bottom-up" procedure, only peptide masses are measured. Therefore, protein identification cannot be performed directly from raw data, but is a multi-step process:

1. Raw data preparations
2. Peptide-to-Spectrum matching
3. Peptide inference
4. Protein inference

A plethora of different software solutions exists for each step. In this tutorial, we will show how to
use the [ProteoWizard](http://proteowizard.sourceforge.net/) tool MSconvert and the [OpenMS](https://openms.de) tool [PeakPickerHiRes](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_PeakPickerHiRes.html) for step 1, and the [Compomics](https://compomics.com/) tools [SearchGUI](https://compomics.github.io/projects/searchgui.html) and [PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html), for the steps 2-4.

For an alternative ID pipeline using only tools provided by the [OpenMS software suite](https://openms.de), please consult [this tutorial]({{site.baseurl}}/topics/proteomics/tutorials/protein-id-oms/tutorial.html).

# Input data
{:.no_toc}

As an example dataset, we will use an LC-MS/MS analysis of HeLa cell lysate published
in [Vaudel et al., 2014, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/24678044). Detailed information
about the dataset can be found on [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD000674).
For step 2 we will use a validated human Uniprot FASTA database without appended decoys.
If you already completed the tutorial on [Database Handling]({{site.baseurl}}/topics/proteomics/tutorials/database-handling/tutorial.html)
you can use the constructed database before the **DecoyDatabase** {% icon tool %} step. You can find a prepared database, as well as the input proteomics data in different file formats on [Zenodo](https://zenodo.org/record/796184).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Preparing Raw Data

Raw data conversion is the first step of any proteomic data analysis. The most common converter is MSConvert from the [ProteoWizard software suite](http://proteowizard.sourceforge.net/), the format to convert to is mzML. Search GUI needs `mgf` format as input, but as we need the `mzML` format for several other tasks, we will convert to `mzML` first. Due to licensing reasons, MSConvert runs only on windows systems and will not work on most Galaxy servers.

Depending on your machine settings, raw data will be generated either in profile mode or centroid mode. For most peptide search engines, the MS2 data have to be converted to centroid mode, a process called "peak picking" or "centroiding". 
Machine vendors offer algorithms to extract peaks from profile raw data. This is implemented in ***msconvert*** {% icon tool %} and can be run in parallel to the mzML conversion. However, the OpenMS tool ***PeakPickerHiRes*** {% icon tool %} is reported to generate slightly better results ([Lange et al., 2006, Pac Symp Biocomput](https://www.ncbi.nlm.nih.gov/pubmed/17094243)) and is therefore recommended for quantitative studies ([Vaudel et al., 2010, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/19953549)).
If your data were generated on a low resolution mass spectrometer, use ***PeakPickerWavelet*** {% icon tool %} instead.

> ### {% icon hands_on %} Hands-On: File Conversion and Peak Picking
>
> We provide the [input data](https://zenodo.org/record/796184) in the original `raw` format and also already converted to `mgf` and `mzML` file formats. If ***msconvert*** {% icon tool %} does not run on your Galaxy instance, please download the preconverted `mzML` as an input.
>
> 1. Create a new history for this Peptide and Protein ID exercise.
> 2. Load the example dataset into your history from Zenodo: [raw](https://zenodo.org/record/892005/files/qExactive01819.raw) [mzML](https://zenodo.org/record/892005/files/qExactive01819_profile.mzml)
> 3. Rename the dataset to something meaningful.
> 4. (*optional*) Run ***msconvert*** {% icon tool %} on the test data to convert to the `mzML` format.
> 5. Run ***PeakPickerHiRes*** {% icon tool %} on the resulting file. Click `+ Insert param.algorithm_ms_levels` and change the entry to "2". Thus, peak picking will only be performed on MS2 level.
> 6. Run ***FileConverter*** {% icon tool %} on the picked mzML. In the **Advanced Options** set the **Output file type** to `mgf`.
>
>   > ### {% icon comment %} Comment: Local Use of MSConvert
>   > The vendor libraries used by MSConvert are only licensed for Windows systems and are therefore rarely implemented in Galaxy instances. If ***msconvert*** {% icon tool %} is not available in your Galaxy instance, please install the software on a Windows computer and run the conversion locally. You can find a detailed description of the necessary steps [here](https://compomics.com/bioinformatics-for-proteomics/identification/) ("Peak List Generation"). Afterwards, upload the resulting mzML file to your Galaxy history.
>  {: .comment}
{: .hands_on}

# Peptide and Protein Identification
MS/MS experiments identify peptides by isolating them and subsequently colliding them with a gas for fragmentation. This method generates a spectrum of peptide fragment masses for each isolated peptide - an MS2 spectrum. To find out the peptide sequences, the MS2 spectrum is compared to a theoretical spectrum generated from a protein database. This step is called peptide-to-spectrum (also: spectrum-to-sequence) matching. Accordingly, a peptide that is successfully matched to a sequence is termed PSM (Peptide-Spectrum-Match). There can be multiple PSMs per peptide, if the peptide was fragmented several times. Different peptide search engines have been developed to fulfill the matching procedure.

It is generally recommended to use more than one peptide search engine and use the combined results for the final peptide inference ([Shteynberg et al., 2013, Mol. Cell. Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23720762)). Again, there are several software solutions for this, e.g. iProphet (TPP) or ConsensusID (OpenMS). In this tutorial we will use ***Search GUI*** {% icon tool %}, as it can automatically search the data using several search engines. Its partner tool ***Peptide Shaker*** {% icon tool %} is then used to combine and evaluate the search engine results.

In bottom-up proteomics, it is necessary to combine the identified peptides to proteins. This is not a trivial task, as proteins are redundant to some degree. Thus, not every peptide can be assigned to only one protein. Luckily, the ***Peptide Shaker*** {% icon tool %} already takes care of protein inference and even gives us some information on validity of the protein IDs. We will discuss validation in a [later step](#evaluation-of-peptide-and-protein-ids) of this tutorial.

> ### {% icon hands_on %} Hands-On: Peptide and Protein Identification
>
> 1. Copy the prepared protein database from the tutorial [Database Handling](../database-handling/tutorial.html) into your current history by using the multiple history view or upload the ready-made database from this [link](https://zenodo.org/record/892005/files/Human_database_%28cRAP_and_Mycoplasma_added%29.fasta).
> 2. Open ***Search GUI*** {% icon tool %} to search the mgf file against the protein database. In the **`Search Engine Options`** select `X!Tandem` and `MS-GF+`. In the **`Protein Modification Options`** add the **`Fixed Modifications`**: `Carbamidomethylation of C` and the **`Variable Modifications`**: `Oxidation of M`.
> 3. Run ***Peptide Shaker*** {% icon tool %} on the Search GUI output. Enable the following outputs: `Zip File for import to Desktop App`, `mzidentML File`, `PSM Report`, `Peptide Report`, `Protein Report`.
>
>   > ### {% icon comment %} Comment: Search GUI Parameters
>   > We ran ***Search GUI*** {% icon tool %} with default settings. When you are processing files of a different experiment, you may need to adjust some of the parameters.
>   > **Search GUI** bundles numerous peptide search engines for matching MS/MS to peptide sequences within a database. In practice, using 2-3 different search engines offers high confidence while keeping analysis time reasonable.
>   > The **`Precursor Options`** have to be adjusted to the mass spectrometer which was used to generate the files. The default settings fit a  high resolution Orbitrap instrument.
>   > In the **`Advanced Options`** you may set much more detailed settings for each of the used search engines. When using X!Tandem, we recommend to switch off the advanced X!Tandem options **`Noise suppression`**, **`Quick Pyrolidone`** and **`Quick Acetyl`**. When using MSGF, we recommend to select the correct **`Instrument type`**.
>   {: .comment}
>
>   > ### {% icon comment %} Comment: PeptideShaker Outputs
>   > Peptide Shaker offers a variety of outputs. 
>   > The `Zip File for import to Desktop App` can be downloaded to view and evaluate the search results in the Peptide Shaker viewer ([Download](https://compomics.github.io/projects/peptide-shaker.html)). 
>   > The several `Reports` contain tabular, human-readable information. 
>   > Also, an `mzidentML` (= `mzid`) file can be created that contains all peptide sequence matching information and can be utilized by compatible downstream software.
>   > The `Certificate of Analysis` provides details on all parameters settings of both Search GUI and Peptide Shaker used for the analysis.
>   {: .comment}
>
>   > ### {% icon question %} Questions:
>   > 1. How many peptides were identified? How many proteins?
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

# Analysis of Contaminants
The FASTA database used for the peptide to spectrum matching contained some entries that were not expected to stem from the HeLa cell lysate, but are common contaminations in LC-MS/MS samples. The main reason to add those is to avoid false assignment of the spectra to other proteins. However, it also enables you to check for contaminations in your samples. **CAVE:** in human samples, many proteins that are common contaminants may also stem from the real sample. The real source can never be determined for sure.

> ### {% icon hands_on %} Hands-On: Analysis of Contaminants
>
> 1. Run ***Select*** {% icon tool %} on the Peptide Shaker Protein Report to select all lines that match the pattern "CONTAMINANT".
> 2. Remove all contaminants from your protein list by running ***Select*** {% icon tool %} on the Peptide Shaker Protein Report. Select only those lines that **do not** match the pattern "CONTAMINANT".
>
>   > ### {% icon question %} Questions
>   > 1. Which contaminants did you identify? Where do these contaminations come from?
>   > 2. What other sources of contaminants exist?
>   > 3. How many mycoplasma proteins did you identify? Does this mean that the analyzed HeLa cells were infected with mycoplasma?
>   > 4. How many false positives do we expect in our list? How many of these are expected to match mycoplasma proteins?
>   >
>   >  <details>
>   >  <summary>Click to view answers</summary>
>   >   <ol type="1">
>   >     <li> TRY_BOVIN is bovine trypsin. It was used to degrade the proteins to peptides. ALBU_BOVIN is bovine serum albumin. It is added to cell culture medium in high amounts.</li>
>   >     <li> Contaminants often stem from the experimenter, these are typically keratins or other high-abundant human proteins. Basically any protein present in the room of the mass spectrometer might get into the ion source, if it is airborne. As an example, sheep keratins are sometimes found in proteomic samples, stemming from clothing made of sheep wool.</li>
>   >     <li> There should be five mycoplasma proteins in your protein list. However, all of them stem from different mycoplasma species. Also, every protein was identified by one peptide only. You can see this in column 17-19 of your output. These observations make it very likely that we are facing false positives here.</li>
>   >     <li> As we were allowing for a false discovery rate of 1 %, we would expect 12 false positive proteins in our list.
>   >     False positives are expected to be randomly assigned to peptides in the FASTA database. Our database consists of about 20,000 human proteins and 4,000 mycoplasma proteins. Therefore, we would expect 17 % (= 2) of all false positives matching to mycoplasma proteins.</li>
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
{:.no_toc}

A premade workflow for this tutorial can be found [here](workflows/wf_proteinID_SG_PS.ga)

# Further Reading
{:.no_toc}

- [Search GUI and Peptide Shaker tutorials at Compomics](https://compomics.com/bioinformatics-for-proteomics/)
- [Using Search GUI and Peptide Shaker in Galaxy](https://drive.google.com/file/d/0B6bIeOvjBkbWVnhMLWxXdGVUY3M/view)