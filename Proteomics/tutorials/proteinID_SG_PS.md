---
layout: tutorial_hands_on
topic_name: Peptide and Protein ID
tutorial_name: proteinID_SG_PS
---

# Introduction

Identifying the proteins contained in a sample is an important step in any proteomic experiment. However, in most settings, proteins are digested to peptides before the LC-MS/MS analysis. In this so-called "bottom-up" procedure, only peptide masses are measured. Therefore, protein identification cannot be performed directly from raw data, but is a multi-step process: 

1. Raw data preparations
2. Peptide-to-Spectrum matching
3. Peptide inference
4. Protein inference

A plethora of different software solutions exists for each step. In this tutorial, we will show how to
use ***msconvert*** :wrench: and ***PeakPickerHiRes*** :wrench: for step 1,
***Search GUI*** :wrench: and ***Peptide Shaker*** :wrench: for the steps 2-4.

As an example dataset, we will use an LC-MS/MS analysis of HeLa cell lysate published
in [Vaudel et al., 2014, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/24678044). Detailed information
about the dataset can be found on [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD000674).
For step 2 we will use a validated human Uniprot FASTA database without appended decoys.
If you already completed the tutorial on [Database Handling](./database-handling.md)
you can use the constructed database before the **DecoyDatabase** :wrench: step.


> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Preparing raw data](#preparing-raw-data)
> 2. [Peptide-to-Spectrum matching](#peptide-and-protein-id)
> 3. [Peptide and Protein Inference](#peptide-and-protein-id) 
> 4. [Analysis of Contaminants](#analysis-of-contaminants)
> 5. [Peptide and Protein Evaluation](#evaluation-of-peptide-and-protein-ids)


# Preparing raw data

Raw data conversion is the first step of any proteomic data analysis. The most common converter is MSConvert, the format to convert to is mzML. SearchGUI takes only mgf format as input, but as we need the mzML format for several other tasks, we will convert to mzML first.

> ### :pencil2: Hands-on: Preparing raw data
>
> 1. Create a new history for this Peptide and Protein ID exercise.
> 2. Load the example dataset into your history from this [link](ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/01/PXD000674/qExactive01819.raw).
> 3. Rename the dataset to "Test data".
> 4. Run ***msconvert*** :wrench: on the test data to convert to the mzML format
>
>	> ### :nut_and_bolt: Comment: Local Use of MSConvert
>	> The vendor libraries used by MSConvert need a Windows system and is therefore rarely implemented in Galaxy instances. If ***msconvert*** :wrench: is not available in your Galaxy instance, please install the software on a Windows computer and run the conversion locally. You can find a detailed description of the necessary steps [here](http://genesis.ugent.be/files/costore/practicals/bioinformatics-for-proteomics/1-Peptide-and-Protein-Identification/1.2-Peak-List-Generation/1.2-Peak-List-Generation.pdf). Afterwards, upload the resulting mzML file to your Galaxy history.
>
> 5. Run ***PeakPickerHiRes*** :wrench: on the resulting mzML file.
>	> ### :nut_and_bolt: Comment: Peak Picking
>	> Depending on your machine settings, raw data will be generated either in profile mode or centroid mode. For most peptide search engines, the data have to be converted to centroid mode, a process called "peak picking". Machine vendors offer algorithms to extract peaks from profile raw data. This is implemented in ***msconvert*** :wrench: and can be run in parallel to the mzML conversion. However, the OpenMS tool ***PeakPickerHiRes*** :wrench: is reported to generate better results ([Lange et al., 2006, Pac Symp Biocomput](https://www.ncbi.nlm.nih.gov/pubmed/17094243)) and is therefore recommended for quantitative studies ([Vaudel et al., 2010, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/19953549)).
>	> If your data were generated on a low resolution mass spectrometer, use ***PeakPickerWavelet*** :wrench: instead.
> 6. Run ***FileConverter*** :wrench: on the picked mzML to convert to mgf format.
> 7. Change the ***Datatype*** of the ***FileConverter*** :wrench: output to mgf by clicking the pencil :pencil: icon.

# Peptide and Protein ID
MS/MS experiments identify peptides by isolating them and subsequently colliding them with a gas for fragmentation. This method generates a spectrum of peptide fragment masses for each isolated peptide - an MS2 spectrum. To find out the peptide sequences, the MS2 spectrum is compared to a theoretical spectrum generated from a protein database. This step is called peptide-to-spectrum (also: spectrum-to-sequence) matching. Accoringly, a peptide that is successfully matched to a sequence is termed PSM (Peptide-Spectrum-Match). There can be multiple PSMs per peptide, if the peptide was fragmented several times. Different peptide search engines have been developed to fulfill the matching procedure. 

It is generally recommended to use more than one peptide search engine and use the combined results for the final peptide inference ([Shteynberg et al., 2013, Mol. Cell. Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23720762)). Again, there are several software solutions for this, e.g. iProphet (TPP) or ConsensusID (OpenMS). In this tutorial we will use ***Search GUI*** :wrench:, as it can automatically search the data using several search engines. Its partner tool ***Peptide Shaker*** :wrench: is then used to combine and evaluate the search engine results. 

In bottom-up proteomics, it is necessary to combine the identified peptides to proteins. This is not a trivial task, as proteins are redundant to some degree. Thus, not every peptide can be assigned to only one protein. Luckily, the ***Peptide Shaker*** :wrench: already takes care of protein inference and even gives us some information on validity of the protein IDs. We will discuss validation in a later [step](#peptide-and-protein-validation) of this tutorial.

> ### :pencil2: Hands-on: Spectrum-to-Sequence matching
>
> 1. Copy the prepared protein database from the tutorial "Database handling" into your current history by using the multiple history view or upload the ready-made database from this [link]().
> 2. Run ***Search GUI*** :wrench: to search the mgf file against the protein database.
> 3. Run ***Peptide Shaker*** :wrench: on the Search GUI output. Enable the following outputs: `Zip File for import to Desktop App`, `mzidentML File`, `PSM Report`, `Peptide Report`, `Protein Report`. You can find a detailed description of Peptide Shaker output [here]().
>
>	> ### :question: Questions: 
>	> 1. How many peptides were identified? How many proteins?
>	> 2. How many peptides with oxidized methionine were identified?
>	>
>	>    <details>
>	>    <summary>Click to view answers</summary>
>	>    	<ol type="1">
>	>			<li> You should have identified 3,325 peptides and 1,170 proteins.</li>
				<li> 328 peptides contain an oxidized methionine (MeO). To get to this number, you can use ***Select*** :wrench: on the Peptide Report and search for either "Oxidation of M" or "M\<ox\>".</li>
>	>   	 </ol>
>	>    </details>

# Analysis of Contaminants
The FASTA database used for the peptide to spectrum matching contained some entries that were not expected to stem from the HeLa cell lysate, but are common contaminations in LC-MS/MS samples. The main reason to add those is to avoid false assignment of the spectra to other proteins. However, it also enables you to check for contaminations in your samples. 

> ### :pencil2: Hands-on: Analysis of Contaminants
>
> 1. Run ***Select*** :wrench: on the Peptide Shaker Protein Report to select all lines that match the pattern "CONTAMINANT".
> 2. Remove all contaminants from your protein list by running ***Select*** :wrench: on the Peptide Shaker Protein Report. Select only those lines that DO NOT match the pattern "CONTAMINANT".
>
>	> ### :question: Questions: 
>	> 1. Which contaminants did you identify? Where do these contaminations come from?
>	> 2. How many mycoplasma proteins did you identify? Does this mean that the analyzed HeLa cells were infected with mycoplasma?
>	>
>	>    <details>
>	>    <summary>Click to view answers</summary>
>	>    	<ol type="1">
>	>			<li> TRY_BOVIN is bovine trypsin. It was used to degrade the proteins to peptides. ALBU_BOVIN is bovine serum albumin. It is added to cell culture medium in high amounts.</li>
				<li> There should be five mycoplasma proteins in your protein list. However, all of them stem from different mycoplasma species. Also, every protein was identified by one peptide only. You can see this is column 17-19 of your output. These observations makes it very likely that we are facing false positives here. As we were allowing for a false discovery rate of 1 %, we would expect 12 false positive proteins in our list. False positives are distributed to random peptides in the FASTA database. Our database consists of about 20,000 human proteins and 4,000 mycoplasma proteins. Therefore, we would expect 20 % of all false positives to match to mycoplasma proteins.</li>
>	>   	 </ol>
>	>    </details>

# Evaluation of Peptide and Protein IDs
***Peptide Shaker*** :wrench: provides us even with validation results for the identified PSM, peptides and proteins. It classifies all these IDs in either "Confident" or "Doubtful". On each level, the meaning differs somewhat. PSMs are marked as "Doubtful" when the measured MS2 spectrum did not fit perfectly to the theoretical spectrum. Peptides have a combined scoring of their PSMs. They are marked as "Doubtful", when the score is below a set threshold. The threshold is defined by the false discovery rate (FDR). At last, proteins are marked as "Doubtful", when they were identified by only a single peptide or when only identified by "Doubtful" peptides.