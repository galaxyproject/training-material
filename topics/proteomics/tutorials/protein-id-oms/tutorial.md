---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: protein-id-oms
---

# Introduction
{:.no_toc}

Identifying the proteins contained in a sample is an important step in any proteomic experiment. However, in most settings, proteins are digested to peptides prior to LC-MS/MS analysis. In this so-called "bottom-up" procedure, only peptide masses are measured. Therefore, protein identification cannot be performed directly from raw data, but is a multi-step process:

1. Raw data preparations
2. Peptide-to-Spectrum matching
3. Peptide inference
4. Protein inference

A plethora of different software solutions exists for each step. In this tutorial, we will use ***msconvert*** {% icon tool %}  for raw data conversion and tools from the [OpenMS software suite](https://openms.de) for all other steps. We will use one peptide search engine at first and later on show how to expand the workflow for using multiple search engines. Protein inference will be performed with the Fido algorithm ([Serang et al, JPR, (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20712337)).
This tutorial covers peptide and protein **identification** only, but you may use the output of this tutorial for the [tutorial on protein quantitation]({{site.baseurl}}/topics/proteomics/tutorials/protein-quant-sil/tutorial.html).

For an alternative protein ID workflow using the [Compomics](https://compomics.com/) tools [SearchGUI](https://compomics.github.io/projects/searchgui.html) and [PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html), please consult [this tutorial]({{site.baseurl}}/topics/proteomics/tutorials/protein-id-sg-ps/tutorial.html).
The latter tutorial does not allow to continue with the tutorial on protein quantitation.

# Input data
{:.no_toc}

As an example dataset, we will use an LC-MS/MS analysis of HeLa cell lysate published
in [Vaudel et al., 2014, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/24678044). Detailed information
about the dataset can be found on [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD000674).
For step 2 we will use a validated human Uniprot FASTA database with appended decoys.

If you already completed the tutorial on [Database Handling]({{site.baseurl}}/topics/proteomics/tutorials/database-handling/tutorial.html)
you can use the constructed database including decoys.
You can find a prepared database, as well as the input LC-MS/MS data in different file formats on [Zenodo](https://zenodo.org/record/796184).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preparing Raw Data

Raw data conversion is the first step of any proteomic data analysis. The most common converter is MSConvert from the [ProteoWizard software suite](http://proteowizard.sourceforge.net/), the format to convert to is mzML.
Due to licensing reasons, MSConvert runs only on windows systems and will not work on most Galaxy servers.

Depending on your machine settings, raw data will be generated either in profile mode or centroid mode. For most peptide search engines, the MS2 data have to be converted to centroid mode, a process called "peak picking" or "centroiding".
Machine vendors offer algorithms to extract peaks from profile raw data. Those are integrated in ***msconvert*** {% icon tool %} and can be run in parallel to the mzML conversion.
However, the OpenMS tool ***PeakPickerHiRes*** {% icon tool %} is reported to generate slightly better results ([Lange et al., 2006, Pac Symp Biocomput](https://www.ncbi.nlm.nih.gov/pubmed/17094243)) and is therefore recommended for quantitative studies ([Vaudel et al., 2010, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/19953549)).
If your data were generated on a low resolution mass spectrometer, use ***PeakPickerWavelet*** {% icon tool %} instead.

> ### {% icon hands_on %} Hands-On: File Conversion and Peak Picking
>
> We provide the [input data](https://zenodo.org/record/796184) in the original `raw` format and also already converted to `mzML`. If ***msconvert*** {% icon tool %} does not run on your Galaxy instance, please download the preconverted `mzML` as an input.
>
> 1. Create a new history for this Peptide and Protein ID exercise.
> 2. Load the example dataset into your history from Zenodo: [raw](https://zenodo.org/record/892005/files/qExactive01819.raw) [mzML](https://zenodo.org/record/892005/files/qExactive01819_profile.mzml)
> 3. Rename the dataset to something meaningful.
> 4. (*optional*) Run ***msconvert*** {% icon tool %} on the test data to convert to the `mzML` format.
> 5. Run ***PeakPickerHiRes*** {% icon tool %} on the resulting file. Click `+ Insert param.algorithm_ms_levels` and change the entry to "2". Thus, peak picking will only be performed on MS2 level.
>
>   > ### {% icon comment %} Comment: Local Use of MSConvert
>   > The vendor libraries used by MSConvert are only licensed for Windows systems and are therefore rarely implemented in Galaxy instances. If ***msconvert*** {% icon tool %} is not available in your Galaxy instance, please install the software on a Windows computer and run the conversion locally. You can find a detailed description of the necessary steps [here](https://compomics.com/bioinformatics-for-proteomics/identification/) ("Peak List Generation"). Afterwards, upload the resulting mzML file to your Galaxy history.
>  {: .comment}
>
>   > ### {% icon comment %} Comment: MS2 peak picking during data acquisition
>   > MS2 peaks are often acquired in centroided mode in first place. The profile data are converted to centroided mode already during data acquisition, resulting in MS2-centroided `raw` files. If your MS2 data are already centroided, simply omit the peak picking step.
>  {: .comment}

{: .hands_on}

# Peptide Identification
MS/MS experiments identify peptides by isolating them and subsequently colliding them with a gas for fragmentation. This method generates a spectrum of peptide fragment masses for each isolated peptide - an MS2 spectrum.
To find out the peptide sequences, the MS2 spectrum is compared to a theoretical spectrum generated from a protein database. This step is called peptide-to-spectrum (also: spectrum-to-sequence) matching. Accordingly, a peptide that is successfully matched to a sequence is termed PSM (Peptide-Spectrum-Match). There can be multiple PSMs per peptide, if the peptide was fragmented several times.

Different peptide search engines have been developed to fulfill the matching procedure. Here, we will use the search engine [X!Tandem](https://www.ncbi.nlm.nih.gov/pubmed/14976030). OpenMS provides "adapters" (wrappers) for several other peptide search engines, like MSGF+ or OMSSA. You may replace the XTandemAdapter by another search engine of your choice.

> ### {% icon hands_on %} Hands-On: Peptide Identification
>
> 1. Copy the prepared protein database from the tutorial [Database Handling](../database-handling/tutorial.html) into your current history by using the multiple history view or upload the ready-made database from this [link](https://zenodo.org/record/892005/files/Human_database_including_decoys_%28cRAP_and_Mycoplasma_added%29.fasta).
> 2. Run the tool ***XTandemAdapter*** {% icon tool %} with:
    - the MS2-centroided mzML as **Input file containing MS2 spectra** and
    - the FASTA protein database as **FASTA file or pro file**.
    - Click `+ Insert param_fixed_modifications` and choose `Carbamidomethyl (C)`.
    - Click `+ Insert param_variable_modifications` and choose `Oxidation (M)`.
> 3. Run the tool ***FileInfo*** {% icon tool %} on the XTandem output.
>
>   > ### {% icon comment %} Comment: Advanced Search Engine Parameters
>   > The OpenMS adapters do not always allow to set every option of the underlying search engine. If an option is missing, you may also run the search engine locally or by using a Galaxy wrapper. Afterwards, convert the search engine output to the OpenMS format `idXML` by running ***IDFileConverter*** {% icon tool %}.
>   >
>   > The search engine X!Tandem features some more advanced options than the ones reflected in the ***XTandemAdapter*** {% icon tool %}. If you need those advanced options, the ***XTandemAdapter*** {% icon tool %} allows for the optional input of a classic X!Tandem parameter file. Upload your parameter file to the history and use it as an input in the field `Default X!Tandem configuration file`. You may also set the option `-ignore_adapter_param` to `Yes` to overwrite all options set by the GUI.
>   {: .comment}
{: .hands_on}

# Peptide FDR filtering
The next step of peptide identification is to decide which PSMs will be used for protein inference. Measured MS2 spectra never perfectly fit the theoretical spectra. Therefore, peptide search engines calculate a score which indicates how well the measured MS2 spectrum was fitting the theoretical spectrum. How do we decide which PSMs are likely true and which are false?

In proteomics, this decision is typically done by calculating false discovery rates (FDRs). Remember that the database we were using for peptide-to-spectrum matching consisted not only of true proteins, but also the same number of "fake entries", the so-called decoys. Those decoys can now be used to estimate the number of false identifications in the list of PSMs.
The calculation is based on a simple assumption: for every decoy protein identified with a given score, we expect one false positive with at least the same score.
The false discovery rate is therefore defined as the number of false discoveries (decoy hits) divided by the number of false and correct discoveries (both target and decoy hits) at a given score threshold.

To calculate FDRs, we first have to annotate the identified peptides to determine which of them are decoys. This is done with the tool ***PeptideIndexer*** {% icon tool %}. Additionally, we will calculate peptide posterior error probabilities (PEPs), because they are needed for the protein inference algorithm used by OpenMS. We will then filter for 1 % FDR and set the score back to PEP.

> ### {% icon hands_on %} Hands-On: Peptide FDR filtering
>
> 2. Run ***IDPosteriorErrorProbability*** {% icon tool %} with
>   - `-prob_correct` set to `Yes`.
> 1. Run ***PeptideIndexer*** {% icon tool %} with
>   - the FASTA protein database as **Input sequence database in FASTA format**, and
>   - **Specificity of the enzyme** set to `none`.
> 3. Run ***FalseDiscoveryRate*** {% icon tool %} with
>   - **Perform FDR calculation on protein level** set to `false`,
>   - **Filter PSMs based on q-value** set to `0.01`, and
>   - `-add_decoy_peptides` set to `Yes`.
> 4. Run ***IDScoreSwitcher*** {% icon tool %} with
>   - **Name of the meta value to use as the new score** set to "Posterior Probability_score", and
>   - **Orientation of the new score`** set to `higher_better`.
> 5. Run ***FileInfo*** {% icon tool %} to get basic information about the identified peptides.
>
>   > ### {% icon question %} Questions:
>   > 1. How many peptides were identified?
>   > 2. How many peptides with oxidized methionine were identified?
>   >
>   >  <details>
>   >  <summary>Click to view answers</summary>
>   >   <ol type="1">
>   >     <li> You should have identified 2,616 unique stripped peptides.</li>
>   >     <li> 503 peptides contain an oxidized methionine (MeO).</li>
>   >   </ol>
>   >  </details>
>   {: .question}
{: .hands_on}

# Protein Inference
In bottom-up proteomics, it is necessary to combine the identified peptides to proteins. This is not a trivial task, as proteins are redundant to some degree. Thus, not every peptide can be assigned to only one protein.
The OpenMS suite implemented the [Fido](https://www.ncbi.nlm.nih.gov/pubmed/20712337) algorithm for protein inference. Fido uses a Bayesian probabilistic model to group and score proteins based on peptide-spectrum matches.

> ### {% icon hands_on %} Hands-On: Protein inference
>
> 1. Run ***FidoAdapter*** {% icon tool %}. Set `-greedy_group_resolution` = `Yes`.
> 2. Run ***FalseDiscoveryRate*** {% icon tool %}. Set the option **`Perform FDR calculation on PSM level`** to `false`.
> 3. Run ***IDFilter*** {% icon tool %}. Set`-prot` to `0.01`.
> 4. Run ***FileInfo*** {% icon tool %} to get basic information about the identified proteins.
>
>   > ### {% icon comment %} Comment: "Greedy" Group Resolution
>   > Protein groups are reported, when an identified peptide maps to multiple proteins in the used database [Nesvizhskii and Aebersold (2005)](https://www.ncbi.nlm.nih.gov/pubmed/16009968). Some peptides may map to different protein groups and can therefore not be used for protein quantitation. The option `-greedy_group_resolution` solves this problem by assigning peptides only to the one most probable protein group, thus enabling to quantify proteins based not only on unique, but also on shared peptides. This usually leads to a much higher number of quantified proteins. However it will introduce noise in the FCs when a peptide was indeed shared by different proteins and the quantity of this peptide was a weighted sum of contributions. The greedy group resolution is similar to Occam's razor.
>   {: .comment}
{: .hands_on}

# Analysis of Contaminants
The FASTA database used for the peptide to spectrum matching contained some entries that were not expected to stem from the HeLa cell lysate, but are common contaminations in LC-MS/MS samples. The main reason to add those is to avoid false assignment of contaminant spectra to other proteins.
It also enables you to check for contaminations in your samples.

**CAVE:** When analyzing human samples, many proteins that are common contaminants may also stem from the sample. Therefore, human contaminants do not have to be excluded from further analysis, but you should keep in mind that the source of these proteins is unclear.

> ### {% icon hands_on %} Hands-On: Analysis of Contaminants
>
> 1. Run ***TextExporter*** {% icon tool %} to convert the idXML output to a human-readable tabular file.
> 1. Run ***Select*** {% icon tool %} to select all lines **Matching** the pattern "CONTAMINANT".
> 3. Run ***Select*** {% icon tool %} to select all lines that **NOT Matching** the pattern "HUMAN".
> 2. Remove all bovine and mycoplasma proteins from your list by running ***Select*** {% icon tool %}. Select only those lines **NOT Matching** match the pattern "HUMAN".
>
>   > ### {% icon question %} Questions
>   > 1. Which contaminants did you identify? Where do these contaminations likely come from?
>   > 2. What other sources of contaminants exist?
>   > 3. How many mycoplasma proteins did you identify? Does this mean that the analyzed HeLa cells were infected with mycoplasma?
>   > 4. How many false positives do we expect in our list?
>   >
>   >  <details>
>   >  <summary>Click to view answers</summary>
>   >   <ol type="1">
>   >     <li> TRY_BOVIN is bovine trypsin. It was used to degrade the proteins to peptides. ALBU_BOVIN is bovine serum albumin. It is added to cell culture medium in high amounts. Also, five human proteins are listed, these are commonly introduced during sample preparation. As we were analyzing a human sample, it is not neccessary to remove these proteins, as they may as well originate from the HeLa cells.</li>
>   >     <li> Contaminants often stem from the experimenter, these are typically keratins or other high-abundant human proteins. Basically any protein present in the room of the mass spectrometer might get into the ion source, if it is airborne. As an example, sheep keratins are sometimes found in proteomic samples, stemming from clothing made of sheep wool.</li>
>   >     <li> One protein stemming from *Acholeplasma laidlawii* (ACHLI) was identified. If you again filter the protein list for "ACHLI", you will see that it was identified by a single peptide. Thus, it is likely a false positive and does not indicate contamination.</li>
>   >     <li> As we were allowing for a false discovery rate of 1 %, we would expect 12 false positive proteins in our list.</li>
>   >   </ol>
>   >  </details>
>   {: .question}
{: .hands_on}

# Using multiple search engines

It is generally recommended to use more than one peptide search engine and use the combined results for peptide inference ([Shteynberg et al., 2013, Mol. Cell. Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23720762)).
By comparing results of multiple search engines, you may improve the *sensitivity* (when accepting peptides that were found by only one of the engines), the *specificity* (when accepting only peptides that were found by all of the search engines) or *both* (when using n>2 search engines and accept peptides found by a fraction of the (e.g. n-1) search engines).

Here, we will use the OpenMS tool [ConsensusID](http://ftp.mi.fu-berlin.de/pub/OpenMS/release-documentation/html/TOPP_ConsensusID.html) to combine the search engine results.

> ### {% icon hands_on %} Hands-On: Multiple search engines
>
> 1. Run ***MSGFPlusAdapter*** {% icon tool %} with
    - the MS2-centroided mzML as **Input file** and
    - the FASTA protein database as **Protein sequence database**.
    - **Precursor monoisotopic mass tolerance** set to `10.0`.
    - **Instrument that generated the data** set to `Q_Exactive`.
    - Click `+ Insert param_fixed_modifications` and choose `Carbamidomethyl (C)`.
    - Click `+ Insert param_variable_modifications` and choose `Oxidation (M)`.
> 3. Run the tool ***FileInfo*** {% icon tool %} on the MSGFPlusAdapter output.
> 2. Run ***IDPosteriorErrorProbability*** {% icon tool %} with
>   - `-prob_correct` set to `No`.
> 5. Run ***IDMerger*** {% icon tool %} with two **Input files [...]**:
>   - the output of **IDScoreSwitcher** based on **XTandemAdapter**
>   - the output of **IDScoreSwitcher** based on **MSGFPlusAdapter**
> 6. Run ***ConsensusID*** {% icon tool %}.
> 1. Run ***PeptideIndexer*** {% icon tool %} with
>   - the FASTA protein database as **Input sequence database in FASTA format**, and
>   - **Specificity of the enzyme** set to `none`.
> 3. Run ***FalseDiscoveryRate*** {% icon tool %} with
>   - **Perform FDR calculation on protein level** set to `false`,
>   - **Filter PSMs based on q-value** set to `0.01`, and
>   - `-add_decoy_peptides` set to `Yes`.
> 4. Run ***IDScoreSwitcher*** {% icon tool %} with
>   - **Name of the meta value to use as the new score** set to "Posterior Error Probability_score",
>   - **Orientation of the new score`** set to `lower_better`, and
>   - **Name to use as the type of the new score** set to "Posterior Error Probability".
> 5. Run ***FileInfo*** {% icon tool %} to get basic information about the identified peptides.
> 9. Proceed with the protein inference as described [above](#protein-inference)
>
>   > ### {% icon question %} Questions:
>   > 1. How many PSMs could be matched with XTandem and MSGFPlus alone? How many peptides were identified?
>   > 2. How many PSMs could be matched after combining the results with ConsensusID? How many peptides were identified?
>   >
>   >  <details>
>   >  <summary>Click to view answers</summary>
>   >   <ol type="1">
>   >     <li> After FDR-filtering, XTandem matched 3,552 PSMs (2,616 unique peptides) and MSGFPlus matched 4,292 PSMs (2,991 peptides).</li>
>   >     <li> Combining the results with ConsensusID leads to matching of 4,299 PSMs (3,041 unique peptides).</li>
>   >   </ol>
>   >  </details>
>   {: .question}
{: .hands_on}

# Premade Workflow

A premade workflow for this tutorial can be found [here](workflows/workflow.ga).

A premade workflow using the search engines XTandem and MSGF+ can be found [here](workflows/workflow_two-search-engines.ga).

# Further Reading

- [Protein inference](https://www.ncbi.nlm.nih.gov/pubmed/16009968)
- [Fido publication](https://www.ncbi.nlm.nih.gov/pubmed/20712337)
- [Evaluation of protein inference algorithms](https://www.ncbi.nlm.nih.gov/pubmed/27498275)
