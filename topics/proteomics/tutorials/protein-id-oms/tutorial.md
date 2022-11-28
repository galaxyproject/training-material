---
layout: tutorial_hands_on

title: "Peptide and Protein ID using OpenMS tools"
zenodo_link: "https://zenodo.org/record/546301"
level: Advanced
questions:
  - "How to convert LC-MS/MS raw files?"
  - "How to identify peptides?"
  - "How to identify proteins?"
  - "How to evaluate the results?"
objectives:
  - "Protein identification from LC-MS/MS raw files."
time_estimation: "45m"
requirements:
  -
    type: "internal"
    topic_name: proteomics
    tutorials:
    - database-handling
key_points:
  - "LC-MS/MS raw files have to be converted to mzML before using GalaxyP on most GalaxyP servers."
  - "OpenMS provides many tools for proteomic analysis and guarantees compatibility by using open file formats."
  - "OpenMS provides several thirdparty search engines and Fido for protein inference."
follow_up_training:
-
  type: "internal"
  topic_name: proteomics
  tutorials:
    - protein-quant-sil
contributors:
  - stortebecker
  - bgruening
subtopic: id-quant
tags: [DDA]
---

# Introduction


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


As an example dataset, we will use an LC-MS/MS analysis of HeLa cell lysate published
in [Vaudel et al., 2014, Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/24678044). Detailed information
about the dataset can be found on [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD000674).
For step 2 we will use a validated human Uniprot FASTA database with appended decoys.

If you already completed the tutorial on [Database Handling]({{site.baseurl}}/topics/proteomics/tutorials/database-handling/tutorial.html)
you can use the constructed database including decoys.
You can find a prepared database, as well as the input LC-MS/MS data in different file formats on [Zenodo](https://zenodo.org/record/796184).

> <agenda-title></agenda-title>
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
If your data were generated on a low resolution mass spectrometer, use **PeakPickerWavelet** {% icon tool %} instead.

> <hands-on-title>File Conversion and Peak Picking</hands-on-title>
>
> We provide the [input data](https://zenodo.org/record/796184) in the original `raw` format and also already converted to `mzML`. If **msconvert** {% icon tool %} does not run on your Galaxy instance, please download the preconverted `mzML` as an input and continue with step 5 of the following hands-on training.
>
> 1. Create a new history for this Peptide and Protein ID exercise.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Load one of the example datasets into your history from Zenodo
>
>    ```
>    https://zenodo.org/record/892005/files/qExactive01819.raw
>    https://zenodo.org/record/892005/files/qExactive01819_profile.mzml
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the dataset to something meaningful
>
> 4. Run {% tool [msconvert](toolshed.g2.bx.psu.edu/repos/galaxyp/msconvert/msconvert/3.0.19052.1) %} on the training `raw` file to convert to the `mzML` format
>    - {% icon param-file %} *"Input unrefined MS data"*: imported `raw` file
>    - *"Do you agree to the vendor licenses"*: set to `Yes`
>    - *"Output Type"*: set to `mzML`
>
> 5. Run {% tool [PeakPickerHiRes](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_peakpickerhires/PeakPickerHiRes/2.3.0) %} on the resulting file
>    - {% icon param-file %} *"input profile data file": output of **msconvert** or `mzML` file
>    - In *"param_algorithm_ms_levels"*
>      - {% icon param-repeat %} Click on *"Insert param_algorithm_ms_levels"*
>      - In *"1: param_algorithm_ms_levels"*
>         - *"List of MS levels for which the peak picking is applied"*: `2`
>
>           Peak picking will only be performed on MS2 level.
>
>   > <comment-title>Local Use of MSConvert</comment-title>
>   > The vendor libraries used by MSConvert are only licensed for Windows systems and are therefore rarely implemented in Galaxy instances. If **msconvert** {% icon tool %} is not available in your Galaxy instance, please install the software on a Windows computer and run the conversion locally. You can find a detailed description of the necessary steps [here](https://compomics.com/bioinformatics-for-proteomics/identification/) ("Peak List Generation"). Afterwards, upload the resulting mzML file to your Galaxy history.
>   {: .comment}
>
{: .hands_on}

> <comment-title>MS2 peak picking during data acquisition</comment-title>
> MS2 peaks are often acquired in centroided mode in first place. The profile data are converted to centroided mode already during data acquisition, resulting in MS2-centroided `raw` files. If your MS2 data are already centroided, simply omit the peak picking step.
{: .comment}

# Peptide Identification

MS/MS experiments identify peptides by isolating them and subsequently colliding them with a gas for fragmentation. This method generates a spectrum of peptide fragment masses for each isolated peptide - an MS2 spectrum.
To find out the peptide sequences, the MS2 spectrum is compared to a theoretical spectrum generated from a protein database. This step is called peptide-to-spectrum (also: spectrum-to-sequence) matching. Accordingly, a peptide that is successfully matched to a sequence is termed PSM (Peptide-Spectrum-Match). There can be multiple PSMs per peptide, if the peptide was fragmented several times.

Different peptide search engines have been developed to fulfill the matching procedure. Here, we will use the search engine [X!Tandem](https://www.ncbi.nlm.nih.gov/pubmed/14976030). OpenMS provides "adapters" (wrappers) for several other peptide search engines, like MSGF+ or OMSSA. You may replace the XTandemAdapter by another search engine of your choice.

> <hands-on-title>Hands-On: Peptide Identification</hands-on-title>
>
> 1. Copy the prepared protein database (Human database including cRAP contaminants and decoys) from the tutorial [Database Handling]({% link topics/proteomics/tutorials/database-handling/tutorial.md %}) into your current history by using the multiple history view
>
>   > <comment-title>You did not run the Database Handling first?</comment-title>
>   > You can upload the ready-made database from Zenodo
>   >
>   > ```
>   > https://zenodo.org/record/892005/files/Human_database_including_decoys_%28cRAP_added%29.fasta
>   > ```
>   {: .comment}
>
> 2. Run {% tool [XTandemAdapter](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_xtandemadapter/XTandemAdapter/2.6+galaxy0) %} with:
>   - {% icon param-file %} *"Input file containing MS2 spectra"*: MS2-centroided mzML file
>   - {% icon param-file %} *"FASTA file"*: FASTA protein database
>   - *"Fragment mass error"*: `10`
>   - *"Fragment monoisotopic mass error units"*: `ppm`
>   - *"Optional Outputs"* select `out (Output file containing search results)`
>
> 3. Run {% tool [FileInfo](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_fileinfo/FileInfo/2.6+galaxy0) %}
>    - {% icon param-file %} *"input file"*: **XTandem** output
{: .hands_on}

> <comment-title>Settings for labelled data</comment-title>
> Several common quantitation methods are based on labels e.g. SILAC, TMT, iTRAQ, Dimethyl. Those labels must be specified in the search engine as a (variable) modification. See also [Peptide and Protein Quantification via Stable Isotope Labelling]({% link topics/proteomics/tutorials/protein-quant-sil/tutorial.md %}).
{: .comment}

> <comment-title>Advanced Search Engine Parameters</comment-title>
> The OpenMS adapters do not always allow to set every option of the underlying search engine. If an option is missing, you may also run the search engine locally or by using a Galaxy wrapper. Afterwards, convert the search engine output to the OpenMS format `idXML` by running **IDFileConverter** {% icon tool %}.
>
> The search engine X!Tandem features some more advanced options than the ones reflected in the **XTandemAdapter** {% icon tool %}. If you need those advanced options, the **XTandemAdapter** {% icon tool %} allows for the optional input of a classic X!Tandem parameter file. Upload your parameter file to the history and use it as an input in the field `Default X!Tandem configuration file`. You may also set the option *"ignore_adapter_param"* to `Yes` to overwrite all options set by the GUI.
{: .comment}

# Peptide FDR filtering

The next step of peptide identification is to decide which PSMs will be used for protein inference. Measured MS2 spectra never perfectly fit the theoretical spectra. Therefore, peptide search engines calculate a score which indicates how well the measured MS2 spectrum was fitting the theoretical spectrum. How do we decide which PSMs are likely true and which are false?

In proteomics, this decision is typically done by calculating false discovery rates (FDRs). Remember that the database we were using for peptide-to-spectrum matching consisted not only of true proteins, but also the same number of "fake entries", the so-called decoys. Those decoys can now be used to estimate the number of false identifications in the list of PSMs.
The calculation is based on a simple assumption: for every decoy peptide identified with a given score, we expect one false positive with at least the same score.
The false discovery rate is therefore defined as the number of false discoveries (decoy hits) divided by the number of false and correct discoveries (both target and decoy hits) at a given score threshold.

To calculate FDRs, we first have to annotate the identified peptides to determine which of them are decoys. This is done with the tool **PeptideIndexer** {% icon tool %}. This tool allows the annotation of protein names and calculates the protein coverage based on the identified peptides. Then we calculate peptide posterior error probabilities (PEPs), because they are needed for the protein inference algorithm Fido, which is used by OpenMS. We will then filter PSMs for 1 % FDR with **FalseDiscoveryRate** {% icon tool %}. and set the score back to PEP with **IDScoreSwitcher** {% icon tool %}.

> <hands-on-title>Hands-On: Peptide FDR filtering</hands-on-title>
>
> 1. Run {% tool [PeptideIndexer](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_peptideindexer/PeptideIndexer/2.6+galaxy0) %} with
>    - {% icon param-file %} *"Input idXML file containing the identifications"*: output of **XTandemAdapter**
>    - {% icon param-file %} *"Input sequence database in FASTA format"*: FASTA protein database
>    - *"If set, the protein sequences are stored as well"*: `Yes`
>    - *"If set, the protein description is stored as well"*: `Yes`
>    - In the *"enzyme"* section set *"Specificity of the enzyme"*: `none`
>
> 2. Run {% tool [IDPosteriorErrorProbability](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_idposteriorerrorprobability/IDPosteriorErrorProbability/2.6+galaxy) %} with
>    - {% icon param-file %} *"input file"*: **PeptideIndexer** output
>    - *"If set scores will be calculated as '1 - ErrorProbabilities' and can be interpreted as probabilities for correct identifications"*: `Yes`
>
> 3. Run {% tool [FalseDiscoveryRate](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_falsediscoveryrate/FalseDiscoveryRate/2.6+galaxy0) %} with
>    - {% icon param-file %} *"Identifications from searching a target-decoy database"*: output of **IDPosteriorErrorProbability**
>    - *"Perform FDR calculation on protein level"*: `false`
>    - In the *"FDR control"* section select *"Filter PSMs based on q-value"*: `0.01`
>    - In the "Parameter section for the FDR calculation algorithm"* section select *"If 'true' decoy peptides will be written to output file, too"*: `Yes`
>
> 4. Run {% tool [IDScoreSwitcher](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_idscoreswitcher/IDScoreSwitcher/2.6+galaxy0) %} with
>    - {% icon param-file %} *"Input file"*: output of **FalseDiscoveryRate**
>    - *"Name of the meta value to use as the new score"*: `Posterior Probability_score`
>    - *"Orientation of the new score"*: `higher_better`
>
> 5. Run {% tool [FileInfo](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_fileinfo/FileInfo/2.6+galaxy0) %} to get basic information about the identified peptides
>    - {% icon param-file %} *"input file"*: **IDScoreSwitcher** output
>
{: .hands_on}

> <question-title></question-title>
> 1. How many peptides were identified?
> 2. How many peptides with oxidized methionine were identified?
>
> > <solution-title></solution-title>
> > 1. You should have identified 3378 non-redundant peptide hits.
> > 2. 804 peptides contain an oxidized methionine.
> > Numbers may slightly vary depending on the versions of the tools and the used FASTA file.
> {: .solution }
{: .question}


# Protein Inference
In bottom-up proteomics, it is necessary to combine the identified peptides to proteins. This is not a trivial task, as proteins are redundant to some degree. Thus, not every peptide can be assigned to only one protein.
The OpenMS suite implemented the [Fido](https://www.ncbi.nlm.nih.gov/pubmed/20712337) algorithm for protein inference. Fido uses a Bayesian probabilistic model to group and score proteins based on peptide-spectrum matches. Afterwards, we keep only proteins with 1% FDR.

> <hands-on-title>Hands-On: Protein inference</hands-on-title>
>
> 1. Run {% tool [FidoAdapter](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_fidoadapter/FidoAdapter/2.6+galaxy0) %}
>    - {% icon param-file %} *"Input: identification results"*: output of **IDScoreSwitcher**
>    - *"Post-process Fido output with greedy resolution of shared peptides based on the protein probabilities"*: `Yes`
>
> 2. Run {% tool [FalseDiscoveryRate](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_falsediscoveryrate/FalseDiscoveryRate/2.6+galaxy0) %}
>    - {% icon param-file %} *"Identifications from searching a target-decoy database"*: output of **FidoAdapter**
>    - *"Perform FDR calculation on PSM level"*: `false`
>    - In the *"FDR control"* section select *"Filter proteins based on q-value"*: `0.01`
>
> 3. Run {% tool [FileInfo](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_fileinfo/FileInfo/2.6+galaxy0) %} to get basic information about the identified proteins
>    - {% icon param-file %} *"input file"*: output of **FalseDiscoveryRate**
{: .hands_on}

> <comment-title>"Greedy" Group Resolution</comment-title>
> Protein groups are reported, when an identified peptide maps to multiple proteins in the used database [Nesvizhskii and Aebersold (2005)](https://www.ncbi.nlm.nih.gov/pubmed/16009968). Some peptides may map to different protein groups and can therefore not be used for protein quantitation. The option `-greedy_group_resolution` solves this problem by assigning peptides only to the one most probable protein group, thus enabling to quantify proteins based not only on unique, but also on shared peptides. This usually leads to a much higher number of quantified proteins. However it will introduce noise in the FCs when a peptide was indeed shared by different proteins and the quantity of this peptide was a weighted sum of contributions. The greedy group resolution is similar to Occam's razor.
{: .comment}

> <question-title></question-title>
> How many proteins were finally identified?
>
> > <solution-title></solution-title>
> > You should have identified 1252 proteins.
> > Numbers may slightly vary depending on the versions of the tools and the used FASTA file.
> {: .solution }
{: .question}

# Analysis of Contaminants

The FASTA database used for the peptide to spectrum matching contained some entries that were not expected to stem from the HeLa cell lysate, but are common contaminations in LC-MS/MS samples. The main reason to add those is to avoid false assignment of contaminant spectra to other proteins.
It also enables you to check for contaminations in your samples.

**CAVE:** When analyzing human samples, many proteins that are common contaminants may also stem from the sample. Therefore, human contaminants do not have to be excluded from further analysis, but you should keep in mind that the source of these proteins is unclear.

> <hands-on-title>Hands-On: Analysis of Contaminants</hands-on-title>
>
> 1. Run {% tool [TextExporter](toolshed.g2.bx.psu.edu/repos/galaxyp/openms_textexporter/TextExporter/2.6+galaxy0) %} to convert the idXML output to a human-readable tabular file.
>    - {% icon param-file %} *"Input file"*: **FalseDiscoveryRate** output
>
> 2. Run {% tool [Select lines that match an expression](Grep1) %}
>    - {% icon param-file %} *"Select lines from"*: **TextExporter**
>    - *"that"*: `Matching`
>    - *"the pattern"*: `CONTAMINANT`
>
> 2. Run {% tool [Select lines that match an expression](Grep1) %} to remove all non human proteins (e.g. bovine)
>    - {% icon param-file %} *"Select lines from"*: **TextExporter**
>    - *"that"*: `Matching`
>    - *"the pattern"*: `HUMAN`
>
{: .hands_on}

> <question-title></question-title>
> 1. Which contaminants did you identify? Where do these contaminations likely come from?
> 2. What other sources of contaminants exist?
> 3. How many false positives do we expect in our list?
>
> > <solution-title></solution-title>
> > 1.  TRY1_BOVIN is bovine trypsin. It was used to degrade the proteins to peptides. ALBU_BOVIN is bovine serum albumin. It is added to cell culture medium in high amounts. Also, eight human proteins are listed, these are commonly introduced during sample preparation. As we were analyzing a human sample, it is not neccessary to remove these proteins, as they may as well originate from the HeLa cells.
> > 2.  Contaminants often stem from the experimenter, these are typically keratins or other high-abundant human proteins. Basically any protein present in the room of the mass spectrometer might get into the ion source, if it is airborne. As an example, sheep keratins are sometimes found in proteomic samples, stemming from clothing made of sheep wool.
> > 3.  As we were allowing for a false discovery rate of 1 %, we would expect 13 (1255/0.01) false positive proteins in our list.
> {: .solution }
{: .question}


# Using multiple search engines

It is generally recommended to use more than one peptide search engine and use the combined results for peptide inference ([Shteynberg et al., 2013, Mol. Cell. Proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23720762)).
By comparing results of multiple search engines, you may improve the *sensitivity* (when accepting peptides that were found by only one of the engines), the *specificity* (when accepting only peptides that were found by all of the search engines) or *both* (when using n>2 search engines and accept peptides found by a fraction of the (e.g. n-1) search engines).

Here, we will use the OpenMS tool [ConsensusID](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_ConsensusID.html) to combine the search engine results.

> <hands-on-title>Hands-On: Multiple search engines</hands-on-title>
>
> 1. Run **MSGFPlusAdapter** {% icon tool %} with
>    - {% icon param-file %} *"Input file"*: the MS2-centroided mzML
>    - {% icon param-file %} *"Protein sequence database*": the FASTA protein database
>    - *"Precursor monoisotopic mass tolerance"*: `10.0`
>    - *"Instrument that generated the data"*: `Q_Exactive`
>    - In *"param_fixed_modifications"*
>      - {% icon param-repeat %} Click on *"Insert param_fixed_modifications"*
>      - In *"1: param_fixed_modifications"*
>         - *"Fixed modifications"*: `Carbamidomethyl (C)`
>    - In *"param_variable_modifications"*
>      - {% icon param-repeat %} Click on *"Insert param_variable_modifications"*
>      - In *"1: param_variable_modifications"*
>         - *"Variable modifications"*: `Oxidation (M)`
>
> 2. Run **FileInfo** {% icon tool %} to get basic information about the identified proteins
>    - {% icon param-file %} *"input file"*: **MSGFPlusAdapter** output
>
> 3. Run **IDPosteriorErrorProbability** {% icon tool %} with
>    - {% icon param-file %} *"input file"*: **MSGFPlusAdapter** output
>    - *"If set scores will be calculated as '1 - ErrorProbabilities' and can be interpreted as probabilities for correct identifications"*: `No`
>
> 4. Run **IDMerger** {% icon tool %}
>    - *"Reduce collectons"*: `reduce collections by aggregating single files of multiple collections`
>      - In *"inputs"*
>        - In *"1: Repeat"*
>          - {% icon param-file %} *"Input files separated by blanks"*: output of **IDScoreSwitcher** based on **XTandemAdapter**
>        - In *"2: Repeat"*
>          - {% icon param-file %} *"Input files separated by blanks"*: output of **IDScoreSwitcher** based on **MSGFPlusAdapter**
>
> 5. Run **ConsensusID** {% icon tool %}
>    - {% icon param-file %} *"input file"*: **IDMerger** output
>
> 6. Run **PeptideIndexer** {% icon tool %} with
>    - {% icon param-file %} *"Input idXML file containing the identifications"*: **ConsensusID** output
>    - {% icon param-file %} *"Input sequence database in FASTA format"*: FASTA protein database
>    - *"Specificity of the enzyme"*: `none`
>
> 7. Run **FalseDiscoveryRate** {% icon tool %} with
>    - {% icon param-file %} *"Identifications from searching a target-decoy database"*: output of **PeptideIndexer**
>    - *"Perform FDR calculation on protein level"*: `false`
>    - *"Filter PSMs based on q-value"*: `0.01`
>    - *"If 'true' decoy peptides will be written to output file, too"*: `Yes`
>
> 8. Run **IDScoreSwitcher** {% icon tool %} with
>    - {% icon param-file %} *"Input file"*: **FalseDiscoveryRate** output
>    - *"Name of the meta value to use as the new score"*: `Posterior Error Probability_score`
>    - *"Orientation of the new score"*: `lower_better`
>    - *"Name to use as the type of the new score"*: `Posterior Error Probability`
>
> 9. Run **FileInfo** {% icon tool %} to get basic information about the identified peptides
>    - {% icon param-file %} *"input file"*: **IDScoreSwitcher** output
>
> 10. Proceed with the protein inference as described [above](#protein-inference)
>
{: .hands_on}

> <question-title></question-title>
> 1. How many PSMs could be matched with XTandem and MSGFPlus alone? How many peptides were identified?
> 2. How many PSMs could be matched after combining the results with ConsensusID? How many peptides were identified?
>
> > <solution-title></solution-title>
> > 1.  After FDR-filtering, XTandem matched 3,552 PSMs (2,616 unique peptides) and MSGFPlus matched 4,292 PSMs (2,991 peptides).
> > 2.  Combining the results with ConsensusID leads to matching of 4,299 PSMs (3,041 unique peptides).
> {: .solution }
{: .question}


# Premade Workflow

- [A premade workflow for this tutorial]({% link topics/proteomics/tutorials/protein-id-oms/workflows/workflow.ga %}).
- [A premade workflow using the search engines XTandem and MSGF+]({% link topics/proteomics/tutorials/protein-id-oms/workflows/workflow_two-search-engines.ga %}).

# Further Reading

- [Protein inference](https://www.ncbi.nlm.nih.gov/pubmed/16009968)
- [Fido publication](https://www.ncbi.nlm.nih.gov/pubmed/20712337)
- [Evaluation of protein inference algorithms](https://www.ncbi.nlm.nih.gov/pubmed/27498275)
