---
layout: tutorial_hands_on

title: MaxQuant and MSstats for the analysis of label-free data
zenodo_link: ''
level: Intermediate
questions:
- How to perform label-free shotgun (DDA) data analysis in MaxQuant and MSstats?
- Which proteins are differentially abundant in the two types of cutaneous squamous cell carcinomas?
objectives:
- Learn how to use MaxQuant and MSstats for the analysis of label-free shotgun (DDA) data
time_estimation: 2H
key_points:
- MaxQuant offers a single tool solution for protein identification and quantification.
- Label-free quantitation reveals the most abundant proteins in serum samples.
contributors:
- foellmelanie
- matthias313
requirements:
  -
    type: "internal"
    topic_name: proteomics
    tutorials:
      - maxquant-label-free
subtopic: id-quant
tags: [label-free]
---


# Introduction


Modern mass spectrometry-based proteomics enables the identification and quantification of thousands of proteins. Therefore, quantitative mass spectrometry represents an indispensable technology for biological and clinical research. Statistical analyses are required for the unbiased answering of scientific questions and to uncover all important information in the proteomic data. Classical statistical approaches and methods from other omics technologies are not ideal because they do not take into account the speciality of mass spectrometry data that include several thousands of proteins but often only a few dozens of samples (referred to as ‘curse of dimensionality’) and stochastic data properties that reflect sample preparation and spectral acquisition (Choi 2014).

In this training we will cover the full analysis workflow from label-free, data dependent acquisition (DDA) raw data to statistical results. We’ll use two popular quantitative proteomics software: MaxQuant and MSstats. MaxQuant allows protein identification and quantification for many different kinds of proteomics data (Cox and Mann 2008).  In case you have no previous experience with MaxQuant, we recommend to go through the [MaxQuant beginners tutorial]({{site.baseurl}}/topics/proteomics/tutorials/maxquant-label-free/tutorial.html) before. MSstats provides statistical functionalities to find differentially abundant peptides or proteins from data dependent acquisition (DDA), data independent acquisition (DIA) or single reaction monitoring (SRM) proteomic experiments.
The training dataset consists of a skin cancer cohort of 19 patients, which is a subset of a [published study](https://doi.org/10.1016/j.matbio.2017.11.004). One fifth of all non melanoma skin cancers are cutaneous squamous cell carcinomas (cSCC) that mainly derive from exposure to ultraviolet light. Most cSCC have a good prognosis but the few metastasizing cSCC have dramatically increased mortality. Here, we compare these metastasizing cSCC to cSCC in patients with the genetic disease recessive dystrophic epidermolysis bullosa (RDEB). RDEB is a genetic skin blistering and extracellular matrix disease caused by collagen VII deficiency. To investigate molecular differences between these two aggressive cSCCs with different origin, we used global proteomic analysis of formalin-fixed paraffin-embedded human cSCC tissues.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data


The annotation file, group comparison file and FASTA file for this training is deposited at [Zenodo](https://zenodo.org/record/4896554). It is of course possible to use another FASTA file with human proteome sequences, but to ensure that the results are compatible we recommend to use the provided FASTA file. MaxQuant not only adds known contaminants to the FASTA file, but also generates the “decoy” hits for false discovery rate estimation itself, therefore the FASTA file is not allowed to have decoy entries. To learn more about FASTA files, have a look at [Protein FASTA Database Handling tutorial]({{site.baseurl}}/topics/proteomics/tutorials/database-handling/tutorial.html). The raw data is available via the [PRIDE repository](https://www.ebi.ac.uk/pride/archive/projects/PXD006914). As this is a real life study, the raw data sizes are large and computation time in MaxQuant is long. To save time and storage capacity, you can skip downloading the raw data and the MaxQuant run and instead continue with the MaxQuant outputs which we provide later on. In this case skip the data upload steps 5-8 which are only necessary for the MaxQuant run.


> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this tutorial and give it a meaningful name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the FASTA database, annotation file and comparison matrix from [Zenodo](https://zenodo.org/record/4896554)
>
>    ```
>    https://zenodo.org/record/4896554/files/input_protein_database.fasta
>    https://zenodo.org/record/4896554/files/input_annotation_file.tabular
>    https://zenodo.org/record/4896554/files/input_comparison_matrix.tabular
>    ```
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Once the files are green, rename the fasta file into 'protein database', the annotation file into 'annotation file' and the comparison matrix file into 'comparison matrix'.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Steps 4 to 7 can be skipped to save time and storage capacity by not running MaxQuant. To run MaxQuant, import the raw data from [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD006914).
>
>    ```
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment105_metast_cSCC1.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment106_metast_cSCC2.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment107_metast_cSCC3.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment109_metast_cSCC4.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment116_metast_cSCC5.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment117_metast_cSCC6.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment118_metast_cSCC7.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment119_metast_cSCC8.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment120_metast_cSCC9.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment121_metast_cSCC10.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment122_metast_cSCC11.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment123_metast_cSCC12.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment124_metast_cSCC13.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment110_RDEB_cSCC1.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment111_RDEB_cSCC2.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment112_RDEB_cSCC3.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment113_RDEB_cSCC4.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment114_RDEB_cSCC5.raw
>    ftp://ftp.pride.ebi.ac.uk/pride-archive/2017/11/PXD006914/Experiment115_RDEB_cSCC6.raw
>
>    ```
> 5. Rename the raw datasets into 'metast_cSCC1.raw', 'metast_cSCC2.raw', etc.. The naming for the raw files have to be exactly this way to later match the file names provided in the MSstats annotation file. 
> 
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 6. Control that the data type of the raw files is 'thermo.raw' otherwise change the datatype into 'thermo.raw'
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="thermo.raw" %}
>
> 7. Generate a collection for all raw files and name it 'raw_files', hide the individual raw files
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
{: .hands_on}

# MaxQuant analysis

The run time of **MaxQuant** {% icon tool %} depends on the number and size of the input files and on the chosen parameters. The run of the training datasets will take a few hours, but the training can be directly continued with the MaxQuant result files from Zenodo. We start the MaxQuant run with the default parameters, with a few adjustments. Protein level quantification parameters do not really matter here, because MSstats will use feature quantifications and perform protein summarization based on them. A quality control report is generated with the [PTXQC functionality](https://pubs.acs.org/doi/10.1021/acs.jproteome.5b00780) that is directly implemented in the MaxQuant Galaxy tool. To continue with statistical analysis in MSstats, the Protein Groups and the Evidence files are needed from MaxQuant.

> <hands-on-title>Optional: MaxQuant analysis</hands-on-title>
>
> 1. {% tool [MaxQuant](toolshed.g2.bx.psu.edu/repos/galaxyp/maxquant/maxquant/1.6.10.43+galaxy3) %} with the following parameters:
>    - In *"Input Options"*:
>        - {% icon param-file %} *"FASTA files"*: `protein database`
>    - In *"Search Options"*:
>    - *"minimum unique peptides"*: `1`
>    - *"Match between runs"*: `yes`
>    - In *"Parameter Group"*:
>        - {% icon param-collection %} *"Infiles"*: `raw_files`
>    - *"Generate PTXQC (proteomics quality control pipeline) report?"*: `Yes`
>    - In *"Output Options"*:
>        - *"Select the desired outputs."*: `Protein Groups` `Evidence`
>
{: .hands_on}

Because the MaxQuant run takes really long, we recommend to download the MaxQuant results from Zenodo and continue with the tutorial. 

> <hands-on-title>Load MaxQuant results from Zenodo</hands-on-title>
>
> 1. Import the files from [Zenodo](https://zenodo.org/record/4896554)
>
>    ```
>    https://zenodo.org/record/4896554/files/MaxQuant_Evidence.tabular
>    https://zenodo.org/record/4896554/files/MaxQuant_proteingroups.tabular
>    https://zenodo.org/record/4896554/files/PTXQC_report.pdf
>    ```
{: .hands_on}


> <question-title></question-title>
>
> 1. How many proteins and features were identified in total?
> 2. In which columns (number) are the potential contaminants in the protein group and evidence file respectively? 
> 3. How large is the proportion of potential contaminants?
>
> > <solution-title></solution-title>
> >
> > 1. 2622 protein groups and ~240000 features were found in total (number of lines of protein group and evidence files)
> > 2. They are in column 118 (protein groups) and 54 (evidence)
> > 3. Up to 60% of the samples intensities derive from potential contaminants (PTXQC plots page 7)
> >
> {: .solution}
>
{: .question}


# MSstats analysis

The protein groups and evidence files of MaxQuant can directly be input into MSstats. MSstats automatically removes all proteins that are labelled as contaminants ('+' sign in the column 'potential contaminant' of both MaxQuant outputs). However, in this skin dataset we expect that the skin proteins are part of the sample and not a contamination. Therefore, we keep all human contaminants by first removing non human proteins (by selecting only lines that contain the word 'HUMAN' or a word from the header line) and then for the human potential contaminants remove the '+' (replacing '+' with ''(empty field) in the 'potential contaminant' column) to keep them for the analysis.

We use the modified MaxQuant protein groups and evidence files as input in MSstats. In addition, an annotation file that describes the experimental design and a comparison matrix is needed. Please start the MSstats run first and while it is running you can find more details on its parameters below.

> <hands-on-title>MSstats Analysis</hands-on-title>
>
> 1. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `Protein Groups` (output of **MaxQuant** {% icon tool %})
>    - *"the pattern"*: `(HUMAN)|(Majority)`
> 2. {% tool [Select](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `Evidence` (output of **MaxQuant** {% icon tool %})
>    - *"the pattern"*: `(HUMAN)|(Sequence)`
> 3. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `select protein groups` (output of **Select** {% icon tool %})
>    - In *"Replacement"*:
>            - *"in column"*: `c118`
>            - *"Find pattern"*: `+`
>    - Once finished, rename the file into `protein groups input for MSstats`
> 4. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `select evidence` (output of **Select** {% icon tool %})
>    - In *"Replacement"*:
>            - *"in column"*: `c54`
>            - *"Find pattern"*: `+`
>    - Once finished, rename the file into `evidence input for MSstats`
> 5. {% tool [MSstats](toolshed.g2.bx.psu.edu/repos/galaxyp/msstats/msstats/4.0.0.0) %}  with the following parameters:
>    - *"input source"*: `MaxQuant`
>        - {% icon param-file %} *"evidence.txt - feature-level data"*: `evidence input for MSstats` (output of **Replace** {% icon tool %})
>        - {% icon param-file %} *"proteinGroups.txt - protein-level data"*: `protein groups input for MSstats` (output of **Replace** {% icon tool %})
>        - {% icon param-file %} *"annotation file"*: `annotation file`
>        - *"Select Protein ID in evidence.txt"*: `Leading razor protein column`
>        - In *"MaxQtoMSstatsFormat Options"*:
>            - *"Remove the proteins which have only 1 peptide and charge"*: `Yes`
>    - In *"dataProcess Options"*:
>        - *"Select outputs"*: `MSstats log` `MSstats FeatureLevelData` `ProteinLevelData` `Group Quantification Matrix Table` `Sample Quantification Matrix Table`
>        - In *"dataProcess Plot Options"*:
>            - *"Select visualization outputs"*: `QCPlot`
>            - *"Select protein IDs to draw plots"*: `Option for QC plot: "allonly" will generate one QC plot with all proteins`
>    	     - In *"Advanced visualization parameters"*:
>                - *"Angle of labels represented each condition at the top of graph"*: `0`
>    - *"Compare Groups"*: `Yes`
>        - {% icon param-file %} *"Comparison Matrix"*: `comparison_matrix`
>        - *"Select outputs"*: `MSstats ComparisonResult`
>        - In *"Comparison Visualization Options"*:
>            - *"Select visualization outputs"*: `MSstats VolcanoPlot`
>    	     - In *"Advanced visualization parameters"*:
>                - *"Involve fold change cutoff or not for volcano plot or heatmap"*: `1.5`
>                - *"Display protein names in Volcano Plot"*: `No`
>
{: .hands_on}


> <question-title></question-title>
>
> 1. How many proteins were removed as potential non human contaminants?
> 2. How many proteins were included into the statistical analysis?
>
> > <solution-title></solution-title>
> >
> > 1. 28 (2622 lines in protein group file minus 2594 lines after select)
> > 2. 1764 (MSstats log)
> >
> {: .solution}
>
{: .question}



## More details on MSstats

MSstats  is designed for statistical modelling of mass spectrometry based proteomic data [Choi 2014](https://doi.org/10.1093/bioinformatics/btu305 ).
Proteomic data analysis requires statistical approaches that reduce bias and inefficiencies and distinguish systematic variation from random artifacts [Käll and Vitek 2011]( https://doi.org/10.1371/journal.pcbi.1002277).

MSstats is directly compatible with the output of several quantitative proteomics software. In addition to the results of the proteomics software an annotation file is needed as input. The annotation file describes the experimental design such as the conditions, biological and technical replicates. To be compatible with MaxQuant results, an additional column with the label type is needed, which only contains L (light) in DDA experiments. A wrong setup of the annotation file is the most common source of errors in MSstats, thus we collected more information in the box below to allow you to adjust the annotation file when analyzing your own experiments. 

> <tip-title>Generating the MSstats annotation file</tip-title>
>
> For label-free MaxQuant data, the annotation file should have 5 columns with exactly these headers: Raw.file, Condition, BioReplicate, Run; IsotopeLabelType
> 1. Raw.file: The names must match exactly to the file names in the MaxQuant evidence.txt "Raw file" column. (e.g. "file1.raw.thermo"). 
> 2. Condition: The conditions which will be compared in the statistical modelling. They are not allowed to start with a number or contain any special characters except for '_'.
> 3. BioReplicate: This column should contain a unique identifier for each biological replicate in the experiment. For example, in a clinical proteomic investigation this should be a unique patient id. If technical replicates are present, all samples from the same biological replicate should have the same id but different run ids. MSstats automatically detects the presence of technical replicates and accounts for them in the model-based analysis.
> 4. Run: This column contains the identifier of a mass spectrometry run. Each mass spectrometry run should have a unique identifier (number or name), regardless of the origin of the biological sample. 
> 5. IsotopeLabelType: This is L (light) for all MaxQuant DDA experiments.
>
{: .tip}

MSstats will compare all conditions that are indicated in the comparison matrix. The comparison matrix has to be setup correctly to avoid errors and wrong statistical modelling. It contains a first column that gives the comparisons a name and one column per condition. This matrix is filled with 1 and -1 to specify the conditions that are compared in each comparison and with 0 for conditions that are not part of the comparison.

> <tip-title>Generating the MSstats comparison matrix</tip-title>
>
> 1. The first column of the comparison matrix contains the names of the comparisons and should have 'names' as header. These names will be used in all MSstats output files, therefore it is important that the names are meaningful and reflect the actual comparison (see below)
> 2. An additional column for each condition that is present in the data. This means each condition present in the annotation file has to be a separate column even when the condition will not be used for any comparison. The header should contain the condition name exactly as written in the annotation file. 
> 3. Fill the matrix: Use 1 and -1 to indicate the conditions to compare and 0 for conditions that are not compared. Multiple groups can be combined by using 0.5. 
> 4. Example: to compare condition1 with condition2: write 1 into the column of condition1 and write -1 in the column of condition2. In the first column of the matrix name this comparison condition1-condition2. The naming of the comparison should reflect the direction of the comparison and thus always have the condition that is set to 1 first and the condition that is set to -1 second (condition1-condition2 and NOT condition2-condition1). 
>
{: .tip}

The first analysis step in MSstats is the conversion of the input data into an MSstats compatible table. For this step several parameters to filter and adjust the input data can be selected. We keep the default parameters and only change one parameter in order to remove proteins which have only a single peptide measurement.
Next, data processing optimizes the data for statistical modelling via log-transformation, and normalization of intensities, feature selection, missing value imputation, and run-level summarization.
Log- transformation is performed to transform multiplicative signals to additive signals which are compatible with linear statistical models and bring the intensity distribution close to a normal distribution. Furthermore, it changes the dependence of variances from the intensity values: in the raw data larger intensities have larger variances but after log transformation lower intensities have larger variances.

Normalization aims to make the intensities of different runs more comparable to each other. The default normalization method, equalize medians, assumes that the majority of proteins do not change across runs and shifts all intensities of a run by a constant to obtain equal median intensities across runs.
A feature in label-free DDA data corresponds to a peptide at a given charge state (m/z value), resulting from the identification of MS2 spectra combined with the quantitative information from the MS1 scans. Feature selection allows the use of either all, only the most abundant features or only high quality peptides for protein summarization.

Missing values and noisy features with outliers are typical in label-free DDA datasets but influence protein summarization. Therefore, it is recommended to perform missing value imputation. Missing values are reported differently in different Softwares. MaxQuant reports them as NA and MSstats assumes that missing intensity values from MaxQuant mean that the intensity was below the limit of quantification. This means the values are not missing for random but for the reason of low abundance. Therefore, the values are only partially known and  called “censored”. This may also be the case for very low intensity values, which might not be reliable. The percentile that is not trusted and should be considered a censored value is defined via the “Maximum quantile for deciding censored missing values” parameter. Censored values are replaced by an intensity that is generated via an accelerated failure time model (AFT). Alternatively censored values may be replaced by the minimum value of the features, runs or both as defined in the “Cutoff value for censoring”. Runs with no intensity measurement for a protein will be removed for any further calculation on this protein.

Protein summarization is by default performed via Tukey’s median polish for robust parameter estimation with median across rows and columns. Run-level summaries are later used for statistical group comparison.

Any two groups can be compared to find differentially abundant proteins between them. MSstats uses a family of linear mixed models that are automatically adjusted for the comparison type according to the information in the annotation file, such as conditions, biological and technical replicates and runs. This allows comparison of groups with different sizes; comparison of the mean of some groups, paired designs and time course experiments.


# Follow up on MSstats results

We obtain several output files from MSstats. MSstats log file contains the MSstats report with warnings and information about the analysis steps.
The MSstats QCPlot visualizes the log transformed intensities for all proteins and runs.
The volcano plot plots the statistical result as transformed p-values vs. the log2 fold change. A fold change of 1.5 means that a protein is 50% more abundand in one condition than the other. The log2 fold change is 0.58.

![qc plot](../../images/maxquant-msstats-lfq/qc_plot.png "QC plot of all proteins")

![volcano plot](../../images/maxquant-msstats-lfq/volcano_plot.png "Volcano plot showing p-values and log2 fold changes for all proteins. Dashed line indicates p-value of 0.05 and log2 fold change of ± 0.58")

The FeatureLevelData file contains the transformed, normalized and imputed intensities for each peptide in each run. ProteinLevelData data summarizes intensities per run and protein.
We’ll count and visualize the number of features per run and calculate the distribution of proteins per sample.

> <hands-on-title>Follow up on MSstats results</hands-on-title>
>
> 1. {% tool [Summary Statistics](Summary_Statistics1) %} with the following parameters:
>    - {% icon param-file %} *"Summary statistics on"*: `ProteinLevelData` (output of **MSstats** {% icon tool %})
>    - *"Column or expression"*: `c8`
> 2. {% tool [Datamash](toolshed.g2.bx.psu.edu/repos/iuc/datamash_ops/datamash_ops/datamash_ops/1.0.6) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `ProteinLevelData` (output of **MSstats** {% icon tool %})
>    - *"Group by fields"*: `4`
>    - *"Input file has a header line"*: `Yes`
>    - *"Print header line"*: `Yes`
>    - *"Sort input"*: `Yes`
>    - In *"Operation to perform on each group"*:
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"On column"*: `c1`
> 3. Click on {% icon galaxy-barchart %} “Visualize this data” on the **Datamash** {% icon tool %} result.
>   - Select `Bar diagram (NVD3)`
>   - *"Provide a title"*: `Number of features per sample`
>   - Click `Select data` {% icon galaxy-chart-select-data %}
>   - *"Data point labels"*: `Column: 1`
>   - Save {% icon galaxy-save %} (file is saved under "User" --> "Visualizations")
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Which sample has the lowest amount of proteins after protein summarization?
> 2. In the complete experiment, how many features has a protein on average?
>
> > <solution-title></solution-title>
> >
> > 1. RDEB cSCC4 
> > ![Number of proteins per sample](../../images/maxquant-msstats-lfq/features_sample.png "Number of proteins per sample (run)")
> > 2. Around 6 features per protein (mean in summary statistics). 
> >
> {: .solution}
>
{: .question}


# Filtering MSstats results

The comparison result table summarizes the statistical results per protein and comparison. First, we keep only the Uniprot ID in column 1 to make the ID less cluttered. This is done by deleting everything before the first pipe '|' and everything after the second pipe '|'. 
Then we keep only statistically significant proteins that means they have an adjusted p-value below 0.05.
Next, we separate up- and down-regulated proteins by filtering for a positive and negative log2FC.
The Sample Quantification Matrix Table contains the summarized intensities per protein and sample.
In order to make its IDs compatible with the ones from the comparison result at a later step, we keep only the Uniprot ID as well.


> <hands-on-title>Filtering MSstats results</hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Comparison Result` (output of **MSstats** {% icon tool %})
>    - In *"Replacement"*:
>            - *"in column"*: `c1`
>            - *"Find pattern"*: `sp\|`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"in column"*: `c1`
>            - *"Find pattern"*: `\|.*`
> 2. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `replaced comparison result` (output of **Replace Text** {% icon tool %})
>    - *"With following condition"*: `c8<0.05`
>    - *"Number of header lines to skip"*: `1`
>    - Once finished, rename the file into `significant proteins`
> 3. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `significant proteins` (output of **Filter** {% icon tool %})
>    - *"With following condition"*: `c3>0.58`
>    - *"Number of header lines to skip"*: `1`
> 4. Add a tag `#metastasized` to the filtered file and rename it into `metastasized filtered`
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
> 5. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `significant proteins` (output of the first **Filter** {% icon tool %})
>    - *"With following condition"*: `c3<-0.58`
>    - *"Number of header lines to skip"*: `1`
> 6. Add a tag `#rdeb` to the filtered file and rename it into `rdeb filtered`
>
{: .hands_on}


> <question-title></question-title>
>
> 1. Why do we filter for the adjusted p-value?
> 2. How many proteins have a p-value below 0.05?
>
> > <solution-title></solution-title>
> >
> > 1. Adjusted p-values control for the multiplicity of testing. Since we fit a separate model, and conduct a separate comparison for each protein, the number of tests equals the number of comparisons. A 0.05 cutoff of an adjusted p-value controls the False Discovery Rate in the collection of tests over all the proteins at 5%. Since they account for the multiplicity, adjusted p-values are more conservative (i.e. it is more difficult to detect a change).
> > 2. 148 in total (first filtering step); 133 are upregulated in metastasized cSCC (metastasized filtered) and 13 are upregulated in RDEB cSCC (rdeb filtered).
> >
> {: .solution}
>
{: .question}

# Finding differentially abundant proteins

For each condition we select only the significant proteins, which are proteins with a p-value above 0 and below 0.05. Proteins with a p-value of 0 are missing in one condition and are therfore discarded in the next steps. We'll keep only the column with the Uniprot ID and extract the average protein intensities per sample from the sample quantification matrix file and vizualize them at heatmap. We do the exact same steps for both conditions, therefore, each time you start a tool (except for the replace step) you can use the multiple input file to start the step for the metastasized and rdeb files at the same time. 

> <hands-on-title>filter differentially abundant proteins</hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `metastasized filtered` (output of **Filter** {% icon tool %})
>    - *"With following condition"*: `c8>0`
>    - *"Number of header lines to skip"*: `1`
>    - Rename the file into `significant metastasized`
> 2. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `significant metastasized ` (output of last **Filter** {% icon tool %})
> 3. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `rdeb filtered` (output of **Filter** {% icon tool %})
>    - *"With following condition"*: `c8>0`
>    - *"Number of header lines to skip"*: `1`
>    - Rename the file into `significant rdeb`
> 4. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1`
>    - {% icon param-file %} *"From"*: `significant rdeb` (output of last **Filter** {% icon tool %})
> 5. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_column/1.1.3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Sample Quantification Matrix Table` (output of **MSstats** {% icon tool %})
>    - In *"Replacement"*:
>            - *"in column"*: `c1`
>            - *"Find pattern"*: `sp\|`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"in column"*: `c1`
>            - *"Find pattern"*: `\|.*`
> 6. {% tool [Join](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_easyjoin_tool/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"1st file"*: `replaced sample quantification matrix` (output of **Replace text** {% icon tool %})
>    - *"Column to use from 1st file"*: `c1`
>    - {% icon param-file %} *"2nd File"*: `metastasized cut` (output of **Cut** {% icon tool %})
>    - *"Column to use from 2nd file"*: `c1`
>    - *"First line is a header line"*: `Yes`
> 7. {% tool [Join](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_easyjoin_tool/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"1st file"*: `replaced sample quantification matrix` (output of **Replace Text** {% icon tool %})
>    - *"Column to use from 1st file"*: `c1`
>    - {% icon param-file %} *"2nd File"*: `rdeb cut` (output of **Cut** {% icon tool %})
>    - *"Column to use from 2nd file"*: `c1`
>    - *"First line is a header line"*: `Yes`
> 8. {% tool [Heatmap2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_heatmap2/ggplot2_heatmap2/3.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Input should have column headers - these will be the columns that are plotted"*: `metastasized join` (output of **Join** {% icon tool %})
>    - *"Plot title"*: `Upregulated proteins in metastasized cSCC`
>    - *"Enable data clustering"*: `No`
>    - *"Data scaling"*: `Scale my data by row`
> 9. {% tool [Heatmap2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_heatmap2/ggplot2_heatmap2/3.0.1) %} with the following parameters:
>    - {% icon param-file %} *"Input should have column headers - these will be the columns that are plotted"*: `rdeb join` (output of **Join** {% icon tool %})
>    - *"Plot title"*: `Upregulated proteins in RDEB cSCC`
>    - *"Enable data clustering"*: `No`
>    - *"Data scaling"*: `Scale my data by row`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many proteins are differentially abundant?
>
> > <solution-title></solution-title>
> >
> > 1. 85 are upregulated in metastasized cSCC and 12 are upregulated in RDEB cSCC (number of lines minus 1 in the first filtering step per condition).
> >
> {: .solution}
>
{: .question}

# Follow up on differentially abundant proteins

In addition we retrieve for each Uniprot ID the corresponding protein names from uniprot to allow an easier interpretation.

> <hands-on-title>MSstats visualizations</hands-on-title>
>
> 1. {% tool [UniProt ID mapping and retrieval](toolshed.g2.bx.psu.edu/repos/bgruening/uniprot_rest_interface/uniprot/0.2) %} with the following parameters:
>    - {% icon param-file %} *"Input file with IDs"*: `metastasized join` (output of **Join** {% icon tool %})
>    - *"ID column"*: `c1`
>    - *"Do you want to map IDs or retrieve data from UniProt"*: `Retrieve: request entries by uniprot accession using batch retrieval`
> 2. {% tool [UniProt ID mapping and retrieval](toolshed.g2.bx.psu.edu/repos/bgruening/uniprot_rest_interface/uniprot/0.2) %} with the following parameters:
>    - {% icon param-file %} *"Input file with IDs"*: `rdeb join` (output of **Join** {% icon tool %})
>    - *"ID column"*: `c1`
>    - *"Do you want to map IDs or retrieve data from UniProt"*: `Retrieve: request entries by uniprot accession using batch retrieval`
> 3. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `metastasized uniprot` (output of **UniProt** {% icon tool %})
>    - *"How many columns to divide title string into?"*: `2`
> 4. {% tool [FASTA-to-Tabular](toolshed.g2.bx.psu.edu/repos/devteam/fasta_to_tabular/fasta2tab/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Convert these sequences"*: `rdeb uniprot` (output of **UniProt** {% icon tool %})
>    - *"How many columns to divide title string into?"*: `2`
>
{: .hands_on}

Three of the differentially abundant proteins found here were also found and stained with antibodies in the original publication: Collagen XIV which is higher in RDEB cSCC than in metastasizing cSCC and Serum amyloid P-component as well as X-ray repair cross-complementing protein 6, which are both higher in metastasizing cSCC than in RDEB cSCC. Collagen XIV is a fibril associated collagen which may have tissue stabilizing function in the dermis. The upregulation of collagen XIV as well as other collagens in RDEB cSCC could be a compensation effort for the impaired collagen VII in RDEB tissues. Collagen VII is indeed only found in some RDEB samples and with such low intensities that it appears as upregulated in metastasized cSCC. 

![col14 staining](../../images/maxquant-msstats-lfq/col14_ihc.png "Immunoflourescence staining of collagen XIV in RDEB and metastasizing cSCC skin tissues")
