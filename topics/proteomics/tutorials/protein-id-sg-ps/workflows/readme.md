# GalaxyP Workflow: Protein Identification (using Search GUI and Peptide Shaker)

This workflow supplements the GalaxyP tutorial on Protein ID ([Link]({{site.url}}/topics/proteomics/tutorials/protein-id-sg-ps/tutorial.html)).
You can find two versions of the same workflow: one is designed for a single MS run as an input. It can also be used for parallel analysis of multiple MS runs. Each MS run will result in a separate output file.
 The second workflow is designed for multiple MS runs to be combined into a single output file.

An overview of the workflow is given below:

![Protein ID Workflow](../../images/wf_proteinID_SG_PS.png)

## Inputs

Two inputs are needed:

1. A protein FASTA database to be searched against. Using the current settings, a database **without decoys** is needed. The decoys will be automatically added by ***Search GUI*** :wrench: . To learn more about databases, please consider the tutorial on [database handling]({{site.url}}/topics/proteomics/tutorials/database-handling/tutorial.html). For creating a database, you can also use a [ready-made workflow](../database-handling/).

2. At least one mass spectrometry data file in the mzML format. With the current settings, you will need a non-centroided mzML (raw data, no prior peak-picking). If you have data in another format or already centroided data, please consider the section [below](#customizing-the-workflow).

## Outputs

The workflow provides the identified proteins, peptides and PSMs as an output. For details on the ***Peptide Shaker*** :wrench: outputs, please consider the tutorial on [database handling]({{site.url}}/topics/proteomics/tutorials/database-handling/tutorial.html).

## Customizing the Workflow

You can customize the workflow after importing it to your Galaxy instance. Click on `Workflows`, choose this workflow and click on `Edit`.

- *Using a centroided mzML file*: Delete the tool ***PeakPickerHiRes*** :wrench: , use the mzML directly as an input for ***FileConverter*** :wrench: .
- *Using another MS file format than mzML*:

  - For `*.mgf`: Delete the tool ***PeakPickerHiRes*** :wrench: and ***FileConverter*** :wrench: , use the mzML directly as an input for ***Search GUI*** :wrench: .
  - For `*.raw`: Convert to mzML before running the workflow. For details, please consider the tutorial on [database handling]({{site.url}}/topics/proteomics/tutorials/database-handling/tutorial.html).
  - For other formats: Use the ***FileConverter*** :wrench: to convert to mzML (profile data) or directly to mgf (centroided data).
