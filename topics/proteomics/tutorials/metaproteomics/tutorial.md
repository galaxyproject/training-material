---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: metaproteomics
---

# Introduction

In this tutorial MS/MS data will be matched to peptide sequences provided through a FASTA file.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Placeholder](#pretreatments)
> 2. [Mapping](#mapping)
> 3. [Analysis of the differential expression](#analysis-of-the-differential-expression)
> 4. [Inference of the differential exon usage](#inference-of-the-differential-exon-usage)
{: .agenda}

# Analysis

## Data upload

There are a many ways how you can upload your data. Three among these are:

*   Upload the files from your computer.
*   Using a direct weblink.
*   Import from the data library if your instance provides the files.

In this tutorial, we will get the data from Zenodo. 

> ### :pencil2: Hands-on: Data upload and organization
>
> 1. Create a new history and name it something meaningful (e.g. *Metaproteomics tutorial*)
> 2. Import the three MGF MS/MS files and the FASTA sequence file from Zenodo
>
>    > ### :bulb: Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**    
>    {: .tip}
>
>    As default, Galaxy takes the link as name.
>
>    > ### :nut_and_bolt: Comments
>    > - Rename the datasets to a more descriptive name
>    {: .comment}
>
> 3. Build a **Dataset list** for the three MGF files
>    - Click the **Operations on multiple datasets** check box at the top of the history panel
>       ![](../../images/operations_icon.png)
>    - Check the three boxes next to the MGF files
>    - Click **For all selected...** and choose **Build dataset list**
>    - Ensure the three control samples are the only ones selected, and enter a name for the new collection (e.g. *MGF files*)
>    - Click **Create list** and exit by clicking again the dataset operations icon
>
{: .hands_on}

## Match MS/MS to peptide sequences

The search database labelled `FASTA_Bering_Strait_Trimmed_metapeptides_cRAP.FASTA` is the input database that
will be used to match MS/MS to peptide sequences via a sequence database search.
For this, the sequence database-searching program called **SearchGUI** will be used.
The created dataset collection of the three *MGF files* in the history is used as the MS/MS input.

#### SearchGUI

> ### :pencil2: Hands-on: SearchGUI
>
> 1. **SearchGUI** :wrench:: Run **SearchGUI** with:
>    - **Protein Database**: `FASTA_Bering_Strait_Trimmed_metapeptides_cRAP.FASTA`(or however you named the `FASTA` file)
>    - **Input Peak lists (mgf)**: `MGF files` dataset collection. 
>
>    > ### :bulb: Tip: Select dataset collections as input
>    >
>    > * Click the **Dataset collection** icon on the left of the input field:
>    >  
>    >      ![](../../images/dataset_button.png)
>    > * Select the appropriate dataset collection from the list
>    {: .tip}    
>    
>    Section **Search Engine Options**:
>
>    - **B-Search Engines**: `X!Tandem`
>
>    > ### :nut_and_bolt: Comment
>    >
>    > The section **Search Engine Options** contains a selection of sequence database searching
>    > programs that are available in SearchGUI. Any combination of these programs can be used for
>    > generating PSMs from MS/MS data. For the purpose of this tutorial, **X!Tandem** we will be used.
>    {: .comment}
>
>    Section **Protein Digestion Options**:
>
>    - **Digestion**: `Trypsin`
>    - **Maximum Missed Cleavages**: `2`
>
>    Section **Precursor Options**:
>
>    - **Precursor Ion Tolerance Units**: `Parts per million (ppm)`
>    - **Precursor Ion Tolerance:**: `10`
>    - **Fragment Tolerance (Daltons)**: `0.02`- this is high resolution MS/MS data
>    - **Minimum charge**: `2`
>    - **Maximum charge**: `6`
>    - **Forward Ion**: `b`
>    - **Reverse Ion**: `y`
>    - **Minimum Precursor Isotope**: `0`
>    - **Maximum Precursor Isotope**: `1`
>
>    Section **Protein Modification Options**:
>
>    - **Fixed Modifications**: `Carbamidomethylation of C`
>    - **Variable modifications**: `Oxidation of M`
>
>    > ### :bulb: Tip: Search for options
>    >
>    > * For selection lists, typing the first few letters in the window will filter the available options.
>    {: .tip}
>
>    Section **Advanced Options**:
>    - **X!Tandem Options**: `Advanced`
>    - **X!Tandem: Total Peaks**: `50`
>    - **X!Tandem: Min Peaks**: `15`
>    - **X!Tandem: Min Frag m/z**: `200`
>    - **X!Tandem: Min Precursor Mass**: `200`
>    - **X!Tandem: Noise Suppression**: `Yes`
>    - **X!Tandem: Dynamic Range**: `100`
>    - **X!Tandem: Quick Acetyl**: `No`
>    - **X!Tandem: Quick Pyrolidone**: `No`
>    - **X!Tandem: Protein stP Bias**: `No`
>    - **X!Tandem: Maximum Valid Expectation Value**: `100`
>    - **X!Tandem: Output Proteins**: `No`
>    - **X!Tandem: Output Sequences**: `No`
>    - **X!Tandem: Output Spectra**: `Yes`
>    - **X!Tandem peptide model refinement**: `Don't refine`
>    - leave everything else as default
>
> 2. Click **Execute**.
>
{: .hands_on}

Once the database search is completed, the SearchGUI tool will output a file (called a
SearchGUI archive file) that will serve as an input for the next section, PeptideShaker.

#### PeptideShaker

There are a number of choices for different data files that can be generated using
PeptideShaker. A compressed file can be made containing all information needed to view the
results in the standalone PeptideShaker viewer. A mzidentML file can be created that contains
all peptide sequence matching information and can be utilized by compatible downstream
software. Other outputs are focused on the inferred proteins identified from the PSMs, as well
as phosphorylation reports, relevant if a phosphoproteomics experiment has been undertaken.
The Certificate of Analysis (selected in our workflow) provides details on all the parameters
used by both SearchGUI and PeptideShaker in the analysis. This can be downloaded from the
Galaxy instance to your local computer in a text file if desired.
Most relevant for this tutorial is the PSM report that was selected in the workflow. Item #8 in
the History contains the PSM report generated for our workflow. Clicking on the item name in
the History will expand the item. Clicking on the View Data “eye” will open up the file for
viewing in the main Viewing Pane.