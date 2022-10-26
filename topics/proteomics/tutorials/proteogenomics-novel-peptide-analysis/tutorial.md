---
layout: tutorial_hands_on

title: "Proteogenomics 3: Novel peptide analysis"
zenodo_link: "https://doi.org/10.5281/zenodo.1489208"
level: Intermediate
questions:
  - "How to verify the spectra of novel proteoforms?"
  - "How to assign genomic allocation to these novel proteoforms?"
objectives:
  - "How to assign and visualize the genomic localization of these identified novel proteoforms?"
time_estimation: "30m"
key_points:
  - "Learning how to visualize proteomic data and to perform its genomic allocation"
follow_up_training:
  -
    type: "internal"
    topic_name: proteomics
    tutorials:
      - proteogenomics-dbcreation
      - proteogenomics-dbsearch

contributors:
  - subinamehta
  - timothygriffin
  - pratikdjagtap
  - jraysajulga
  - jj-umn
  - pravs3683
subtopic: multi-omics
tags: [proteogenomics]
---

# Introduction
{: .no_toc}

The third and the last workflow in the proteogenomics tutorial is to identifying the "**Novel peptides**" using BlastP and to localize the peptides to its genomic coordinates. Inputs from both workflow 1 and 2 will be used in this workflow.
Please look at the following tutorials in this proteogenomics series before starting this tutorial:
1. [Proteogenomics database creation]({% link topics/proteomics/tutorials/proteogenomics-dbcreation/tutorial.md %})
2. [Proteogenomics database search]({% link topics/proteomics/tutorials/proteogenomics-dbsearch/tutorial.md %})

![Workflow](../../images/Third_workflow.png)

> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
> 1. TOC
> {:toc}
>
{: .agenda}

# Pretreatments

{: .no_toc}

All the files to run this workflow can be obtained from the [second tutorial]({% link topics/proteomics/tutorials/proteogenomics-dbsearch/tutorial.md %}) output. Once the tabular output is generated, we convert this tabular report into a FASTA file. This can be achieved by using the Tabular to FASTA convertion tool.

> <hands-on-title>data organization</hands-on-title>
>
> 1. The inputs for this workflow are:
>    - **Tabular file** – “**Peptides for BlastP analysis**”
>    - **Tabular file** – “**PeptideShaker_PSM**”
>    - **Mz to sqlite**
>    - **Genomic mapping sqlite**
>
> If you do not have these files from the previous tutorials in this series, you can import them from Zenodo:
> ```
> https://zenodo.org/record/1489208/files/Peptides_for_Blast-P_analysis.tabular
> https://zenodo.org/record/1489208/files/PeptideShaker_PSM.tabular
> https://zenodo.org/record/1489208/files/mz_to_sqlite.mz.sqlite
> https://zenodo.org/record/1489208/files/genomic_mapping_sqlite.sqlite
> ```
{: .hands_on}



# Peptide Selection

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is a web based tool used to compare biological sequences. BlastP, matches protein sequences against a protein database. More specifically, it looks at the amino acid sequence of proteins and can detect and evaluate the amount of differences between say, an experimentally derived sequence and all known amino acid sequences from a database. It can then find the most similar sequences and allow for identification of known proteins or for identification of potential peptides associated with novel proteoforms.

The first step in this tutorial is to perfrom BLAST-P analysis using the NCBI-NR database. The output from BLASTP will determine the identification of the novel peptides. The result is a tabular file with 25 columns containing all the information regarding the alignment of these peptides with the sequences in the NCBI-NR database.

> <hands-on-title>NCBI BLAST+ blastp</hands-on-title>
>
> 1. {% tool [NCBI BLAST+ blastp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/0.3.3) %} with the following parameters:
>    - {% icon param-file %} **Protein query sequence(s)** - `Peptides for Blast-P analysis.tabular`
>    - {% icon param-select %} **Subject database/sequences** - `Locally installed BLAST database`
>      - {% icon param-select %} **Protein BLAST database** - `NCBI-NR(dated)`
>    - {% icon param-select %} **Type of BLAST** - `blast-p short`
>    - {% icon param-text %} **Set expectation value cutoff** - `200000.0`
>    - {% icon param-select %} **Output format** - `Tabular format (25 columns)`
>    - {% icon param-select %} **Advanced Options** - `show advanced options`
>      - {% icon param-check %} **Filter out low complexity regions (with SEG)** - `No`
>      - {% icon param-select %} **Scoring matrix and gap costs** - `PAM30`
>         - {% icon param-select %} **Gap Costs** - `Extension:9 Extension:1`
>      - {% icon param-text %} **Maximum hits to consider/show** - `1`
>      - {% icon param-text %} **Maximum number of HSPs (alignments) to keep for any single query-subject pair** - `1`
>      - {% icon param-text %} **Word size for wordfinder algorithm** - `2`
>      - {% icon param-text %} **Multiple hits window size: use 0 to specify 1-hit algorithm, leave blank for default** - `40`
>      - {% icon param-text %} **Minimum score to add a word to the BLAST lookup table** - `11`
>      - {% icon param-select %} **Composition-based statistics** - `0: no composition-based statistics`
>      - {% icon param-check %} **Should the query and subject defline(s) be parsed?** - `No`
>      - {% icon param-select %} **Restrict search of database to a given set of ID's** - `No restriction, search the entire database`
>      - {% icon param-check %} **Minimum query coverage per hsp (percentage, 0 to 100)?** - `0`
>      - {% icon param-check %} **Compute locally optimal Smith-Waterman alignments** - `No`
>  2. Click **Execute** and inspect the query results file after it turned green.
>
{: .hands_on}

Once Blast-P search is performed, it provides a tabular output containing “**Novel peptides**”. Now this output is further processed by comparing the Novel Peptide output with the PSM report for selecting only distinct peptides which meet the criteria.


> <hands-on-title>Query Tabular</hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.0.0) %} with the following parameters:
>    - {% icon param-repeat %} **Insert Database Table**
>      - Section **Table Options**
>        - *"Specify Name for Table"*: `blast`
>        - *"Use first line as column names"* : `No`
>        - *"Specify Column Names (comma-separated list)"*:
>`qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,sallseqid,score,nident,positive,gaps,ppos,qframe,sframe,qseq,sseq,qlen,slen,salltitles`
>        - *"Only load the columns you have named into database"*: `Yes`
>        - {% icon param-repeat %} **Insert Table Index**
>          - *"Table Index"*: `No`
>          - *"Index on Columns"*: `qseqid`
>
>    - {% icon param-repeat %} **Insert Database Table**
>      - Section **Filter Dataset Input**
>        - {% icon param-repeat %} **Filter Tabular Input Lines**
>          - *"Filter by"*:  `skip leading lines`
>          - *"Skip lines"*: `1`
>      - Section **Table Options**
>        - *"Specify Name for Table"*: `psm`
>        - *"Use first line as column names"* : `No`
>        - *"Specify Column Names (comma-separated list)"*: `ID,Proteins,Sequence,AAs_Before,AAs_After,Position,Modified_Sequence,Variable_Modifications,Fixed_Modifications,Spectrum_File,Spectrum_Title,Spectrum_Scan_Number,RT,mz,Measured_Charge,Identification_Charge,Theoretical_Mass,Isotope_Number,Precursor_mz_Error_ppm,Localization_Confidence,Probabilistic_PTM_score,Dscore,Confidence,Validation`
>        - *"Only load the columns you have named into database"*: `Yes`
>
>
>    - *"Save the sqlite database in your history"*: `No`
>
>       > <comment-title>Querying an SQLite Database</comment-title>
>       >
>       > **Query Tabular** can also use an existing SQLite database. Activating `Save the sqlite database in your history`
>       > will store the created database in the history, allowing to reuse it directly.
>       >
>       {: .comment}
>
>    - *"SQL Query to generate tabular output"*:
>      ```
>      SELECT distinct psm.*
>      FROM psm join blast ON psm.Sequence = blast.qseqid
>      WHERE blast.pident < 100 OR blast.gapopen >= 1 OR blast.length < blast.qlen
>      ORDER BY psm.Sequence, psm.ID
>      ```
>
>       > <comment-title>Query information</comment-title>
>       >
>       > The query wants a tabular list of peptides in which the lenght of the PSM sequence is equal to the length of the Blast sequence, where in the pident (percentage identity) is less that 100 i.e. Peptide cannot be a 100% identical to the NCBI-nr reference database. Or it should fulfill the criteria that there should be atleast 1 gap present (blast.gapopen >= 1) or the length of the peptide in NCBI-nr should be less than the length of the query length. If the peptide follows all this then it is accepted as a "Novel" proteoform.
>       >
>       {: .comment}
>    - *"Include query result column headers"*: `Yes`
>
>  2. Click **Execute** and inspect the query results file after it turned green.
>
{: .hands_on}

Once this step is completed, a tabular output containing novel proteoforms are displayed. These novel proteforms fulfill our criteria of not being present in the existing NCBI repository. The next step is to remove any duplicate sequences. For this, we use the Query tabular tool again to select distinct sequences from the tabular output.

> <hands-on-title>Query Tabular</hands-on-title>
>
> 1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.0.0) %}
>    - {% icon param-repeat %} **Insert Database Table**
>      - Section **Filter Dataset Input**
>        - {% icon param-repeat %} **Insert Filter Tabular Input Lines**
>          - *"Filter by"*:  `skip leading lines`
>          - *"Skip lines"*: `1`
>      - Section **Table Options**
>        - *"Specify Name for Table"*: `psm`
>        - *"Use first line as column names"* : `No`
>        - *"Specify Column Names (comma-separated list)"*:`ID,Proteins,Sequence`
>        - *"Only load the columns you have named into database"*: `Yes`
>
>    - *"SQL Query to generate tabular output"*:
>      ```
>      SELECT distinct Sequence from psm
>      ```
>    - *"include query result column headers"*: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green.
>
{: .hands_on}


# Multiomics Visualization Platform (MVP)

The Multiomics Visualization Platform is a Galaxy visualization plugin that allows the user to browse the selected proteomics data. It uses the SQlite database which allows the data to be filtered and aggregated in a user defined manner. It allows various features such as; the PSM can be displayed with a lorikeet spectral view, the selected peptide can be displayed in a protein view and an IGV browser is also available for the selected protein. The step by step guide shown below will provide a walkthrough on how to use this plugin (NOTE: the example shown below is a representative peptide which is subjected to change, so while you are running this tool please take a look at the "Novel Peptide" output from the previous steps).

> <hands-on-title>Guide to MVP</hands-on-title>
>
> The spectra belonging to these "Novel peptides" can be viewed using MVP,this can be achieved by selecting the output from the `mz to sqlite tool` (Generated in the second workflow).
> Here is a step by step guide to obtain the proteogenomic view of the "Novel peptides".
>
>
> 1) Click on the **Visualize in MVP application**, it will open up options for visualization application in the center pane, Select **MVP Application** from the options (or Right click to open in a new window).
>
> ![mz to sqlite](../../images/Visualize.png){:width="20%"}
>
> This will open in the Center Pane. Select the MVP application.
>
>![MVP](../../images/Open_MVP.png){:width="20%"}
>
> This is how it will look.
>
> ![Overview](../../images/Peptide_overview.png){:width="60%"}
>
> 2) Click on **Load from Galaxy** to open the list of peptides you would like to view.
>
> ![load from galaxy](../../images/load_from_Galaxy.png){:width="80%"}
>
>This will open a dropdown list. Select the novel peptides from there.
>
>
> ![dropdown](../../images/DropDown_novel.png){:width="50%"}
>
>
> 3) Select **Novel Peptides** from the right hand side.
>
> ![novel peptides view](../../images/novel_peptides_view.png){:width="70%"}
>
> 4) Select any peptide, For eg: `DGDLENPVLYSGAV`, and then click on **Selected Peptide PSMs**.
>
> ![select_pep](../../images/Filtering_novelpeptides.png){:width="70%"}
>
> 5) If you scroll down, the PSM associated with the peptide will be displayed. By clicking on the PSM, the *Lorikeet*
> values will be shown. The Lorikeet visualization is interactive, i.e the user can change the values or select any
> parameter and click on Update button to view these changes.
>
> ![PSM](../../images/Peptide-Protein_Viewer.png){:width="90%"}
>
> ![lorikeet](../../images/Lorikeet_Viewer.png){:width="70%"}
>
>
> 6) For a Protein centric view, click on **View in Protein** , it will open up all the proteins associate with the
> peptides. For eg: Select the `DGDLENPVLYSGAV` peptide and click on the first protein. The chromosome location
> of the peptide will be displayed.
>
> ![view in protein](../../images/Selecting_peptide.png){:width="70%"}
>
> Once you click on protein it will show the list of proteins the belongs to the peptides.
>
> ![select protein](../../images/opening_IGV.png){:width="60%"}
>
> Once you select the protein that you want to visualize you can click on the protein view.
>
> ![protein view](../../images/ProteinViewer.png){:width="80%"}
>
>
> 7) Clicking on the arrow marks will open up the IGV(js) visualization tool, where-in the genomic localization of the
> peptide will be displayed.
>
>![IGV](../../images/IGV_1.png){:width="90%"}
>
>
> 8) To add tracks to your IGV viewer, click on **Add Track**. This will open up a list of tracks that are compatible
> to view in your IGV viewer. For eg. Select the `Pep_gen_coordinate.bed` file and then click on **Load Track**.
> This will open up the bed will below the nucleotide sequence.
>
> ![tracks align](../../images/Options_for_Track.png){:width="90%"}
>
> 9) By clicking the wheel, you can select the **three frame translate** which will show the three frame translated
> region of your sequence.
>
> ![IGV viewer](../../images/IGV_viewer.png){:width="90%"}
>
>
> 10) The IGV is inbuilt in the MVP viewer and is very interactive, you could also load more tracks such as the aligned
> proBAM file (from HISAT) or the identified probam file (one of the input file).
> MVP has many useful features beyond those covered in this workshop and is under active development.
>
>
>![IGV2](../../images/IGV_2.png){:width="60%"}
>
{: .hands_on}

The next tool in the workflow is the Peptide genomic coordinate tool which takes the "novel peptides" as the input along with the mztosqlite file and the genomic mapping sqlite file (obtained during creation of the database). This tool helps create a bed file with the genomic coordinate information of the peptides based on the sqlite files.




# Obtain Peptide genomic Coordinates

Gets genomic coordinate of peptides based on the information in mzsqlite and genomic mapping sqlite files. This program
loads two sqlite databases (mzsqlite and genomic mapping sqlite files) and calculates the genomic coordinates of the
peptides provided as input. This outputs bed file for peptides.

> <hands-on-title>Peptide genomic Coordinate</hands-on-title>
>
> 1. Run {% tool [Peptide genomic Coordinate](toolshed.g2.bx.psu.edu/repos/galaxyp/peptide_genomic_coordinate/peptide_genomic_coordinate/0.1.1) %} with the following parameters:
>    - *"Input"*: `Peptide list file`, `mzsqlite sqlite DB file`, and `genomic mapping sqlite DB file`
>    - *"Output"*: `Tabular BED file with all the columns`
>
>      ![Output PGC](../../images/Output_PGC.png){:width="50%"}
>
>
> 2. Click **Execute** and inspect the resulting files
>
{: .hands_on}


# Classify Peptides

Given chromosomal locations of peptides in a BED file, PepPointer classifies them as CDS, UTR, exon, intron, or intergene.

> <hands-on-title>Peppointer</hands-on-title>
>
> 1. {% tool [Peppointer](toolshed.g2.bx.psu.edu/repos/galaxyp/pep_pointer/pep_pointer/0.1.3) %} with the following parameters:
>   - {% icon param-select %} *"Choose the source of the GTF file"* - `From History`
>   - {% icon param-file %} *"GTF file with the genome of interest"* - `edited_Mus_Musculus_GRCm38.90_Ensembl_GTF`
>   - {% icon param-file %} *"BED file with chromosomal coordinates of peptides"*: `Bed file from Peptide genomic coordinate tool`
>
> 2. Click **Execute** and inspect the query results file after it turned green.
>
> This tool provides a bed output with the classification of the genomic location of the peptides.The Mus-musculus GTF file will be in your history if you have completed the proteogenomics 1 tutorial.
>
> ![Output PP](../../images/Output_PP.png){:width="50%"}
>
{: .hands_on}

The final tool for this workflow generates a tabular output that summarizes the information after running these workflows. The final summary output consists of the Peptide sequence, the spectra associated with the peptides, the protein accession number, chromosome number, Start and Stop of the genomic coordinate, the annotation, the genomic coordinate entry for viewing in Integrative Genomics Viewer (IGV), MVP or UCSC genome browser and the URL for viewing it on UCSC genome browser. This summary is created with the help of the query tabular tool.


# Final Summary Output

> <hands-on-title>Query Tabular</hands-on-title>
>
>  1. {% tool [Query Tabular](toolshed.g2.bx.psu.edu/repos/iuc/query_tabular/query_tabular/3.0.0) %}
>     - {% icon param-repeat %}  **Insert Database Table**
>       - Section **Table Options**:
>         - *"Specify Name for Table"*: `bed_pep_pointer`
>         - *"Use first line as column names"* : `No`
>         - *"Specify Column Names (comma-separated list)"*:`chrom,start,end,peptide,score,strand,annot`
>         - *"Only load the columns you have named into database"*: `No`
>
>     - {% icon param-repeat %} **Insert Database Table**
>       - Section **Filter Dataset Input**
>         - {% icon param-repeat %} **Insert Filter Tabular Input Lines**
>           - *"Filter by"*:  `skip leading lines`
>           - *"Skip lines"*: `1`
>       - Section **Table Options**
>         - *"Specify Name for Table"*: `psm`
>         - *"Use first line as column names"* : `No`
>         - *"Specify Column Names (comma-separated list)"*: `ID,Proteins,Sequence,AAs_Before,AAs_After,Position,Modified_Sequence,Variable_Modifications,Fixed_Modifications,Spectrum_File,Spectrum_Title,Spectrum_Scan_Number,RT,mz,Measured_Charge,Identification_Charge,Theoretical_Mass,Isotope_Number,Precursor_mz_Error_ppm,Localization_Confidence,Probabilistic_PTM_score,Dscore,Confidence,Validation`
>         - *"Only load the columns you have named into database"*: `No`
>     - *"SQL Query to generate tabular output"*:
>       ```
>       SELECT psm.Sequence AS PeptideSequence,
>       count(psm.Sequence) AS SpectralCount,
>       psm.Proteins AS Proteins,
>       bed_pep_pointer.chrom AS Chromosome,
>       bed_pep_pointer.start AS Start,
>       bed_pep_pointer.end AS End,
>       bed_pep_pointer.strand AS Strand,
>       bed_pep_pointer.annot AS Annotation,
>       bed_pep_pointer.chrom||':'||bed_pep_pointer.start||'-'||bed_pep_pointer.end AS GenomeCoordinate,
>       'https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&position='||bed_pep_pointer.chrom||'%3A'||bed_pep_pointer.start||'-'||bed_pep_pointer.end AS UCSC_Genome_Browser
>       FROM psm
>       INNER JOIN bed_pep_pointer on bed_pep_pointer.peptide = psm.Sequence
>       GROUP BY psm.Sequence
>       ```
>
>     - *"include query result column headers"*: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>
>    ![Final Summary](../../images/final_summary.png){:width="100%"}
>
> The Final summary displays a tabular output containing the list of novel peptides and its corresponding protein. It also provides the users with the chromosomal location of the novel proteoform along with the peptide's start and end position.	The output also features the strand information, gene	annotation and the genomic coordinates in a specific format that could be used on IGV or UCSC browser. It also provides the user with a	UCSC Genome Browser link, which the user can directly copy and paste it on a web browser to learn more about the novel proteoform. Here we are demonstrating the use of proteogenomics workflow on an example trimmed mouse dataset. This study explores the possibilities for downstream biological /functional analysis of peptides corresponding to novel proteoforms.
>
{: .hands_on}

### Conclusion

This completes the proteogenomics workflow analysis. This training workflow uses mouse data. For any other organism the data, tool paramters and the workflow will need to be modified accordingly.This workflow is also available at [usegalaxy.eu](https://usegalaxy.eu/).

This workflow was developed by the Galaxy-P team at the University of Minnesota.
For more information about Galaxy-P or our ongoing work, please visit us at [galaxyp.org](http://galaxyp.org)
