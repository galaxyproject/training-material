---
layout: tutorial_hands_on
topic_name: sequence-analysis
tutorial_name: "annotation-with-prokka"
---

# Introduction
{:.no_toc}

In this section we will use a software tool called Prokka to annotate a draft genome sequence. Prokka is a “wrapper”; it collects together several pieces of software (from various authors), and so avoids “re-inventing the wheel”.

Prokka finds and annotates features (both protein coding regions and RNA genes, i.e. tRNA, rRNA) present on on a sequence. Note, Prokka uses a two-step process for the annotation of protein coding regions: first, protein coding regions on the genome are identified using [Prodigal](http://prodigal.ornl.gov/); second, the *function* of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases. Prokka is a software tool that can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats. More information about Prokka can be found [here](https://github.com/tseemann/prokka).

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Import the data

Prokka requires assembled contigs.

> ### {% icon hands_on %} Hands-on: Obtaining our data
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    > ### {% icon tip %} Starting a new history
>    >
>    > * Click the **gear icon** at the top of the history panel
>    > * Select the option **Create New** from the menu
>    {: .tip}
>
> 2. **Import Sample Data.**
>   - Obtain data directly from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1156405.svg)](https://doi.org/10.5281/zenodo.1156405)
>   - Download `contigs.fasta`
>   - Upload the file to your history.
> <br><br>
>
{: .hands_on}

## Annotate the genome

Now we will run the tool called Prokka.

> ### {% icon hands_on %} Hands-on: Annotate genome
>
> 1. **Prokka** {% icon tool %} with the following parameters (leave everything else unchanged)
>    - contigs to annotate: `contigs.fasta`
>    - Locus tag prefix (--locustag): P
>    - Force GenBank/ENA/DDJB compliance (--compliant): No
>    - Sequencing Centre ID (--centre): V
>    - Genus Name: Staphylococcus  
>    - Species Name: aureus  
>    - Use genus-specific BLAST database: No  
>    - Your tool interface should look like this:
>    - ![prokka interface](images/interface.png)
>    - Click Execute
> <br><br>
{: .hands_on}

## Examine the output

Once Prokka has finished, examine each of its output files.

 - The GFF and GBK files contain all of the information about the features annotated (in different formats.)
 - The .txt file contains a summary of the number of features annotated.
 - The .faa file contains the protein sequences of the genes annotated.
 - The .ffn file contains the nucleotide sequences of the genes annotated.
 <br><br>


## View annotated features in JBrowse

Now that we have annotated the draft genome sequence, we would like to view the sequence in the JBrowse genome viewer. First, we have to make a JBrowse file. Then, we can view it within Galaxy.

> ### {% icon hands_on %} Hands-on: Visualize the annotation
>
> 1. Search for **JBrowse** {% icon tool %} and run it with the following parameters
>    - "Reference genome to display" to `Use a genome from history`
>    - "Select the reference genome" to `Prokka on data XX.fna`.
>       
>       This sequence will be the reference against which annotations are displayed
>
>    - "Produce Standalone Instance" to `Yes`
>    - "Genetic Code" to `11: The Bacterial, Archaeal and Plant Plastid Code`
>    - "JBrowse-in-Galaxy Action" to `New JBrowse Instance`
>    - "Track Group"
>    - We will now set up one track - each track is a dataset displayed underneath the reference sequence (which is displayed as nucleotides in FASTA format).
>    - We will choose to display the annotations (the Prokka.gff file).
>
>       - **Track 1 - sequence reads**: Click on `Insert Track Group` and fill it with
>           - "Track Cateogry" to `gene annotations`
>           - Click on `Insert Annotation Track` and fill it with
>               - "Track Type" to `GFF/GFF3/BED/GBK Features`
>               - "GFF/GFF3/BED Track Data" to `Prokka on data XX:gff`
>               - "Track Visibility" to `On for new users`
>               - "JBrowse Track Type [Advanced]" to `Canvas Features`
>               - Click on "JBrowse Styling Options [Advanced]"
>               - "JBrowse style.label" to `product,name,id`
>               - "Track Visibility" to `On for new users`
>               - Click Execute
>
> A new file will be created in your history, this contains the JBrowse interactive visualisation. We will now view its contents and play with it
> 2. Inspect the `JBrowse on data XX and data XX - Complete` file by clicking on the eye icon
>
>    The JBrowse window will appear in the centre Galaxy panel.
>
> 3. Display all the tracks and practice maneuvering around
>    1. Click on the tick boxes on the left to display the tracks
>    2. Select contig 1 in the drop down box. You can only see one contig displayed at a time.
>    1. Zoom out by clicking on the `minus` button to see sequence reads and their coverage (the grey graph)
>    1. Zoom in by clicking on the `plus` button to see annotations.
>    1. JBrowse displays the sequence and a 6-frame amino acid translation.
>    1. Right click on a gene/feature annotation (the bars on the annotation track), then select View Details to see more information.
>      - gene name
>      - product name
>      - you can download the FASTA sequence by clicking on the disk icon.
> ![JBrowse](images/jbrowse6.png)
{: .hands_on}
