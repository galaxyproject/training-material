---
layout: tutorial_hands_on
topic_name: Genome-Annotation
tutorial_name: microbial-annotation
---

# Introduction

Microbial (bacterial) genomes are small and there are usually no introns in the genes. These genomes can be annotated using the tool "Prokka". 

Prokka finds and annotates features (both protein coding regions and RNA genes, i.e. tRNA, rRNA) present on on a sequence. Note, Prokka uses a two-step process for the annotation of protein coding regions: first, protein coding regions on the genome are identified using [Prodigal](http://prodigal.ornl.gov/); second, the *function* of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases. Prokka is a software tool that can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats. More information about Prokka can be found [here](https://github.com/tseemann/prokka).

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. [Data input](#input)
> 2. [Genome Annotation](#annotation)
> 3. [Visualising the annotation](#JBrowse)
> {: .agenda}

# Input data

Prokka requires a genome assembly; e.g. assembled contigs.

Upload the file called "assembled contigs" from the Zenodo link (to do).

> ### :pencil2: Hands-on: Data upload
>
> 1. Step1
> 2. Step2
>
>    > ### :nut_and_bolt: Comments
>    > A comment
>    {: .comment}
>
>    > ### :bulb: Tip: A tip
>    >
>    > * Step1
>    > * Step2
>    {: .tip}
{: .hands_on}

# Run Prokka

> ### :pencil2: Hands-on: Data upload
> - In Galaxy, go to Tools: NGS Analysis: NGS: Annotation: Prokka  
> - Set the following parameters (leave everything else unchanged):
>    - Contigs to annotate: your input assembly  
>    - Locus tag prefix (--locustag): P
>    - Force GenBank/ENA/DDJB compliance (--compliant): *No*
>    - Sequencing Centre ID (--centre): V
>    - Use genus-specific BLAST database *No* 
{: .hands_on}

# Examine the output

> ### :pencil2: Hands-on:

> Once Prokka has finished, examine each of its output files.
{: .hands_on}

- The GFF and GBK files contain all of the information about the features annotated (in different formats.)
- The .txt file contains a summary of the number of features annotated.
- The .faa file contains the protein sequences of the genes annotated.
- The .ffn file contains the nucleotide sequences of the genes annotated.

## View annotated features in JBrowse

Now that we have annotated the draft genome sequence, we would like to view the sequence in the JBrowse genome viewer.

- Go to Statistics and Visualisation: Graph/Display Data: JBrowse

- Under JBrowse-in-Galaxy Action choose *New JBrowse Instance*.

- Under Reference genome to display choose *Use a genome from history*.

- Under Fasta sequences choose Prokka on data XX:fna. This .fna sequence is the fasta nucleotide sequence, and will be the reference against which annotations are displayed.

- For Produce a Standalone Instance select *Yes*.

- For Genetic Code choose *11: The Bacterial, Archaeal and Plant Plastid Code*.

- Click Insert Track Group

- Under Track Category type in *gene annotations*.

- Click Insert Annotation Track

- For Track Type choose *GFF/GFF3/BED/GBK Features*

- For GFF/GFF3/BED Track Data select Prokka on data XX:gff  [Note: not wildtype.gff]


- Under JBrowse Track Type[Advanced] select *Canvas Features*.

- Click on JBrowse Styling Options <Advanced]

- Under JBrowse style.label correct the word "prodcut" to "product".

- Under Track Visibility choose *On for new users*.

Your tool interface should look like this:

![JBrowse interface](images/jbrowse_interface.png)

- Click Execute

- A new file will be created, called JBrowse on data XX and data XX - Complete. Click on the eye icon next to the file name. The JBrowse window will appear in the centre Galaxy panel.

- Under Available Tracks on the left, tick the box for Prokka on data XX:gff.

- Select contig 6 in the drop down box. You can only see one contig displayed at a time.

![JBrowse](images/jbrowse_01.png)

- Use the plus and minus buttons to zoom in and out, and the arrows to move left or right (or click and drag within the window to move left or right).

- Zoom in to see the reference sequence at the top. JBrowse displays the sequence and a 6-frame amino acid translation.

Zoomed in view:

![JBrowse](images/jbrowse4.png)

- Right click on a gene/feature annotation (the bars on the annotation track), then select View Details to see more information.
    - gene name
    - product name
    - you can download the FASTA sequence by clicking on the disk icon.

> ### :nut_and_bolt: Comment
>
> Do you want to learn more about the principles behind mapping? Follow our [training](../../NGS-mapping)
> {: .comment}

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
