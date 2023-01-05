---
layout: tutorial_hands_on

title: "NCBI BLAST+ against the MAdLand"
zenodo_link: 'https://zenodo.org/record/4710649'
questions:
- "How to access Galaxy server?"
- "How can we perform Blast analysis on Galaxy?"
- "What is MAdLand DB?"

objectives:
- "Load fasta sequence into Galaxy"
- "Perform NCBI-Blast+ analysis on Galaxy"

time_estimation: "20m"
key_points:
- Blast tool searches a database of sequences for similar sequences to a query sequence.
- MAdLand is a database of fully sequenced plant and algal genomes, with an emphasis on non-seed plants and streptophyte algae that can be use for sequence similarity search.

contributors:
- deeptivarshney

---

# Introduction

<!-- This is a comment. -->

MAdLand is a collection of fully sequenced plant and algal genomes, with a focus on non-seed plants and streptophyte algae. For comparison, it includes genomes from fungi, animals, the SAR group, bacteria, and archaea. It is developed and maintained by the [Rensing lab](http://plantco.de). The species are abbreviated by a 5 letter code, which consists of the first three letters of the genus and the first two of the species name, e.g. CHABR for Chara braunii. We add the gene ID to that and additional shortcuts, like whether it is plastome encoded (pt) or transcriptome-based (tr, in cases when no genome is available yet).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Access the Galaxy server

> <hands-on-title>Log in or register</hands-on-title>
>
> 1. Open your favorite browser (Chrome/Chromium, Safari, or Firefox, but not Internet Explorer/Edge!)
> 2. Browse to [Galaxy](https://usegalaxy.eu) 
> 3. Choose *Login or Register* from the navigation bar at the top of the page
> 4. If you have previously registered an account with this particular instance of Galaxy (user accounts are *not* shared between public servers!), proceed by logging in with your registered *public name*, or email address, and your password.
>
>    If you need to create a new account, click on *Register here* instead.
>

# Upload data on Galaxy

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial and give it a proper name
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>    {% snippet faqs/galaxy/histories_rename.md %}


# Perform NCBI Blast+ on Galaxy 

> After successfully logging in to the Galaxy server, Go to the [NCBI-Blast+](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) tool.  
> Since MAdLandDB is the collection of protein sequences, You can perform [BLASTp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) and [BLASTx](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastx_wrapper/2.10.1+galaxy2) tools.

> <hands-on-title> Similarity search against MAdLand Database </hands-on-title>
>
> 1. {% tool [NCBI BLAST+ blastp](toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastp_wrapper/2.10.1+galaxy2) %} with the following parameters:
>    -  *"Protein query sequence(s)"*:  `Aminoacid input sequence`
>    - *"Subject database/sequences"*:  `Locally installed BLAST database`
>    - *"Protein BLAST database"*: `MadLandDB (Genome zoo) plant and algal genomes with a focus on non-seed plants and streptophyte algae (22 Dec 2022)`
>    - *"Type of BLAST"*:  `blastp - Traditional BLASTP to compare a protein query to a protein database`
>    - *"Set expectation value cutoff"*: `0.001`
>    - *"Output format"*: 
>    - In *"Output Options"*: `Tabular (extended 25 columns)` 

> <img src="../../images/ncbi-blast-against-the-madland/blast-example.png" alt="blast against madland" width="80%">

# Blast output 

>{% icon tool %} The following 12 columns are predicted in the tabular formatted output file of the Blast result (you can select desired output format options):  

<table>
  <tr>
    <th>Column</th>
    <th>NCBI name</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>1</td>
    <td>qseqid</td>
    <td>Query Seq-id (ID of your sequence)</td>
  </tr>
  <tr>
    <td>2</td>
    <td>sseqid</td>
    <td>Subject Seq-id (ID of the database hit)</td>
  </tr>
   <tr>
    <td>3</td>
    <td>pident</td>
    <td>Percentage of identical matches</td>
  </tr>
   <tr>
    <td>4</td>
    <td>length</td>
    <td>Alignment length</td>
  </tr>
   <tr>
    <td>5</td>
    <td>mismatch</td>
    <td>Number of mismatches</td>
  </tr>
   <tr>
    <td>6</td>
    <td>gapopen</td>
    <td>Number of gap openings</td>
  </tr>
   <tr>
    <td>7</td>
    <td>qstart</td>
    <td>Start of alignment in query</td>
  </tr>
   <tr>
    <td>8</td>
    <td>qend</td>
    <td>End of alignment in query</td>
  </tr>
   <tr>
    <td>9</td>
    <td>sstart</td>
    <td>Start of alignment in subject (database hit)</td>
  </tr>
   <tr>
    <td>10</td>
    <td>send</td>
    <td>End of alignment in subject (database hit)</td>
  </tr>
    <tr>
    <td>11</td>
    <td>evalue</td>
    <td>Expectation value (E-value)</td>
  </tr>
    <tr>
    <td>12</td>
    <td>bitscore</td>
    <td>Bit score</td>
  </tr>
</table>

For more details for BLAST analysis and output, we recommand you to follow the [Similarity-searches-blast](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/genome-annotation/tutorial.html#similarity-searches-blast) tutorial.

> <details-title>Further Reading about BLAST Tools in Galaxy</details-title>
>
> Cock et al. (2015): [NCBI BLAST+ integrated into Galaxy](http://biorxiv.org/content/early/2015/05/04/014043.full-text.pdf+html)
>
> Cock et al. (2013): [Galaxy tools and workflows for sequence analysis with applications in molecular plant pathology](https://peerj.com/articles/167/)
{: .details}