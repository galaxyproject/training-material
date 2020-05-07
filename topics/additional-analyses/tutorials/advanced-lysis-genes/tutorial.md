---
layout: tutorial_hands_on
topic_name: additional-analyses
tutorial_name: advanced-lysis-genes
title: "Advanced Lysis Gene Search"
---
> ### Agenda
>
> 1. Tools for Finding Spanins
> 2. Tools for Finding Holins
> 3. Tools for Finding Endolysins
> 4. Tools for Checking the Proximity of Potential Lysis Genes
> {:toc}
>
{: .agenda}

This tutorial expands on what has been discussed in the [Finding and Annotating Lysis Genes Tutorial](https://cpt.tamu.edu/training-material/topics/additional-analyses/tutorials/finding-lysis-genes/tutorial.html). When lysis genes were not identified by inspection of Functional Workflow outputs, try some of the tools described below.

> ### {% icon tip %} Reminder!
> Bioinformatic tools cannot always be used to identify novel lysis genes. A conservative approach is preferred to calling something a lysis gene when it isn't, but researchers should still do their due diligence before concluding that their phage's lysis genes are not identifiable with the tools/databases in their current state.
{: .tip}


# Tools for Finding Spanins
*Brief para on spanin features, and what functional workflow looks for.*

There are three tools one can use to help find spanins. The [ISP Candidates tool](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt2.spanin.generate-putative-isp) constructs a putative list of potential i-spanin from an input genomic FASTA file.After this tool successfully runs it generates FASTA, gff3, and txt files consisting of potential i-spanins for your phage. An example output is below.


[Image of ISP output here]


A second useful spanin tool is the [OSP candidates tool](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt2.spanin.generate-putative-osp). Similar to the ISP Candidate tool, this tool constructs a putative list of potential o-spanin from an input genomic FASTA file. An example output is below. 


[Image of OSP output here]


Finally, the [Find Spanin tool](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt2.spanin.findSpanin) can be run to narrow down the putative spanins lists to hopefully obtain more accurate candidate i-spanin and o-spanin pairs. To run this tool, use the FASTA output files from the ISP candidates and OSP candidates tools, select the preferred distance between each spanin gene, choose the strand, and click "Execute". Upon successful completion, the tool will output a file for each potential type of candidates (eg. overlap_results.txt) and  a basic summary statistics file (findSpanin_summary.txt).

[Image of Spanin Tool output here]


> ### {% icon tip %} Note:
> Alternatively, you can import the [link worflow that chains all 3 together here]()to run all three tools together. 
{: .tip}


*Describe how gff3 can be piped to Apollo.*


> ### {% icon tip %} Note:
> Some additional tools that may prove helpful include the following
> * [LipoP tool](https://cpt.tamu.edu/galaxy-pub/root?tool_id=geiger.tamu.edu/toolshed/repos/esr/cpt_external_programs/LipoP/1.0.0)
>    > * prediction of lipoproteins and for discriminating between lipoprotein signalpetides, other signal peptides and n-terminal membrane helices in Gram negative bacteria
> * [LipoP to GFF3 Tool](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt.gff3.lipoP_to_gff3)
>    > * adds LipoP results to GFF3
> * [Identify Lipoboxes Tool](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt.fasta.lipory)
>    > * identifies possible LipoBoxes from an input GFF3 and FASTA
{: .tip}


# Tools for Finding Holins
Brief description of holin types.

Link to TMHMM tools (there are two, one makes graphs and the other the gff3 that usually goes to Apollo), and the search file tool (tell them to check the holin/antiholin boxes, search through the blast results and interpro results).


# Tools for Finding Endolysins
Search file tool (tell them to check the endolysin boxes, search through the blast results and interpro results).


# Tools for Checking the Proximity of Potential Lysis Genes
Discuss cassettes, which doesn't apply to some larger phages. Leave space to describe prox to lysis concept. Link to AMH new tools. In reality this can be a Coming Soon section for now, since these tools aren't all finished.
