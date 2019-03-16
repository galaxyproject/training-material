---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: genome-close-reopen
---

# Genome Closure from a Sequencing Run & Re-opening

> ### Agenda
>
> * Genome Closure
>    > * Confirmation PCR
>    > * Closure PCR
>    >    > * Background
>    > * Trimming
>    > * Renaming
>
> * Contig Re-Opening
>    > * Deciding Where to Re-Open
>    >    > * BLASTn Analysis
>    >    > * PhageTerm Analysis
>    > * Re-Opening in Galaxy
> {:toc}
>
{: .agenda}

# Genome Closure

# Contig Re-Opening

In the final deposited phage genome, base 1 should be in a logical place. The convention is to follow the accepted standard in the field. We try to make it mostly syntenic with genomes already in the database. If your phage is similar to another genome in the database, best-practice will open your genome in the same place. It should not be in the middle of a feature (especially genes). Sometimes, there is no precedent, or the precedent might not make sense in light of published data, in which case, a different course can be taken.

For your phage genome, re-opening will likely need to take place. Many researchers will save this for the last step, or penultimate step, after gene calling and similarity analyses have been determined through the annotation process. In some cases, such as the undergraduate teaching setting, it may be desirable to re-open the genome prior to annotation. Below are described the considerations that might be taken into account when:

1) deciding where to re-open an unannotated genome, and
2) understanding the mechanics of re-opening a genome in Galaxy.

## Deciding Where to Re-Open

After contigs are closed, the closed contigs should be re-opened properly according to their genome types. Without annotating and analyzing the entire genome, the genome type of a given contig could be predicted by running:

> * BLASTn against the NCBI nr database
> * PhageTerm, a software package that uses raw reads and its genomic reference sequence to predict the termini positions. *It is worth noting that such prediction is not always successful or accurate*.

### BLASTn Analysis

Running BLASTn using your genome sequence against the nr database serves two purposes: to get a quick idea on the **orientation** or your genome sequence, and to get an idea on the  **phage type** by looking at the related phages in GenBank. For immediate results, use the public NCBI web-based BLAST page. To save results for the future reference, do the Galaxy procedure.

**_Procuedure_** for running [BLASTn at the NCBI public site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=blastp), and on Galaxy with the [NCBI BLAST + blastn tool](https://cpt.tamu.edu/galaxy/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastn_wrapper/0.1.01). If this was already done in the [ASSEMBLY TUTORIAL LINK], then this section can be skipped.

**_NCBI BLASTn_**
> 1. Copy the entire closed contig sequence into the query sequence box. Usually default parameters for megablast against nr is sufficient.
> 2. Run BLAST
> 3. Sort through hits and their alignments for the below info.

**_GALAXY BLASTn_**
> 1. In the appropriate history with the closed contig fast file, run [NCBI BLAST + blastn tool](https://cpt.tamu.edu/galaxy/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/ncbi_blast_plus/ncbi_blastn_wrapper/0.1.01) with the latest nt database.
> 2. Copy the rest of the directions from the Assembly tutorial.

**Orientation:**
After running BLASTn, check the alignments to the top closest hits. If the top genomes are correctly oriented (most genes are codes on the plus strand), and your genome shows the same orientation, you do not need to reverse complement your sequence. Otherwise, you need to reverse complement your FATSA sequence before re-opening it.
> * For DTR: when PhageTerm gave us the boundaries and we don't have to do RC, then open at the first base in the boundaries (it is an inclusive number).
> * If you have to RC: 1) do that, then re-run PhageTerm (or, do math) and then check the boundary number and re-open OR 2) re-open one base after the left end of the boundary and then RC (re-run PhageTerm to double check).
>    > 1. Check the first 2-3 top BLASTn hits for orientation of all our alignments with it.
>    > 2. If it is ambiguous, or a poorly annotated genome that may have been deposited incorrectly, in that case we will *not* have a convention to follow.
>    > 3. <!-- ADD IN PROCEDURE HERE OR SOMEWHERE IN PROTOCOL. -->
**Phage type:**
Check the top BLAST hits list. Record the accession #, identity, and coverage information of the top BLAST hit for easy reference. Keep a sharp eye out for phage names that will have papers published about them, the canonical phages. Having a genome like a known phage should inform both your re-opening and entire annotation process. If a well-known/studied phage is not among the BLAST hit list, it may help to check the closest type phage in the NCBI taxonomy or on ICTV. **Use this strategy with caution.**

If your genome's phage type can be determined, re-open your genome to make it syntenic to the canonical phage in that phage type. See below for the different phage type scenarios.

### PhageTerm Analysis

[PhageTerm](https://www.nature.com/articles/s41598-017-07910-5) [PMID:28811656](https://www.ncbi.nlm.nih.gov/pubmed/?term=28811656) predicts termini and packaging mechanisms using the raw reads of a phage sequenced with technologies that rely on random fragmentation and its genomic reference sequence. While not fully verified, the tool provides a good guide for genome with well-described end types. Sometimes this prediction is informative when closing a genome (see assembly protocol); it can also be useful for deciding where to re-open a genomic sequence. After BLASTn, run PhageTerm in Galaxy as detailed below.

<!-- COPY IN PROCEDURE FROM THE ASSEMBLY PROTOCOL -->

## Re-Opening genomes in Galaxy After BLASTn and PhageTerm Analysis

After integrating the information from BLASTn and the PhageTerm analysis, the following re-opening guidelines can be followed to re-open your genome with the [Galaxy re-opening tool](https://cpt.tamu.edu/galaxy/root?tool_id=edu.tamu.cpt.fasta.reopen).
> * Run the [Re-open FASTA sequence]() tool