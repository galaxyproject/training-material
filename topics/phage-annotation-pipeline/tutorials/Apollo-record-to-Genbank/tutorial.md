---
layout: tutorial_hands_on
topic_name: phage-annotation-pipeline
tutorial_name: Apollo-record-to-Genbank
---
This guide  covers the final steps to prepare an annotated Apollo phage genome record for deposit into Genbank. 

> ### Agenda
>
> 1. Polishing Annotation in Apollo
>    > *  Reviewing all features and naming in Apollo
> 2. Retrieve Apollo Annotation Record into Galaxy
> 3. Fix Special Features not Supported by Apollo
>    > * Locus Tags
>    > * Intron-containing Genes
>    > * Frameshifts
>    > * Merged CDS on the Minus Strand
>    > * Terminal Repeats
>    > * Cos End Sequence
> 4. Final Formatting for Genbank Deposit
>    > * 5-column Table Text File
>    > * DNA Sequence FASTA File
>    > * Source Information
>
{: .agenda}

# Step 1: Polishing Annotation in Apollo
Finalize the genome annotation in Apollo. Check the following before proceeding to the next step.

> * Ensure names conform to the [International Protein Nomenclature Guidelines](https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/), which have been unified between all major organizations that store sequence databases including NCBI and UniProt.
> * As you review the genome, delete unnecessary comments and clean up comments with special characters (such as “) that may break other softwares. 
> * Ensure the frameshift proteins are annotated correctly. See [this guideline](https://cpt.tamu.edu/training-material/topics/phage-annotation-pipeline/tutorials/annotating-tmp-chaperone-frameshifts/tutorial.html) for properly annotated tape measure chaperone frameshift products.
> * Review other special annotation cases such as those covered in [the interrupted genes tutorial](https://cpt.tamu.edu/training-material/topics/additional-analyses/tutorials/finding-interrupted-genes/tutorial.html). 


# Step 2: Retrieve Apollo Annotation Record into Galaxy

> * In Galaxy, run the [Retrieve Data](https://cpt.tamu.edu/galaxy-pub/root?tool_id=export) tool (use this [tool link](https://cpt.tamu.edu/galaxy/root?tool_id=export) if you are CPT internal user).
> * If the genome needs to be re-opened, then do that now following this [tutorial](https://cpt.tamu.edu/training-material/topics/additional-analyses/tutorials/reopening-apollo-with-annotations/tutorial.html). If not, proceed with the gff3 to Genbank step below.
> * Run the [GFF3 to GenBank](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt.gff.gff2gb) conversion tool (use this [tool link](https://cpt.tamu.edu/galaxy/root?tool_id=edu.tamu.cpt.gff.gff2gb) if you are CPT internal user).
> * Run the [GenBank to 5 Column Table](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt.genbank.GBKtoFiveCol) tool (use this [tool link](https://cpt.tamu.edu/galaxy/root?tool_id=edu.tamu.cpt.genbank.GBKtoFiveCol) if you are CPT internal user) to generate the 5 column table text file, which will be verified and edited before depositing to GenBank.


# Step 3: Fix Special Features not Supported by Apollo
In the 5 column table generated from **Step 2**, examine if locus tags assigned to all gene features are correct and in order.  Also examine the genome features that are not reliably supported by Apollo, and fix manually in the 5-column table if necessary.  This includes making sure frameshift proteins are fused properly at the correct position, special features such as terminal repeats or cos site sequences are added, and special qualifiers required by NCBI are added.  The [NCBI rules](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/) are the ultimate guidelines on the format that needs to be used in the polished version to deposit. It includes the current standard for features like introns, inteins, etc. Read it to make sure you are still up to date.

> ### {% icon comment %} Note:
>  An alternative to fixing the 5 column table is to use a third party Genbank file editing software, such as [Artemis](http://sanger-pathogens.github.io/Artemis/) to edit and fix problematic features.  Open your Genbank file in Artemis to proceed with the editing and generating a 5 column table (sequin table format) required for Genbank submission.  This guide does not cover how to edit a Genbank file in Artemis.  
{: .tip}

Typically we verify the following:
>    > #### I. Locus Tags
>    > 
>    >Verify if the locus tags are assigned properly for the tricky features per [NCBI locus tag rules](https://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf).  If not correct, you need to run the genbank file through [Renumbering](https://cpt.tamu.edu/galaxy-pub/root?tool_id=edu.tamu.cpt.genbank.RelabelTags) tool (use this [tool link](https://cpt.tamu.edu/galaxy/root?tool_id=edu.tamu.cpt.genbank.RelabelTags) if you are CPT internal user) and re-generatet 5 column table to verify again.
>    >
>    > Basically, locus_tags should be assigned to all protein coding and non-coding genes such as structural RNAs.  Locus_tag should appear on gene, mRNA, CDS, 5'UTR, 3'UTR, intron, exon, tRNA, rRNA, misc_RNA, etc within a genome project submission. Repeat_regions do not have locus_tag qualifiers. The same locus_tag should be used for all components of a single gene. 
>    >
>    > For example, all of the exons, CDS, mRNA and gene features for a particular gene would have the same locus_tag. There should only be one locus_tag associated with one gene.  For tRNA, the locus tag assigned usually follows the consecutive order from CDS, rather than assigning a specific prefix.  There should be one gene and one product assigned to each tRNA under the same locus tags.  See below for an example 5 column table.
>    >

> ### {% icon comment %} Note:
>  In theory we should check each feature for its locus tag. In practice, check all tRNAs, a couple genes on the plus strand and a couple genes on the minus strand, then spot check through the rest of the way down your genome.
{: .tip}

>    >
>    > ![](../../images/Apollo-record-to-Genbank-screenshots/1-feature_table.PNG)
>    >
>    > #### II. Intron-containing Genes
>    >
>    > Verify the base locations of the intron-containing genes are correct.  Gene feature span should be a single span covering all exons and introns.  The actual CDS feature should be annotated with sets of nucleotide spans showing how the exons are joined to create the correct product. Informative notes such as “Intron splice site predicted by Blast homology to XXXXX” or “intron contains VSR homing endonuclease” can be added.  When exon boundaries can not be identified, only the gene span covering all exons and introns is reported. See [the interrupted genes tutorial](https://cpt.tamu.edu/training-material/topics/additional-analyses/tutorials/finding-interrupted-genes/tutorial.html) for more info. 
>    >
>    > #### III. Frameshifts
>    >
>    > > > * Make sure the base coordinates of the two joint CDS, such as “join(55078..55479,55478..56029)”, are correct. 
>    > 
>    > > > * There should be a note “exception=ribosomal slippage” added to the frameshift product.
>    >
>    > > > * Make sure the previously added qualifier “Frameshift protein a” and any other qualifier added in Apollo is stripped.
>    >
>    > > > * Each protein (non-frameshift version and frameshift version) should have their own RBS (which will have the exact same coordinates) and their own locus tags.  
>    >
>    > See below for an example with properly formatted frameshifts.
>    > ![](../../images/Apollo-record-to-Genbank-screenshots/2-formatted-frameshifts.PNG)
>    >
>    > #### IV. Merged CDS on Minus Strand
>    >
>    > For cases with joint exons or joint frameshift proteins located on minus strand, make sure the order of the two CDS in the 5 column table is correct in order to avoid an error message after NCBI reviews the submission.  Here is an example for the proper order of two merged CDS in a 5 column table (*note the fragment from right side of the genome is listed first in the table*):           
>    >
>    > ![](../../images/Apollo-record-to-Genbank-screenshots/3-two-merged-cds.PNG)
>    >
>    > #### V. Terminal Repeats
>    >
>    > If applicable, add "repeat_region" for terminal repeat with the defined coordinates. The CPT convention is to put the TR at the 5’ end, not repeating the sequence at the 3’ end of the genome.  Add notes such as "direct terminal repeat predicted by PhageTerm", and "right end of genome sequence not duplicated in this record" to indicate that the bases are not repeated in the genome deposited, but that that sequence (and its genes) are present in the viral genome.  See below for an example.
>    > ![](../../images/Apollo-record-to-Genbank-screenshots/4-repeat-region.PNG)
>    >
>    > #### VI. Cos End Sequence
>    >
>    > If applicable, add "misc_feature" for cos end sequence with the defined coordinates to indicate the cos end sequence . See below for an example.
>    > ![](../../images/Apollo-record-to-Genbank-screenshots/5-misc-feature.PNG)
>    >

# Step 4: Final Formatting for GenBank Deposit
For Genbank submission through [BankIt](https://www.ncbi.nlm.nih.gov/WebSub/?tool=genbank), you need to provide a 5 column table text file, a DNA sequence FASTA file, and fill in genome source information as required. 

>    > #### I. 5-column Table Text File
>    >
>    > > > * Check the feature table header, it should be >Feature Phagename (such as >Feature Minorna).  All the locus tags in the table can be “CPT_phagename_00X” (such as CPT_ Minorna_00X).
>    > >
>    > > > * Change phage name if needed.  In the 5 column table file, execute **Find and Replace** in a text editor.  Be careful if editing in Excel or Word because it can insert characters upon saving that **break** the NCBI uploader.
>    >
>    > #### II. DNA Sequence FASTA File
>    >
>    > > > * Check the header, it should be >Phagename (such as >Minorna).  **This header has to match the name after “Feature” in 5 column table header.**
>    >
>    > #### III. Source Information
>    >
>    > Gather all source information (source info example is shown below), which are needed when you follow [BankIt](https://www.ncbi.nlm.nih.gov/WebSub/?tool=genbank) to deposit the 5 column table file and the DNA FASTA file to Genbank. 
>    >
>    > ![](../../images/Apollo-record-to-Genbank-screenshots/6-source-info.PNG)
>    >






