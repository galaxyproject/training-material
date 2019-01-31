---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: annotation-in-apollo
---

> ### Agenda
>
> * Prerequisites
> * Making an Annotation
> * Making the Best Prediction
>    > * Gene Calls
>    > * BLAST
>    >    > 1. NT (Nucleotide) database
>    >    > 2. NR (non-redundant) protein database
>    > * Phage Analyses
>    >    > 1. Candidate ISPs/OSPs
>    >    > 2. Possible Intron Locations
>    >    > 3. Possible Frameshifts
>    > Sequence Analyses
>    >    > 1. InterProScan
>    >    > 2. TMHMM (Transmembrane hidden Markov model)
>    >    > 3. Terminators
>    >    > 4. tRNA and tmRNA
>
{: .agenda}

# Prerequisites

Before beginning annotation within Galaxy ([CPT Public Galaxy](https://cpt.tamu.edu/galaxy-pub), [CPT TAMU Galaxy](https://cpt.tamu.edu/galaxy)), it is necessary that there is a genome loaded into Apollo. Additionally, the [structural annotation]({{site.baseurl }}//topics/genome-annotation/tutorials/structural-annotation-workflow/tutorial.html) *must* be complete, and the [functional annotation workflow]({{site.baseurl }}//topics/genome-annotation/tutorials/functional-annotation-workflow/tutorial.html) has **already been run**. The functional annotation workflow opens up the necessary evidence tracks for annotation.

# Making an Annotation

Generally, the annotation process is a synthesis between the understanding of phage genomics and the available evidence tracks. The [Center for Phage Technology](https://cpt.tamu.edu) encourages **phage** annotators on the CPT’s Apollo instance to follow some specific conventions (Field -> *Recommended Input*):

> * Name -> *Gene name* (Could be something like **terminase small subunit** or **hypothetical protein**.) Follow the universal naming conventions at [NCBI](https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/) and [UniProt](https://www.uniprot.org/docs/International_Protein_Nomenclature_Guidelines.pdf). If you are basing your prediction off a BLAST hit (Canonical phage database, SwissProt, TrEMBL, or nr), the hit name *may* be the one you should use as well.
> * Symbol -> *Do not use.*
> * Description -> *Do not use.*
> * DBXRefs -> *Only use if the annotator is experienced; please ensure formatting is correct.*
> * Attributes -> *Do not use,* except in the case of frame shifted proteins.
> * PubMed IDs -> *Do not use.*
> * Gene Ontology IDs -> *Do not use.*
> * Comments -> *Apply any free-text comments here.* (Could be something like **the e-value(s)** between the annotated gene and homologs or notes to one’s self.)

> ### {% icon tip %} Note that…
> Calling genes is [covered in another tutorial.]({{ site.baseurl }}//topics/genome-annotation/tutorials/structural-annotation-workflow/tutorial.html)
{: .tip}

To annotate a gene that has been called, right click on the gene in the pale yellow  User-Created Annotations track, and select “Edit information (alt-click).”

![](../../images/annotation-in-apollo-screenshots/0_right_click.png)

A screen with various fillable fields will appear where information about the gene can be manually entered. These are free text fields.

![](../../images/annotation-in-apollo-screenshots/1_gene_information.png)

Reference the list above to see how the [CPT](https://cpt.tamu.edu) would prefer to have genes annotated.

> ### {% icon comment %} Naming Guidelines
> It is imperative to follow suit with the [UniProt](https://www.uniprot.org/docs/International_Protein_Nomenclature_Guidelines.pdf) and [NCBI](https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/) international naming conventions. It allows for standardization and consistency in naming proteins, subsequently aiding data retrieval and improving communication. Follow the convention for capitalization and hypothetical protein naming.
{: .comment}

> ### {% icon tip %} Ensuring Changes in Gene Information are Saved
> There are occasional small bumps on the road when annotating in Apollo, many of which are encountered when editing information for a gene. It helps to be aware of how to avoid them, and where to fix issues when they arise. For example, the information being shown is that of the gene highlighted in red; this is a gene in the phage P22 genome. The name has been changed from “gene name” (as is seen in the User-Created Annotations track in the background) to “gtrB” in the Information Editor window. There has been no other action outside of typing in “gtrB”.
>
> ![](../../images/annotation-in-apollo-screenshots/6_gene_name_before.png)
>
> Clicking anywhere outside of the “Name” field (or the most recently adjusted field) OR hitting tab on the keyboard should automatically save the change. If it has successfully been saved, the change will be immediately noticeable in the User-Created Annotations track, even without closing the popup window.
>
> ![](../../images/annotation-in-apollo-screenshots/7_gene_name_saved.png)
>
> If the changes are not being saved, refresh the page and try again. If this continues to happen, try opening Galaxy ([CPT Public Galaxy](https://cpt.tamu.edu/galaxy-pub), [CPT TAMU Galaxy](https://cpt.tamu.edu/galaxy)) in an **incognito window**.
{: .tip}


# Making the Best Prediction

The [Center for Phage Technology](https://cpt.tamu.edu/) integrates many data sources when making predictions about gene function. Please contact the [CPT](https://cpt.tamu.edu) IT staff (cpt@tamu.edu) to suggest adding a data source not currently available in Apollo. The [CPT](https://cpt.tamu.edu) staff will assess your recommendation and may add it to the Phage Annotation Pipeline (PAP).

> ### {% icon tip %} A word on genome annotation
> While we use all the best bioinformatic tools available to complete these analyses, judgments based off sequence alone are still predictions. Hypothetical. Every prediction is subject to being wrong and getting corrected with new information. This is the nature of genome annotation, and science in general. **Due diligence and thorough work is expected, but it is a *waste of time to agonize over any single gene prediction*.** To be useful, these predictions need to specific and accurate. Where necessary, it is better to trade specificity for accuracy (be more general when unclear of the specific function).
{: .tip}

Many of the protein prediction tracks yield multiple homologs for the same gene. To learn more about each homolog, hover over the homolog displayed on the track, or right-click on it and select “View details.”

![](../../images/annotation-in-apollo-screenshots/12_track_right_click.png)

The details presented will vary between evidence types (e.g.BLAST vs. InterProScan). It is often necessary to further investigate the evidence for that entry on the orginating database website. 

### Gene Calls

The [CPT’s](https://cpt.tamu.edu) PAP integrates gene calls from numerous sources, specifically *MetaGeneAnnotator* and *Glimmer3*. These gene callers are generally very accurate and should always be the preferred source of gene features where they are available. However, should those fail to find a gene, *SixPack* is a backup gene caller.

![](../../images/annotation-in-apollo-screenshots/2_gene_calls.png)

> ### {% icon tip %} Note that…
> Gene calling is part of structurally annotating a genome in Apollo. For more information on structural annotation and gene calling, please look at [this tutorial]({{ site.baseurl }}//topics/genome-annotation/tutorials/structural-annotation-workflow/tutorial.html).
{: .tip}

### BLAST

BLAST, an acronym standing for Basic Local Alignment Search Tool, is a widely used tool suite available from the National Center for Biotechnology Information website, also known as [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi). Several variations of BLAST exist to search different kinds of sequences. BLAST breaks down the input sequence (called the query) into *k*-mers and compares these to sequences in the database; default *k*-mers are currently 6 for BLASTp (protein database) and 28 for BLASTn (nucleotide sequence database). Matching *k*-mers are then extended to the left and right until the score drops below a threshold, *T*. Sometimes the extended coverage encompasses the entire sequence, and other times the alignment is broken up into homologous segments that are separated by low similarity regions. Alignments with high similarity are referred to as a *high-scoring pair* (HSP); this may be the entire sequence pair, or only part of the query.

Although BLAST is accessible through the [NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi), the databases detailed below are those that are available to an annotator in Apollo. In Apollo, right-clicking on a hit in a BLAST evidence track and selecting ‘View Details’ will pop up a summary of the details about the protein.

![](../../images/annotation-in-apollo-screenshots/13_blast_homolog_details.png)

> * Score = *E-value*; the lower the score, the better the alignment with the query. This is also reflected in the color highlight intensity of the feature displayed in the evidence track.
> * Description = may be informative, but this depends on the quality of the annotation in the protein record. In the above example, the accession number is underlined; this same information can be found in the Description portion. Searching the accession number in the [NCBI protein database](https://www.ncbi.nlm.nih.gov/protein/) will yield more information about the protein, including the paper in which this protein was originally reported. If you are predicting your gene function based off a strong BLAST hit, you *may* want to name it the same way, if it conforms to the proper [UniProt](https://www.uniprot.org/docs/International_Protein_Nomenclature_Guidelines.pdf) and [NCBI](https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/) international naming conventions.

> ### {% icon tip %} Note that…
> In practice, an E-value of **less than** 1e-3 or 1e-5 are considered relevant, **if that hit covers most or all of the protein!** This means that the smaller an E-value, the more confident we can about the alignment between those sequences.
{: .tip}

##### 1. NT (Nucleotide) database

Megablast is run against a copy of NCBI’s NT database. Hovering over a hit segment will show where in the target genome the region aligns.

![](../../images/annotation-in-apollo-screenshots/3_blast_nt.png)

Where you have hits in this track, your phage is similar at the nucleotide level to another sequence in the database. Pay attention to this track as it will give you hints of any mosaicism that may be present.

##### 2. NR (non-redundant) protein database

BLASTp is run against three databases (in the 2018 PAP iteration):

> * [CPT’s](https://cpt.tamu.edu) Canonical Phage database, a select collection of high-quality and well-studied representative phage proteomes
> * SwissProt (curated from [UniProt](https://www.uniprot.org/))
> * TrEMBL (from [UniProt](https://www.uniprot.org/))
> * nr (from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/))

Apollo represents similarity score between the query/novel protein and the homologs by the saturation of the color of the evidence tracks. The more saturated the color of the track, the greater the similarity between the two proteins. An example from TrEMBL is shown below.

![](../../images/annotation-in-apollo-screenshots/4_blast_trembl.png)

Hits in these databases can be strong clues as to the function of your phage protein. Investigate them closely. In order of priority, consider SwissProt>TrEMBL/nr because the SwissProt entries have been manually curated by a human being. Usually by going directly to the [UniProt](https://www.uniprot.org/) or [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/) websites and searching for the accession number displayed in Apollo, you can find more detailed information about the function, location, source, publications, possible structures, *etc*... 

If you are predicting your gene function based off a strong nr hit, you *may* want to name it the same way, if it conforms to the proper [UniProt](https://www.uniprot.org/docs/International_Protein_Nomenclature_Guidelines.pdf) and [NCBI](https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/) international naming conventions.

### Phage Analyses

The [CPT](https://cpt.tamu.edu) has developed a number of phage analysis tools specific to *phage* annotation. These are supplementary bits of information to consider in your analysis, but they must be looked at critically. In an effort to cover every reasonable option, there is a high rate of noise in the tools, unfortunately giving high false positive rates.

![](../../images/annotation-in-apollo-screenshots/5_phage_analyses.png)

##### 1. Candidate ISPs/OSPs

Phage lysis genes are [notoriously poorly annotated](https://www.ncbi.nlm.nih.gov/pubmed/30219026). Often they are missed or completely misattributed. To combat this problem for phage spanin proteins, lysis proteins specific to disrupting the outer membrane of gram-negative bacterial hosts, the [CPT](https://cpt.tamu.edu) utilizes the candidate ISP (i-spanin) and OSP (o-spanin) tool output.

> ### {% icon tip %} Note that…
> These tracks will feature a *huge* number of false positives. Be sure that the data occurs somewhere around the phage’s lysis cluster (where applicable). Additionally, know what to look for in a lipobox in these potential spanin genes.
{: .tip}

The ISP track naïvely searches the genome for every possible CDS, and then analyzes them with TMHMM. This happens even in the case of a mis-called or entirely missed i-spanin. The OSP track searches through every possible CDS which contains a lipobox as defined by the [CPT](https://cpt.tamu.edu). A lipobox (L/V)X(G/A/S)-C is required at the N-terminus of a u-spanin or OSP. *Both* of these datasets are filtered for proximity. Co-incidence of a possible ISP gene and a possible OSP gene is a good sign, but the genomic context information will need to be taken into account to complete the functionality inference.

##### 2. Possible Intron Locations

This track analyzes BLASTp against NR data for locations where two or more called, disjointed CDSs, match separate locations on the same target protein. Below is an example alignment from phage K.

![](../../images/annotation-in-apollo-screenshots/8_intron_phage_k.png)

Both 195a and 195b align to distinct regions of the same protein, based on BLAST data. It can be theorized that these are actually *one* protein with *one* intron and *two* exons; however, **this evidence should not be taken as 100% correct**. Similar results may happen for other reasons, such as separation of domains from a single protein due to evolution, sequencing errors, and a myriad of other possibilities.

##### 3. Possible Frameshifts

Like the Possible Intron Locations track, the Possible Frameshifts track is very optimistic in what it considers a possible frameshift; it is searching for a frameshift that may indicate phage tape measure protein chaperones. It searches for an XXXYYYZ nucleotide pattern (allowing for some wobble) wherein a frameshift would **not** change both codons. This is based on evidence shown in [this research paper.](https://www.sciencedirect.com/science/article/pii/S1097276504005398?via%3Dihub) Further information on annotating the tape measure chaperone proteins can be found in [this tutorial]({{ site.baseurl }}//topics/genome-annotation/tutorials/annotating-tmp-chaperone-frameshifts/tutorial.html).

### Sequence Analyses

Additional analyses run in the PAP are listed in the Annotations Track on the left under both the structural and functional sections. 

![](../../images/annotation-in-apollo-screenshots/9_structural_annotation_tracks.png)

![](../../images/annotation-in-apollo-screenshots/10_functional_annotation_tracks.png)

#### 1. InterProScan

InterProScan is an extremely useful domain finder. It is hosted by [EMBL-EBI](https://www.ebi.ac.uk/) (European Molecular Biology Laboratory - European Bioinformatics Institute) and integrated into [UniProt](https://www.uniprot.org/) (a freely accessible database of protein sequence and functional information, of which many entries are derived from genome sequencing projects). InterProScan searches a protein sequence against the member databases and detects similarity to conserved domains. It integrates 14 other conserved domain databases to assign a single [InterPro](https://www.ebi.ac.uk/interpro/) ID to related domains. Hits from InterProScan predict protein function based on conserved domains (**beware of domain swapping!**)

> ### {% icon details %} Domain Swapping 
> * It is not uncommon for proteins to contain more than one domain that function independent of each other. These domains can be found alone. Sometimes you may see a BLAST hit to a protein that has two domains, but the name of the protein only reflects one domain's function. **When you see a hit, make sure that the domain the name is based off actually aligns to your query.**
{: .details}

These InterPro conserved domain hits will be more sparse than BLAST hits and often only align to a portion of the protein. Right clicking on a feature from the InterProScan evidence track and selecting ‘View Details’ will show a summary of that feature.

![](../../images/annotation-in-apollo-screenshots/14_interproscan_homolog_details.png)

If the domain hit is part of InterPro, there will be a “Dbxref” entry. Searching either the name or Dbxref fields in Google (or use the IPR number at [InterPro's website](http://www.ebi.ac.uk/interpro/)) will often find the domain entry in either the member database or InterPro. The Score for these conserved domain searches is a probability or a quality score, and its calculation varies between databases. This tool, unlike all other tools found in Apollo, has been pre-calibrated to show only results the authors think are significant, so the annotator does not need to evaluate this value.

##### 2. TMHMM (Transmembrane hidden Markov model)

The [TMHMM](http://www.cbs.dtu.dk/services/TMHMM/) tool is used to predict transmembrane helices in proeteins. Here, TMHMM is run over the genome to pick out genes that contain likely transmembrane domains (TMDs). TMHMM data is used in a number of other tracks and analyses as well.

There are three variations of the TMHMM track - *Inside Probability, Membrane Probability,* and *Outside Probability.* The Membrane Probability track yields the actual scores from the TMHMM track plotted to the genome. The Inside Probability track shows the probability of a region of a protein being within the cytoplasm. The Outside Probability shows the probability of a region of a protein being within the periplasm or extracellular.

![](../../images/annotation-in-apollo-screenshots/11_tmhmm_tracks.png)

This tool provides useful evidence when trying to determine if any given protein has the characteristics expected for a potential homolog. Use these predictions only when looking for transmembrane phage proteins, like the holins.

##### 3. Terminators

Terminators are produced from [TransTermHP.](http://transterm.ccb.jhu.edu/) TransTermHP finds rho-independent transcription terminators in bacterial genomes. Each terminator found by the program is assigned a confidence value that estimates its probability of being a true terminator.

> ### {% icon tip %} Be conservative…
> This track will feature a large number of false positives. Only call terminators that have: a score greater than 90, a stem of at least a 5 bp without mismatches, and a polyT of at least 4 in length. Additionally the terminator should be within a logical context, *e.g.* in the correct orientation downstream of a gene.
>
{: .tip}

> ### {% icon tip %} Also note that…
> This track can be found underneath the “Sequence Analysis” section of the Structural Annotation portion all of the tracks.
>
> ![](../../images/annotation-in-apollo-screenshots/9_structural_annotation_tracks.png)
>
{: .tip}

##### 4. tRNA and tmRNA

[ARAGORN](https://www.ncbi.nlm.nih.gov/pubmed/14704338) predicts highly conserved secondary structure found in tRNAs. An annotator should feel confident annotating tRNA and tmRNAs in Apollo using this track as the tool provides high-quality annotation results. Recall that tRNAs are **not** likely to be embedded within genes.

> ### {% icon tip %} Note that…
> This track can be found underneath the “Sequence Analysis” section of the Structural Annotation portion all of the tracks.
>
> ![](../../images/annotation-in-apollo-screenshots/9_structural_annotation_tracks.png)
>
{: .tip}
