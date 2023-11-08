---
layout: tutorial_hands_on

title: Analyse unaligned ncRNAs
zenodo_link: 'https://zenodo.org/record/3482616'
questions:
- Which biological questions are addressed by the tutorial?
  - Is my set of RNAs of interest homologous?
  - Are the RNAs in the set coding for a protein?
- Which bioinformatics techniques are important to know for this type of data?
  - Alignment
  - RNA secondary structure folding
objectives:
- Knowing for which cases the "Analyze non coding RNAs" workflow can be used
- Being able to apply the analyze non coding RNAs workflow
- Learn about the tools used in the workflow
time_estimation: 1H
key_points:
- Homology
- RNA
- Secondary structure
- Coding potential
contributions:
  authorship:
    - eggzilla
    - bebatut
  editing:
    - pavanvidem
---

# Introduction
{:.no_toc}

Non coding RNAs (ncRNA) that are homologous, meaning that they share a common ancestor RNA, can be identified
by sharing a common sequence or structure. Since the biological function depends on the
spatial structure of a molecule, it is often strongly conserved for ncRNAs.

This workflow investigates a set of sequences for their sequence and structure conservation and for
their alignability. We will inspect a set of dual function RNAs and a set of epidermal growth factor sequences with the workflow and compare the different outcomes.

The first step in any analysis of unaligned ncRNAs is to run a multi alignment tool on the sequences of the ncRNAs. Using the multi alignment, we can extract some representative sequences, analyze them for protein coding potential, compute consensus structure and predict stable secondary structure.

For the multi-alignment, two types of tools are used: one generic multi-aligner and one RNA-secondary structure aware.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

> <hands-on-title> Data upload </hands-on-title>
>
> 1. Create a new history for this tutorial
>
> 2. Import the files from [Zenodo](https://zenodo.org/record/3482616) or from the shared data library
>
>    We will retrieve two files one file containing homolog sequences from the
>    SgrS family, which has been discovered to be a dual-function RNA.
>    This allows us to show the ability of the workflow to test for protein
>    coding potential, as well as for predicting functional non-coding RNAs.
>
>    Search for Upload File in the Tool search field and select it from the search results.
>    The popup for file upload that opens offers the Paste/Fetch data option.
>    Select it and paste following URLs into the URL field and press Start:
>    https://zenodo.org/record/3482616/files/GCF_000005845.2_ASM584v2_genomic.fna?download=1
>    https://zenodo.org/record/3482616/files/RF00534.fasta?download=1
>
> 3. Rename the datasets
>
>    The two datasets uploaded from Zenodo will now be renamed with less verbose names.
>    The history on the left shows two items which symbolize the uploaded datasets.
>    We rename the first from
>    https://zenodo.org/record/3482616/files/GCF_000005845.2_ASM584v2_genomic.fna?download=1
>    to  genome.fa
>    and we also rename
>    https://zenodo.org/record/3482616/files/RF00534.fasta?download=1
>    to sgrs.fa
>
> 4. Check that the datatype
>
>    Both uploaded and renamed files should automatically be detected as fasta format files.
>
{: .hands_on}

# Analysis using generic multi alignment tool

## Align sequences and extract representative sequences

> <hands-on-title> Run multi-alignment tool and extract representative sequences from the alignment </hands-on-title>
>
> 1. **MAFFT** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequences to align"*: `output` (Input dataset)
>    - *"MAFFT flavour"*: `fftns`
>    - *"Output format"*: `ClustalW`
>
> 2. **Select Sequences** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input clustal alignment"*: `outputAlignment` (output of **MAFFT** {% icon tool %})
>
{: .hands_on}

MAFFT is a general purpose alignment tool that outputs a multiple seqence alignment. We will use this
alignment to test the sequences for coding potential and to test if there is a conserved secondary
structure, without using secondary structure information to guide the alignment.

> <question-title></question-title>
>
> 1. What additional characters can you observe in the multiple sequence alignment constrast to the individual sequences?
>
> > <solution-title></solution-title>
> >
> > 1. We can observe gap characters that are representing insertions and deletions
> >    in case the individual sequences are of different length.
> {: .solution}
>
{: .question}

We can now analyze more in depth the representative sequences of the alignment: analyze for potential protein coding regions, compute the consensus structure and predict the stable RNA secondary structure.

## Analysis for potential protein coding regions

We would like first to know if there is any potential protein coding regions left in the alignment. We use **RNAcode** {% icon tool %} which predicts protein coding regions in an alignment of homologous nucleotide sequences. The prediction is based on evolutionary signatures typical for protein genes, i.e. the presence of synonyomous/conservative nucleotide mutations, conservation of the reading frame and absence of stop codons.

**RNAcode** {% icon tool %} does not rely on any species specific sequence characteristics whatsoever and does not use any machine learning techniques.

> <hands-on-title> Identify potential protein coding regions </hands-on-title>
>
> 1. **RNAcode** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Multiple Alignment"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Scoring parameters"*: `Default`
>    - *"Create colored plots in EPS format"*: `Create Plots`
>
{: .hands_on}

**RNAcode** {% icon tool %} reports local regions of unusual high coding potential together with an associated p-value.


> <question-title></question-title>
>
> 1. What are features that RNAcode uses to detect reagions with coding potential?
>
> > <solution-title></solution-title>
> >
> > 1. RNAcode uses the conservation of a reading frame, synonymous mutations and the absence of stop codons
> >    in the individual sequences of the different species to detect coding potential.
> >
> {: .solution}
>
{: .question}

## Consensus structure prediction

We now calculate secondary structures for a set of aligned RNAs using **RNAalifold** {% icon tool %}. This tool reads aligned RNA sequences and calculate their minimum free energy (mfe) structure, partition function (pf) and base pairing probability matrix.

> <hands-on-title> Compute consensus structure </hands-on-title>
>
> 1. **RNAalifold** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Clustal file"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - In *"Algorithm Options"*:
>        - *"Calculate partition function"*: `1: Calculate the partition function and base pairing probability matrix`
>    - In *"General Options"*:
>        - *"Colored secondary structure image"*: `Yes`
>        - *"Colored and annotated alignment image"*: `Yes`
>    - In *"Structure constraints"*:
>        - *"Constraints"*: `Don't use constraints`
>        - *"Shape reactivity data"*: `Don't use shape reactivity data`
>    - In *"Model Options"*:
>        - *"Weight of the covariance term"*: `0.5`
>        - *"Penalty for non-compatible sequences in the covariance term"*: `0.6`
>        - *"Use ribosum scoring matrix"*: `Yes`
>
{: .hands_on}

**RNAalifold** returns the mfe structure in bracket notation, its energy, the free energy of the thermodynamic ensemble and the frequency of the mfe structure in the ensemble to stdout. It also produces Postscript files with plots of the resulting secondary structure graph and a "dot plot" of the base pairing matrix.

> <question-title></question-title>
>
> 1. What is covariation of structure in terms of homologous RNA sequences?
>
> > <solution-title></solution-title>
> >
> > 1. Covariation is the conservation of functionally important base-pairs after mutation events.
> >    If a base-pair is disrupted by a mutation that relevant for obtaining the correct spatial structure
> >    a second mutation of the pairing nucleotide can restore the base-pair. The more of these mutations
> >    that can be observed over multiple specices the higher the likelyhood that this base-pair is of
> >    biological relevance.
> >
> {: .solution}
>
{: .question}

## Stable RNA secondary structure prediction

To predict structurally conserved and thermodynamically stable RNA secondary structures from our data, we use RNAz. It can be used in genome wide screens to detect functional RNA structures, as found in noncoding RNAs and cis-acting regulatory elements of mRNAs.

> <hands-on-title> Predict the stable RNA secondary structure </hands-on-title>
>
> 1. **RNAz** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input Alignment File"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Use mononucleotide shuffled z-scores"*: `Yes`
{: .hands_on}

RNAz screens the input alignment for conserved secondary structures which are classified by a Support Vector Machine
to be structured RNA.


> <question-title></question-title>
>
> 1. Why is a conserved secondary structure an indicator that the region detected by RNAz possessed a biologial function?
>
> > <solution-title></solution-title>
> >
> > 1. Biologically functional RNA molecules need to fold into the necessary spatial structure to perform their function.
> >    If the secondary structure is conserved over larger phylogenetic distances this is a strong indicator that
> >    the function is preserved in multiple species and therefore has been maintained despite mutation events.
> >
> {: .solution}
>
{: .question}

# Analysis using RNA-aware multi alignment tool

We will now run a similar analysis using a RNA-aware multi alignment tool and will expand the analysis to the RNA family study.

## Align sequences and extract representative sequences

> <hands-on-title> Run multi-alignment tool and extract representative sequences from the alignment </hands-on-title>
>
> 1. **LocARNA Multiple Aligner** {% icon tool %} with the following parameters:
>    - *"Input type"*: `Fasta input (strict)`
>        - {% icon param-file %} *"Sequence input"*: `output` (Input dataset)
>    - *"Alignment mode"*: `Global alignment
                (LocARNA)`
>    - *"Output options"*: ``
>    - In *"Scoring parameters"*:
>        - *"Indel opening score"*: `-900`
>        - *"Indel score"*: `-250`
>        - *"Type of sequence score contribution"*: `Use ribofit`
>    - In *"Heuristic parameters"*:
>        - *"Restrict alignable positions by maximum difference"*: `Maximal difference of aligned positions`
>    - In *"Constraint parameters"*:
>        - *"Anchor constraints"*: `Don't load anchor constraints from bed file`
>
> 2. **Select Sequences** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input clustal alignment"*: `clustal_strict` (output of **LocARNA Multiple Aligner** {% icon tool %})
>
{: .hands_on}

In contrast to the alignment with MAFFT we now aligned the individual sequences considering the secondary
structure. This will make it easier to test for a conserved structure with RNAz, but it will
make it more difficult to detect common open reading frames and coding potential.


> <question-title></question-title>
>
> 1. What is the difference of aliging the individual sequnces with MAFFT or with LocaRNA?
>
> > <solution-title></solution-title>
> >
> > 1. MAFFT is a general purpose alignment tool based on sequence information only, while LocaRNA uses
> >    secondary structure information. In the MAFFT alignment we expect to see nucleotides to be in the
> >    same column based on reading frame information. In LocaRNA output there is a preference to align
> >    nucleotides in the same column if they are part of the same base-pair, even when the sequence is
> >    not similar.
> >
> {: .solution}
>
{: .question}

## Analysis of the representative sequences from alignment

As with output of MAFFT, we can now analyze more in depth the representative sequences of the alignment: analyze for potential protein coding regions, compute the consensus structure and predict the stable RNA secondary structure.

### Analysis for potential protein coding regions

> <hands-on-title> Identify potential protein coding regions </hands-on-title>
>
> 1. **RNAcode** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Multiple Alignment"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Break long alignment blocks"*: `Process original alignment`
>    - *"Scoring parameters"*: `Default`
>    - *"Create colored plots in EPS format"*: `Create Plots`
{: .hands_on}


> <question-title></question-title>
>
> 1. Is it easier to detect coding potential with the alignment output of MAFFT or of LocaRNA?
>
> > <solution-title></solution-title>
> >
> > 1. It is easier to detect with the MAFFT output, since LocaRNA gives preference to maintain base-pairs over
> >    the reading frames.
> >
> {: .solution}
>
{: .question}

### Consensus structure prediction

> <hands-on-title> Predict the consensus structure </hands-on-title>
>
> 1. **RNAalifold** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Clustal file"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - In *"Algorithm Options"*:
>        - *"Calculate partition function"*: `1: Calculate the partition function and base pairing probability matrix`
>        - *"Most Informative Sequence"*: `Yes`
>    - In *"General Options"*:
>        - *"Colored secondary structure image"*: `Yes`
>        - *"Colored and annotated alignment image"*: `Yes`
>    - In *"Structure constraints"*:
>        - *"Constraints"*: `Don't use constraints`
>        - *"Shape reactivity data"*: `Don't use shape reactivity data`
>    - In *"Model Options"*:
>        - *"Weight of the covariance term"*: `0.6`
>        - *"Penalty for non-compatible sequences in the covariance term"*: `0.5`
>        - *"Use ribosum scoring matrix"*: `Yes`
{: .hands_on}


> <question-title></question-title>
>
> 1. Compare the predicted secondary structure with the one predicted with the MAFFT alignment?
>
> > <solution-title></solution-title>
> >
> > 1. TODO add answer
> >
> {: .solution}
>
{: .question}

### Stable RNA secondary structure prediction

> <hands-on-title> Predict the stable RNA secondary structure </hands-on-title>
>
> 1. **RNAz** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input Alignment File"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Use mononucleotide shuffled z-scores"*: `Yes`
>    - *"Use decision model for structural alignments"*: `Yes`
{: .hands_on}


> <question-title></question-title>
>
> 1. Compare the result of the prediction vs the one performed with the MAFFT alignment?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. TODO add answer
> >
> {: .solution}
>
{: .question}

## Analysis of RNA family

We would like now to construct a probabilistic model representing structure and sequence of the RNA family and extract information about this family

### RNA family building

**cmbuild** {% icon tool %} is a tool which makes consensus RNA secondary structure profiles, and uses them to search nucleic acid sequence databases for homologous RNAs, or to create new structure-based multiple sequence alignments. It builds a covariance model of an RNA multiple alignment and uses the consensus structure to determine the architecture of the CM.

> <hands-on-title> Create RNA family model </hands-on-title>
>
> 1. **Build covariance models** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequence database"*: `stockholm` (output of **LocARNA Multiple Aligner** {% icon tool %})
>    - *"These options control how consensus columns are defined in an alignment"*: `automatic (--fast)`
>    - *"Options controlling relative weights"*: `Henikoff (--wgb)`
>    - *"Options controlling effective sequence number"*: `Turn off the entropy weighting strategy (--enone)`
>
{: .hands_on}

The output of **cmbuild** contains information about the size of the input alignment (in aligned columns and # of sequences), and about the size of the resulting model.

In addition to writing CM(s) to the output file, **cmbuild** also outputs a single line for each model created to stdout. Each line has the following fields:
- `aln`: the index of the alignment used to build the CM
- `idx`: the index of the CM in the output file
- `name`: the name of the CM
- `nseq`: the number of sequences in the alignment used to build the CM
- `eff nseq`: the effective number of sequences used to build the model
- `alen`: the length of the alignment used to build the CM
- `clen`: the number of columns from the alignment defined as consensus (match) columns
- `bps`: the number of basepairs in the CM
- `bifs`: the number of bifurcations in the CM
- `rel entropy`: CM: the total relative entropy of the model divided by the number of consensus columns
- `rel entropy`: HMM: the total relative entropy of the model ignoring secondary structure divided by the number of consensus columns
- `description`: description of the model/alignment.


> <question-title></question-title>
>
> 1. What does a RNA family model represent?
> 2. What is the practical application of a RNA family model?
>
> > <solution-title></solution-title>
> >
> > 1. RNA family model capture the sequence and structure information of the aligned individual sequences of a RNA family with a probabilistic model, also called covariance model.
> > 2. RNA family models can be used for homology search that considers the secondary structure of the RNA family and enables to identify novel members of the RNA family.
> {: .solution}
>
{: .question}

### RNA family statistics and visualization

To get a better idea about the RNA family, we can extract summary statistics for the covariance model generated by **cmbuild**, e.g. number of nodes, information content of sequence and structure.

> <hands-on-title> Generate RNA family statistics </hands-on-title>
>
> 1. **Summary statistics** {% icon tool %} with the following parameters:
>    - *"Subject covariance models"*: `Covariance model from your history`
>        - {% icon param-file %} *"Covariance models file from the history."*: `cmfile_outfile` (output of **Build covariance models** {% icon tool %})
>
{: .hands_on}

CMstat outputs different statistics about the generated RNA family model. Most useful is the information content of the HMM and the CM part of the model.
Since the HMM part models the sequence information and the CM part models also the structure information we expect a rather low information content
for protein coding sequences. However for functional ncRNAs we expect that there is a relevant secondary structure encoded by the model and therefore
a higher CM information content.


> <question-title></question-title>
>
> 1. How is the information content for HMM and CM components of the RNA family model different for ncRNA and protein coding genes?
>
> > <solution-title></solution-title>
> >
> > 1. Protein coding genes often posses no or only limited conserved secondary structure, we therefore expect that the CM information content that models the secondary structure will be low.
> >    For non-coding RNAs a higher CM information content is expected.
> >
> {: .solution}
>
{: .question}

We can also visualize the RNA family model by showing nodes and probabilites of the model.

> <hands-on-title> Visualize the RNA family model </hands-on-title>
>
> 1. **cmv** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input model"*: `cmfile_outfile` (output of **Build covariance models** {% icon tool %})
>    - In *"Common parameters"*:
>        - {% icon param-file %} *"Input stockholm alignment"*: `stockholm` (output of **LocARNA Multiple Aligner** {% icon tool %})
{: .hands_on}

The text output of RNA family models is very verbose and therefor clunky to inspect in detail. If you are interested how the probabilies and
state contained in the model the visualization simplyfies to investigate the model by displaying the different nodes of the model
with additional information.


> <question-title></question-title>
>
> 1. How many nodes does the newly constructed model posses?
>
> > <solution-title></solution-title>
> >
> > 1. TODO answer question
> >
> {: .solution}
>
{: .question}

## Homologous RNA sequence search

The RNA family model can also be used to detect further potentially homologous RNA sequences using **cmsearch** from Infernal.

Infernal is used to search sequence databases for homologs of structural RNA sequences, and to make sequence- and structure-based RNA sequence alignments. Infernal needs a profile from a structurally annotated multiple sequence alignment of an RNA family with a position-specific scoring system for substitutions, insertions, and deletions. Positions in the profile that are basepaired in the consensus secondary structure of the alignment are modeled as dependent on one another, allowing Infernal’s scoring system to consider the secondary structure, in addition to the primary sequence, of the family being modeled. Infernal profiles are probabilistic models called “covariance models”, a specialized type of stochastic context-free grammar (SCFG)

> <hands-on-title> Search for homologous RNA sequences </hands-on-title>
> 1. Import the sequence database file from [Zenodo]() or from the shared data library
> 2. **Search covariance model(s)** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequence database"*: `output` (Input dataset)
>    - *"Subject covariance models"*: `Covariance model from your history`
>        - {% icon param-file %} *"Covariance models file from the history."*: `cmfile_outfile` (output of **Build covariance models** {% icon tool %})
{: .hands_on}

Finally we want to search other genomes with the RNA family model constructed for our input alignment.

> <question-title></question-title>
>
> 1. How many potentially homolog sequences were detected by searching with the model?
>
> > <solution-title></solution-title>
> >
> > 1. TODO answer question
> >
> {: .solution}
>
{: .question}

# Conclusion
{:.no_toc}

In summary we applied to different alignment strategies to test if the input alignment is encoding a protein and if there is a conserved RNA secondary structure indicating
a biological function as non-coding RNA. We learned that using RNA secondary structure information during alignment can be beneficial to detect conserved RNA secondary
structures and disrupt open-reading frame information. Therefore to test for RNA-secondary structure it is useful to apply a structure aware alignment tool, while
testing for proteins coding potential should be done with a generic alignment tool. We used a dual-function RNA for the examples to simplify the examples
