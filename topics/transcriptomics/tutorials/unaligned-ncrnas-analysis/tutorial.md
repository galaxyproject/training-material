---
layout: tutorial_hands_on

title: Analyse unaligned ncRNAs
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- bebatut
- pavanvidem
---

# Introduction
{:.no_toc}

<!-- This is a comment. -->

***TODO***: *Add a general introduction to the topic*

***TODO***: *Add a description of the dataset and what you want to do with it*

**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

The first step in any analysis of unaligned ncRNAs is to run a multi alignment tool on the sequences of the ncRNAs. Using the multi alignment, we can extract some representative sequences, analyze them for protein coding potential, compute consensus structure and predict stable secondary structure. 

For the multi-alignment, two types of tools are used: one generic multi-aligner and one RNA-aware. In this tutorial, we will see how to deal with both types of aligners to analyze unaligned ncRNAs.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]() or from the shared data library
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% include snippets/change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% include snippets/add_tag.md %}
>
{: .hands_on}

# Analysis using generic multi alignment tool

## Align sequences and extract representative sequences

> ### {% icon hands_on %} Hands-on: Run multi-alignment tool and extract representative sequences from the alignment
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

***TODO***: *Comment quickly the generated outputs*

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

We can now analyze more in depth the representative sequences of the alignment: analyze for potential protein coding regions, compute the consensus structure and predict the stable RNA secondary structure.

## Analysis for potential protein coding regions

We would like first to know if there is any potential protein coding regions left in the alignment. We use **RNAcode** {% icon tool %} which predicts protein coding regions in an alignment of homologous nucleotide sequences. The prediction is based on evolutionary signatures typical for protein genese, i.e. the presence of synonyomous/conservative nucleotide mutations, conservation of the reading frame and absence of stop codons.

**RNAcode** {% icon tool %} does not rely on any species specific sequence characteristics whatsoever and does not use any machine learning techniques. 

> ### {% icon hands_on %} Hands-on: Identify potential protein coding regions
>
> 1. **RNAcode** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Multiple Alignment"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Scoring parameters"*: `Default`
>    - *"Create colored plots in EPS format"*: `Create Plots`
>
{: .hands_on}

**RNAcode** {% icon tool %} reports local regions of unusual high coding potential together with an associated p-value.

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Consensus structure prediction

We now calculate secondary structures for a set of aligned RNAs using **RNAalifold** {% icon tool %}. This tool reads aligned RNA sequences and calculate their minimum free energy (mfe) structure, partition function (pf) and base pairing probability matrix. 

> ### {% icon hands_on %} Hands-on: Compute consensus structure
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

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Stable RNA secondary structure prediction

To predict structurally conserved and thermodynamically stable RNA secondary structures from our data, we use RNAz. It can be used in genome wide screens to detect functional RNA structures, as found in noncoding RNAs and cis-acting regulatory elements of mRNAs.

> ### {% icon hands_on %} Hands-on: Predict the stable RNA secondary structure
>
> 1. **RNAz** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input Alignment File"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Use mononucleotide shuffled z-scores"*: `Yes`
{: .hands_on}

***TODO***: *Comment quickly the generated outputs*

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

# Analysis using RNA-aware multi alignment tool

We will now run a similar analysis using a RNA-aware multi alignment tool and will expand the analysis to the RNA family study.

## Align sequences and extract representative sequences

> ### {% icon hands_on %} Hands-on: Run multi-alignment tool and extract representative sequences from the alignment
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

***TODO***: *Comment quickly the generated outputs*

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Analysis of the representative sequences from alignment

As with output of MAFFT, we can now analyze more in depth the representative sequences of the alignment: analyze for potential protein coding regions, compute the consensus structure and predict the stable RNA secondary structure.

### Analysis for potential protein coding regions

> ### {% icon hands_on %} Hands-on: Identify potential protein coding regions
>
> 1. **RNAcode** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Multiple Alignment"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Break long alignment blocks"*: `Process original alignment`
>    - *"Scoring parameters"*: `Default`
>    - *"Create colored plots in EPS format"*: `Create Plots`
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding and compare with results from MAFFT*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

### Consensus structure prediction

> ### {% icon hands_on %} Hands-on: Predict the consensus structure
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

***TODO***: *Consider adding a question to test the learners understanding and compare with results from MAFFT*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

### Stable RNA secondary structure prediction

> ### {% icon hands_on %} Hands-on: Predict the stable RNA secondary structure
>
> 1. **RNAz** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input Alignment File"*: `clustal` (output of **Select Sequences** {% icon tool %})
>    - *"Use mononucleotide shuffled z-scores"*: `Yes`
>    - *"Use decision model for structural alignments"*: `Yes`
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding and compare with results from MAFFT*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Analysis of RNA family

We would like now to construct a probabilistic model representing structure and sequence of the RNA family and extract information about this family

### RNA family building

**cmbuild** {% icon tool %} is a tool which makes consensus RNA secondary structure profiles, and uses them to search nucleic acid sequence databases for homologous RNAs, or to create new structure-based multiple sequence alignments. It builds a covariance model of an RNA multiple alignment and uses the consensus structure to determine the architecture of the CM.

> ### {% icon hands_on %} Hands-on: Create RNA family model
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

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

### RNA family statistics and visualization

To get a better idea about the RNA family, we can extract summary statistics for the covariance model generated by **cmbuild**, e.g. number of nodes, information content of sequence and structure.

> ### {% icon hands_on %} Hands-on: Generate RNA family statistics
>
> 1. **Summary statistics** {% icon tool %} with the following parameters:
>    - *"Subject covariance models"*: `Covariance model from your history`
>        - {% icon param-file %} *"Covariance models file from the history."*: `cmfile_outfile` (output of **Build covariance models** {% icon tool %})
>
{: .hands_on}

***TODO***: *Comment quickly the generated outputs*

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

We can also visualize the RNA family model by showing nodes and probabilites of the model.

> ### {% icon hands_on %} Hands-on: Visualize the RNA family model
>
> 1. **cmv** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input model"*: `cmfile_outfile` (output of **Build covariance models** {% icon tool %})
>    - In *"Common parameters"*:
>        - {% icon param-file %} *"Input stockholm alignment"*: `stockholm` (output of **LocARNA Multiple Aligner** {% icon tool %})
{: .hands_on}

***TODO***: *Comment quickly the generated outputs*

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Homologous RNA sequence search

The RNA family model can also be used to detect further potentially homologous RNA sequences using **cmsearch** from Infernal.

Infernal is used to search sequence databases for homologs of structural RNA sequences, and to make sequence- and structure-based RNA sequence alignments. Infernal needs a profile from a structurally annotated multiple sequence alignment of an RNA family with a position-specific scoring system for substitutions, insertions, and deletions. Positions in the profile that are basepaired in the consensus secondary structure of the alignment are modeled as dependent on one another, allowing Infernal’s scoring system to consider the secondary structure, in addition to the primary sequence, of the family being modeled. Infernal profiles are probabilistic models called “covariance models”, a specialized type of stochastic context-free grammar (SCFG)

> ### {% icon hands_on %} Hands-on: Search for homologous RNA sequences
> 1. Import the sequence database file from [Zenodo]() or from the shared data library
> 2. **Search covariance model(s)** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Sequence database"*: `output` (Input dataset)
>    - *"Subject covariance models"*: `Covariance model from your history`
>        - {% icon param-file %} *"Covariance models file from the history."*: `cmfile_outfile` (output of **Build covariance models** {% icon tool %})
{: .hands_on}

***TODO***: *Comment quickly the generated outputs*

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.