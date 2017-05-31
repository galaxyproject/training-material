# Introduction

Extraction of useful information from raw microbiota sequences is a complex process
with numerous steps:

.. _data_processing:

.. figure:: /assets/images/framework/workflows/asaim_workflow.png

These steps can be get together in 4 main steps (corresponding to the 4 colors), with numerous sub-steps.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. [Pretreatments](#pretreatments)
> 2. [Taxonomic analyses](#taxonomic_analyses)
> 3. [Functional analyses](#functional_analyses)
> 4. [Combination of taxonomic and functional results](#taxonomic_functional_analyses)
> {: .agenda}

# Pretreatments

Before any analyses (taxonomic or functional), raw sequences have to be pre-processed 
with quality control and sequence sorting. 

In this tutorial, we will only focus on steps for single-end input sequences. But 
several notes will teach how to do with paired-end sequences.

Four main files are created during pretreatments:

- Non rRNA sequences
- 16S rRNA sequences
- 18S rRNA sequences
- Other rRNA sequences

## Quality control and treatment

As described in :ref:`description of quality estimation <framework-tools-available-pretreatments-control-quality-estimation>`, several sequence parameters must be checked to ensure that raw data looks good and then reduce bias in data analysis.

In this tutorial, quality control and treatments are made using PRINSEQ, as described in :ref:`quality treatments <framework-tools-available-pretreatments-control-quality-treatment>`:

- Elimination of sequences:
    - with length inferior to 60 bp (to eliminate sequences with too few information)
    - with a mean quality score inferior to 15 (to eliminate bad sequences)
    - with more than 2% of N bases (to eliminate sequences with too few usefull information)
- Trimming of sequences on right end when the mean quality score over a window of 5 bp is inferior to 20 (to improve conserved part of sequences)

To apply quality treatment with PRINSEQ on raw sequences, you click on `Control quality` on `Pretreatment` section on left pannel. Two tools will be proposed and you choose `PRINSEQ`. The central panel will be then filled with possible options to execute PRINSEQ

- The type of library (single-end, here)
- The FastQ file (input sequence file, it will be automatically proposed the downloaded file) 
- Parameters
    - Filtering
        - Filtering of sequences based on their length
            - Filtering of too smal sequences with a minimum length of 60 bp
            - No filtering of too big sequences
        - Filtering of sequences based on quality score
            - No filtering of sequences based on their minimum score
            - No filtering of sequences based on their maximum score
            - Filtering of sequences based on their mean score
                - Filtering of sequences with too small mean score with a minimum mean score of 15
                - No filtering of sequences with too high mean score
        - Filtering of sequences based on their base content
            - No filtering of sequences based on their GC percentage
            - No filtering of sequences based on their number of N bases
            - Filtering of sequences based on their percentage of N bases with a maximal N percentage of 2%
            - No filtering of sequences with characters other than A, T, C, G and N
        - Filtering of sequences based on their complexity with a maximum DUST score of 7
    - Trimming
        - No trimming from 3'-end
        - No trimming from the ends
        - No tail trimming
        - Trimming by quality score
            - No trimming by quality score from the 5'-end
            - Trimming by quality score from the 3'-end with a minimum mean quality threshold of 20, computed on a sliding window of 5bp moving by 5bp


For paired-end sequences, both sequence files are quality treated together to limit issue during paired-end assembly. This is done by selected `paired-end` in library type

Once the boxes are green, quality treatments are done. With chosen parameters on our datasets (217,386 input sequences with a mean length of bp), we expect conservation of 215,444 (99.09%) sequences with a mean length of 239.29 bp. Sequences are filtered mainly because of length.

In quality control of EBI metagenomic workflow, 88.43% of sequences (192,248) are conserved. 

> ### :nut_and_bolt: Comments
> For downstream tools such as SortMeRNA, one sequence file is expected. Paired-end sequences have then to be assembled. You can use FastQJoin with default parameters:
>
>    - Minimum of 6 bp overlap is required to join pairs
>    - Maximum 8% differences within region of overlap
{: .comment}


# Dereplication

Dereplication corresponds to identification of unique sequences in a dataset to conserve only one copy of each sequence in the dataset and then reduce the dataset size without loosing information.

In ASaiM, this task can be done with `VSearch dereplication` of `VSEARCH` suite :cite:`rognes_vsearch:_2015`. 

Sequence file with good quality sequences (PRINSEQ) are in FASTQ format. `VSearch` tools require FASTA file. So, the file has to be formatted using `Extract` in  `Manipulate sequence files` section (`Common tools`):

3 files are generated:

- A file with sequences in FASTA format
- A file with quality sequences
- A file with a report

To dereplicate, you execute `Vsearch dereplication` on the sequence file:

In the dataset, 3 sequences (<1%) are removed using dereplication. In EBI metagenomic workflow, ~ 4.74% of sequences are removed during this step of dereplication.

# rRNA (rDNA) sorting

Metagenomic and metatranscriptomic data are constitued of different types of sequences: sequences corresponding to CDS, sequences corresponding to ribosomal sequences (rDNA or rRNA), ... 

Useful functional information are present in sequences corresponding to CDS, and taxonomic information in sequences corresponding to ribosomomal sequences. To enhance downstream analysis such as extraction of functional or taxonomic information, it is important to :ref:`sort sequences into rRNA and non rRNA <framework-tools-available-pretreatments-manipulate-rna>`.

For this task, we use SortMeRNA :cite:`kopylova_sortmerna:_2012`. This tool filter RNA sequences based on local sequence alignment (BLAST) against rRNA databases. With SortMeRNA, 8 rRNA databases are proposed:

- A Rfam database for 5.8S eukarya sequences
- A Rfam database for 5S archea/bacteria sequences
- A SILVA database for 16S archea sequences
- A SILVA database for 16S bacteria sequences
- A SILVA database for 18S eukarya sequences
- A SILVA database for 23S archea sequences
- A SILVA database for 23S bacteria sequences
- A SILVA database for 28S eukarya sequences

16S and 18S sequences are mainly used in taxonomic analyses. So to limitate bias due to numerous sequences, it is interesting to extract these sequences from the other rRNA sequences.

So, the step of sequence sorting is split into 3 sub-step:

![](../../images/sequence_sorting.png)

SortMeRNA has to be executed 3 times, with different databases. 

For this sequence sorting, you click on `Manipulate RNA` in `Pretreatment` section on left pannel. You click then on `Filter with SortMeRNA`. Such as for PRINSEQ, parameters for SortMeRNA can be chosed in central panel:

2 output files are generated:

- A sequence file with `aligned reads` (sequences similar to rRNA databases)
- A sequence file with `rejected reads` (sequences non similar to rRNA databases)


In this tutorial, SortMeRNA has to be executed 3 times:

1. First execution of SortMeRNA to split sequences between rRNA and non rRNA sequences
    - Query sequences: output of dereplication
    - rRNA databases: All
2. Second execution of SortMeRNA to split rRNA sequences between 16S rRNA and non 16S rRNA sequences
    - Query sequences: `Aligned reads` of first SortMeRNA execution (rRNA sequences)
    - rRNA databases: SILVA 16S archea and SILVA 16S bacteria
3. Third execution of SortMeRNA to split non 16S rRNA sequences between 18S rRNA and other rRNA sequences
    - Query sequences: `Rejected reads` of second SortMeRNA execution (non 16S rRNA sequences)
    - rRNA databases: SILVA 18S eukarya

For the dataset, we obtain 1,550 sequences (0.72%) predicted as rRNA sequences. In EBI metagenomics, the percentage is similar with 0.50%.

With these sequences sorted as rRNA and non rRNA, we can run:

- Taxonomic analyses
- Functional analyses 

# Taxonomic analyses

To identify micro-organisms populating a sample and their proportion, we use :ref:`taxonomic and phylogenetic approaches <framework-tools-available-taxonomic-assignation>`. Indeed, in these approaches, each reads or sequences are assigned to the most plausible microbial lineage.

Most tools used to estimate sample composition use 16S rRNA genes as marker for bacteria and archea and 18S rRNA genes for eukaryota. Other tools, such as MetaPhlAn :cite:`segata_metagenomic_2012,truong_metaphlan2_2015`, propose alternatives based on more general clade-specific marker genes.

We use such approaches to analyze taxonomy of sequences with MetaPhlAn2 :cite:`truong_metaphlan2_2015` on all rRNA sequences. This tool infer the presence and read coverage of clade-specific markers to detect taxonomic clades and their relative abundance.

## Taxonomic assignation

This tool is available on left panel in `Assign taxonomy on non rRNA sequences` (`STRUCTURAL AND FUNCTIONAL ANALYSIS TOOLS`). In this tutorial, you execute it on all sequences before :ref:`SortMeRNA execution <framework-tutorial-pretreatments-rna-sorting>`.

In the dataset, we obtain a text file with 59 lines, each line corresponding to a taxonomic assignation (represented at different taxonomic level) with its relative abundance.

.. _metaphlan_2_output:

.. figure:: /assets/images/framework/tutorial/metaphlan_2_output.png

## Formatting

The output in plain text is not easy to interpret: all taxonomic levels are mixed together.


## Visualization

To extract information, we can also use visualization tools to get graphical representations of MetaPhlAn 2 output. Two solutions can be used: GraPhlAn or KRONA.

### Interactive visualization

Krona :cite:`ondov_interactive_2011` is a visualization tool for intuitive exploration of relative abundances of taxonomic classifications. Krona requires a formatted input file.

MetaPhlAn2 output has then to be formatted using `Format MetaPhlAn2 output for Krona` (in `Assign taxonomy on non rRNA sequences`, `Taxonomic assignation`):

.. _format_metaphlan2:

.. figure:: /assets/images/framework/tutorial/format_metaphlan2.png

Krona can then be called (`Krona pie chart from taxonomic profile`, in `Visualize data`, in `Post-treatments`)

.. _krona_call:

.. figure:: /assets/images/framework/tutorial/krona_call.png

Krona produces an interactive HTML file. The content can be visualized inside Galaxy environment by clicking on `View data` on top right of Krona output in right panel.

.. _krona_output:

.. figure:: /assets/images/framework/tutorial/krona_output.png

This visualization is similar to the one on EBI metagenomic.

## Static, easy to export visualization

Alternatively, `GraPhlAn <https://bitbucket.org/nsegata/graphlan/wiki/Home>`_ is a tool for producing circular representation of taxonomic analyses, easily exportable. This tool requires 2 files: a tree and an annotation file.

However, MetaPhlAn produces only a :ref:`text file <metaphlan_2_output>`. We need to use a tool to extract tree and annotations from MetaPhlAn output. We use `export2graphlan <https://bitbucket.org/CibioCM/export2graphlan>`_, available in section `Visualize data` (in `Post-treatments`). Numerous parameters modulates informations in annotation file. For our dataset, we fix :

- Levels to annotate in the tree: 5
- Levels to annotate in the external legend: 6,7
- Title font size: 15
- Default size for clades not found as biomarkers: 10
- Minimum value of biomarker clades: 0
- Maximum value of biomarker clades: 250
- Font size: 10
- Minimum font size: 8
- Maximum font size: 12
- Font size for the annotation legend: 11
- Minimum abundance value for a clade to be annotated: 0
- Number of clades to highlight: 100
- Row number contaning the names of the features: 0
- Row number containing the names of the samples: 0

We decide to display the maximum of clade (100, here). If you want more or less, you can modulate the number of clades to highlight. And if you want to change displayed annotations, you can change levels to annotate.

This tool will generate two outputs (a tree and an annotation files). These two outputs have to be combined in first GraPhlAn script (`Modify an input tree for GraPhlAn`, in `Visualize data`):

.. _graphlan_annotate_parameters:

.. figure:: /assets/images/framework/tutorial/graphlan_annotate_parameters.png

This tool generates a PhyloXML file, input file for GraPhlAn.

GraPhlAn is available in `Visualize data` section (`Post-treatments`). It generates an output file (an image) corresponding to circular representation of MetaPhlAn outputs. Available parameters have impact on output file format, size, ...

.. _graphlan_parameters:

.. figure:: /assets/images/framework/tutorial/graphlan_parameters.png

With our dataset, we obtain a nice graphical representation of taxonomic diversity inside our sample, with circle radius being proportional to relative abundance of the corresponding clade.

.. _graphlan_metaphlan_output:

.. figure:: /assets/images/framework/tutorial/graphlan_metaphlan_output.svg

After these taxonomic analyses, we can then run :ref:`functional analyses <framework-tutorial-functional-analysis>`.

# Functional analyses

Investigation of sample composition give an insight on "What organisms are present in our sample". We now want to know "What are they doing in that sample ?", with metabolic analyses.

## Metabolic analysis

For this investigation, we need to affiliate sequences to a protein database. We choose for this task to use `HUMAnN2 <http://huttenhower.sph.harvard.edu/humann2>`_. HUMAnN profiles the presence/absence and abundance of gene families and microbial pathways in a community from metagenomic or metatranscriptomic sequencing data.

HUMAnN2 is available in `Analyze metabolism` (`Functional assignation` section). We execute HUMAnN2 on non rRNA sequences:

.. _humann2_param:

.. figure:: /assets/images/framework/tutorial/humann2_param.png

3 output files are generated:

- A file with abundance of found UniRef50 gene families
- A file with coverage of found Metacyc pathways
- A file with abundance of found Metacyc pathways

These 3 files give detailed insights into gene families and pathways. This is interesting when we want to look to a particular pathway or to check abundance of a given gene families. However, when we want a broad overview of metabolic processes in a community, we need tools to regroup gene families or pathways into global categories.

## Broad overview of metabolic analysis

To get global categories from HUMAnN2 outputs, we decide to use `Gene Ontology <http://geneontology.org/>`_. Gene ontology project is a collaborative project to get a 3 structures ontologies to describe gene products in terms of their associated biological processes, cellular components and molecular functions. 

HUMAnN2 gives opportunity to regroup UniRef50 gene family abundances into GO abundances. However, these GO terms are still too precise to get a good overview of metabolic processes. 

Gene Ontology Consortium proposes `GO slim <http://geneontology.org/page/go-slim-and-subset-guide>`_, which are cut-down versions of the GO ontologies to give a broad overview of the ontology content. In our case, we use metagenomic GO slim terms developed by Jane Lomax and the InterPro group. 

To regroup HUMAnN2 output containing UniRef50 gene family abundances into abundances of metagenomic GO slim term, we use `Group humann2 uniref50 abundances to Gene Ontology (GO) slim terms <https://github.com/ASaiM/group_humann2_uniref_abundances_to_GO>`_ :cite:`batut_group_humann2_uniref_abundances_to_go:_2016`. This tool uses `GoaTools <https://github.com/tanghaibao/goatools>`_ :cite:`tang_goatools:_2016` to map GO terms to GO slim terms, HUMAnN2 to regroup abundances of UniRref50 gene families into abundances of metagenomc GO slim terms and custom Python scripts.

Tool to group HUMAnN2 UniRef50 abundances to Gene Ontology (GO) slim terms is available in `Analyse metabolism` (`Functional assignation` section). We execute it on HUMAnN2 output containing UniRef50 gene family abundance:

.. _group_humann2_uniref_abundances_to_go_param:

.. figure:: /assets/images/framework/tutorial/group_humann2_uniref_abundances_to_go_param.png

This tool generates 3 tabular outputs:

- A file with abundances of GO terms corresponding to molecular functions
- A file with abundances of GO terms corresponding to biological processes
- A file with abundances of GO terms corresponding to cellular components

## Visualization of metabolic analysis

The 3 previously generated tabular files with relative abundances of Gene Ontology slim terms can be visualised with barplots. A tool `Plot barplot with R` is available in `Visualize data` (`Post treatments` section):

.. _plot_barplot:

.. figure:: /assets/images/framework/tutorial/plot_barplot.png

Several graphical options are available such as margins, labels, bar color, ...

For our dataset, we execute this tool 3 times to obtain the following 3 graphics 

.. _cellular_component_abundance:

.. figure:: /assets/images/framework/tutorial/cellular_component_abundance.png
    :scale: 50 %

    Relative abundance of GO slim terms corresponding to cellular components
    

.. _biological_process_abundance:

.. figure:: /assets/images/framework/tutorial/biological_process_abundance.png
    :scale: 50 %

    Relative abundance of GO slim terms corresponding to biological processes


.. _molecular_function_abundance:

.. figure:: /assets/images/framework/tutorial/molecular_function_abundance.png
    :scale: 50 %

    Relative abundance of GO slim terms corresponding to molecular fonctions

Our analyses are done, you may want to :ref:`download the generated files <framework-tutorial-download>` and to :ref:`extract these numerous step into a workflow to reproduce it <framework-tutorial-transform>`. When you are done, do not forgot to :ref:`stop the Galaxy instance and clean the environment <framework-tutorial-clean>`. 

.. bibliography:: /assets/references.bib
   :cited:
   :style: plain
   :filter: docname in docnames

# Combination of taxonomic and functional results

.. _framework-tutorial-taxonomic-functional-analysis:

# Combination of taxonomic and functional analyses

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.
