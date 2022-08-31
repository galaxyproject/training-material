---
layout: tutorial_hands_on

title: Comparative gene analysis
zenodo_link: https://zenodo.org/record/7034885
questions:
- I have several genomes assemblies that are not annoatated (or I do not trust annotations)
- I am interested to compare structure of a particular gene across these genome assemblies
- How do I do that?
objectives:
- Provide a quick method for identifying genes of interest in unannotated or newly assembled genomes
time_estimation: 2H
key_points:
- You can easily 
contributors:
- nekrut
tags:
- annotation
- vgp
- cookbook
requirements:
  -
    type: "internal"
    topic_name: sequence-analysis
    tutorials:
      - quality-control
  -
    type: "internal"
    topic_name: galaxy-interface
    tutorials:
      - collections
abbreviations:
  ORF: Open Reading Frame
subtopic: eukaryote
priority: 8
---

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Introduction
{:.no_toc}

Despite the rapidly increasing number of fully assembled genomes few genomes are well annotated. This is especially true for large eukaryotic genomes with their complex gene structure and abundance of pseudogenes. And of course do not forget about the [Murthy's law](https://en.wikipedia.org/wiki/Murphy%27s_law): if you are intersted in a particular gene the chances are that it will not be annotated in your genome of interest. In this tutorial we will demonstrate how to compare gene structures across a set of vertebrate genomes. So ...

> > ### {% icon question %} What I want:
> > - I work with a gene _X_
> > - I would like to compare the structure of gene _X_ across _N_ genomes
> {: .code-in}
>
> > ### {% icon galaxy-chart-select-data %} What I have:
> > - I know the gene's name
> > - I know which species I'm interested in
> > - I know where to find genomes of these species
> {: .code-out}
{: .code-2col}

> ### {% icon interactive_tour %} What I will get:
> - Interactive graphs showing location of the gene across your species of choice. These will allow you to see the absence/presence of the genes across genomes, to detect potential duplications, predogenization events, re-arrangements etc.
> - Phylogenetic trees for individual exons of the gene. The trees will give you an idea of potential unusual evolutionary dynamics for the gene.
{: .warning}

------

# Logic

The analysis follow the following logic (also see the following figure):

1. Pick a gene and select genomes to analyze
1. Get amino acid translations for all exons of my gene of interest
1. Identify genomes of interest
1. Extract amino acid sequences and genome coordinates for all possible open reading frames (ORFs) from all genomes in my set

    The following steps 5, 6, and 7 are performed using a single Galaxy workflow (grey area in the following figure):

1. Align - Find matches between exons of the gene of interest and ORFs
1. Intersect - Compute genome coordinates of matches
1. Build phylogenetic trees for individual exons

    The final step is performed using Jupyter notebook

1. Create a comparative plot showing distribution of matches across genomes of interest

![Logic](../../images/gene-centric/gene_analysis.svg)

------

# Example history

This [example history](https://usegalaxy.org/u/cartman/h/xbp1vgpsample) contains results of the analysis described in this tutorial. 

{% snippet faqs/galaxy/histories_import.md %}

You can use this history to inderstand the input datasets as well as outputs of the entire analysis. The key items in the history are labelled with <kbd>tags</kbd>:

> ### {% icon code-in %} Input dataset in the example history
> - <kbd>EXONS</kbd> - amino acid translation of exons of the gene of interest (*XBP-1*)
> - <kbd>ORF_BED</kbd> - coordinates of predicted ORFs in the genomes of interest
> - <kbd>DiamondDB</kbd> - database and amino acid translations of predicted ORFs in the genomes of interest
{: .code-in}

> ### {% icon code-out %} Outputs in the sample history
> - <kbd>PlottingData</kbd> - summary necessary for plotting comparative genome graphs
> - <kbd>Trees</kbd> - phylogenetic trees for each exon 
{: .code-out}

> ### {% icon warning %} A suggestion!
> Importing and looking around this history is very helpful for understanding how this analysis works!
{: .warning}

-------

# Practice

## Step 1: Pick a gene and select genomes to analyze

In this example we will compare structure of X-box protein 1 gene (*XBP-1*) across a set of five vertebrate genomes sequenced by [VGP](https://vertebrategenomesproject.org/) consortium.

## Step 2: Get amino acid translations for all exons of my gene of interest

This step is manual and is performed outside Galaxy.

The best annotated vertebrate is ... human. To obtain amino acid translation of individual exons of *XBP-1* you can use [UCSC Genome browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A28794555%2D28800597&hgsid=1433476153_N1Iet0RvAaPQIpPuVl2qbi74BrC3). The "UCSC annotations of RefSeq RNAs (NM_\* and NR_\*)" track shows two alternative transcripts produced by *XBP-1* : spliced and unspliced (here 'spliced' and 'unspliced' do not refer to the "normal' splicing process. To learn more about this fascinating gene read here).  

![XBP-1](../../images/gene-centric/xbp-1.svg)

Clicking on the names of the transcripts open a UCSC page specific to that transcript. From there you can obtain genomic sequences corresponding to all exons. You can then translate these sequences using any available online tool (e.g., [NCBI ORFinder](https://www.ncbi.nlm.nih.gov/orffinder/)). 

In this particular case we did not create translation of all exons. Instead, we created translation of exon 1  and two terminal exons. The first exon is shared between the two transcripts. It will allow us to anchor the beginning of the gene. The two terminal exons are different between the two transcripts. Because of this we create two alternative tralations per exon: `s-p2` and `s-p12` for the spliced version and `u-p1` and `u-p2` for unspliced transcript. A FASTA file containing all translation is shown below:

```
>xbp-1u-p1
GNEVRPVAGSAESAALRLRAPLQQVQAQLSPLQNISPWILAVLTLQIQ
>xbp-1u-p2
LISCWAFWTTWTQSCSSNALPQSLPAWRSSQRSTQKDPVPYQPPFLCQWG
RHQPSWKPLMN
>xbp-1s-p12
GNEVRPVAGSAESAAGAGPVVTPPEHLPMDSGGIDSSDSE
>xbp-1s-p3
SDILLGILDNLDPVMFFKCPSPEPASLEELPEVYPEGPSSLPASLSLSVG
TSSAKLEAINELIRFDHIYTKPLVLEIPSETESQANVVVKIEEAPLSPSE
NDHPEFIVSVKEEPVEDDLVPELGISNLLSSSHCPKPSSCLLDAYSDCGY
GGSLSPFSDMSSLLGVNHSWEDTFANELFPQLISV
>xbp1_ex1
MVVVAAAPNPADGTPKVLLLSGQPASAAGAPAGQALPLMVPAQRGASPEA
ASGGLPQARKRQRLTHLSPEEKALR
```

## Step 3: Identify and upload genomes of interest

This tutorial has been initially designed for the analysis of data produced by the Veretebrate Genome Consortium [VGP](https://vertebrategenomesproject.org/). However, it is equally suitable for any genomic sequences (including prokaryotic ones).

In this section we first show how to upload sample datasets. These datasets were intentionally made small. All subsequent steps of this tutoirial are performed using these sample data. We then demonstrate how to upload full size genomes from the VGP data repository called [GenomeArk](https://vgp.github.io/genomeark/).

### Uploading sample data

> ### {% icon hands_on %} Sample Data upload
>
> - Create a new history for this tutorial
> - Import the files from [Zenodo]({{ page.zenodo_link }}) using the following URLs:
>
> ```
> https://zenodo.org/record/7034885/files/aGasCar1.fa
> https://zenodo.org/record/7034885/files/bTaeGut2.fa
> https://zenodo.org/record/7034885/files/fScoJap1.fa
> https://zenodo.org/record/7034885/files/mCynVol1.fa
> https://zenodo.org/record/7034885/files/mHomSapT2T.fa
> ```
> {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

### Uploading VGP data from GenomeArk

Galaxy provides a direct connection to [GenomeArk](https://vgp.github.io/genomeark/) from its "Upload Data" tool. To access GenomeArk data you need to:

> ### {% icon hands_on %} Upload VGP data from GenomeArk
>
> - Create a new history for this tutorial
> - Import genome assembly FASTA files from GenomeArk:
>
>    - Open the file {% icon galaxy-upload %} __Upload Data__ menu
>    - Click on __Choose remote files__ button at the bottom of the Upload interface
>    - Choose __Genome Ark__
>    - Pick assembles you would like to include in your analysis (use Search box on top)
>
>
> Generally for a given species you want the assembly with the highest version number. For example, Human (_Homo sapiens_) had three assemblies at the time of writing: `mHomSap1`, `mHomSap2`, and `mHomSap3`. Pick `mHomSap3` in this situation. 
>
> Inside each assembly folder you will see typically see different subfolders. The one you should generally be interested in will be called `assembly_curated`. Inside that folder you will typically see separate assemblies for maternal (indicated with `mat`) and paternal (indicated with `pat`) haplotypes. In the analyses here we are choosing as assembly representing heterogametic haplotype (`pat` in the case of Human).
>
> - Repeat this process for all assemblies you are interested in.
>
{: .hands_on}

### Combine uploaded sequences into a dataset collection

After uploading data from GenomeArk you will end up with a number of FASTA datasets in your history. To proceed with the analysis you will need to combine these into a single dataset collection:

> ### {% icon hands_on %} Combine genome assemblies into a dataset collection
>
> To create a collection from genome assemblies you just uploded into GenomeArk follow the steps in the video below. Obviously, you only want to click checkboxes next to FASTA files (assemblies) you want to include and leave everything else out.
>
> {% snippet faqs/galaxy/collections_build_list.md %}
>
> <p align="center"><iframe width="560" height="315" src="https://www.youtube.com/embed/6ZU9hFjnRDo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></p>
>
{: .hands_on}

## Step 4: Extract amino acid sequences and genome coordinates for all ORFs

### Extract ORFs from genomic sequences

Now we can go ahead and find all ORFs in all genome assemblies we have bundled into a collection. For this we will use tool called [ORFiPy](https://github.com/urmi-21/orfipy):

> ### {% icon hands_on %} Find ORFs witrh ORFiPy
>
> Run {% tool [ORFiPy](toolshed.g2.bx.psu.edu/repos/iuc/orfipy/orfipy/0.0.4+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Find ORFs in"*: The collection containing genome assemblies we just created (click "Dataset collection" {% icon param-collection %} button on left-side of this input field)
>    - *"Select options"*: check the box next to `BED` and `Peptides`
>    - *"Minimum length of ORFs"*: type the number `99`. We do not want ORFs to be too short as their number will be too high
>    - *"Output ORFs bound by Stop codons?"*: set this option to `Yes`. 
>
{: .hands_on}

Because this operation is performed on a dataset collection. it will produce a dataset collection as an output. Actually, it will produce two because we selected `BED` and `Pepetides` options. One collection will contain amino acid sequences for all ORFs (it will typically be called `ORFs on collection ... (FASTA format)`) and another their coordinates (called `ORFs on collection ... (BED format)`)
)

This will produce two new dataset collections in your histiory: one containing coordinates of ORFs and the other containing their amino acid translations. Because genomes are large these will likely be tens of millions of ORFs identified this way.

### Creating __Diamond__ database

Because we will be using the [Diamond](https://github.com/bbuchfink/diamond) tool to find matches between our gene of interest and ORF translations we need to convert FASTA files into Diamond database using __Diamond makedb__ tool:

> ### {% icon hands_on %} Create Diamond database
>
> Run {% tool [Diamond makedb](toolshed.g2.bx.psu.edu/repos/bgruening/diamond/bg_diamond_makedb/2.0.15+galaxy0) %} on the collection containing amino acid translations of ORFs generated using __ORFiPy__:
>    - {% icon param-collection %} "Input reference file in FASTA format"*: Select collections containing amino acid (FASTA Protein) output of __ORFiPy__ 
>
{: .hands_on}

At this point we have three input datasets that would allow us to find and visualize matches between the gene of interest and genome asseblies:

1. Diamond database of ORF translations from genome assemblies
2. Coordinates and frame information about ORFs in BED format
3. Amino acid translation of exons from the gene of interest

## Steps 5, 6, and 7: Finding matches and building trees

To find location of genes, we will use the following [workflow](/training-material/topics/genome-annotation/tutorials/gene-centric/workflows/) that is available as a part of this tutorial. To use this workflow you need to import it into your Galaxy instance as described [here](/training-material/topics/genome-annotation/tutorials/gene-centric/workflows/).

![WF](../../images/gene-centric/wf.png)

The workflow takes three inputs: 

> ### {% icon code-in %} Workflow inputs
>
> 1. <kbd>EXONS</kbd> - Amino acid translation of exons from the gene of interest ([Step 2](#step-2-get-amino-acid-translations-for-all-exons-of-my-gene-of-interest) of this tutorial)
> 1. <kbd>DiamondDB</kbd> - Diamond database of ORF translations from genome assemblies ([Step 4](#step-4-extract-amino-acid-sequences-and-genome-coordinates-for-all-orfs) of this tutorial)
> 1. <kbd>ORF BED</kbd> - Coordinates and frame information about ORFs in BED format ([Step 4](#step-4-extract-amino-acid-sequences-and-genome-coordinates-for-all-orfs)  of this tutorial)
>
{: .code-in}

It produces two primary outputs:

> ### {% icon code-out %} Results
>
> 1. <kbd>Trees</kbd> - Phylogenetic trees for each input exon as [Newick](https://en.wikipedia.org/wiki/Newick_format) file
> 1. <kbd>PlottingData</kbd> - A summary table of exon matches for each genome
>
{: .code-out}

The overall logic of the workflow is as follows:

1. Perform [__Diamond__](https://github.com/bbuchfink/diamond) search to identify matches between the gene of interest (<kbd>EXONS</kbd>) and ORF sets from each genomes (<kbd>DiamondDB</kbd>).
2. Intersect information about the matches with BED file containing ORF coordinates (<kbd>ORF BED</kbd>). This allows us to know genomic position of the ORFs and their frames (1, 2, 3 or -1, -2, -3).
3. Extract matching parts of the sequences and generated multiple alignments using [__MAFFT__](https://mafft.cbrc.jp/alignment/software/). This is done for each amino acid sequence in <kbd>EXONS</kbd>. Thus in the case of *XBP-1* there will be five sets of alignemnts. 
4. Build phylogenetic tree for each alignment produced in the previous step using the [Neighbor Joining method](https://en.wikipedia.org/wiki/Neighbor_joining) - a simple and quick way to obtain a rough phylogeny.
5. Use information from step 1 to compute genome coordinates of matches. This is done by combining th einformation about genomic positions of ORFs from the <kbd>ORF BED</kbd> file and local alignment information generated during step 1.

## Step 7: Looking at the trees

After running the workflow phylogenetic trees will be saved into a collection named `Join neighbors on ...`. This collection will also be labelled with tag <kbd>Trees</kbd>. To visualize the trees:

> ### {% icon hands_on %} Visualize the trees
>
> - Expand the collection by clicking on it
> - Click on any dataset
> - Click the {% icon galaxy-barchart %} icon:
>     ![Graph icon](../../images/gene-centric/display_tree.png)
> - A selection of tree visualization tools will appear in the center pane of the interface
> - Select `Phylogenetic Tre Visualization` and play with it:
>     ![The tree](../../images/gene-centric/tree.png)
{: .hands_on}

## Step 8: Generating comparative genome graph

Another workflow output will represent a single file summarzing genomic location of matches between each of the genomes in our dataset and amino acid translation fo exons from the gene of interest. It will be called `Mapping report` and will have tag <kbd>PlottingData</kbd> associated with it. To plot the data contained in this file we will use external Jupyter notebook (note that Jupyter can be run directly from Galaxy, but to make this tutorial runnable on any Galaxy instance we will use an internal notebook server). 

### Starting notebook

The notebook can be accessed [from here](https://colab.research.google.com/drive/1smTpRejBb7c02LiIxMPNVDjOucLAsD81?usp=sharing). You do need to have a Google Account to be able to use the notebook.

### Proving data

The data is provided to the notebook by setting the `dataset_url` parameter in this cell:

```python
# Paste link to the dataset here
dataset_url = 'https://usegalaxy.org/api/datasets/f9cad7b01a472135a8abd43f91f8d3cf/display?to_ext=tabular'
```

This is a URL pointing to one of the workflow outputs: `Mapping report` with the <kbd>PlottingData</kbd> tag. To copy the URL click on the {% icon galaxy-link %} icon adjacent to the dataset:

![Copying the link](../../images/gene-centric/copy_link.png)

Running the notebook will generate two graphs explained in the next section.

### Interpreting the graph




# About the gene

XBP-1 is another example of a gene with overlapping reading frames. In this case the switch between two reading frames occurs not because of alternative splicing as was the case for CDKN2a. Instead, a highly specialized RNA endonuclease, IRE1, excises a 26 nucleotide spacer from XBP-1 mRNA. This converts the so-called “unspliced” form of the transcript (XBP-1U) to the “spliced” form XBP-1S (Note that the term “spliced” is misleading here. The XBP-1U is already processed by splicing machinery and contains no canonical introns. The “spliced” is simply used to indicate the removal of the 26 nucleotide spacer). Because the 26 is not divisible by three the XBP-1S transcript is translated in a different frame from the point of cleavage. The activation of IRE1 and resulting removal of the spacer is triggered by presence of unfolded proteins in the endoplasmic reticulum and the IRE1-XBP1 pathway is one of the three major unfolded protein response systems in higher eukaryotes that is conserved all the way to yeast (with XBP-1 homologue HAC-1). The XBP-1 is present in all assembled genomes described here. Yet the most interesting aspect is the duplication history of this gene revealed by our analysis. XBP-1 duplications have occurred at several points within vertebrate evolution. For example, humans contain an XBP-1 pseudogene on chromosome 5 and our analysis indicates the presence of another pseudogene in Tammar wallaby (Fig. SX). However, in philippine flying lemur (Cynocephalus volans) there appear to be two functional copies (Fig. XBP Top panel) located within the same linkage group (scaffold) and separated by over 160 Mb of sequence. The copies have distinct arrangements. The XBP-1-Right is the “original” copy retaining the gene structure expected from comparison with the human orthologue. The XBP-1-Left copy appears to be an insertion of a processed mRNA -- a hallmark of retrotransposition (REF). Exon 1 and “unspliced” fragments p1 and p2 (see Methods) lay within the same reading frame and in +1 phase relative to the “spliced” fragments as would be in the processed transcript. Both reading frames appear to be intact and sharing XX% identity with the XBP-1-Right copy and IRE1 recognition site is undisturbed. 


> ### {% icon tip %} Tip: Getting help
>
> For questions about using Galaxy, you can ask in the [Galaxy help forum](https://help.galaxyproject.org/).
>
{: .tip}

# Conclusion
{:.no_toc}

{

# Acknowledgements
{:.no_toc}

Thank to VGP consortium