---
layout: tutorial_hands_on

title: "Unicycler Assembly"
zenodo_link: "https://doi.org/10.5281/zenodo.940733"
level: Introductory
tags:
  - prokaryote
questions:
  - "I have short reads and long reads. How do I assemble a genome?"
objectives:
  - "Perform Quality Control on your reads"
  - "Perform a Small genome Assembly with Unicycler"
  - "Evaluate the Quality of the Assembly with Quast"
  - "Annotate the assembly with Prokka"
follow_up_training:
  -
    type: "internal"
    topic_name: assembly
    tutorials:
      - ecoli_comparison
time_estimation: "4h"
key_points:
  - "We learned about the strategies used by assemblers for hybrid assemblies"
  - "We performed an hybrid assembly of a bacterial genome and its annotation"
  - "Unicycler is a pipeline bases on Spades and Pilon dedicated to hybrid assembly of Small genomes"
  - "Combination of short and long reads helped us produce an almost perfect assembly"
contributors:
  - nekrut
  - delphine-l
  - slugger70
---

# The goal: *E. coli* C-1 assembly

In this tutorial we assemble and annotate the genome of *E. coli* strain [C-1](http://cgsc2.biology.yale.edu/Strain.php?ID=8232). This strain is routinely used in experimental evolution studies involving bacteriophages. For instance, now classic works by Holly Wichman and Jim Bull ({% cite Bull1997 %}, {% cite Bull1998 %}, {% cite Wichman1999 %}) have been performed using this strain and bacteriophage phiX174.

To sequence the genome we have obtained the strain from the [Yale E. coli Stock Center](http://cgsc2.biology.yale.edu/). The stock center sent us a filter paper disk infused with cells. The disk was placed in the center of an LB-agar plate. A single colony was picked and resuspended in a liquid LB medium, grown overnight, and genomic DNA was isolated. The DNA was then sequenced using two methods. To obtain high coverage, high accuracy data we used Illumina miSEQ to generated 250-bp paired end reads. To generate high length reads we used the Oxford Nanopore MinION machine.

Our goal is to reconstruct and annotate the full genome of *E. coli* C-1. As you will see in this tutorial a combination of many short, high accuracy reads with long, error-prone reads helps us produce an almost perfect assembly.

> ### Outline step-by-step
>
> In this tutorial, we will deal with:
>
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Background on data and tools


## The data

In this tutorial we will assemble a genome using two types of input data: (1) Illumina 250 bp paired-end reads and (2) Oxford Nanopore reads.

### Illumina data

We generated 9,345,897 250 bp read pairs (library preparation performed on genomic DNA fragmented to mean size of 600 bp). However, to make sure that you can complete this tutorial in a finite amount of time we have downsampled (reduced in size) to 1,000,000 paired end reads - just enough to produce an accurate assembly.

### Oxford Nanopore Data

There are 12,738 [2d-reads](http://www.nature.com/nmeth/journal/v12/n4/fig_tab/nmeth.3290_SF13.html). Maximum read length is 27,518 bp. The distribution of reads lengths looks like this:

![Nanopore read length distribution](../../images/ont_length.png "Distribution of nanopore read lengths.")

You can see that there many reads under the second peak with median of approximately 7.5 kb.

> ### {% icon warning %} Oxford Nanopore Data Format
> Oxford Nanopore machines output
 data in [fast5](http://bioinformatics.cvr.ac.uk/blog/exploring-the-fast5-format/) format that contains additional information besides sequence data. In this tutorial we assume that this data is *already* converted into [fastq](https://en.wikipedia.org/wiki/FASTQ_format). An additional tutorial dedicated to handling fast5 datasets will be developed shortly.
{: .warning}


## The tools

In this analysis we will perform two tasks: (1) assembly and (2) annotation. Below we will briefly outline the main ideas behind these two procedures and will describe the tools we will be using.

### Assembly

> ### {% icon comment %} Knowing your assembly
>
> Here we assume that you know a thing or two about assembly process. If you don't: look at the slides accompanying this tutorial as well as other tutorials is this section.
{: .comment}

![Logo unicycler](https://github.com/rrwick/Unicycler/raw/master/misc/logo.png)

For assembly we will be using [Unicycler](https://github.com/rrwick/Unicycler) (also see publication {% cite Wick2017 %}). Unicycler is designed specifically for *hybrid assembly* (that is, using both short- and long-read sequencing data) of small (e.g., bacterial, viral, organellar) genomes. In our hands it has produced complete high quality assemblies. Unicycler employs a multi-step process that utilizes a number of software tools:

![Unicycler process](../../images/unicycler.png "Simplified view of the Unicycler assembly process (From {% cite Wick2017 %}) In short, Unicycler uses SPAdes (see below) to produce an assembly graph, which is then bridged (simplified) using long reads to produce the longest possible set of contigs. These are then polished by aligning the original short reads against contigs and feeding these alignments to Pilon - an assembly improvement tool. {% cite Wick2017 %}")



As you can see Unicycler relies heavily on SPAdes ({% cite Bankevich2012 %}) and [Pilon](https://github.com/broadinstitute/pilon/wiki). We will briefly describe these two tools.

#### Spades

##### Multisized deBruijn graph

Assemblers usually construct graphs for *k*-mers of a fixed size. We have noted that when *k* is small it is difficult to resolve the repeats. If *k* is too large a corresponding graph may be fragmented (especially if read coverage is low). SPAdes uses several values for *k* (that are either manually set or inferred automatically) to create a *multisized* graph that minimized tangledness and fragmentation by combining various *k*-mers ({% cite Bankevich2012 %})):

![Multigraph approach implemented in SPAdes](../../images/multiGraph.jpg "Multisized de Bruijn graph. A circular Genome CATCAGATAGGA is covered by a set of Reads consisting of nine 4-mers, {ACAT, CATC, ATCA, TCAG, CAGA, AGAT, GATA, TAGG, GGAC}. Three out of 12 possible 4-mers from Genome are missing from Reads (namely {ATAG,AGGA,GACA}), but all 3-mers from the Genome are present in the Reads. (A) The outside circle shows a separate black edge for each 3-mer from Reads. Dotted red lines indicate vertices that will be glued. The inner circle shows the result of applying some of the glues. (B) The graph DB(Reads, 3) resulting from all the glues is tangled. The three h-paths of length 2 in this graph (shown in blue) correspond to h-reads ATAG, AGGA, and GACA. Thus Reads<sub>3,4</sub> contains all 4-mers from Genome. (C) The outside circle shows a separate edge for each of the nine 4-mer reads. The next inner circle shows the graph DB(Reads, 4), and the innermost circle represents the Genome. The graph DB(Reads, 4) is fragmented into 3 connected components. (D) The multisized de Bruijn graph DB (Reads, 3, 4). Figure and text from {% cite Bankevich2012 %}.")

##### Read pair utilization

While the use of paired reads and mate pairs is key to genome assembly, and not new, SPAdes utilizes so called paired DeBruin graphs to take the advantage of the paired end data. One of the key issues with paired DeBruin graphs is that the resulting genome assemblies do not tolerate variability in insert sizes: The initial formulation of paired DeBruijn graphs assumed constant distance between pairs of reads. In practice this distance is always variable. SPAdes performs *k*-bimer (these are *k*-mers derived from *paired* reads) adjustment to identify exact or nearly-exact distances for each *k*-bimer pair.

##### Error correction

Sequencing data contains a substantial number of sequencing errors that manifest themselves as deviations (bulges and non-connected components) within the assembly graph. One way to improve the graph before assembly it is to minimize the number of sequencing errors by performing error correction. SPAdes uses BayesHammer ({% cite Nikolenko2013 %}) to correct the reads. Here is a brief summary of what it does:

1. SPAdes (or rather BayesHammer) counts *k*-mers in reads and computes *k*-mer statistics that take into account base quality values.
2. A [Hamming graph](https://en.wikipedia.org/wiki/Hamming_graph) is constructed in which *k*-mers are nodes. In this graph edges connect nodes (*k*-mers) if they differ from each other by a number of nucleotides up to a certain threshold (the [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance)). The graph is central to the error correction algorithm.
3. Then Bayesian subclustering is done on the graph from the previous step. For each *k*-mer we now know the center of its subcluster.
4. **Solid** *k*-mers are derived from cluster centers and are assumed to be *error free*.
5. Solid *k*-mers are mapped back to the reads.
6. Reads are corrected using solid *k*-mers:

![Read correction with BayesHammer](../../images/readCorrection.jpg "Read correction. Black <em>k</em>-mers are solid. Grey <em>k</em>-mers are non-solid. Red <em>k</em>-mers are the centers of the corresponding clusters (two grey <em>k</em>-mers striked through on the right are non-solid singletons). As a result, one nucleotide is changed based on majority rule. (From {% cite Nikolenko2013 %})")

In the case of the full dataset, SPAdes error correction changed 14,013,757 bases in 3,382,337 reads - a substantial fraction of the full ~18 million read dataset.

#### Pilon

Pilon improves draft assemblies by using the information from the original reads aligned to the draft assembly. The following image from a publication by {% cite Walker2014 %} highlights the steps of this process:

![Pilon workflow](../../images/pilon.png "The left column depicts the conceptual steps of the Pilon process, and the center and right columns describe what Pilon does at each step while in assembly improvement and variant detection modes, respectively. During the first step (top row), Pilon scans the read alignments for evidence where the sequencing data disagree with the input genome and makes corrections to small errors and detects small variants. During the second step (second row), Pilon looks for coverage and alignment discrepancies to identify potential mis-assemblies and larger variants. Finally (bottom row), Pilon uses reads and mate pairs which are anchored to the flanks of discrepant regions and gaps in the input genome to reassemble the area, attempting to fill in the true sequence including large insertions. The resulting output is an improved assembly and/or a VCF file of variants. (From {% cite Walker2014 %})")

### Annotation

For annotation we are using [Prokka](https://vicbioinformatics.com/software.prokka.shtml) (also see {% cite Seemann2014 %}). It scans the assembly generated with Unicycler with a set of feature prediction tools and compiles a list of genome annotation. It predicts the following features (Table from {% cite Seemann2014 %}):

| Feature | Tool used by Prokka |
|---------|----------------------|
| Protein-coding sequences (CDS) | [Prodigal](https://github.com/hyattpd/Prodigal) |
| Ribosomal RNA genes | [RNAmmer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1888812) |
| Transfer RNA genes | [Aragorn](https://www.ncbi.nlm.nih.gov/pubmed/14704338) |
| Signal leader peptides | [SignalP](https://www.ncbi.nlm.nih.gov/pubmed/21959131) |
| Non-coding RNA genes | [Infernal](http://eddylab.org/infernal/) |

Prokka predicts protein-coding regions using a two step process. It first identifies coordinates of putative genes using [Prodigal](https://github.com/hyattpd/Prodigal) and then compares the gene sequence against databases of known sequences at protein level using [Blast+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) and [HMMer](http://hmmer.org/).

# Let's try it


## Load data and assess quality

In this example we will use a downsampled version of *E. coli* C-1 Illumina and ONT sequencing data. These include 3 files: forward and reverse reads for Illumina, and Long read file produced by ONT. All data are in [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format.

> ### {% icon hands_on %} Hands-on: Obtaining our data
>
> 1. Make sure you have an empty analysis history. Give it a name.
>
>    {% include snippets/create_new_history.md %}
>
> 2. Import the following file from [Zenodo](https://zenodo.org/record/940733)
>
>    ```
>    https://zenodo.org/record/940733/files/illumina_f.fq
>    https://zenodo.org/record/940733/files/illumina_r.fq
>    https://zenodo.org/record/940733/files/minion_2d.fq
>    ```
>
>    {% include snippets/import_via_link.md %}
>    {% include snippets/import_from_data_library.md %}
>
{: .hands_on}

If all goes well you will see datasets uploading and changing states from gray to green as shown below. The figure below also shows how datasets can be tagged.

![Datasets in History](../../images/starting_data.png  "Sequencing data loaded into Galaxy history. The full progression from gray (scheduling) to green (all OK) state is shown. To make it easier to identify datasets as we progress through the analysis we use so-called <em>Hashtags</em>. To tag a dataset: click on dataset to expand it (as shown in panel four); click the tag icon (<i class='fa fa-tags' aria-hidden='true'></i>) and a text field will appear. Add a tag (in this case <b>F</b>) pre-pended with hash (#). Hit enter. Do this for all three datasets and it will appear as in panel five.")

### Assess Read Quality

To assess quality we will use two tools: FastQC ({% cite FastQC %}) to generate quality statistics and multiQC ({% cite Ewels2016 %}) to summarize these statistics.

> ### {% icon hands_on %} Hands-on: Quality Control
>
> 1. **FastQC** {% icon tool %}:
>
>    - {% icon param-files %} *"Short read data from your current history"*: Select all three FastQ datasets simultaneously
>
> 2. **MultiQC** {% icon tool %}: to generate a summary of the FastQC reports with
>   - *"Which tool was used generate logs?"*: `FastQC`
>   - *"Type of FastQC output"*: `Raw data`
>   - {% icon param-files %} *"FastQC Output"*: `RawData` outputs of FastQC
>
{: .hands_on}

A quick look at quality score distribution will show a confusing picture:

![QC reported zoomed out](../../images/multiqc1.png "Because Illumina reads (green) are <b>much</b> shorter that ONT reads (red) the plot looks strange. ONT reads generally have low quality scores and so they are not really meaningful in the context of this technology. However, in the case of Illumina data they mean a lot...")

So let's zoom in into Illumina data:

![QC reported zoomed in](../../images/multiqc2.png "Zooming in shows quality distribution for Illumina reads. This is excellent data with mean base qualities above 30 across all reads.")

## Assembly with Unicycler

Now it is time to perform assembly.

> ### {% icon hands_on %} Hands-on: Unicycler Assembly
>
> 1. **Unicycler** {% icon tool %} with the following parameters :
>   - "Paired or Single end data?" to `Paired`
>   - "First Set of reads" to the forward reads file `f`
>   - "Second Set of reads" to the reverse reads file `r`
>   - "Long reads" to the MinION file
>   - Use default parameters
>
{: .hands_on}


> ### <i class="fa fa-cutlery" aria-hidden="true"></i> <i class="fa fa-coffee" aria-hidden="true"></i> Assembly takes time!
>
> There is no such thing as Assembly in real time. It takes time so it is a good time to have lunch or at least coffee. This Unicycler run will take anywhere between 90 minutes and two hours.
{: .warning}

## Assess Assembly quality with Quast

Quast ({% cite Gurevich2013 %}) is a tool providing quality metrics for assemblies, and can also be used to compare multiple assemblies. The tool can also take an optional reference file as input, and will provide complementary metrics.

> ### {% icon hands_on %} Hands-on: Assembly Quality
>
> 1. **Quast** {% icon tool %}: with the following parameters
>   - *"Contigs/scaffolds output file"*: Select the fasta file resulting from the Unicycler assembly.
>
{: .hands_on}

The Quast tool outputs assembly metrics as an html file with metrics and graphs. The image below looks exceptionally boring. This is a **good** thing:

![Quast Interface](../../images/quast_output.png  "Quast Output: Quast provides different statistics such as the number of contigs or scaffolds, the N50 and N75, and the total length of the assembly. You can also access 3 plots, the cumulative length of the contigs, the Nx, or the GC content.")

One can see that there are two (!) contigs. The largest contig is 4,576,290 bp (for comparison *E. coli* K12 MG1655 strain genome length is [4,656,144 bp](https://www.ncbi.nlm.nih.gov/nuccore/NZ_APIN00000000.1)) and the smallest is 4,581,676 (total length) - 4,576,290 (length of the largest) = 5,386 bp. When we analyzed this dataset for the first time we were initially puzzled by this second contig. But we quickly realized that this is simply the genome of bacteriophage [phiX174](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1) which is routinely used as a spike-in in Illumina sequencing. Thus we have two genomes: the one of *E.coli* C-1 and phiX174! We can now use Prokka to annotate our two genomes.

## Annotation with Prokka


> ### {% icon hands_on %} Hands-on: Annotation
>
> 1. **Prokka** {% icon tool %}:
>   - *"Contigs to annotate"*: Select the assembly ouput of `Unicycler`
>   - *"Genus name"*: `Escherichia`
>   - *"Species name"*: `coli`
>   - *"Strain name"*: `C-1`
>   - *"Use genus-specific BLAST database"*: `yes`
>
{: .hands_on}

Prokka outputs 10 datasets (including two log files). These are in various formats:

* **txt** : Provides Statistics on the annotation : number of CDS predicted, number of rRNA etc.
* **tbl** : Provides a tabulated list of annotated features.
* **fsa** : Nucleotide fasta file of the input contig sequence.
* **sqn** : [ASN.1](https://www.ncbi.nlm.nih.gov/Structure/asn1.html) format file for submission to GenBank.
* **ffn** : Nucleotide FASTA file of all the prediction transcripts.
* **faa** : Protein FASTA file of the translated CDS sequences.
* **fna** : Nucleotide fasta file of the input contig sequence.
* **gbk** : [GenBank](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) file.
* **gff** : [gff3](http://gmod.org/wiki/GFF3) file.

## Visualize the results in IGV

Let's look at the entire assembly and its annotation in the genome browser. We can do this using Integrated Genome Browser ([IGV](http://software.broadinstitute.org/software/igv/)).

Visualization requires a local installation of IGV. If you have IGV installed - just start it. If you don't - read on.

#### Starting IGV

Go to IGV [download page](http://software.broadinstitute.org/software/igv/download) and select one of the options. The one I would try first would be **Java Web Start**. Simply click the **Launch** button for 10 GB distribution.


> ### {% icon hands_on %} Hands-on: Visualization in IGV
>
>
> 1. Start IGV. It will look something like this:
>
>      ![IGV just started](../../images/igv1.png "IGV started by using Java Web Start. Note that it is currently showing Human Genome (hg38) build. This is obviously <b>not</b> what we want.")
>
> 2. Locate the output of Unicycler and expand it :
>
>      ![Unicycler result Visualization](../../images/unicycler_igv.png  "Expanding Unicycler assembly reveals IGV link.")
>
> 3. Click on the *local* link highlighted with orange outline. The browser will change:
>
>      ![IGV with Unicycler assembly](../../images/igv2.png "Now you can see the major contig shown in the browser window.")
>
> 4. Let's add Prokka annotations to the browser image. For this simply expand Prokka's GFF3 dataset and click on the *local* link:
>
>      ![Expanded GFF3 dataset representing Prokka annotations](../../images/prokka_item.png "Expanded GFF dataset generated with Prokka. Click on the <em>local</em> link (highlighted with orange outline) to display this dataset within IGV.")
>
{: .hands_on}

You will now see the annotations within the browser window:


![Prokka result Visualization](../../images/igv3.png  "Viewing annotated assembly withing IGV.")
