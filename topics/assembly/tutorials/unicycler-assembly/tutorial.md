---
layout: tutorial_hands_on
topic_name: assembly
tutorial_name: unicycler-assembly
---

# The goal: *E. coli* C-1 assembly

In this tutorial we assemble and annotate genome of *E. coli* strain [C-1](http://cgsc2.biology.yale.edu/Strain.php?ID=8232). This strain is routinely used in experimental evolution studies involving bateriophages. For instance, now classical works by Holly Wichman and Jim Bull ([Bull et al. 1997](https://www.ncbi.nlm.nih.gov/pubmed/9409816), [Bull & Wichman 1998](https://www.ncbi.nlm.nih.gov/pubmed/9767038), [Wichman et al. 1999](https://www.ncbi.nlm.nih.gov/pubmed/10411508)) have been performed using this strain and bacteriophage phiX174. 

To sequence the genome we have obtained the strain from the [Yale E. coli Stock Center](http://cgsc2.biology.yale.edu/). The stock center sent us a filter paper disk infused with cells. The disk was placed in the center of an LB-agar plate. A single colony was picked and resuspended in a liquid LB medium, grown overnight, and genomic DNA was isolated. The DNA was then sequenced using two methods. To obtain high coverage, high accuracy data we used Illumina miSEQ to generated 250-bp paired end reads. To generate high length reads we used the Oxford Nanopore MinIon machine. 

Our goal is to reconstruct and annotate the full genome of *E. coli* C-1. As you will see in this tutorial a combination of many short, high accuracy reads with long, error-prone reads helps us produce an almost perfect assembly.

# Background on data and tools

## The data

In this tutorial we assemble genome using two types of input data: (1) Illumina 250 bp paired-end reads and (2) Oxford Nanopore reads. 

### Illumina data

We generated 9,345,897 250 bp read pairs (library preparation performed on genomic DNA fragmented to mean size of 600 bp). However, to make sure that you can complete this tutorial in a finite amount of time we have downsampled (reduced in size) this to 50,000 paired end reads - just enough to produce a reasonable assembly.

### Oxford Nanopore Data

There are 12,738 [2d-reads](http://www.nature.com/nmeth/journal/v12/n4/fig_tab/nmeth.3290_SF13.html). Maximum read length is 27,518. The distribution of reads lengths looks like this:

![Nanopore read length distribution](../../images/ont_length.png "Distribution of nanopore read lengths.")

You can see that there many reads under the second peak with median of approximately 7.5 kb. 

> ### <i class="fa fa-warning" aria-hidden="true"></i> Oxford Nanopore Data Format
> Oxford Nanopore machines output
 data in [fast5](http://bioinformatics.cvr.ac.uk/blog/exploring-the-fast5-format/) format that contains additional information besides sequence data. In this tutorial we assume that these data *already* converted into [fastq](https://en.wikipedia.org/wiki/FASTQ_format). An additional tutorial dedicated to handling of fast5 datasets will be developed shortly. 
{: .warning-box}

## The tools

In this analysis we will perform two tasks: (1) assembly and (2) annotation. Below we will briefly outline main ideas behind these two procedures and will describe the tools we will be using.

### Assembly

> ### <i class="fa fa-lightbulb-o" aria-hidden="true"></i> Knowing your assembly
>
> Here we assume that you know a thing or two about assembly process. If you don't: look at the [slides](./slides) accompanying this tutorial as well as other tutorials is this section.
{: .info-box}

For assembly we will be using [Unicycler](https://github.com/rrwick/Unicycler) (also see publication by Wick:[2017](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595)). Unicycler is designed specifically for *hybrid assembly* (the one combines short and long read sequencing data) of small (e.g., bacterial, viral, organellar) genomes. In out hard it gas produced complete high quality assemblies. Unicycler employs a multi-step process that utilizes a number of software tools:

![Unicycler process](../../images/unicycler.png "Simplified view of the Unicycler assembly process (From Wick:2017). In short, Unicycler uses SPAdes (see below) to produce an assembly graph, which is then bridged (simplified) using long reads to produce longest possible set pf contigs. These are then polished by aligning original short reads against produced contigs and feeding these alignment to Pilon - an assembly improvement tool.")

As you can see Unicycler relies heavily on [SPAdes](http://cab.spbu.ru/software/spades/) and [Pilon](https://github.com/broadinstitute/pilon/wiki). We will briefly describe these two tools.

#### Spades

##### Multisized deBruijn graph

Assemblers usually constructing graphs for *k*-mers of a fixed size. We have noted that when *k* is small it is difficult to resolve the repeats. If *k* is too large a corresponding graph may become fragments (especially if read coverage is low). SPAdes uses several values for *k* (that are either manually set or inferred automatically) to create a *multisized* graph that minimized tangledness and fragmentation by combining various *k*-mers (see [Bankevich:2012](http://online.liebertpub.com/doi/full/10.1089/cmb.2012.0021)):

![Multigraph approach implemented in SPAdes](../../images/multiGraph.jpg "Multisized de Bruijn graph. A circular Genome CATCAGATAGGA is covered by a set Reads consisting of nine 4-mers, {ACAT, CATC, ATCA, TCAG, CAGA, AGAT, GATA, TAGG, GGAC}. Three out of 12 possible 4-mers from Genome are missing from Reads (namely {ATAG,AGGA,GACA}), but all 3-mers from Genome are present in Reads. (A) The outside circle shows a separate black edge for each 3-mer from Reads. Dotted red lines indicate vertices that will be glued. The inner circle shows the result of applying some of the glues. (B) The graph DB(Reads, 3) resulting from all the glues is tangled. The three h-paths of length 2 in this graph (shown in blue) correspond to h-reads ATAG, AGGA, and GACA. Thus Reads3,4 contains all 4-mers from Genome. (C) The outside circle shows a separate edge for each of the nine 4-mer reads. The next inner circle shows the graph DB(Reads, 4), and the innermost circle represents the Genome. The graph DB(Reads, 4) is fragmented into 3 connected components. (D) The multisized de Bruijn graph DB (Reads, 3, 4). Figure from [Bankevich:2012]") 

##### Read pair utilization

While the use of paired reads and mate pairs is not new (and key) to genome assembly, SPAdes utilizes so called paired DeBruin graphs to take the advantage of the paired end data. One of the key issues with paired DeBruin graphs is the resulting genome assemblies do not tolerate variability in insert sizes (the initial formulation of paired DeBruijn graphs assumed constant distance between pairs of reads). In practice this distance is always variable. SPAdes performs *k*-bimer (these are *k*-mers derived from *paired* reads) adjustment to identify exact of nearly-exact distances for each *k*-bimer pair.

##### Error correction

Sequencing data contains a substantial number of sequencing errors that manifest themselves as deviations (bulges and non-connected components) within the assembly graph. One of the ways to improve the graph even constructing it is to minimize the amount sequencing errors by performing error correction. SPAdes uses [BayesHammer](https://goo.gl/1iGkMe) to correct the reads. Here is a brief summary of what it does (see [Nikolenko:2013](https://goo.gl/1iGkMe)):

1. SPAdes (or rather BayesHammer) counts *k*-mers in reads and computed *k*-mer statistics that takes into account base quality values. 
2. [Hamming graph](https://en.wikipedia.org/wiki/Hamming_graph) is constructed for *k*-mers is which *k*-mers are nodes. In this graph edges connect nodes (*k*-mers) is they differ from each other by a number of nucleotides up to a certain threshold (the [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance)). The graph is central to the error correction algorithm.
3. At this step Bayesian subclustering of the graph produced in the previous step. For each *k*-mer we now know the center of its subcluster. 
4. **Solid** *k*-mers are derived from cluster centers and are assumed to be *error free*.
5. Solid *k*-mers are mapped back to the reads.
6. Reads are corrected using solid *k*-mers:

![Read correction with BayesHammer](../../images/readCorrection.jpg "Read correction. Black <em>k</em>-mers are solid. Grey <em>k</em>-mers are non-solid. Red <em>k</em>-mers are the centers of the corresponding clusters (two grey <em>k</em>-mers striked through on the right are non-solid singletons). As a result, one nucleotide is changed based on majority rule. (From [Nikolenko:2013])")

In the case of the full dataset SPAdes error correction changed 14,013,757 bases in 3,382,337 reads - a substantial fraction of the full ~18 million read dataset.

#### Pilon

Pilon improves draft assemblies by using the information from the original reads aligned to the draft assembly. The following image from a publication by [Walker:2014](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963) highlights the steps of this process:

![Pilon workflow](../../images/pilon.png "The left column depicts the conceptual steps of the Pilon process, and the center and right columns describe what Pilon does at each step while in assembly improvement and variant detection modes, respectively. During the first step (top row), Pilon scans the read alignments for evidence where the sequencing data disagree with the input genome and makes corrections to small errors and detects small variants. During the second step (second row), Pilon looks for coverage and alignment discrepancies to identify potential mis-assemblies and larger variants. Finally (bottom row), Pilon uses reads and mate pairs which are anchored to the flanks of discrepant regions and gaps in the input genome to reassemble the area, attempting to fill in the true sequence including large insertions. The resulting output is an improved assembly and/or a VCF file of variants. (From Walker:2014)")

### Annotation

For annotation we are using [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml) (also see [Seeman:2014](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu153)). It scans the assembly generated with Unicycler with a set of feature prediction tools and compiles a list of genome annotation. It predicts the following features (Table from [Seeman:2014](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu153)):

| Feature | Tool used by Prokka |
|---------|----------------------|
| Protein-coding sequences (CDS) | [Prodigal](https://github.com/hyattpd/Prodigal) |
| Ribosomal RNA genes | [RNAmmer](http://www.cbs.dtu.dk/cgi-bin/nph-runsafe?man=rnammer) |
| Transfer RNA genes | [Aragorn](https://www.ncbi.nlm.nih.gov/pubmed/14704338) |
| Signal leader peptides | [SignalP](https://www.ncbi.nlm.nih.gov/pubmed/21959131) |
| Non-coding RNA genes | [Infernal](http://eddylab.org/infernal/) | 

Prokka predicts protein-coding regions using a two step process. It first identifies coordinates of putative genes using [Prodigal](https://github.com/hyattpd/Prodigal) and then compares the gene sequence against databases of known sequences at protein level using [Blast+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) and [HMMer](http://hmmer.org/).

# Let's try it

> ### Outline step-by-step
>
> In this tutorial, we will deal with:
>
> 1. [Get the data](#get-the-data)
> 2. [Assess reads quality](#assess-read-quality)
> 3. [Assembly with Unicycler](#assemble-with-unicycler)
> 4. [Assess Assembly quality with Quast](#quast)
> 5. [Annotate with Prokka](#annotate-with-prokka)
> 6. [Visualize the results](#visualize-the-result)
{: .agenda}

## <a name="get-the-data">Load data and assess quality

In this example we will use a downsampled version of *E. coli* C-1 Illumina and ONT sequencing data. These include 3 files: forward and reverse reads for Illumina, and Long read file produced by ONT. All data are in [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format. 

### Load data into History

To load data into your Galaxy instance log in into Galaxy and create new history (if you are new to Galaxy see [Galaxy 101 tutorial](/topics/introduction/tutorials/galaxy-intro-101/tutorial.html) first). 

Click **Get data** icon as shown below (see [these slides](/topics/introduction/tutorials/galaxy-intro-get-data/slides.html) for an introduction on how to load data into Galaxy):

<hr>
![Get Data](../../images/get_data.png "Getting data into history starts with clicking <b>Get data</b> button")
<hr>

Open [Zenodo](https://zenodo.org/record/842795#.WcKdQtOGPOZ) link in a **new browser window** and right-click on dataset links:

<hr>
![Get Data](../../images/zenodo.png "Right click on links to copy them into clipboard")
<hr>

And paste them into the **Galaxy upload**:

<hr>
![Upload file](../../images/upload_file.png  "Uploading data into Galaxy. First (1) click <b>Paste/Fetch data</b> link. Next, paste URL copied from Zenodo. Finally (2), set type of all datasets to <tt>fastqsanger</tt>. Click <b>Start</b>.")
<hr>

If all goes well you will see this:

-----
![Datasets in History](../../images/starting_data.png  "Sequencing data loaded into Galaxy history.")
-----

### <a name="assess-read-quality">Assess Read Quality

To assess quality we will use two tools: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate quality statistics and [multiQC](http://multiqc.info/) to summarize these statistics.

First, run **FastQC** on all three datasets simultaneously:

-----
![Screenshots of FastQC interface](../../images/fastc_interface.png  "Using <b>fastqc</b> to compute the quality statistics for all reads. Note that multiple dataset selection button (<i class='fa fa-files-o' aria-hidden='true'></i>) is pressed and all three datasets are selected at the same time. ")
-----

Although FastQC generated graphical reports for each dataset we can look at everything at once using multiQC by selecting "raw" output of FastQC:

-----
![Screenshot of multiQC interface](../../images/multiqc_interface.png "Running <b>multiQC</b> requires selecting <em>RawData</em> output of <b>FastQC</b>. Again, note that multiple dataset selection button (<i class='fa fa-files-o' aria-hidden='true'></i>) is pressed and all <em>RawData</em> inputs are selected.") 
-----

A quick look at quality score distribution will show a confusing picture:

-----
![QC reported zoomed out](../../images/multiqc1.png "Because Illumina reads (green) are <b>much</b> shorted that ONT reads (red) the plot looks strange. ONT reads generally have low quality scores and so they are not really meaningful in the context of this technology. However, in case of Illumina data they mean a lot...")
-----

So let's zoom in into Illumina data:

-----
![QC reported zoomed in](../../images/multiqc2.png "Zooming in shows quality distribution for Illumina reads. This is excellent data with mean base qualities above 30 across all reads.")
-----

Our data is great and we can jump directly to assembly process. 

## <a name="assemble-with-unicycler"></a>Assembly with Unicycler 

The Unicycler tool takes fastqsanger files as inputs. If your files are identified as generic fastq files you will need to change the type of your files.

![Edit dataset attributes](../../images/edit_attribute.png  "Edit dataset attributes. Click on the pen bouton of a dataset in your history to edit its attributes.")

![Dataset attribute interface](../../images/change_type.png  "Change datatype. Click on the Datatype tab and select the appropriate type in the list, here fastqsanger.")

Repeat the process for the three datasets.

You can now run Unicycler to perform the assembly with the following parameters: 

* **Paired or Single end data?** : Select the appropriate option to describe you data. In this example we are using Paired end Data.
* **Select first set of reads** : Specify the dataset containing the forward reads, often specified by a "-1" in the file name, but specified here by the "R1".
* **Select second set of reads** : Specify the dataset containing the forward reads, often specified by a "-2" in the file name, but specified here by the "R2".
* **Select long reads** : Optional, here specify you Oxford Nanopore dataset.

![Unicyler Interface](../../images/Unicycler_interface.png  "Unicycler interface. Run Unicycler with your sequencing dataset in fastqsanger format.")

Unicycler returns two output files: a fasta file containing the result of the assembly, and a graph file.  You can then evaluate the quality of the resulting alignments by using the Quast tool on the fasta file.

## <a name="quast">Assess Assembly quality with Quast

[Quast](http://bioinf.spbau.ru/quast) is a tool providing quality metrics for assemblies, and can also be used to compare multiple assemblies. The tool can also take an optional reference file as input, and will provide complementary metrics.
For this tutorial we will simply use quast on the fasta file resulting from the Unicycler assembly.

![Quast Interface](../../images/Quast_Interface.png  "Quast Interface")

The Quast tool outputs assembly metrics as an html file with metrics and graphs.

![Quast Interface](../../images/quast_output.png  "Quast Output: Quast provides different statistics such as the number of contigs or scaffolds, the N50 and N75, and the total length of the assembly. You can also access 3 plots, the cumulative length of the contigs, the Nx, or the GC content.")

We can now use Prokka to annotate our genome.

### <a name="annotate-with-prokka">Annotation with Prokka

Run Prokka with the following paramters:

* **Contigs to annotate** : Specify the fasta file resulting from your assembly with Unicycler.
* **Locus tag prefix** : Specify the format you desire for your locus tags. By Default PROKKA.
* **Locus tag counter increment** : By default 1, but a 10 increment facilitate the insertion of new genes when manually correcting the annotation.
* **Force GenBank/ENA/DDJB compliance** : Select "yes" if you desire to force the GenFank Locus tag formatting. If you do so be aware of the length limitation. Here we select no for more convenience.
* **Add 'gene' features for each 'CDS' feature** : Select yes to get the gene feature in addition to the CDS feature in the gff3.
* **Genus name** : Specify the Genus of your organism. Here "Escherichia".
* **Species name** : Specify the species of your organism. Here "Coli".
* **Strain name** : Specify the strain of your organism. Here "C".
* **Kingdom** : Select the kingdom to which your organism belong. Here "Bacteria".
* **Use genus-specific BLAST database** : Select "yes" to use the genus-specific Blast.

![Prokka Interface](../../images/Prokka_Interface.png  "Prokka Interface")

Prokka outputs 10 datasets. One of the is the Prokka log, another is the error repport,  but 8 are diverse result files : 
* **txt file** : Provides Statistics on the annotation : number of CDS predicted, number of rRNA etc.
* **tbl file** : Provides a tabulated list of annotated features.
* **fsa file** : Nucleotide fasta file of the input contig sequence.
* **sqn file** : ASN1 format file for submission to GenBank.
* **ffn file** : Nucleotide FASTA file of all the prediction transcripts.
* **faa file** : Protein FASTA file of the translated CDS sequences.
* **fna file** : Nucleotide fasta file of the input contig sequence.
* **gbk file** : GenBank file.
* **gff file** : gff3 file.

### <a name="visualize-the-result">Visualization

You can visualize Prokka annotations using Integrative Genome Browser (IGV).
First, download and install [IGV](http://software.broadinstitute.org/software/igv/) Open an instance of IGV on you computer, and then import the genome file from galaxy by clicking on the "display with IGV local" 

![Unicycler result Visualisation](../../images/unicycler_result.png  "Unicycler output can be sent to a local instance of IGV")

You can then send the gff file resulting from the annotation with Prokka.

![Prokka result Visualisation](../../images/prokka_result.png  "Prokka output can then be sent to a local instance of IGV")

You can then visualize the result of your analysis in IGV 


![IGV whole assembly visualization](../../images/IGV.png  "IGV Whole assembly visualization. You can visualize the whole assembly and gene density by selecting the 'all' view.")


![IGV zoomed visualization](../../images/zoomed_igv.png  "IGV zoomed visualization. For more detail view, select a contig and zoom on a region. You can see more information about a feature by passing over it with you cursor")





