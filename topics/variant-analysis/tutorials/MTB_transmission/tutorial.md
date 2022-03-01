---
layout: tutorial_hands_on

title: 'Identifying tuberculosis transmission links: from SNPs to transmission clusters'
zenodo_link: ''
objectives:
- Create a SNP alignment
- Calculate pairwise SNP distances between MTB samples
- Identify transmission clusters based on SNP distances
- Study the emergence and spread of drug resistance based on transmission analysis.
time_estimation: 2H
enable: true
level : Intermediate
key_points:
- Clustering is a useful tool to detect transmission links between patients and oubreak investigation.
- Clustering can be used to investigate the transmission of certain traits, like drug resistance.
- Clustering does not provide information about particular transmission events nor their directionality (who infected whom).
- Clustering is very much influenced by sampling. Lower sampling proportions and shorter sampling timeframes lead to lower clustering rates that shoud not be confounded with lack of transmission.
contributors:
- Galo A. Goig
- Daniela Brites
- Christoph Stritt
---



# Introduction
{:.no_toc}


Now you are familiar with the process of genome sequencing, quality control of sequencing data, mapping and
variant calling. Additionally, you know how to run `TB-profiler` in order to detect drug resistance-conferring
mutations and predict drug-resistant phenotypes. However you have only analyzed one sample unrelated to our study for the purpose of
learning. Before starting the transmission analysis, we would need to analyze the 20 samples that we
have been asked to analyze. The process to follow is exactly the same for any sample, meaning that
we would need to analyze the 20 samples following the same steps... twenty times??? Of course not.

The way you can perform the same analysis tasks on multiple samples at a a time in Galaxy (being it
a simple step like running Trimmomatic, or a complete analysis pipeline composed of many different steps)
is by using `Dataset collections`. To know more about dataset collections,
you can [have a look at the following material](https://galaxyproject.org/tutorials/collections/), however the idea is very simple:
all samples to be analyzed are merged into a "dataset collection", and then you analyze that "collection"
as you would analyze a single sample, with the difference that all steps will be performed for each
sample individually. Do not worry if you do not fully understand this right now, we will be creating
a dataset collection at the beginning of this tutorial and we will analyze it.   

To save you some time, (and server load) we have run the pipeline that you used in the previous tutorial for
the 20 samples that we have been asked to analyze. Thus, we now have 20 VCF files that describe the
mutations found for each of the samples. **These 20 VCFs files will be the starting point of this tutorial.**
 If you want to perform the mapping and variant calling for
all of the samples, feel free to do it. You can find the respective [FASTQ files here](https://zenodo.org/record/5911437), and
the [Galaxy workflow that was used to analyze the samples here](https://usegalaxy.eu/u/galo_a_goig/w/from-fastqs-to-vcfs-and-bams).
However this is completely optional, and we would suggest to do it after you have finished all the
tutorials of this workshop.

Before starting, bear in mind that this tutorial assumes that you watched the respective webinars of
this lesson (links) and therefore you understand 1) How genotypic drug susceptibility is determined
based on WGS analysis 2) The concept of clustering.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get the data
As mentioned in the introduction, we have performed mapping and variant calling for the 20 samples
that we need to analyze. The result are the respective 20 VCF files that describe the mutations found
for each of the samples. Before starting the analysis of such mutations, we will need to import them
into Galaxy:
> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
> Import the VCF files containing the variants of each sample
>    ```
> https://zenodo.org/record/6010176/files/ERR1203059.vcf
> https://zenodo.org/record/6010176/files/ERR181435.vcf
> https://zenodo.org/record/6010176/files/ERR2659153.vcf
> https://zenodo.org/record/6010176/files/ERR2704678.vcf
> https://zenodo.org/record/6010176/files/ERR2704679.vcf
> https://zenodo.org/record/6010176/files/ERR2704687.vcf
> https://zenodo.org/record/6010176/files/ERR313115.vcf
> https://zenodo.org/record/6010176/files/ERR551620.vcf
> https://zenodo.org/record/6010176/files/ERR5987300.vcf
> https://zenodo.org/record/6010176/files/ERR5987352.vcf
> https://zenodo.org/record/6010176/files/ERR6362078.vcf
> https://zenodo.org/record/6010176/files/ERR6362138.vcf
> https://zenodo.org/record/6010176/files/ERR6362139.vcf
> https://zenodo.org/record/6010176/files/ERR6362156.vcf
> https://zenodo.org/record/6010176/files/ERR6362253.vcf
> https://zenodo.org/record/6010176/files/ERR6362333.vcf
> https://zenodo.org/record/6010176/files/ERR6362484.vcf
> https://zenodo.org/record/6010176/files/ERR6362653.vcf
> https://zenodo.org/record/6010176/files/SRR13046689.vcf
> https://zenodo.org/record/6010176/files/SRR998584.vcf
>    ```
> We will also need the reference genome that was used for SNP calling
>    ```
> https://zenodo.org/record/3497110/files/MTB_ancestor_reference.fasta
>    ```
>
> Finally  Create a **Dataset List (Collection)** for all the VCFs.
>
> Use a meaningful name, for example **MTB VCFs**.
{: .hands_on}

> ### {% icon tip %} Tip
> To create a dataset collection:
> - Click on the box under the history name (bottom right) that says "Operations on multiple datasets"
> - Tick all 20 VCFs files that we just imported.
> - In the "For all selected..." box, select "Build Dataset List"
> - You can also [watch this one minute video](https://www.youtube.com/watch?v=F3qNs1_675g) about how to create collections
{: .tip}


# Generate a SNP alignment

In this tutorial we aim to calculate the genetic differences (in SNPs) between pairs of MTB genomes.
To do so, we need to compare, between each pair of genomes (thus pairwise), the nucleotides that are
observed at each position. Each time we find a different nucleotide at a given position, we will sum
1 SNP of genetic distance. For example, if two strains have a genetic distance of 5 SNPs, that would
mean that their respective genomes are almost identical, except for 5 positions along all the genome
in which they have different nucleotides.

To do such calculation we need to first build an alignment of all the genomes (multiple-sequence
  alignment, or MSA). Afterwards, we will use specific software to analyze this MSA, count SNPs,
  and thus calculate the genetic distance between each pair of samples.

![MSA of 5 MTB genomes](./images/MSA.png "A SNP is highlighted in a MSA of five MTB genomes")

## Generate complete genomes
The first step to generate the genomes MSA will be... to get the complete genomes of our samples!
In the [MTB Variant Analysis tutorial](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/tb-variant-analysis/tutorial.html) we have analyzed short-read high-throughput sequencing data (Illumina) to
obtain the respective VCF files that describe the mutations found in each of our samples, as compared
to the reference genome. We can now use these VCF files to build the complete genome of each of our
samples.

### Filter VCF files for epidemiological/phylogenetic investigation
Interpreting **mixed calls** or **indels** in phylogenetic/epidemiological applications can be very
complicated. That is the reason why we tipically use alignments that only contain **fixed SNPs**. (REFS?)

Thus, the first step in this tutorial will be to filter the VCFs so we are sure that
they only contain **fixed SNPs**. As it was introduced in *part X of lesson Y ofthe workshop (@Daniela)*, we will consider fixed
those variants at a frequency equal or greater than 90%. We will be using here the tool
**TB Variant Filter**

> ### {% icon hands_on %} Hands-on: Filter VCF files for epidemiological investigation
>
> 1. {% tool [TB Variant Filter](toolshed.g2.bx.psu.edu/repos/iuc/tb_variant_filter/tb_variant_filter/0.1.3+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"VCF file to be filter"*: `MTB VCFs` (Select Dataset Collection instead of Single Dataset)
>    - *"Filters to apply"*: `Only accept SNVs`, `Filter variants by percentage alt allele`
>    - *"Show options for the filters"*: `Yes`
>    - *"Minimum alternate allele percentage to accept"*: `90.0`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1.**`TB Variant Filter`** reads the VCF and output only SNVs that have, at leat, 90% frequency.
> How can this sofware extract such information from the VCF files?
>
> > ### {% icon solution %} Solution
> >
> > 1. That information is contained, for each mutation, in the VCF:
> > - The `TYPE` field within the INFO string will tell us if the mutation is a SNP (TYPE=snp)
> > - You can look for other types of mutations like insertions (TYPE=ins)
> > - The `AF` field within the INFO string describes the estimated **A**llele **F**requency
> > - An alterntive way to calculate it, would be to divide the number of observations of the
> > alternate allele (`AO`) by the total depth at that position (`DP`)
> {: .solution}
>
{: .question}

### Reconstruct the complete genome of each sample with **bcftools consensus**
Our VCF files now only contain **fixed SNPs** that were found in the genome of the respective strains.
Genomic positions not in the VCF mean that, at that particular position, the strain has the
same nucleotide than the reference genome. Knowing this information, one could reconstruct the
complete genome of each strain pretty easily. From the first to the last position in the genome,
one would put the same nucleotide than in the reference if that position is not in the VCF, or the
SNP described in the VCF otherwise. This is exactly what `bcftools consensus` will do for us,
given the reference genome and the VCF of the strain we want to reconstruct the genome for.

> ### {% icon hands_on %} Hands-on: Reconstruct the complete genome of each sample
>
> 1. {% tool [bcftools consensus](toolshed.g2.bx.psu.edu/repos/iuc/bcftools_consensus/bcftools_consensus/1.10+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"VCF/BCF Data"*: `MTB VCFs filtered` (output of **TB Variant Filter** {% icon tool %})
>    - *"Choose the source for the reference genome"*: `Use a genome from the history` (The reference genome that was used for SNP-calling)
>        - {% icon param-file %} *"Reference genome"*: `(MTB_ancestor_reference.fasta)`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. Imagine that we forgot to filter the VCFs to contain only fixed variants, and there are also
> SNPs with frequencies, of 15%, 30%, or 56.78%. Which allele do you think bcftools consensus would
> insert in the genome?
>
> > ### {% icon solution %} Solution
> > 1. The behaviour of bcftools consensus in this case can be specified with the option `--haplotype`
> > For example, we can set `haplotype=2` so the second allele will be used... wait... what?
> {: .solution}
{: .question}

> ### {% icon question %} Second allelle!?
>
> 1. What do you think that things like "second allele" or "*The* alterntive allele" mean here?
>
> > ### {% icon solution %} Solution
> >  1. Many of the bioinformatic programs are developed to analyze eukaryotic genomes, particularly
> >  human genomes. That means that these programs have in mind that the genomes
> >  are diploid and thus each posible position in the genome has two possible alleles. In
> >  bacterial genomics, in contrast, we are **always sequencing a population** of cells with
> >  potential genetic diversity (with the exception of single-cell sequencing).
> >  That does not mean that we cannot use this type of software, we can (and we do!) but it is
> >  good to know what they are ment for, and their possible limitations
> {: .solution}
{: .question}

## Multiple-sequence alignment (MSA) of all genomes
Multiple sequence alignment is a process in which multiple DNA, RNA or protein sequences are arranged
to find regions of similarity that are supposed to reflect the consequences of different evolutionary
processes (Figure 1). MSAs are used to test hypotheses about this evoluitionary processes and infer phylogenetic
relationships, and for these reasons we build MSA for sequences for which we already assume some sort
of evolutionary relationship. More on MSAs here (link here to a webpage or Christoph's tutorial.).
**Galo: maybe this is too much, or too complex explanation or not necessary here**

We can build MSAs of a particular gene, a particular sequence or the complete genome. Building MSAs
of several complete genomes can be a complicated process and computationally demanding. To perform
such task there are many software packages available like `Muscle`, `MAFFT` or `Clustal` just to mention some.
Which one are we going to use in this tutorial? Well, we are going to use a trick. We are going to
just stack one genome on top of each other within a text file. (More on why we can do this below).

Our aim is to generate a **multifasta** file in which the genomes of our samples are aligned.
Something that looks like this:

```
>Sample 1
TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGC...
>Sample 2
ATGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGC...
>Sample 3
TTCGTTCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACG...
```


But given the output of **bcftools consensus**, what we have **for each** sample looks like this:
```
>MTB_anc
TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGC...
```

That is to say, because we had no option to specify the name of the sample, all our samples are now
identified in the fasta header as `">MTB_anc"` (the name of the reference genome that was used as
  template by `bcftools`), so we now have no way of distinguish them!

Performing text transformations is very common when doing
bioinformatics. Sometimes the output of one software/step is not completely compatible with the
downstream analysis, so we have to transform it to make it compatible. In this case, we have
to find a way of renaming the genomes, so they contain the name of the respective sample as in the example box
 above. There is no "correct" way of doing this, and many possible strategies/tools could be used. We
 could do it by hand, for example.
 Because we are using Galaxy, we will use tools available in Galaxy to perform such task.

In Galaxy, because we are
analzying a collection of 20 samples instead of one sample at a time (meaning that renaming manually
each sample is not an option, and never should be!) we can do the renaming in two steps:

First we will use **Add input name as column**. This will add a colum with the name of the input
file that we will be able to use as identifier of the sample. For example, after using it, for the
sample ERR6362078, the file will look like this:

```
>MTB_anc  ERR6362078.vcf
TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTT  ERR6362078.vcf
AACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTG  ERR6362078.vcf
ACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTT  ERR6362078.vcf
GCTCTGTTATCCGTGCCGAGCAGCTTTGTCCAAAACGAAATCGAGCGCCATCTGCGGGCC  ERR6362078.vcf
CCGATTACCGACGCTCTCAGCCGCCGACTCGGACATCAGATCCAACTCGGGGTCCGCATC  ERR6362078.vcf
```

Now we want to remove "MTB_anc" from the fasta header, so instead of ">MTB_anc", it reads ">ERR6362078.vcf".
We  also want to remove "ERR6362078.vcf" from the rest of the lines, as they should only contain the genome
sequence, so it loos like:

```
>ERR6362078.vcf
TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTT
AACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTG...
```

We can easily do this using **Text reformatting with AWK**. AWK is a programming language integrated in the linux
terminal that is very useful and powerful to manipulate text. **Don't worry about this now**. You are
not supposed to know **awk** language at this stage. However if your aim is to become a fully dedicated
bioinformatician, you will need to learn the basics of the linux terminal, and some of its tools like
`AWK`, `sed` or `grep`. You can find plenty of information and tutorial like
[this one.](https://riptutorial.com/awk) in the internet.

### **Add input name as column**

> ### {% icon hands_on %} Hands-on: Add file name as a column
>
> 1. {% tool [Add input name as column](toolshed.g2.bx.psu.edu/repos/mvdbeek/add_input_name_as_column/addName/0.2.0) %} with the following parameters:
>    - {% icon param-file %} *"to Dataset Collection"*: `Consensus fasta` (output of **bcftools consensus** {% icon tool %})
>    - *"input contains a header line?"*: `No`
>
>    > ### {% icon comment %} Comment
>    >
>    > We need to specify that input does **not** contain a header so the name is also added to the
>    > first line (the fasta header). Otherwise this tool will only add the name starting from the second line.
>    {: .comment}
>
{: .hands_on}

### **Text reformatting with AWK**

> ### {% icon hands_on %} Hands-on: Reformat text to generate a correct fasta file
>
> 1. {% tool [Text reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Dataset collection` (output of **Add input name as column** {% icon tool %})
>    - *"AWK Program"*: copy-paste the following text: `$1 == ">MTB_anc" { print ">"$2 } $1 != ">MTB_anc"  { print $1 }`
>
{: .hands_on}


> ### {% icon question %} Questions
>
> 1. As a game in which you have to decipher a rebus puzzle: can you decipher what the AWK command used above means?
>
> > ### {% icon solution %} Solution
> >
> > AWK reads files line by line, and performs a set of actions given a set of conditions that we
> > programmed. In our case, for each line:
> >
> > * `$1 == >MTB_anc { print >$2 }` This part reads: **if** the first field of the line ($1) is **equal to** ">MTB_anc", then print ">" followed by the second field ($2) of the line (ERR6362078).
> > * `$1 != >MTB_anc  { print $1 }` This part reads: **if** the first field of the line ($1) is **NOT equal to**
> > ">MTB_anc", then print (only) the first field ($1) of the line (the sequence).
> {: .solution}
>
{: .question}

### Build a multiple-sequence alignment from complete genomes with "**Concatenate datasets**"
We have now the genome of each sample with a proper identifier in the fasta header. As we already
mentioned, we are going to build a MSA by stacking one genome on top of each other in a single text
file. **A multifasta file**.

```
>Sample 1
TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGC...
>Sample 2
ATGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGC...
>Sample 3
TTCGTTCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACG...
```

This could be done manually, by copy-pasting all genomes in a single text file.
However we can do the same with a specific command that *concatenates* files.

> ### {% icon hands_on %} Hands-on: Concatenate genomes to build a MSA
>
> 1. {% tool [Concatenate datasets tail-to-head (cat)](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: `Dataset collection` (output of **Text reformatting** with awk {% icon tool %})
> 2. The output of **concatenate datasets** may be of type tabular. Make sure Galaxy sees this file as
> a fasta file by editing its attributes. Click in the pencil icon, select "Datatypes" and then select fasta.
>
{: .hands_on}

Now we have a multifasta file, where each position of each genome corresponds to the same position
of the rest of genomes in the file. This can be seen and used as a multiple-sequence alignment of
all of our genomes! However, it is important that you understand the following question...

> ### {% icon question %} Question
>
> 1. Generating multiple-sequence alignments can be complicated and computationally demanding, and there are
> many software packages to perform such task. How is then possible that we were able to build a MSA by just
> stacking genomes one on top of each other? Can you think about what makes our case special, so we can
> just use this "trick"?
>
> > ### {% icon solution %} Solution
> >
> > We have generated the complete genome of each sample by substituting in the reference genome
> > those **SNPs** that we found in that said sample (described in the VCF). Remember that **we
> > removed indels from the VCF when filtering!** Because the complete genomes we generated do not
> > contain insertions or deletions, ALL the genomes have the same length (the length of the
> > reference genome) and each nucleotide corresponds to the same genomic coordinate (the one
> > also in the reference genome). So we are not aligning genomes *per se* but, knowing this, we
> > can build a MSA by just stacking genomes that **are the same length >>AND<< have the same coordinates**.
> {: .solution}
>
{: .question}

### Remove invariant positions with **Finds SNP sites**
We have generated a MSA that is the basis for the transmission (clustering) and phylogenetic
analysis. Although we could already use this MSA for such analysis, it is common practice to remove
the invariant sites from the alignment. Think that our file now contains 20 genomes of 4.4 Mb each.
MTB genomes have a very low genetic diversity, meaning that in reality, there are only some hundreds or
thousands of SNPs *in total*, because the genomes are >99% identical between them. Identical positions
in a MSA provide no information **(link here about invariant positions?)**, so we can remove those and
generate a **SNP alignment** that only contain variant positions with phylogenetic information. By
doing this, we will generate a much smaller file, that will be easier to handle by downstream applications.

We can exemplify this with a couple of pictures:

In the following picture, a polymorphic position in the alignment (SNP)
is highlighted in blue and green. Invariant positions are not highlighted and have an asterisc `*`
on the upper part of the alignment. Most of the MSA is composed of these invariant positions, making
our file larger than necessary.

![MSA of 5 MTB genomes](./images/MSA.png "A SNP is highlighted in a MSA of five MTB genomes")

After removing invariant positions, we end up with a SNP alignment like the following.

![MSA of 5 MTB genomes](./images/SNP_MSA.png "A SNP alignment where all positions are polymorphic")

> ### {% icon hands_on %} Hands-on: Removing invariant sites from a MSA
>
> 1. {% tool [Finds SNP sites](toolshed.g2.bx.psu.edu/repos/iuc/snp_sites/snp_sites/2.5.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"FASTA file"*: `Single dataset` (output of **Concatenate datasets** {% icon tool %})
>    - *"Output"*: `Sequence alignment / VCF`
>        - *"Output formats"*: `Multi-FASTA alignment file`
>
{: .hands_on}

In SNP alignments you have to **bear in mind that positions do not longer correspond to genomic
coordinates**, meaning that two contiguous nucleotides may correspond to coordinates thousands of
positions apart.

# Identify transmission clusters

## Calculate pairwise SNP distances
Now we are all set to calculate pairwise SNP distances between samples and decide whether two
patients are within the same transmission cluster or not. Having a SNP alignment, this is fairly
easy. We will use **SNP distance matrix**, that will generate a matrix with pairwise SNP distances.

> ### {% icon hands_on %} Hands-on: Distance matrix from SNP alignment.
>
> 1. {% tool [SNP distance matrix](toolshed.g2.bx.psu.edu/repos/iuc/snp_dists/snp_dists/0.6.3+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"FASTA multiple sequence alignment"*: `Single dataset`
> (output of **Finds SNP sites** {% icon tool %})
>
>    > ### {% icon comment %} Comment
>    >
>    > Have a look at the distance matrix to make sure you understand the whole process. Given that
>    > we only have 20 samples, you could already spot some samples that are involved in the same
>    > transmission chain (samples with a small number of SNPs between them as explained below).
>    {: .comment}
>
{: .hands_on}

## Determine transmission clusters based on a SNP threshold
Now that we have a distance matrix that describes the SNP distance between each pair of samples, we
could already describe the transmission clusters based on a **SNP threshold**, as explained in the
respective webinar(link). If two samples are at a distance below that threshold, we will say that they
belong to the same transmission cluster, because they are close enough genetically speaking. Again,
we could do this manually, but we are doing bioinformatics, and we want to be able to do the same
analysis regardless of whether we are analyzing two or two million samples. Also, **note that two samples
that are dozens of SNPs apart may belong to the same transmission cluster** if there are other samples
linking them in between as exemplified in the picture below.

![MSA of 5 MTB genomes](./images/Cluster.png "Example of transmission cluster using 10 SNPs threshold")



> ### {% icon question %} Very Important Question
>
> 1. In the image above exemplifying a transmission cluster, the distance between samples A and E is
> 17 SNPs. Being the other pairwise distances in the figure same,
> would it be possible that the distance between A and E is different?
>
> > ### {% icon solution %} Solution
> >
> > 1. The figure used above as an example is a **flagrant oversimplification**. In the figure not all pairwise
> > distances are represented (for example between sample A and C).
> > Most importantly you have to remember that **transmission clusters do not reflect transmission
> > events**.
> > In fact it may happen that transmission does not happen within the cluster we are analyzing!
> > For example, when a patient that is not sampled is the source of infection of all the cases in
> > the cluster (may be a superspreading event). You have to consider that, taking into account the
> > same SNP distances, another possible (yet still oversimplified) scenario could be...
> > ![MSA of 5 MTB genomes](./images/Cluster2.png "Example of transmission cluster using 10 SNPs threshold")
> {: .solution}
{: .question}




## Determine transmission clusters using Rscript
Currently there is not tool in Galaxy to perform the exact task that we need, but we can use the
statistics package `R` and the library `cluster` to perform such task. Again, don't worry about this,
programming in `R` is beyong the scope of this workshop. We have prepared a simple `R` script that
identifies transmission clusters based on the matrix of SNP distances. This script is available
within the Galaxy workflow we are following, so feel free to use it for your own analysis!

> ### {% icon question %} Questions
>
> 1. How many transmission clusters did we find? How many samples are linked to recent transmission in our dataset?
>
> > ### {% icon solution %} Answer
> > 1. We have found two transmission clusters with respective IDs 10 and 12. Transmission cluster 10
> > is composed by two samples linked by recent transmission and transmission cluster 12 by three samples
> > linked by recent transmission. For example samples ERR6362484 and ERR5987352 are linked by
> > recent transmission.
> {: .solution}
{: .question}

The output of the R script is a table containing, for samples that were found within a cluster,
their respective names and the cluster id (an arbitrary number) they belong to:

| Sample                  | cluster_id |
|-------------------------|------------|
| ERR6362484.vcf    | 10         |
| ERR6362138.vcf | 12         |
| ERR6362156.vcf | 12         |
| ERR6362253.vcf | 12         |
| ERR5987352.vcf    | 10         |


> ### {% icon question %} Question
>
> 1. Let's assume that we have the isolation dates of samples ERR6362484 and ERR5987352, which
> belong to the same transmission cluster. Sample ERR6362484 was isolated on January 2021, while sample
>  ERR5987352 was isolated on June 2021. Would you be able to determine who was the infector and who the infectee?
>
> > ### {% icon solution %} Solution
> > 1. **NO**
> >
> >  Isolation dates have been used traditionally to define **index cases** within transmission clusters
> > under the assumption that the most likely scenario is the first isolated sample to be the
> > source of transmission. Today we know that this assumption often leads
> > to misidentification of index cases. **Remember:** we cannot rule out the possibility
> > that patients within the cluster were infected by an index case that was not sampled.
> >
> >  Read [Xu et al., 2019](https://doi.org/10.1371/journal.pmed.1002961) for more information on this topic.  
> {: .solution}
{: .question}

# Using clustering to investigate the emergence of drug resistance

 Although we have stressed the fact that **clustering cannot be used to delinate
transmission events**, clustering is very useful to investigate outbreaks and determine which cases are
involved in the same transmission chain. We can leverage this information to investigate the relationship
between tuberculosis transmission and particular biological or clinical traits.

In this part of the tutorial, we will investigate the emergence and spread of drug resistance based
on our clustering analysis.

## Get the data

In the [MTB Variant Analysis tutorial](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/tb-variant-analysis/tutorial.html)
you have used **TB-profiler** to generate a report with determinants of drug resistance of a
particular MTB strain, and predict its genotypic drug susceptibility. We have done **exactly the same**
for the 20 samples that we used in the clustering analysis, so we have now the TB-profiler report for
all of them.

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
> Import the TB profiler report of each sample
>    ```
> https://zenodo.org/record/6010176/files/ERR1203059.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR181435.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR2659153.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR2704678.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR2704679.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR2704687.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR313115.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR551620.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR5987300.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR5987352.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362078.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362138.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362139.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362156.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362253.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362333.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362484.TBprof.txt
> https://zenodo.org/record/6010176/files/ERR6362653.TBprof.txt
> https://zenodo.org/record/6010176/files/SRR13046689.TBprof.txt
> https://zenodo.org/record/6010176/files/SRR998584.TBprof.txt
>    ```
> Create a **Dataset List (Collection)** for all the report files.
>
> Use a meaningful name, for example **TBprofiler reports**.
{: .hands_on}

## Summarize the data
TB-profiler reports are very useful and comprehensive, and we will use them to better investigate
drug resistance in our dataset. However, it is always useful to summarize the data on a per-sample basis,
on a table, so we can quickly check which strains are, for example, MDR, an which are pan-susceptible.
We would like to generate a table like the following:

| Sample                  | DR profile |
|-------------------------|------------|
| Sample A |    Sensitive      |
| Sample B | Sensitive         |
| Sample C | MDR        |
| Sample Z | XDR        |

If we have a look at a TB profiler report we can see that there is one line describing the genotypic drug susceptibility.
Let's have a look at the first part of the TB-profiler report for sample ERR6362653:

```

TBProfiler report
=================

The following report has been generated by TBProfiler.

Summary
-------
ID: tbprofiler
Date: Fri Jan 28 13:14:47 2022
Strain: lineage2.2.1
Drug-resistance: MDR

Lineage report
--------------
Lineage	Estimated Fraction	Family	Spoligotype	Rd
lineage2	1.000	East-Asian	Beijing	RD105
lineage2.2	0.996	East-Asian (Beijing)	Beijing-RD207	RD105;RD207
lineage2.2.1	0.999	East-Asian (Beijing)	Beijing-RD181	RD105;RD207;RD181
```

As you can see, this strain is multi-drug resistant as inficated by (`Drug-resistance: MDR`)
We could then look for this information in each TB-profiler report and generate this table manually,
for example in a spreadsheet. However this is not feasible when analyzing hundreds or thousands of samples
(and very **error-prone!**).

We are here to learn bioinformatics, so let's generate this table using Linux commands.

The process will consist on three steps (of which you already know two of them):
* 1) Select the line containing the drug resistance profile with **grep**:

```
Drug-resistance: MDR
```

* 2) Prepend the name of the sample with **Add input name as column**:

```
ERR6362653.txt Drug-resistance: MDR
```

* 3) Concatenate results from all samples in a single file with **Concatenate datasets**:

```
ERR6362653.txt Drug-resistance: MDR
ERR313115.txt Drug-resistance: Sensitive
ERR5987300.txt Drug-resistance: Pre-XDR
.... etc
```

* 4) As an optional step, we can reformat the table with **sed** to get rid of the `.txt` and `Drug-resistance:`
so the table looks like:

```
ERR6362653 MDR
ERR313115 Sensitive
ERR5987300 Pre-XDR
```

#### 1. Select the line containing the drug resistance profile with `grep`

**grep** is used to search patterns of text within text files. Each time **grep** finds that
pattern, it will print as a result **the complete line** containing such pattern.

If we use **grep** to search for the pattern "`Drug-resistance`", in a TB-profiler file, we will get
as output the complete line, for example `Drug-resistance: MDR`

#### *Search in textfiles (grep)*

> ### {% icon hands_on %} Search for `Drug-resistance` in TB-profiler files
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/1.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `Dataset collection` (Collection of the TB profiler reports we just imported) {% icon tool %})
>    - *"Regular Expression"*: `Drug-Resistance`
{: .hands_on}


#### 2. Prepend the sample name
We have already used this command before. This time we will *prepend* the column with the sample name
so it appears as the first column. This is arbitrary and just a matter of personal taste:

#### *Add input name as column*

> ### {% icon hands_on %} Prepend the sample name to the DR profile
>
> 1. {% tool [Add input name as column](toolshed.g2.bx.psu.edu/repos/mvdbeek/add_input_name_as_column/addName/0.2.0) %} with the following parameters:
>    - {% icon param-file %} *"to Dataset"*: `Dataset collection` (output of **Search in textfiles** {% icon tool %})
>    - *"input contains a header line?"*: `No`
>    - *"Prepend the colum"*: `Yes`
>
{: .hands_on}

#### 3. Concatenate results

#### *Concatenate datasets*

> ### {% icon hands_on %}
>
> 1. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.1) %} with the following parameters:
>    - {% icon param-file %} *"Datasets to concatenate"*: `Dataset collection` (output of **Add input name as column** {% icon tool %})
>
{: .hands_on}

#### 4. Cleanup the table (optional)
In this step we will use a simple tool that searches and replaces text. We want to remove the ".txt"
at the end of sample names, and the string "Drug-resistance:". So we will tell the tool to search
for these two *patterns* and to replace them with "*nothing*"

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/1.1.2) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Single file` (output of **Concatenate datasets** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `.txt`
>            - *"Replace with:" (leave this in blank)*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `Drug-resistance:`
>           - *"Replace with:"(leave this in blank)*:
>
{: .hands_on}


> ### {% icon question %} Question
>
> 1. How many MDR strains did we find in the dataset?
> 2. What does it mean to be Pre-MDR?
>
> > ### {% icon solution %} Solution
> > 1. Eight MDR strains, three of which are pre-XDR because they have additional resistance to fluoroquinolones. (You can look into details by looking into the TB profiler reports).
> > 2. As MDR means to be resistant to INH and RIF, pre-MDR means to be either INH-monoresistant or RIF-monoresistant.
> > If we have a look at the respective TB-profiler reports, we can see that these three strains are RIF-monoresistant.
> {: .solution}
{: .question}

# Put everything together
Now that we have performed a clustering analysis and know which DR mutations carry each strain,
let's try answer a series of questions about how DR may be emergind and spreading in our study
population.
You will see that we will be supporting our findings in the results of our analysis, but general
knowledge in the TB field will also help us to conduct our investigation!

| Sample       | Cluster_id | DR profile | Clustering  |
|--------------|------------|------------|-------------|
| ERR1203059   | -          | Sensitive  | Unclustered |
| ERR181435    | -          | Sensitive  | Unclustered |
| ERR2659153   | -          | Sensitive  | Unclustered |
| ERR2704678   | -          | Sensitive  | Unclustered |
| ERR2704679   | -          | Sensitive  | Unclustered |
| ERR2704687   | -          | Sensitive  | Unclustered |
| ERR313115    | -          | Sensitive  | Unclustered |
| ERR551620    | -          | MDR        | Unclustered |
| ERR5987300   | -          | Pre-XDR    | Unclustered |
| ERR6362078   | -          | MDR        | Unclustered |
| ERR6362139   | -          | Pre-MDR    | Unclustered |
| ERR6362333   | -          | Pre-XDR    | Unclustered |
| ERR6362653   | -          | MDR        | Unclustered |
| SRR13046689  | -          | Other      | Unclustered |
| SRR998584    | -          | Sensitive  | Unclustered |
| ERR5987352   | 10         | Pre-MDR    | Clustered   |
| ERR6362484   | 10         | Pre-MDR    | Clustered   |
| ERR6362138   | 12         | MDR        | Clustered   |
| ERR6362156   | 12         | Pre-XDR    | Clustered   |
| ERR6362253   | 12         | MDR        | Clustered   |



> ### {% icon question %} Question
>
> Assuming that we have a very good sampling of the outbreak. Which strains **may** represent instances
> of *de novo* evolution of drug resistance and which ones instances of *transmitted* (primary) resistance?
> Remember that you can look at the TB-profiler reports of independent samples for detailed information.
>
> > ### {% icon solution %} Solution
> > In a simplistic scenario, we could consider clustered strains as instances of transmission and
> > unclustered strains as instances of de novo evolution of DR. Thus, we see that for example there
> > are three MDR strains (ERR551620, ERR6362078, ERR6362653) that are unclustered and therefore may
> > represent cases in which drug resistance evolved independently as response to treatment within the respective
> > patients. However you need to always bear in mind that, although **within our population** those MDR
> > strains doesn't seem to be linked to transmission, this does NOT rule out the possibility that
> > some of these patients were infected with an MDR strain somewhere else.
> >
> > When looking at clustered strains, distinguishing between transmitted and *de-novo* may be tricky. Note
> > that, for example, for the two RIF-monoresistant strains linked within the same transmission
> > cluster there are, at least, a couple of scenarios possibe: one in which a RIF-monoresistant strain
> > evolved in one patient and was transmitted to the other patient afterwards, and one in which both
> > were infected with the same RIF-monoresistant strain from a third patient that we have not sampled.
> > We need to note that, in the first scenario, drug resistance evolved *de novo* in one patient,
> > and was *transmitted* to the other patient, whereas in the second scenario drug resistance was
> > *transmitted* in both cases.
> {: .solution}
{: .question}



> ### {% icon question %} Question
>
> 1. The same principles than those explained above apply to the three MDR strains that are
> within the same transmission cluster. However in this case there is one strain that shows clear
> evidence of *de-novo* evolution of DR. Do you know which strain and why?
> 2. Are there possible scenarios other than *de-novo* evolution of DR for this strain?
>
> > ### {% icon solution %} Solution
> > 1. Within this cluster of MDR strains, there is one tagged as Pre-XDR by TB-profiler. If we have
> > a look at the TB profiler report, we can see that this strain carries an additional mutation in
> > *gyrA* that confers resistance to fluorioquinolones. This is compatible with an scenario in which
> > fluoroquinolone resistance evolved independently within this patient after being infected with
> > the MDR strain.
> > 2. Remember, although it is reasonable to think that fluorioquinolone resistance evolved in this
> > patient, we cannot rule out that actually it evolved in another patient who we did not sample and
> > was the transmitter of the fluoroquinolone-resistant strain.
> {: .solution}
{: .question}


> ### {% icon question %} Question
>
> 1.  There is one strain with a DR profile "other", because it is only resistant to pyrazinamide. This
> strain is not within a transmission cluster. Therefore, we conclude that pyrazinamide resistance
> most likely evolved *de-novo* in this patient. But we are wrong. Do you know why?  
>
> > ### {% icon solution %} Solution
> > 1. The strain is indeed PZA-resistant. And indeed this is strain is NOT linked to transmission
> > within our population. However, if we have a look at the TB-profiler report, we observe that this
> > is a *M. bovis* strain, wich are known to be intrinsically resistant to PZA (all carry the same
> > mutation conferring resistance to PZA).
> {: .solution}
{: .question}

> ### {% icon question %} Question
>
> 1. Is it possible to find in the same transmission cluster two RIF-monoresistant strains that
> carry different rpoB mutations?
> 2. Is it possible to find in the same transmission cluster strains of different MTB sublineages?
>
> > ### {% icon solution %} Solution
> > 1. Yes, it is **possible**. In that scenario, both patiens were **recently**
> > transmitted with the **same susceptible strain**, and RIF resistance evolved **independently** in both.
> > 2. No, by definition. Remember that clustering is based on a threshold that we set of genetic
> > distance measured in SNPs. We want to cluster samples that are genetically so similar that we
> > can consider them as the same genotype, that is to say, as the same strain. Two different
> > sublineages, by definition, do not belong to the same genotype and will have a distance in SNPs
> > between well beyond any SNP threshold we could use.
> {: .solution}
{: .question}


# Conclusion
{:.no_toc}

You have learned how to perform a clustering analysis to identify patients that are linked by events of
**recent transmission**. Clustering analysis is very useful in outbreak investigation and
can also be used to describe the emergence and spread of drug-resistance within a population. You have
also learned, however, that interpreting clustering results requires careful considerations, given the
limitations of the methodology. Clustering analysis is better complemented with phylogenetic analysis,
which may help overcome some of these limitations.

In the following tutorial [this will be a link] you will perform a phylogenetic analysis of these
same 20 strains.
