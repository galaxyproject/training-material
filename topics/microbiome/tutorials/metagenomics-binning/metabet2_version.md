## Bin contigs using MetaBAT 2

**MetaBAT** stands for "Metagenome Binning based on Abundance and Tetranucleotide frequency". It is:

> Grouping large fragments assembled from shotgun metagenomic sequences to deconvolute complex microbial communities, or metagenome binning, enables the study of individual organisms and their interactions. Here we developed automated metagenome binning software, called MetaBAT, which integrates empirical probabilistic distances of genome abundance and tetranucleotide frequency. On synthetic datasets MetaBAT on average achieves 98percent precision and 90% recall at the strain level with 281 near complete unique genomes. Applying MetaBAT to a human gut microbiome data set we recovered 176 genome bins with 92% precision and 80% recall. Further analyses suggest MetaBAT is able to recover genome fragments missed in reference genomes up to 19%, while 53 genome bins are novel. In summary, we believe MetaBAT is a powerful tool to facilitate comprehensive understanding of complex microbial communities.
{: .quote author="Kang et al, 2019" }

MetaBAT is a popular software tool for metagenomics binning, and there are several reasons why it is often used:
- *High accuracy*: MetaBAT uses a combination of tetranucleotide frequency, coverage depth, and read linkage information to bin contigs, which has been shown to be highly accurate and efficient.
- *Easy to use*: MetaBAT has a user-friendly interface and can be run on a standard desktop computer, making it accessible to a wide range of researchers with varying levels of computational expertise.
- *Flexibility*: MetaBAT can be used with a variety of sequencing technologies, including Illumina, PacBio, and Nanopore, and can be applied to both microbial and viral metagenomes.
- *Scalability*: MetaBAT can handle large-scale datasets, and its performance has been shown to improve with increasing sequencing depth.
- *Compatibility*: MetaBAT outputs MAGs in standard formats that can be easily integrated into downstream analyses and tools, such as taxonomic annotation and functional prediction.

The first step when using tools like MetaBAT or MaxBin2 is to compute contig depths from the raw alignment data. Both tools require per-contig depth tables as input, as their binning algorithms rely on summarized coverage statistics at the contig level. However, standard BAM files store read-level alignment information, which must first be processed to generate the necessary contig-level coverage data. This preprocessing step ensures compatibility with the input requirements of these binning tools.
> <hands-on-title> Calculate contig depths </hands-on-title>
>
> 1. {% tool [Calculate contig depths](toolshed.g2.bx.psu.edu/repos/iuc/metabat2_jgi_summarize_bam_contig_depths/metabat2_jgi_summarize_bam_contig_depths/2.17+galaxy0) %} with the following parameters:
>    - *"Mode to process BAM files"*: `One by one`
>        - {% icon param-file %} *"Sorted bam files"*: output of **Samtools sort** {% icon tool %}
>        - *"Select a reference genome?"*: `No`
>
>
{: .hands_on}
We can now launch the proper binning with MetaBAT 2
> <hands-on-title>Individual binning of short-reads with MetaBAT 2</hands-on-title>
> 1.  {% tool [MetaBAT 2](toolshed.g2.bx.psu.edu/repos/iuc/metabat2/metabat2/2.17+galaxy0) %} with parameters:
>     - *"Fasta file containing contigs"*: `Contigs`
>     - In **Advanced options**, keep all as **default**.
>     - In **Output options:**
>       - *"Save cluster memberships as a matrix format?"*: `"Yes"`
>
{: .hands_on}

MetaBAT 2 generates several output files during its execution, some of which are optional and only produced when explicitly requested by the user. These files include:

1. The final set of genome bins in FASTA format (`.fa`)
2. A summary file with information on each genome bin, including its length, completeness, contamination, and taxonomy classification (`.txt`)
3. A file with the mapping results showing how each contig was assigned to a genome bin (`.bam`)
4. A file containing the abundance estimation of each genome bin (`.txt`)
5. A file with the coverage profile of each genome bin (`.txt`)
6. A file containing the nucleotide composition of each genome bin (`.txt`)
7. A file with the predicted gene sequences of each genome bin (`.faa`)

These output files can be further analyzed and used for downstream applications such as functional annotation, comparative genomics, and phylogenetic analysis.

> <question-title>Binning metrics</question-title>
>
> 1. How many bins where produced by MetaBAT 2 for our sample?
> 2. How many contigs are in the bin with most contigs? 
> > <solution-title></solution-title>
> >
> > 1. There is only one bin for this sample.
> > 2. 52 (these numbers may differ slightly depending on the version of MetaBAT2). So not all contigs were binned into this bin!
> >
> {: .solution}
>
{: .question}