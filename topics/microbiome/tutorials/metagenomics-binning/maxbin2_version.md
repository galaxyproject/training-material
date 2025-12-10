## Bin contigs using MaxBin2

**MaxBin2** ({% cite maxbin2015 %}) is an automated metagenomic binning tool that uses an Expectation-Maximization algorithm to group contigs into genome bins based on abundance, tetranucleotide frequency, and single-copy marker genes.

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
We can now launch the proper binning with MaxBin2
> <hands-on-title> Individual binning of short-reads with MaxBin2 </hands-on-title>
>
> 1. {% tool [MaxBin2](toolshed.g2.bx.psu.edu/repos/mbernt/maxbin2/maxbin2/2.2.7+galaxy6) %} with the following parameters:
>    - {% icon param-collection %} *"Contig file"*: `Contigs` (Input dataset collection)
>    - *"Assembly type used to generate contig(s)"*: `Assembly of sample(s) one by one (individual assembly)`
>        - *"Input type"*: `Abundances`
>            - {% icon param-file %} *"Abundance file"*: `outputDepth` (output of **Calculate contig depths** {% icon tool %})
>    - In *"Outputs"*:
>        - *"Generate visualization of the marker gene presence numbers"*: `Yes`
>        - *"Output marker gene presence for bins table"*: `Yes`
>        - *"Output marker genes for each bin as fasta"*: `Yes`
>        - *"Output log"*: `Yes`
>
>
{: .hands_on}

> <question-title>Binning metrics</question-title>
>
> 1. How many bins where produced by MaxBin2 for our sample?
> 2. How many contigs are in the bin with most contigs? 
> > <solution-title></solution-title>
> >
> > 1. There are two bin for this sample.
> > 2. 35 and 24 in the other bin. 
> >
> {: .solution}
>
{: .question}