## Bin contigs using CONCOCT

**CONCOCT** ({% cite Alneberg2014 %}) is an *unsupervised metagenomic binner* that groups contigs using both **sequence characteristics** and **differential coverage across multiple samples**. In contrast to SemiBin, it does **not** rely on pretrained models or marker-gene constraints; instead, it clusters contig fragments purely based on statistical similarities.

> CONCOCT jointly models contig abundance profiles from multiple samples using a Gaussian mixture model. By taking advantage of differences in coverage across samples, it can separate genomes with similar sequence composition but distinct abundance patterns. CONCOCT also introduced the now-standard technique of splitting contigs into fixed-length fragments, allowing more consistent and accurate clustering.
> {: .quote author="Alneberg et al., 2014" }

CONCOCT is widely used in metagenomic binning due to:

* **Unsupervised probabilistic clustering**: No marker genes, labels, or pretrained models are required.
* **Strong performance with multiple samples**: Differential coverage helps disentangle closely related genomes
* **Reproducible, transparent workflow**: Its stepwise pipeline—fragmentation, coverage estimation, clustering—yields interpretable results.
* **Complementarity to other binners**: Frequently used alongside SemiBin, MetaBAT2, or MaxBin2 in ensemble pipelines (e.g., MetaWRAP, nf-core/mag).


Before initiating the binning process with **CONCOCT**, the input data must be preprocessed to ensure compatibility with its Gaussian mixture model. This model treats each contig fragment as an individual data point, necessitating a critical preprocessing step: dividing contigs into **equal-sized fragments**, usually around 10 kb in length.
Fragmentation serves several essential purposes:

- **Balancing Influence**: It mitigates bias between long and short contigs, ensuring each contributes equally to the analysis.
- **Uniform Data Points**: It creates consistent data points, which are crucial for accurate statistical modeling.
- **Detecting Local Variations**: It helps identify potential misassemblies or variations within long contigs.
- **Enhanced Resolution**: It improves the detection of abundance differences across genomes, leading to more precise binning results.

This fragmentation step is mandatory for CONCOCT to operate effectively and deliver reliable results.

> <hands-on-title> Fragment contigs </hands-on-title>
>
> 1. {% tool [CONCOCT: Cut up contigs](toolshed.g2.bx.psu.edu/repos/iuc/concoct_cut_up_fasta/concoct_cut_up_fasta/1.1.0+galaxy2) %} with the following parameters:
>    * {% icon param-collection %} *"Fasta contigs file"*: `Contigs` (Input dataset collection)
>    * *"Concatenate final part to last contig?"*: `Yes`
>    * *"Output bed file with exact regions of the original contigs corresponding to the newly created contigs?"*: `Yes`
>
{: .hands_on}

After fragmentation, CONCOCT calculates the **coverage of each fragment** across all samples. This coverage data serves as a measure of abundance, which, combined with sequence composition statistics, forms the input for the tool's analysis.

> <hands-on-title> Generate coverage table </hands-on-title>
>
> 1. {% tool [CONCOCT: Generate the input coverage table](toolshed.g2.bx.psu.edu/repos/iuc/concoct_coverage_table/concoct_coverage_table/1.1.0+galaxy2) %} with the following parameters:
>    * {% icon param-file %} *"Contigs BEDFile"*: `output_bed` (output of **CONCOCT: Cut up contigs** {% icon tool %})
>    * *"Type of assembly used to generate the contigs"*: `Individual assembly: 1 run per BAM file`
>      * {% icon param-file %} *"Sorted BAM file"*: `output1` (output of **Samtools sort** {% icon tool %})
{: .hands_on}

These **coverage profiles**, along with fundamental sequence features, are fed into CONCOCT's **Gaussian mixture model clustering algorithm**. This model groups fragments into distinct bins based on their statistical similarities. 

> <hands-on-title> Run CONCOCT </hands-on-title>
>
> 1. {% tool [CONCOCT](toolshed.g2.bx.psu.edu/repos/iuc/concoct/concoct/1.1.0+galaxy2) %} with the following parameters:
>    * {% icon param-file %} *"Coverage file"*: `output` (output of **CONCOCT: Generate the input coverage table** {% icon tool %})
>    * {% icon param-file %} *"Composition file with sequences"*: `output_fasta` (output of **CONCOCT: Cut up contigs** {% icon tool %})
>    * In *"Advanced options"*:
>      * *"Read length for coverage"*: `{'id': 1, 'output_name': 'output'}`
{: .hands_on}

Finally, the clustering results are mapped back to the **original contigs, allowing each contig to be assigned to its corresponding bin. This process ensures accurate and meaningful binning of metagenomic data.

> <hands-on-title> Merge fragment clusters </hands-on-title>
>
> 1. {% tool [CONCOCT: Merge cut clusters](toolshed.g2.bx.psu.edu/repos/iuc/concoct_merge_cut_up_clustering/concoct_merge_cut_up_clustering/1.1.0+galaxy2) %} with the following parameters:
>    * {% icon param-file %} *"Clusters generated by CONCOCT"*: `output_clustering` (output of **CONCOCT** {% icon tool %})
{: .hands_on}


While **CONCOCT** generates a table mapping contigs to their respective bins, it does not automatically produce **FASTA files** for each bin. To obtain these sequences for further analysis, users must employ the **`CONCOCT: Extract a FASTA file`** utility. This tool combines the original contig FASTA file with CONCOCT’s clustering results, extracts contigs assigned to a specific bin, and outputs a **FASTA file representing a single metagenome-assembled genome (MAG)**. This step is crucial for enabling downstream genomic analyses.







> <hands-on-title> Extract MAG FASTA files </hands-on-title>
>
> 1. {% tool [CONCOCT: Extract a fasta file](toolshed.g2.bx.psu.edu/repos/iuc/concoct_extract_fasta_bins/concoct_extract_fasta_bins/1.1.0+galaxy2) %} with the following parameters:
>    * {% icon param-collection %} *"Original contig file"*: `output` (Input dataset collection)
>    * {% icon param-file %} *"CONCOCT clusters"*: `output` (output of **CONCOCT: Merge cut clusters** {% icon tool %})
{: .hands_on}

> <question-title>Binning metrics</question-title>
>
> 1. How many bins where produced by MaxBin2 for our sample?
> 2. How many contigs are in the bin with most contigs? 
> > <solution-title></solution-title>
> >
> > 1. There are 10 bins for this sample.
> > 2. 50 - while all other bins only contain one contig each ! 
> >
> {: .solution}
>
{: .question}