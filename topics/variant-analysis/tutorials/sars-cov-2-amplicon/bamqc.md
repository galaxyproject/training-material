#### Examining mapping results using BamQC

Amplicon sequencing works by PCR amplifying viral sequencing using pools of *tiled amplicons*. As with any PCR reaction,
multiple factors including primer binding affinity, abundance of sub-genomic DNA fragments and annealing temperature parameters can result very widely varying abundance of different amplicons. For more information on the use and optimisation of the ARTIC protocol read this [guide by Josh Quick](https://www.protocols.io/view/ncov-2019-sequencing-protocol-v3-locost-bh42j8ye) and this [technical note from the Sanger Institute](https://s3.amazonaws.com/protocols-files/public/a07b6b20b1986ca1f94d2c5604d9d6b99912950b9f571acee28a457e5264d59b/cmm8bcqrx.pdf). Also note that as SARS-CoV-2 evolves, approaches to sequencing it need to evolve to. See, for example, this [preprint](https://www.medrxiv.org/content/10.1101/2021.06.01.21258181v1.full) on the impact of mutations in emerging variant of concern lineages on SARS-CoV-2 sequencing.

The {% tool [BamQC](toolshed.g2.bx.psu.edu/repos/iuc/qualimap_bamqc/qualimap_bamqc/2.2.2d+galaxy1) %} tool collects both statistics on mapping (i.e. on a BAM file) and also produces plots illustrating the read coverage across the genome. We will use this tool to do a detailed inspection of our filtered BAM output.

> ### {% icon hands_on %} Examining genomic coverage with *BamQC*
> 1. Select the {% tool [BamQC](toolshed.g2.bx.psu.edu/repos/iuc/qualimap_bamqc/qualimap_bamqc/2.2.2d+galaxy1) %} tool and select the {% icon param-collection %} *Filter SAM or BAM, output SAM or BAM on data XX* dataset collection as input.
> 2. In the *Skip duplicate reads* options, ensure that the *Duplicates detected by Qualimap* option is selected
> 3. *Execute* the tool
{: .hands_on}

> ### {% icon question %} Interpreting coverage depth plots
>
> *BamQC* produces two sets of output: reports and raw data. The raw data is designed for other programs (e.g. *MultiQC*) to read.
> We will focus on the reports, and specifically the one for sample ERR4970105. 
>
> In the *BamQC* output, examine the report for , pay special attention to the *Mean Coverage* (in section *Coverage*) and the *Coverage across reference*, *Coverage Histogram* and *Coverage Histogram (0-50X)* plots.
>
> 1. What was the mean coverage across the genome?
> 2. How many reads were mapped?
> 3. How many reads were mapped in a proper pair?
> 4. How many bases mapped against the genome?
> 5. Do you think the mapping was successful?
>
> > ### {% icon solution %} Solution
> >
> > 1. The sample had 1393396 reads
> > 2. 1393232 (almost all) of the reads mapped against the reference
> > 3. 1391622 reads were in [proper pairs](https://www.biostars.org/p/8318/). This means that both ends of the paired end reads mapped and the insertion size in the reads was of an expected size.
> > 4. A total of 191464525 mapped against the genome. The SARS-CoV-2 genome is approximately 30,000 bases long, so the sample contains approximately 6400 genomes worth of reads.
> > 5. Yes, 99% of the reads mapped against the reference genome and almost all of them in proper pairs. This suggests that the reads represent SARS-CoV-2 sequence and give us some confidence that variant calling and genome consensus reconstruction will succeed.
> {: .solutio