#### Mapping statistics with `samtools stats`

As we do our analysis we should keep an eye on results to spot any errors that might occur. For mapping,
we can use {% tool [samtools stats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.2+galaxy2)}
to get a summary of how well the mapping step worked.

To check the quality of the alignment using `samtools stats`
> ### {% icon hands_on %} Evaluating mapping with `samtools stats`
> 1. Find the {% tool [samtools stats](toolshed.g2.bx.psu.edu/repos/devteam/samtools_stats/samtools_stats/2.0.2+galaxy2) %} tool from the tool bar
> 2. For the input `BAM file` you need to select the {% icon param-collection %} *Dataset collection* option. This will allow you to select the "Map with BWA-MEM on collection XX (mapped reads in BAM format)" collection.
> 3. Leave the rest of the parameters at their default values, scroll to the bottom and click the *Execute* button.
> 4. When the tool has finished running, click the name of the *Samtools stats* collection. This will open the collection so that you can see the individual datasets within it.
{: .hands_on}

> ### {% icon question %} Understanding mapping statistics
>
> The samtools stats report is described in detail in the corresponding [man page](http://www.htslib.org/doc/samtools-stats.html).
> The first section summarises all mapped reads. We will focus on the stats for sample ERR4970105.
>
> 1. How many reads are in this sample?
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
> {: .solution}
>
{: .question}