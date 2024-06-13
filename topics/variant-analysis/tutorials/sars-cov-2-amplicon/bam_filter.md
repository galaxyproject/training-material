#### Filtering mapped reads by quality and pairing

Earlier in the tutorial we used `fastp` to filter reads by read quality. Read quality is repored using
a [logarithmic model]({{ site.baseurl }}/topics/sequence-analysis/tutorials/quality-control/slides.html#10) of
anticipated accuracy. The `bwa mem` mapper reports mapping quality using a [similar model](https://bioinformatics.stackexchange.com/questions/2417/meaning-of-bwa-mem-mapq-scores). Thus each mapped read has a map quality (MapQ) score.

In addition to the mapping quality reported for each individual read, we can use reads where both reads in a pair are mapped with more confidence than reads where only one of a pair was mapped.

We can use the {% tool [Filter SAM or BAM](toolshed.g2.bx.psu.edu/repos/devteam/samtool_filter2/samtool_filter2/1.8+galaxy1) %} tool
to filter both based on *MapQ* and on flags in the aligned reads (BAM) that determine if the reads are mapped correctly.

> ### {% icon hands_on %} Filtering aligned reads with
> 1. Select the {% tool [Filter SAM or BAM](toolshed.g2.bx.psu.edu/repos/devteam/samtool_filter2/samtool_filter2/1.8+galaxy1) %} tool and select the {% icon param-collection %} *Map with BWA-MEM on collection XX* dataset collection as input.
> 2. Set the *Minimum MAPQ quality score* to *20*.
> 3. Set *Filter on bitwise flag* to *yes* and in the *Only output alignments with all of these flag bits set* section select the *Read is paired* and *Read is mapped in a proper pair* options
> 4. *Execute* the tool
{: .hands_on}