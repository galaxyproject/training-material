#### Trimming the aligned reads with `ivar trim`

Amplicon sequencing (like that used in the ARTIC protocol) uses primers to enrich the sequencing of viral vs host genome. As you have seen, these reads are mapped to the genome, but the region of the read containing the primer should not be included in the variant calling step (see [this reference](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-1073) for background on the technique). You can trim primers using the {% tool [ivar trim](toolshed.g2.bx.psu.edu/repos/iuc/ivar_trim/ivar_trim/1.3.1+galaxy0) %} tool.

> ### {% icon hands_on %} Read trimming with `ivar trim`
> 1. Find the {% tool [ivar trim](toolshed.g2.bx.psu.edu/repos/iuc/ivar_trim/ivar_trim/1.3.1+galaxy0) %} tool.
> 2. Select the {% icon param-collection %} *Filter SAM or BAM, output SAM or BAM on data XX: bam* dataset collection as input.
> 3. Select **History** for the *Source of primer information* and select the "SARS-CoV-2-ARTICv3.bed" dataset as the *BED file with primer sequences and position*.
> 4. Select **Yes** for *Include reads not ending in any primer binding sites?*
> 5. Enter Quality score (FASTQ quality) of **20** and click on *Execute*
{: .hands_on}