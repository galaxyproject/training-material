## Trimming the aligned reads with `ivar trim`
You need to trim the primers using `ivar trim` from the tool bar.
> 1. The input is the bam dataset
> 2. Select **Built-in** from the *Source of primer information*
> 3. Select **SARS-CoV-2-ARTICv3** from the *Primer scheme name*
> 4. Select **Yes** for *Include reads not ending in any primer binding sites?*
> 5. Enter Quality score (FASTQ quality) of **20**
> 6. Enter Frequency threshold of **0.7**. This threshold is preferred for picking out common variants 
