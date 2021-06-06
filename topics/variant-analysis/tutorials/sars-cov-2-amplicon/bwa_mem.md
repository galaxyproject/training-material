## Mapping FASTQ reads with `BWA MEM`
Now that you have clean reads, you will align also known as map the fastq collection to the reference genome using `BWA MEM` tool from the tool bar. This step tells you which precise location in the genome each base pair in each sequencing read comes from. Mapped reads from here will be used to identify genetic variants in the downstream steps. Here is the documentation of [bwa](http://bio-bwa.sourceforge.net/bwa.shtml) tool. 
> 1. Input is **fastp out**. 
> 2. For **Reference genome** select *Use a genome from history and build index* from the drop down arrow. 
> 3. Select **paired collection** from the *Single or Paired-end reads* drop down arrow.
> 4. Select the first and second set of reads appropriately
