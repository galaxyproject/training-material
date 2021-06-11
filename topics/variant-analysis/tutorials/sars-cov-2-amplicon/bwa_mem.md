
### Mapping FASTQ reads with `BWA MEM`
Now that you have clean reads, you will align, or map, the cleaned reads collection to the reference genome using the [`bwa mem`](https://github.com/lh3/bwa) tool. This step tells you which precise location in the genome each base pair in each sequencing read comes from. Mapped reads from here will be used to identify genetic variants in the downstream steps. If you want to understand more about how mapping works, read the tutorials in the [Sequence Analysis]({{ site.baseurl }}/topics/sequence-analysis) section.

> ### {% icon hands_on %} Read mapping with `bwa mem`
> 0. We will be using the {% tool [bwa mem](toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.17.1) %} tool.
> 1. For the **reference genome** select *Use a genome from history and build index* from the drop down arrow and select the `MN908947_3_Wuhan-Hu-1.fasta` dataset as the reference 
> 3. Select **Paired Collection** from the *Single or Paired-end reads* drop down arrow.
> 4. Select the `list:paired` collection that you created (i.e. the one named `samples`) as the input reads. Galaxy will run `bwa mem` for each sample in the collection.
> 5. Leave the other options on their default values and click *Execute*
{: .hands_on}