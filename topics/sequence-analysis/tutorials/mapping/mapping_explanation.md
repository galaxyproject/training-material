We first need to figure out where the sequenced DNA fragments originated from in the genome. The reads must be aligned to a reference genome to identify the {{ include.to_identify }}.

This is equivalent to solving a jigsaw puzzles, but unfortunately, not all pieces are unique.

In principle, we could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way may, however, take a couple of weeks. And we do not really care about exact base to base correspondence (**alignment**). We are more interesting on "where did my reads come from?", an approach called **mapping**.

Nowadays, there are many read alignment programs for shotgun sequenced DNA. We will use [{{ include.mapper }}]({{ include.mapper_link }}).

We need a reference genome on which we map the reads.

> ### {% icon question %} Questions
>
> 1. What is a reference genome?
> 2. For human, several possible reference genomes are available: hg18, hg19, hg38, etc. What do they correspond to?
> 2. Which reference genome should we use?
>
> > ### {% icon solution %} Solution
> > 1. A reference genome (or reference assembly) is a DNA sequence assembled as a representative example of a species' DNA sequence. As they are often assembled from the sequencing of the DNA from different donors, they do not accurately represent the set of genes of any single person but a mosaic of different DNA sequences from each donor.
> > 2. As the cost of DNA sequencing falls, and new full genome sequencing technologies emerge, more genome sequences continue to be generated. Using these new sequences, new alignments are builts and the reference genomes improved (fewer gaps, fixed misrepresentations in the sequence, etc). The different reference genomes correspond to the different versions released (named build)
> > 3. {{ include.answer_3}}
> {: .solution }
{: .question}