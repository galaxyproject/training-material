We first need to figure out where the sequenced DNA fragments originated from in the genome. The reads must be aligned to a reference genome to identify the {{ include.to_identify }}.

This is equivalent to solving a jigsaw puzzles, but unfortunately, not all pieces are unique.

In principle, we could do a BLAST analysis to figure out where the sequenced pieces fit best in the known genome. Aligning millions of short sequences this way may, however, take a couple of weeks. And we do not really care about the exact base to base correspondence (**alignment**). What we are really interested in is "where these reads came from". This approach is called **mapping**.

Nowadays, there are many mapping programs for shotgun sequenced DNA. We will use [{{ include.mapper }}]({{ include.mapper_link }}).

{% include topics/sequence-analysis/tutorials/mapping/ref_genome_explanation.md %}
