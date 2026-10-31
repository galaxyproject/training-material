> <comment-title>Background: Taxonomic assignment</comment-title>
>
> In this tutorial we will use the RDP classifier and reference taxonomy for classification, but there are several different taxonomic
> assignment algorithms and reference databases available for this purpose.
>
> An overview of different methods is given by {% cite Liu2008 %} and shown below:
>
> ![overview of different methods for taxonomy assignment](../../images/classification_methods.jpg){:width="50%"}
>
> The most commonly used reference taxonomies are:
>
>  - SILVA ({% cite Quast2012 %})
>  - GreenGenes ({% cite DeSantis2006 %})
>  - RDP ({% cite Cole2013 %})
>  - NCBI Taxonomy Database ({% cite Federhen2011 %})
>
> The choice of taxonomic classifier and reference taxonomy can impact downstream results. The figure from {% cite Liu2008 %}
> given below shows the taxonomic composition determined when using different classifiers and reference taxonomies, for different primer sets (16S regions).
>
> ![comparison of reference taxonomies](../../images/reference_taxonomy_comparison.jpg "Compositions at the phylum level for each of the three datasets: (a) Guerrero Negro mat, (b) Human gut and (c) Mouse gut, using a range of different methods (separate subpanels within each group). The x-axis of each graph shows region sequenced. The y-axis shows abundance as a fraction of the total number of sequences in the community. The legend shows colors for phyla (consistent across graphs)."){:width="75%"}
>
> <!-- figure captions not working in includes, hardcode html til fixed -->
> <figure><figcaption><span class="figcaption-prefix">Figure:</span> Compositions at the phylum level for each of the three datasets: (a) Guerrero Negro mat, (b) Human gut and (c) Mouse gut, using a range of different methods (separate subpanels within each group). The x-axis of each graph shows region sequenced. The y-axis shows abundance as a fraction of the total number of sequences in the community. The legend shows colors for phyla (consistent across graphs). </figcaption></figure>
>
> Which reference taxonomy is best for your experiments depends on a number of factors such as the type of sample and variable region sequenced.
>
> Another discussion about how these different databases compare was described by {% cite Balvoit2017 %}.
>
{: .comment}

