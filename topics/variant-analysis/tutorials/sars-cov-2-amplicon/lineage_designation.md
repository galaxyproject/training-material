## Lineage designation

Phylogenetic lineage assignment is important in tracing both existing and emerging SARS-CoV-2 variants.
Mutations occur in the genome of SARS-CoV-2 virus to produce a fingerprint that can be used for to infer ancestral relationships. The specific collection of mutations in the viral genome leads to a particular lineage/clade of SARS-CoV-2 variant. Here variant means genetically distinct SARS-CoV-2 strain from the reference genome. There are two projects that current track lists of variants: [Nextclade](https://clades.nextstrain.org/) (part of the [Nextstrain project](https://nextstrain.org/)) and the [PANGO network](https://cov-lineages.org/). Some variants have been noted a Variants of Concern (VOCs) and Variants of Interest (VOIs). The World Health Organisation (WHO) maintains a [list](https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/) of these VOCs and VOIs with easy to use names (e.g. alpha is PANGO lineage B.1.1.7 and Nextstrain clade 20I (V1), beta is PANGO lineage B.1.351 and Nextstrain claid 20H (V2), etc). Public Health England also maintains a [list](https://github.com/phe-genomics/variant_definitions) of VOCs and VOIs (not necessarily the same as the WHO list) with associated amino acid mutations and finally the US CDC maintains its own VOC and VOI [list](https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html) with associated spike protein mutations.

>### Phylogenetic lineage assignment with `Pangolin` and `Nextclade`
> 1. Using the {% tool [pangolin](toolshed.g2.bx.psu.edu/repos/iuc/pangolin/pangolin/3.0.3+galaxy0) %} tool, select the dataset collection {% icon param-collection %} output by `ivar consensus` as the input.
> 2. Leave the **Download the latest Pangolin from web** option as default for the pangoLEARN source section and click *Execute*
> 3. Find the {% tool [nextclade](toolshed.g2.bx.psu.edu/repos/iuc/nextclade/nextclade/0.14.3+galaxy1) %} tool and again use the `ivar consensus` output as its input, then click *Execute*
> 4. **Optional:** using the {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/0.1.0) %} tool on the dataset collections output by `pangolin` and `nextclade`, concatenate each collection into a single file for easier viewing.
>
{: .hands_on}

>### Examining lineage classification
>
> 1. For sample *ERR4970107*, what was the `pangolin` classification?
> 2. What was the `nextclade` classification for sample *ERR4970107*?
> 3. Pangolin can optionally use a new classification model, *UShER*. Re-run `pangolin` with the UShER model - does anything change?
> 4. Both `pangolin` and `nextclade` classify the diversity of SARS-CoV-2 using a tree model. At the time of this writing there were [18 nextstrain clades](https://clades.nextstrain.org/) used by `nextclade`. Are there more PANGO lineages than nextstrain clades?
>
> > ### {% icon solution %} Solution
> > 1. ERR4970107 was classified as PANGO lineage *B.1*. According to the [PANGO network description](https://cov-lineages.org/lineage_description_list.html), B.1 is a large lineage that arose during the outbreak in Northern Italy in early 2020. Before the emergence and spread of the Variants of Concern (VOCs) many SARS-CoV-2 samples fell into this lineage but in many parts of the world it has now been replaced by VOCs.
> > 2. Nextclade assigns ERR4970107 to 20A, which the Nextstrain team [calls](https://nextstrain.org/blog/2021-01-06-updated-SARS-CoV-2-clade-naming) "a basal lineage bearing S 614G that's globally distributed". Examining the SnpEff annotated VCF and the Nextclade output we can see that ERR4970107 contains the D614G mutation in the S (spike) protein.
> > 3. Using the UShER model, sample ERR4970106 is designated in PANGO lineage B.1.1 whereas previously it was designed as lineage B.1. Both UShER and non-UShER models use alignment and tree based models for lineage assignments (you can read about the original pangolin model [here](https://virological.org/t/pangolin-web-application-release/482)), but they can sometimes come to different conclusions on the details of lineage assignment. The PANGO network defines B.1.1 as a descendant of the B.1 lineage that has 3 specific nucleotide changes in positions 28881, 28882 and 28883 of the genome. Samples ERR4970105 and ERR4970106 have these changes, ERR4970107 does not. Note that if problems with sequences lead to ambiguous nucleotides (especially Ns), some of the lineage-defining characteristics can be lost, leading to weak or incorrect lineage classification.
> > 4. There are many more PANGO lineages than Nextstrain clades. The Nextclade model was built with stability in mind and thus focuses on classifying viruses into a smaller set of stable clades, whereas the PANGO lineages try and capture a larger diversity of viral variants.
> > {: .solution}
{: .question}