### Variant annotation with ivar variants and SnpEff

Variant calling is the process of using sequence mapping to infer evidence about the presence of genomic variants. There are many
variant callers available, for example the `ivar variants` command, [lofreq](https://csb5.github.io/lofreq/) and [freebayes](https://github.com/freebayes/freebayes). Most modern variant callers use Bayesian statistics to infer the likelihood of reads predicting a variant. This means that they use evidence from read mapping together with other assumptions about the likelihood of events (called priors) to determine the likelihood of a variant in the genome being present in a particular form. You can read more on the theory of variant calling in [this presentation](https://training.galaxyproject.org/training-material/topics/variant-analysis/slides/introduction.html#1).

We will use the {% tool [ivar variants](toolshed.g2.bx.psu.edu/repos/iuc/ivar_variants/ivar_variants/1.3.1+galaxy2) %} tool for variant calling. This tool is widely used in SARS-CoV-2 genome analysis e.g. in the [ncov2019-artic-nf](https://github.com/connor-lab/ncov2019-artic-nf) workflow and in the [Thiagen Titan workflow](https://github.com/theiagen/public_health_viral_genomics). The Galaxy COVID-19 intrahost allelic variation workflow uses a different, more complex workflow that you can find [here](https://covid19.galaxyproject.org/artic/#a-galaxy-workflow-for-the-analysis-of-illumina-paired-end-sequenced-artic-amplicon-data). Interest in the best tools for SARS-CoV-2 variation analysis, and the evolution of the virus itself, continues to drive further development of tools in this domain, and this tutorial will be updated as new information on approaches becomes available. For those working in a public health context, the [PHA4GE Bioinformatics Pipelines and Visualization working group](https://pha4ge.org/bioinformatics-pipelines-and-visualization/) is maintaining a set of resources describing commonly used workflows and recommended practices at this [GitHub repository](https://github.com/pha4ge/pipeline-resources).

> ### Variant calling and annotation
> 1. Open the {% tool [ivar variants]{toolshed.g2.bx.psu.edu/repos/iuc/ivar_variants/ivar_variants/1.3.1+galaxy2} %} tool.
> 2. Select the collection of BAM files created by {% icon param-collection %} *ivar trim* as input.
> 3. Enter minimum quality score (FASTQ quality) of **20** and minimum frequency threshold of **0.7**. This frequency threshold is preferred for picking out common variants. After setting these parameters *Execute* the tool.
> 4. Select *Both Tabular and VCF* as the Output format.
>
> #### Rename the reference sequence with `Text transformation with sed`
> 1. The *SnpEff* tool requires a genome reference named *NC_045512.2*. To ensure that the VCF has this reference, find the {% tool [Text transformation with sed](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sed_tool/1.1.1) tool %} and...
> 2. Paste this `/^[^#]/s/^[^\t]+\t/NC_045512.2\t/` sed code in the SED Program window and *Execute* the code.
>
> #### Annotate the SARS-CoV-2 variants with the `SnpEff eff` tool. 
> 1. Find the SARS-CoV-2 version of the {% tool [SnpEff eff](toolshed.g2.bx.psu.edu/repos/iuc/snpeff_sars_cov_2/snpeff_sars_cov_2/4.5covid19) %} tool from the tool bar.
> 2. Select the {% icon param-collection %} collection created by the *Text transformation with sed* tool and *Execute* the tool.
> At the end of running these tools you will have a collection of variants in tabular output from the *ivar variants* tool and variants with their effects on the genome annotated (in VCF format) from the *SnpEff* tool.
{: .hands_on}

> ### {% icon question %} Examining variants in SARS-CoV-2
>
> Have a look at the tabular and VCF output of variants for sample *ERR4970105*. The fields in the tabular output are described in [iVar manual page](https://andersen-lab.github.io/ivar/html/manualpage.html).
> 
> 1. How many variants were identified?
> 2. Look at the variant at position 23403. How many reads support this variant?
> 3. The variant at position 23403 involved a change from a *A* to a *G*. What protein was affected by this variant, and how?
> 4. Are there any variants in the list with frequency close to the allele frequency threshold?
>
> > ### {% icon solution %} Solution
> >
> > 1. The sample had 12 variants (the tabular output has 13 lines, but one is a header line).
> > 2. 6966 reads support the variant (more than 99% of reads mapped at this position).
> > 3. The annotation of this variant is `ANN=G|missense_variant|MODERATE|S|GU280_gp02|transcript|GU280_gp02|protein_coding|1/1|c.1841A>G|p.Asp614Gly|`. This format is described [in this PDF](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf). Specifically it means that that the alternative base is a *G* at this position, and is a *missense variant* (i.e. it changes the amino acid encoded by the codon the base is in). It impacts the *S* (spike) protein and changes amino acid 614 in this protein from a Aspartate (Asp or D) to a Glycine (Gly or G). This is the *D614G* mutation which became well known in 2020 and is thought to [confer a selective advantage](https://www.sciencedirect.com/science/article/pii/S0092867420315373).
> > 4. No, all the variants identified have an alternative allele frequency over 99%. The parameters chosen in this tutorial are designed to select for majority variants, and thus we will not see variants with a lower frequency (including those from intra-host viral evolution).
> {: .solution}
>
{: .question}