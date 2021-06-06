## Variant annotation
### Ivar variants
> 1. Enter Quality score (FASTQ quality) of **20**
> 2. Enter Frequency threshold of **0.7**. This threshold is preferred for picking out common variants

### Convert the iVar tabular variants output to a `.vcf` file format with iVar Variants tool
> Find the `Ivar Variants to VCF` tool from tool bar
### Rename the reference sequence with `Text transformation with sed`
> Find the `Text transformation with sed` from tool bar
> Paste this `/^[^#]/s/^[^\t]+\t/NC_045512.2\t/` sed code in the SED Program window. Here you are renaming the Wuhan reference genome
### Annotate the SARS-CoV-2 variants with the `SnpEff eff` tool. 
> Find the `SnpEff eff` tool from the tool bar. Make sure you use the SARS-CoV-2 version
> Make sure you select the reference sequence you have renamed above.
