## Bin contigs using SemiBin

**SemiBin** ({% cite Pan2022 %}) is a *semi-supervised deep learning method* for metagenomic binning. It uses both **must-link** and **cannot-link** constraints derived from single-copy marker genes to guide binning, allowing higher accuracy than purely unsupervised methods.

> Metagenome binning is essential for recovering high-quality metagenome-assembled genomes (MAGs) from environmental samples. SemiBin applies a semi-supervised Siamese neural network that learns from both contig features and automatically generated constraints. It has been shown to recover more high-quality and near-complete genomes than MetaBAT2, MaxBin2, or VAMB across multiple benchmark datasets. SemiBin also supports single-sample, co-assembly, and multi-sample binning workflows, demonstrating excellent scalability and versatility.
> {: .quote author="Pan et al., 2022" }

SemiBin is increasingly popular for metagenomic binning due to:

* **High reconstruction quality**: SemiBin usually recovers **more high-quality and near-complete MAGs** than traditional binners, including MetaBAT2.

* **Semi-supervised learning**: it combines deep learning with automatically generated constraints to better separate similar genomes.

* **Flexible binning modes**: it works with:

  * individual samples
  * co-assemblies
  * multi-sample binning

* **Support for multiple environments**: SemiBin includes trained models for:

  * human gut
  * dog gut
  * ocean
  * soil
  * a *generic pretrained model*.


> <hands-on-title>Individual binning of short reads with SemiBin</hands-on-title>
>
> 1. {% tool [SemiBin](toolshed.g2.bx.psu.edu/repos/iuc/semibin/semibin/2.1.0+galaxy1) %} with the following parameters:
>
>    * *"Binning mode"*: `Single sample binning (each sample is assembled and binned independently)`
>
>      * {% icon param-collection %} *"Contig sequences"*: `Contigs` (Input dataset collection)
>      * {% icon param-file %} *"Read mapping to the contigs"*: output of **Samtools sort** {% icon tool %}
>      * *"Reference database"*: `Use SemiBin ML function`
>      * *"Environment for the built-in model"*: `<choose Global or None>`
>
>    * *"Method to set up the minimal length for contigs in binning"*: `Automatic`
>
>    > <comment-title>Environment for the built-in model</comment-title>
>    >
>    > SemiBin provides several pretrained models. If a model matching your environment is available, selecting it can improve binning performance.
>    >
>    > If no environment-specific model fits your data, you may choose:
>    >
>    > * **Global** — a general-purpose pretrained model trained across many environments.
>    > * **None** — no pretrained model is used. SemiBin then runs in fully unsupervised mode, which is recommended when your environment differs substantially from all available pretrained models.
>    {: .comment}
>
{: .hands_on}

> <question-title>Binning metrics</question-title>
>
> 1. How many bins where produced by SemiBin for our sample?
> 2. How many contigs are in the bin with most contigs? 
> > <solution-title></solution-title>
> >
> > 1. There is only one bin for this sample.
> > 2. 50 (these numbers may differ slightly depending on the version of SemiBin). So not all contigs where binned into this bin !
> >
> {: .solution}
>
{: .question}