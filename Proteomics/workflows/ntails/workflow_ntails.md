# GalaxyP Workflow: N-Tails Data Analysis

N-Tails is a special Proteomics technique to analyze peptide abundancy changes of protein N-termini. Prior to the MS measurement, N-Tails enriches unmodified, as well as acetylated N-termini. Both common and "unusual" N-termini are identified, where "unusual" means that the protein N-terminus was changes. 
This is best explained by an example: directly after translation, a protein has exactly one N-terminus. When a protease is cutting the protein in half, each half has its own N-terminus. While the N-terminus of the first half protein is exactly the same as the one of the full protein precursor, the N-terminus of the second half is different and depends on the amino acid sequence where the protein was cut.

The figure below illustrates the mechanism of N-Tails. It was originally published by [Stefan Tholen (2014)](link zu doktorarbeit in freidok). Further reading on N-Tails and other N-terminal techniques, see [Tholen et al., Springer Vienna, 2013](http://dx.doi.org/10.1007/978-3-7091-0885-7_5).

[!N-Tails technique](../images/WF_ntails_technique.png)

The N-Tails technique was originally designed to research protease biology and has most often been used in this field. It was originally published in [Kleifeld et al., Nat. Biotechnol., 2010](http://www.ncbi.nlm.nih.gov/pubmed/20208520).

> ### :exclamation: Warning: Interpretation of N-Tails results
>
> Do not overinterprete the results of N-Tails experiments. While the technique **is** fit to identify direct protease substrates, it does not discriminate direct from indirect ("downstream") effects. Thus, most of the identified N-termini will **not** be direct protease substrates, even if their change in protein abundance is statistically significant.
> This warning is not specific for the N-Tails technique, but applies to every proteomic N-terminal technique (e.g. COFRADIC).
> {: .warning}

This workflow was originally built in the OpenMS framework "TOPPAS" and published in [Lai et al., Proteomics, 2015](https://www.ncbi.nlm.nih.gov/pubmed/26013158). It was rebuild for the Galaxy framework by Melanie Föll.
(We compared the output of this Galaxy workflow to the original data and confirm that the WF gives the same results as the original one.)

The figure below gives an overview of the used nodes. For further description of the workflow, please consider the [original publication](https://www.ncbi.nlm.nih.gov/pubmed/26013158).

[!N-Tails Galaxy Workflow](../images/WF_ntails.png)