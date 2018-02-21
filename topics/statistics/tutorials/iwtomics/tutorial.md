---
layout: tutorial_hands_on
topic_name: statistics
tutorial_name: iwtomics
---

# Introduction
{:.no_toc}

IWTomics implements the Interval-Wise Testing (IWT) for omics data. This
inferential procedure tests for differences in "Omics" data between two groups
of genomic regions, and does not require fixing location and scale at the
outset.
It is composed of three steps (corresponding to three tools):

 1. loading and pre-processing (IWTomics Load Smooth and Plot);
 2. performing Interval-Wise Testing (IWTomics Test and Plot);
 3. selecting test scale (IWTomics Plot with Threshold on Test Scale).

# Real data example

**Recombination hotspots in fixed ETn vs Control**

> This example contains two region datasets "ETn fixed", "Control" and one feature "Recombination hotspots content". In particular, the region dataset "ETn fixed" contains 1296 genomic regions of 64 kb surrounding fixed ETns elements (32-kb flanking sequences upstream and 32-kb flanking sequences downstream of each element). The region dataset "Control" contains 1142 regions of 64 kb without elements, used as control in the test. The regions are aligned around their center (i.e. around the ETn integration sites).

> Recombination hotspots measurements are associated to each "ETn fixed" and "Control" region. In particular, this feature is measured in 1-kb windows, so that each region is associated to a recombination hotspots curve made of 64 values. The measurement used is the feature content, i.e. the fraction of the 1-kb window that is covered by recombination hotspots

# Loading and pre-processing

> The first tool (IWTomics Load Smooth and Plot) imports a collection of genomic region datasets, and associates to each region multiple genomic feature measurements. It allows to align the regions in multiple ways (center, left, right or scale alignment), to smooth the feature curves (possibly filling gaps in the measurements) and to create a graphical representation of the feature measurements in each region datasets (aligned curves or pointwise quantile curves).

# Performing Interval-Wise Testing

> The second tool (IWTomics Test and Plot) statistically evaluates differences in genomic features between groups of regions along the genome. In particular, it implements the Interval-Wise Testing for omics data, an extended version of the Interval-Wise Testing for functional data presented in Pini and Vantini (2017).

> It allows to perform multiple two sample permutation tests between pairs of region datasets, on several features. It returns the adjusted p-value curves for every test and all possible scales. Moreover, it creates a graphical representation of the Interval-Wise Testing results and a summary plot (optional) with p-values at the maximum scale. The tool IWTomics Plot with Threshold on Test Scale permits to select the scale to be used in the plots.

# Selecting test scale

> The third tool (IWTomics Plot with Threshold on Test Scale) allows to select the scale for the Interval-Wise Testing results. In particular it returns the p-value curves for the different tests performed at the selected scale, and it creates a graphical representation of the Interval-Wise Testing results and a summary plot (optional) at the selected scale.

> ### {% icon tip %} Additional resources:
>
> read more about **IWTomics** [here](https://bioconductor.org/packages/release/bioc/vignettes/IWTomics/inst/doc/IWTomics.pdf).
{: .tip}  
