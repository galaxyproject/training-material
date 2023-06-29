---
layout: tutorial_hands_on
title: Calculating α and β diversity with Krakentools
zenodo_link: xxx
questions:
- How many different species are present in my sample? How do I additionally take their relative abundance into account?
- How similar or how dissimilar are my samples?
- What are the different metrics used to calculate the diversity of my samples?
objectives:
- Explain what diversity is
- Explain different metrics to calculate α and β diversity 
- Apply Krakentools to calculate α and β diversity and understand the output
level: Introductory
key_points:
- There are 2 different types of diversity metrics (α and β diversity)
- Krakentools can be used in Galaxy for calculating the diversity
time_estimation: 20M
contributions:
   authorship:
    - sophia120199
    - bebatut
tags:
- metagenomics
- diversity
---

# Introduction

A **diversity index** is a quantitative measure that is used to assess the level of diversity or variety within a particular system, such as a biological community, a population, or a workplace. It provides a way to capture and quantify the distribution of different types or categories within a system.

In various fields, diversity indexes are employed to understand and compare the composition and richness of various elements. Apart from ecology, fields such as social and cultural science are interested in the diversity within a population or workplace. In these cases, the indexes may consider factors like age, gender, ethnicity, or other relevant characteristics to assess the diversity and inclusiveness of a group or organization.

To study microbiome data, indirect methods like **metagenomics** can be used. Metagenomic samples contain DNA from different organisms at a specific site, where the sample was collected. Metagenomic data can be used to find out which organisms coexist in that niche and which genes are present in the different organisms.

Once we know which species are present in a metagenomic sample ([Tutorial on Taxonomic Profiling and Visualization of Metagenomic Data](https://training.galaxyproject.org/training-material/topics/metagenomics/tutorials/taxonomic-profiling/tutorial.html])), we can do diversity analyses.

Related to ecology, the term **diversity** describes the number of different species present in one particular area and their relative abundance. More specifically, several different metrics of diversity can be calculated. The most common ones are α, β and γ diversity:

- **α diversity** describes the diversity within a community

    It considers the number of different species in an environment (also referred to as species **richness**). Additionally, it can take the abundance of each species into account to measure how evenly individuals are distributed across the sample (also referred to as species **evenness**). 

- **β diversity** measures the distance between two or more separate entities

    It therefore describes the difference between two communities or ecosystems.

- **γ diversity** is a measure of the overall diversity for the different ecosystems within a region.

    ![α, β and γ diversity](../images/diversity_differences.png "α, β and γ diversity")

In this analysis we will use Galaxy for calculating the Shannon's alpha diversity index and the Bray-Curtis dissimilarity index for β diversity. 

# Background on data

The dataset we will use for this tutorial comes from an oasis in the Mexican desert called Cuatro Ciénegas ({% cite Okie.2020 %}). The researchers were interested in genomic traits that affect the rates and costs of biochemical information processing within cells. They performed a whole-ecosystem experiment, thus fertilizing the pond to achieve nutrient enriched conditions.

Here we will use 2 datasets:
- `JP4D`: a microbiome sample collected from the Lagunita Fertilized Pond
- `JC1A`: a **control** samples from a control mesocosm.

The datasets differ in size, but according to the authors this doesn't matter for their analysis of genomic traits. Also, they underline that differences between the two samples reflect trait-mediated ecological dynamics instead of microevolutionary changes as the duration of the experiment was only 32 days. This means that depending on available nutrients, specific lineages within the pond grow more successfully than others because of their genomic traits.

The datafiles are named according to the first four characters of the filenames.

Originally, it was a collection of paired-end data with R1 being the forward reads and R2 being the reverse reads. The samples have than been analysed as explained in the [Taxonomic profiling tutorial]({% link topics/sequence-analysis/tutorials/taxonomic-profiling/tutorial.md %}).

In a nutshell, taxonomic labels have been assigned to the metagenomics data using [Kraken2](toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.1.1+galaxy1) to find out which species are present in the samples. Finally, species abundance was estimated using [Bracken](toolshed.g2.bx.psu.edu/repos/iuc/bracken/est_abundance/2.7+galaxy1). For this tutorial, we will use the output file of Bracken.

Here, to get an overview, you can find a Krona chart visualizing the different species present in the two samples.

<iframe id="krona" src="krona-kraken.html" frameBorder="0" width="100%" height="900px"> ![Krona chart with multi-layered pie chart representing the community profile with in the center the higher taxonomy levels (i.e. domain) and on the exterior the more detailed ones (i.e. species)](./images/krona-kraken.png) </iframe>

The dataset we will work with in this tutorial is the output file of Bracken, which estimates species abundance.

![Output file of Bracken](../images/bracken_output.png "Output file of Bracken")

xxx output file description

# Prepare Galaxy and data

Any analysis should get its own Galaxy history. So let's start by creating a new one:

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

We need now to import the data

> <hands-on-title>Import datasets</hands-on-title>
>
> 1. Import the following samples via link from [Zenodo]({{ page.zenodo_link }}) or Galaxy shared data libraries:
>
>    ```text
>    {{ page.zenodo_link }}/files/xxx
>    {{ page.zenodo_link }}/files/xxx

>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 2. 3. Create a paired collection.
>
>    {% snippet faqs/galaxy/collections_build_list_paired.md %}
>
{: .hands_on}

# Calculating α diversity

**α diversity** describes the diversity within a community. There are several different indexes used to calculate α diversity because different indexes capture different aspects of diversity and have varying sensitivities to different factors. These indexes have been developed to address specific research questions, account for different ecological or population characteristics, or highlight certain aspects of diversity. 

There are various measures of alpha diversity accessible:
- **richness** indexes that estimate the quantity of distinct species within a sample
- **evenness** indexes that evaluate the relative abundances of species rather than their total count
- **diversity** indexes that incorporate both the relative abundances and total count of distinct species

In the table below you can find a list of commonly used indexes to calculate α diversity and their description.

![α diversity](../images/alphadiversity_metrics.png "α diversity")

![richness and evenness](../images/alpha_diversity_richness_evenness.png "richness and evenness")

| Indices for α diversity | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | Class     |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- |
| Shannons                | Calculates the uncertainty in predicting the species identity of an individual that is selected from a community.                                                                                                                                                                                                                                                                                                                                                                                                                 | Diversity |
| Berger-Parker           | Expresses the proportional importance of the most abundant type. Highly biased by sample size and richness.                                                                                                                                                                                                                                                                                                                                                                                                                       | Diversity |
| Simpsons                | Calculates the probability that two individuals selected from a community will be of the same species. Obtains small values in datasets of high diversity and large values in datasets of low diversity.                                                                                                                                                                                                                                                                                                                          | Diversity |
| Inverse Simpons         | Transformation of Simpsons index that increases with increasing diversity.                                                                                                                                                                                                                                                                                                                                                                                                                                                        | Diversity |
| Fishers                 | Describes the relationship between the number of species and the number of individuals in those species. Parametric index of diversity that assumes that the abundance of species follows a log series distribution.                                                                                                                                                                                                                                                                                                              | Diversity |
| Pielou’s evenness       | Quantifies how close the community’s diversity is to the maximum possible diversity. This index is calculated by taking the Shannon Diversity Index (which measures the overall diversity of the community) and dividing it by the maximum possible diversity given the observed species richness.                                                                                                                                                                                                                                | Evenness  |
| Margalef’s richness     | Indicates the estimated species richness, accounting for the community size. This metric takes into account that a larger community size can support a greater number of species.                                                                                                                                                                                                                                                                                                                                                 | Richness  |
| Chao1                   | Estimates the true species richness or diversity of a community, particularly when there might be rare or unobserved species. Chao1 estimates the number of unobserved species based on the number of singletons and doubletons. It assumes that there are additional rare species that are likely to exist but have not been observed. The estimation considers the number of unobserved singletons and doubletons and incorporates them into the observed species richness to provide an estimate of the true species richness. | Richness  |
| ACE                     | ACE (Abundance-based Coverage Estimator) takes into account the abundance distribution of observed species and incorporates the presence of rare or unobserved species. ACE estimates the number of unobserved species based on the abundance distribution and incorporates it into the observed species richness. It takes into account the relative rarity of observed species and uses this information to estimate the true species richness.                                                                                 | Richness  |                                                                                        |

> <details-title> Mathematical expressions for calculating α diversity</details-title>


| Index               | Mathematical Expression                                                                         | Description                                                                                                                                                                                                                                       |
| ------------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Simpson             | D = ∑<sub>i=1</sub><sup>S</sup> (n<sub>i</sub>/N)<sup>2</sup>                                   | ni is the number of individuals in species i, N = total number of individuals of all species, and ni/N = pi (proportion of individuals of species i), and S = species richness.                                                                   |
| Shannon             | H' = -∑<sub>i=1</sub><sup>S</sup> p<sub>i</sub> \* ln(p<sub>i</sub>)                            | pi = proportion of individuals of species i, and ln is the natural logarithm, and  S = species richness.                                                                                                                                          |
| Berger-Parker       | D = n<sub>max</sub>/N                                                                           | n<sub>max</sub> is the abundance of the most dominant species,  and N is the total number of individuals (sum of all abundances).                                                                                                                 |
| Fisher's alpha      | S\=a\*ln(1+n/a)                                                                                 | S is number of taxa, n is number of individuals and a is the Fisher's alpha.                                                                                                                                                                      |
| Pilou's Evenness    | J = H'/ln(S)                                                                                    | H' is Shannon Weiner diversity and S is the total number of species in a sample, across all samples in dataset.                                                                                                                                   |
| Margalef's richness | D = (S - 1) / Log (n)                                                                           | S is the total number of species, and n is the total number of individuals in the sample                                                                                                                                                          |
| Chao1               | S<sub>chao1</sub> = S<sub>obs</sub> + (n<sub>1</sub>(n<sub>1</sub> - 1))/(2(n<sub>2</sub> + 1)) | S<sub>obs</sub> is the observed species richness, n<sub>1</sub> represents the number of species represented by a single individual (singletons), and n<sub>2</sub> represents the number of species represented by two individuals (doubletons). |

{: .details}

Krakentools introduction + metrics available there

> <hands-on-title>Calculate α diversity with Krakentools</hands-on-title>
>
> 1. {% tool [Krakentools: Calculate alpha diversity]([toolshed.g2.bx.psu.edu/view/iuc/krakentools_alpha_diversity/9d0330e23bfd)) %} with the following parameters:
>    - *"Abundance file"*: `Dataset Collection`: uploaded Bracken output file
>      
>    - *"Specify alpha diversity type"*: `Shannon's alpha diversity`
>
>
{: .hands_on}



> <question-title></question-title>
>
> 1. Calculate the 5 different alpha indexes available in Krakentools and compare the results. What do these numbers tell you?
> 2. Are the results consistent among the different indexes?
>
> > <solution-title></solution-title>
> >
> > 1. 

|           | JC1A      | JP4D      | Explanation                                                                                                                                                                                                                                                                                                                                                                                              |
|-----------|-----------|-----------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Shannon   | 5.3441    | 6.4429    | When the Shannon index is given as a value of 5, it indicates a relatively high level of diversity within the community. The index ranges from 0 to a maximum value that depends on the number of species and their relative abundances. The higher the Shannon index value, the greater the diversity within the community.                                                                                                                                           |
| Berger-Parker | 0.2299    | 0.0581    | When the Berger-Parker index is given as a value of 0.23, it suggests that a single species dominates the community, as it represents 23% of the total individuals in the community. This indicates a relatively low level of species evenness, meaning that the abundance of individuals is heavily skewed towards one dominant species. In contrast to the Shannon index, which considers both species richness and evenness, the Berger-Parker index emphasizes the dominance of a particular species. A value of 0.23 indicates that the community is heavily influenced by one species, while the other species in the community are less abundant. In the case of JP4D, the dominant species accounts for only 5% of the total individuals, which implies a more balanced distribution of individuals among different species compared to a higher Berger-Parker index value.                                                                                                                                                                                       |
| Simpson   | 0.9401    | 0.9926    | When the Simpson's index is given as a value of 0.94, it indicates a high level of species diversity and evenness within the community. The index ranges from 0 to 1, with 1 representing maximum diversity. Therefore, a Simpson's index of 0.94 suggests that the community is highly diverse, with a relatively even distribution of individuals among different species. In other words, the value of 0.94 indicates that if you were to randomly select two individuals from the community, there is a 94% probability that they would belong to different species. This implies a rich and balanced community where multiple species coexist in relatively equal abundance. |
| Inverse Simpson | 16.6941   | 136.0287  | When the Inverse Simpson's index is given as a value of 16.69, it suggests a relatively low level of species diversity within the community. The index ranges from 1 to the total number of species in the community, with higher values indicating higher diversity. Therefore, a value of 16.69 indicates a lower diversity compared to a higher index value. An Inverse Simpson's index of 136 suggests a relatively high level of species diversity within the community. The index ranges from 1 to the total number of species in the community, with higher values indicating greater diversity. Therefore, a value of 136 indicates a higher diversity compared to a lower index value. The Inverse Simpson's index is the reciprocal of the Simpson's index, which quantifies species diversity and evenness within a community. A higher Inverse Simpson's index value signifies a community with a greater number of species and a more even distribution of individuals among those species. |
| Fisher    | 3240.0957 | 9163.5027 |                                                                                                           > >
> > 2. The results are consistent as all indexes show JP4D to be the more diverse sample compared to JC1A.
                                                                      |

> {: .solution}
>
{: .question}


> <comment-title></comment-title>

>Apart from Krakentools, there are two more tools available in Galaxy that can be used to calculate diversity indexes, QIIME2 and Vegan.


> QIIME 2 (Quantitative Insights Into Microbial Ecology 2) is a powerful open-source bioinformatics software package that provides a comprehensive suite of tools and methods for processing, analyzing, and visualizing microbiome data. It offers a modular approach to microbiome analysis, allowing researchers to build flexible analysis pipelines tailored to their specific research goals. The software supports a wide range of data types, including 16S rRNA gene sequencing, metagenomics, metatranscriptomics, and others.
> 
>Some of the key features and functionalities of QIIME 2 include:
>1. Data Import and Preprocessing: QIIME 2 supports the import of raw sequencing data and performs quality control and data preprocessing steps, such as demultiplexing, quality filtering, and primer removal.
>2. Taxonomic Assignments: The software enables taxonomic classification of microbial sequences using various algorithms and reference databases.
>3. Diversity Analysis: QIIME 2 allows users to explore and quantify microbial diversity within and between samples. It provides metrics for alpha diversity (within-sample diversity) and beta diversity (between-sample diversity).
>4. Community Analysis: Users can investigate the composition and structure of microbial communities, including taxonomic summaries, abundance profiles, and statistical comparisons between groups.
>5. Phylogenetic Analysis: QIIME 2 supports the construction of phylogenetic trees to infer evolutionary relationships among microbial taxa and perform phylogenetic diversity analysis.
>6. Statistical Analysis: The software offers a wide range of statistical methods for differential abundance analysis, correlation analysis, multivariate analysis, and other types of statistical tests.
>7. Visualization: QIIME 2 provides interactive and customizable visualizations to aid in the exploration and interpretation of microbiome data, including heatmaps, bar plots, PCoA plots, and taxonomic trees.

>The vegan package is a community ecology package in the R programming language. It provides a wide range of tools and methods for analyzing and interpreting ecological data, particularly in the context of community ecology. The package is designed to handle multivariate data and offers various statistical techniques for studying species composition, diversity, and community dynamics.

>The vegan package encompasses several functionalities, including:

>1. Diversity Analysis: vegan offers numerous diversity indices, such as species richness, Shannon diversity index, Simpson index, and many others. These indices allow researchers to quantify the diversity of species within a community and compare diversity between different samples or groups.
>2. Community Similarity: The package provides tools for measuring community similarity or dissimilarity, including popular metrics such as Bray-Curtis dissimilarity and Jaccard index. These metrics allow researchers to assess the degree of similarity between communities and perform clustering or ordination analyses.
>3. Ordination Techniques: vegan includes several ordination methods, such as Principal Component Analysis (PCA), Correspondence Analysis (CA), Non-Metric Multidimensional Scaling (NMDS), and Canonical Correspondence Analysis (CCA). These techniques aid in visualizing and exploring patterns in multivariate ecological data.
>4. Community Classification: The package offers tools for performing community classification and assessing the significance of group differences. It includes methods such as Permutational Multivariate Analysis of Variance (PERMANOVA) and Analysis of Similarities (ANOSIM).
>5. Ecological Network Analysis: vegan provides functions for analyzing ecological networks, including network visualization, calculation of network metrics (e.g., connectance, centrality), and testing network structure.
>6. Ecological Indices: The package includes various ecological indices, such as niche overlap indices, indicator species analysis, and null model analysis for testing community patterns against null hypotheses.
>7. Plotting and Visualization: vegan offers flexible plotting functions to visualize ecological data, including bar plots, scatter plots, biplots, and ordination plots.

{: .comment}

# Calculating β diversity 

**β diversity** measures the distance between two or more separate entities. It therefore describes the difference between two communities or ecosystems. 

There are **multiple indexes** used to calculate β diversity because different indexes emphasize different aspects of compositional dissimilarity between communities or sites.

These indexes have been developed to address specific research questions, accommodate different data types, or provide insights into different dimensions of β diversity. In the table below you can find a list of commonly used indexes to calculate β diversity and their description.



| Indices for β diversity                          | Description                                                                          |
|---------------------------------|--------------------------------------------------------------------------------------|
| Jaccard Index                   | Measures the proportion of shared species between two samples                        |
| Sørensen Index                  | Similar to Jaccard Index, but accounts for species abundance                          |
| Bray-Curtis Dissimilarity       | Measures the dissimilarity of species abundances between two samples                 |
| Kulczynski Dissimilarity        | Measures the dissimilarity in the proportional abundances of shared species          |
|UniFrac|Incorporates information on phylogenetic distances between observed species in the computation. Can be calculated either weighted (accounts for abundances) or unweighted (accounts only for richness).|

> <details-title>More details on calculating β diversity</details-title>

| Index                     | Mathematical Expression                                                      | Description                                                                                                                                                                                                                                                                                                     |
| ------------------------- | ---------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Jaccard Index             | J(X, Y) = \vert X ∩ Y \vert / \vert X ∪ Y\vert                                                  | X ∩ Y represents the intersection of sets X and Y (elements common to both sets), and X ∪ Y represents the union of sets X and Y (all unique elements from both sets combined).                                                                                                                                 |
| Sørensen Index            | DSC = 2 \vert X ∩ Y \vert / \vert X \vert + \vert Y \vert                                                   | X ∩ Y represents the intersection of sets X and Y (elements common to both sets),  and \vert X \vert and \vert Y \vert are the cardinalities of the two sets (i.e. the number of elements in each set)                                                                                                                              |
| Bray-Curtis Dissimilarity | BC<sub>ij</sub> = 1 - (2C<sub>ij</sub> / (S<sub>i</sub> + S<sub>j</sub>))    | C<sub>ij</sub> represents the sum of the absolute differences in abundances between corresponding species in samples i and j, S<sub>i</sub> represents the total abundance or sum of species abundances in sample i, and S<sub>j</sub> represents the total abundance or sum of species abundances in sample j. |
| Kulczynski Dissimilarity  | D = 1 - (S<sub>AB</sub> / (S<sub>A</sub> + S<sub>B</sub> - 2S<sub>AB</sub>)) | S<sub>AB</sub> the number of shared OTUs between communities A and B, S<sub>A</sub> the number of OTUs in community A, and S<sub>B</sub> the number of OTUs in community B                                                                                                                                      |


xxx
![UniFrac](../images/unifrac.png "UniFrac")


{: .details}

## Hands on: Calculate β diversity with Krakentools

> <hands-on-title>Calculate α and β diversity with Krakentools</hands-on-title>
>
>
>  {% tool [Krakentools: Calculate beta diversity (Bray-Curtis dissimilarity)]([https://toolshed.g2.bx.psu.edu/view/iuc/krakentools_beta_diversity/b33f117e9b67]) %} with the following parameters:
>     - *"Taxonomy file"*: `Dataset Collection`: uploaded Bracken output file
>      
>    - *"Specify type of input file"*: `Bracken species abundance file`
{: .hands_on}

> <question-title></question-title>
>
> 1. What is the Bray-Curtis dissimilarity calculated for the two samples?
> 2. What does this number tell you?
>
> > <solution-title></solution-title>
> >
> > 1. xxx
> > 2. The Bray-Curtis dissimilarity  measures the dissimiliraty of two samples. Consequently, an output of 0 represents two samples that are exactly the same, while an output of 1 means they are maximally divergent. In our case, xxx
> {: .solution}
>
{: .question}

# Conclusion

In this tutorial, we look how to calculate α and β  diversity from microbiome data. We apply **Krakentools** to calculate the α and β  diversity of two microbiome sample datasets.
