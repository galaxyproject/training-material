---
layout: tutorial_hands_on

title: "Determining multivariable association between various meta’omic features using MaAslin2"
subtopic: metagenomics
tags: 
    - Heatmap
    - metagenomics
    - associations
    - differential analysis
    - statistics

level: Introductory
zenodo_link: https://zenodo.org/records/12614561
questions:
- How do I find associations between microbial features and specific metadata variables ?
objectives:
- Identify statistically significant associations between microbial features and metadata variables (such as clinical conditions, environmental factors, or demographic information) in microbiome data.
- Uncover potential biomarkers associated with specific disease states.
time_estimation: 15M
contributions:
   authorship:
    - renu-pal
    - paulzierep
   editing:
    - shiltemann



---
# Microbiome Association Detection with MaAsLin2

The importance of identifying associations between microbial features and metadata variables using tools like MaAsLin2 lies in several key areas:

- **Understanding Disease Mechanisms:** These associations can provide insights into how changes in microbial composition may contribute to the development or progression of diseases. This understanding is crucial for advancing knowledge of disease mechanisms.

- **Potential Diagnostic Markers:** Identifying microbial biomarkers associated with specific diseases or conditions can potentially lead to the development of diagnostic tests. These tests could aid in earlier detection, more accurate diagnosis, and monitoring of disease progression.

- **Personalized Medicine:** By identifying associations between gut microbes and individual characteristics (such as diet, medication use, or genetic factors), personalized treatment strategies can be developed. This approach may lead to more effective and tailored interventions for patients.

- **Therapeutic Targets:** Associations can reveal potential targets for therapeutic interventions, such as probiotics, prebiotics, or dietary changes, aimed at modulating the gut microbiota to promote health or mitigate disease risks.

- **Advancing Microbiome Research:** Building a comprehensive understanding of microbial associations with various factors enhances microbiome research. This knowledge can contribute to broader insights into microbial ecology, evolution, and interactions within the human body and the environment.

In addition to MaAslin2, Galaxy offers several other differential analysis tools that are widely used in both transcriptomics and microbiome studies. These tools are designed to handle different types of data (e.g., RNA-seq, microbial count data), with varying strengths in terms of statistical power, handling of sparsity, and treatment of compositional data. Some of them are mentioned below: 


| Tool                | Strengths                                    | Weaknesses                               | Comparison to MaAsLin2                          |
|---------------------|----------------------------------------------|------------------------------------------|------------------------------------------------|
| **ANCOM-BC**        | Compositionality and bias correction          | Computationally intensive for large datasets | ANCOM-BC is good for simpler designs, but MaAsLin2 handles complex metadata better. |
| **LEfSe**           | Easy to interpret, focuses on effect size     | No covariates, may overfit               | LEfSe is simpler but lacks the flexibility and multivariable depth of MaAsLin2. |
| **ALDEx2**          | Robust to sparsity, small sample sizes        | Limited handling of complex metadata     | ALDEx2 is suitable for small datasets, but MaAsLin2 is superior in handling multivariable data and covariates. |
| **MetagenomeSeq**   | Handles zero-inflation, sparse data           | Computationally heavy for large datasets | MetagenomeSeq is great for zero-inflated data, but lacks MaAsLin2's multivariable modeling capacity. |
| **Corncob**         | Models both abundance and variability        | Complex to use, requires R expertise     | Corncob excels at overdispersion analysis, but MaAsLin2 is easier for broader multivariable models. |
| **Phyloseq + DESeq2**| Strong for RNA-seq and transcriptomics; integrates with Phyloseq | Lacks compositionality awareness         | While DESeq2 works for microbiome data, MaAsLin2 offers more suitable options for compositional data and covariate handling. |
| **Limma-Voom**      | Effective for RNA-seq and microarray data, handles low counts | Not tailored for compositional microbiome data | Limma-Voom is well-suited for gene expression, but MaAsLin2 better accounts for the unique characteristics of microbiome data. |

- ANCOM-BC and MaAsLin2, outperform general-purpose tools like DESeq2 and limma-voom when it comes to microbiome data. This is due to their handling of the compositional nature of microbiome data and the sparsity typical of microbial datasets.[PMID: 36617187](https://pubmed.ncbi.nlm.nih.gov/36617187/)
- While general methods like DESeq2 and limma-voom are reliable for gene expression analysis, they do not handle the unique properties of microbiome data as effectively as MaAsLin2. The latter provides a more accurate estimation of differential abundance in the presence of metadata confounders.
[PMID: 1009442](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009442 )

![sensitivity and false discovery rate (FDR) across different tools](https://journals.plos.org/ploscompbiol/article/figure/image?size=large&id=10.1371/journal.pcbi.1009442.g004 "Source: <a href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009442#pcbi-1009442-g004">sensitivity and false discovery rate (FDR) across different tools</a>"){:width="60%"}

- The above figure compares various tools for differential abundance detection (Panel A) and multivariable association detection (Panel B) in microbiome studies, based on sensitivity and false discovery rate (FDR).
- **Sensitivity** measures how well the methods detect true signals ,higher values leads to better performance.
- **False discovery rate (FDR)** measures the proportion of false positives among detected signals (lower FDR is better).
- MaAsLin2 is the clear standout for both differential abundance detection and multivariable association detection, showing high sensitivity and maintaining a low FDR.

> <comment-title></comment-title>
>
> For more information on MaAslin2, [click here](https://huttenhower.sph.harvard.edu/maaslin/).
{: .comment}

MaAsLin2 requires the following input files:

- **Taxonomy (or features) file** : \
        This file is tab-delimited.\
        Formatted with features as columns and samples as rows.\
        The transpose of this format is also okay.\
        Possible features in this file include microbes, genes, pathways, etc.
- **Metadata file** : \
        This file is tab-delimited.\
        Formatted with features as columns and samples as rows.\
        The transpose of this format is also okay.

The Taxonomy file can contain samples not included in the metadata file (or vice versa). For both cases, those samples not included in both files will be removed from the analysis. Also the samples do not need to be in the same order in the two files.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Get the data
In this tutorial, the two input files used are:
-  `HMP2_taxonomy.tsv` or taxonomy file
-  `HMP2_metadata.tsv` or metadata file

The files provided were generated from the HMP2 data. To download [Click here](https://ibdmdb.org/)
 
 **Origin** : \
The **HMP2_taxonomy.tsv** and **HMP2_metadata.tsv** files are part of the **Human Microbiome Project 2 (HMP2)**, which is a key component of the Inflammatory Bowel Disease Multi'omics Database [**(IBDMDB)**](https://ibdmdb.org/). The IBDMDB is a large-scale, multi-omic research initiative aimed at understanding the microbiome's role in IBD progression by integrating various omics data like metagenomics, metabolomics, and host genetics. 

The **HMP2_taxonomy.tsv** file contains microbiome data (species abundances) collected from IBD patients and healthy controls, while the **HMP2_metadata.tsv** file includes clinical and demographic metadata for these samples, such as IBD diagnosis (non-IBD, ulcerative colitis(UC), or Crohn’s disease(CD)), dysbiosis state, and treatments like antibiotics.


> <hands-on-title>Getting the data</hands-on-title>
> 1. Create and name a new history for this tutorial.
> 
>    {% snippet faqs/galaxy/histories_create_new.md %}
> 
> 2. Import the files from [Zenodo](https://zenodo.org/records/12614561) or from the data library:
> 
>    ```
> 
>    https://zenodo.org/records/12614561/files/HMP2_taxonomy.tsv
>    https://zenodo.org/records/12614561/files/HMP2_metadata.tsv
>
>    ```
>
>    > <tip-title>Importing data via links</tip-title>
>    >
>    > * Copy the link location (Right-click on the filename then "Copy Link Address")
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field
>    > * Press **Start**
>    {: .tip}
>
> 3. Change the name of the files to `taxonomy` and `metadata` .
>
>    As a default, Galaxy uses the link as the name of the new dataset. It also does not link the dataset to a database or a reference genome.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
> 4. Inspect the content of a file.
>
>    > <tip-title>Inspecting the content of a dataset</tip-title>
>    >
>    > * Click on the {% icon galaxy-eye %} (eye) icon next to the relevant history entry
>    > * View the content of the file in the central panel
>    {: .tip}
>
>    > <question-title></question-title>
>    >
>    > 1.  What is the main difference between the two files?
>    >
>    >
>    > > <solution-title></solution-title>
>    > > 1. The metadata file describes sample characteristics (e.g., clinical data, demographics) while the taxonomy or feature file contains microbial data (e.g., taxa abundance) used for analysis in microbiome studies.
>    > {: .solution }
>    {: .question}
>
{: .hands_on}


# Find associations between the two files
Now we will find significant associations between microbial features( taxonomy file) and metadata variables (metadata file) using the **MaAslin2** tool

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [MaAsLin 2](toolshed.g2.bx.psu.edu/repos/iuc/maaslin2/maaslin2/1.16.0+galaxy0) %} with the following parameters:
>    - *"Interactions: Fixed effects"*: `c3:age`, `c4:diagnosis`
>    - *"Random effects"*: `c5:subject`
>    - *"Reference"*: `diagnosis,CD`
>
>  Keep rest of the default values as it is.
{: .hands_on}

# Understanding parameters in the tool
Lets now understand the role of each parameter in the tool.

1. **Interactions:Fixed effects** : Fixed effects are the factors in your model that you want to study and draw conclusions about. These are the variables you hypothesize have a direct and consistent influence on the outcome. For example, you are studying how different diets affect gut microbiome composition, then diet would be a fixed effect because you’re specifically interested in understanding how different diets influence the microbiome. You might also include other fixed effects like age and gender to control for their impact.

2. **Random effects** : In some studies, like those following people over time or studying families, samples from the same group can be similar. MaAsLin2 helps handle this by letting researchers choose a grouping factor. This helps make sure the statistical analysis is more accurate. For example, setting random_effects = "Subject_ID" helps control for the correlation between samples that come from the same individual.

3. **Reference** : It allows researchers to establish a baseline or standard category against which other categories are compared, helping to interpret and understand the effects of different variables on microbial features. 

   > <comment-title></comment-title>
   > - In MaAslin2, reference level is must for variables with more than two distinct kind of values.
   > - Reference for a variable with more than two levels is provided as a string of `variable,reference`.
   > - Reference for more than one variable having more than two levels each is provided as a string of `variable1,reference1,variable2,reference2` .
   > - Example, both diagnosis and site variable have more than two levels hence reference can be provided as `diagnosis,CD,site,Cedars-Sinai`.
   {: .comment}

**Additional options** :

4. **min_abundance** [ Default: 0 ] 
- The minimum abundance for each feature within a single sample.
- If a feature's abundance (or presence) is lower than this threshold, it won't be included in the analysis.
- For example, if you set min_abundance to 0.01, only features that make up at least 1% of the total abundance in at least one
   sample will be analyzed. 
- Setting the min_abundance parameter to 0 means that no abundance threshold is applied. In other words, all features (regardless
   of how rare or abundant they are) will be included in the analysis, as long as they are present in the data.
5. **min_prevalence** [ Default: 0.1 ] :
- The minimum proportion of samples in which a feature must be detected to be included in the analysis. 
- It filters features based on their frequency across the entire dataset.
- For example, if min_prevalence is set to 0.1 (10%), only features that appear in at least 10% of the samples are kept,
    regardless of how abundant they are in those samples.
6. **max_significance** [ Default: 0.25 ] :
- When you set a value for max_significance, you are specifying the maximum q-value a feature can have to be
 considered significant.
- A max_significance of 0.25 means you are considering results with a q-value of 0.25 or lower as
 statistically significant.
7. **normalization**:  [ Default: "TSS" ] [ Options: "TSS", "CLR", "CSS", "TMM", "NONE"]
- Ensures that features or samples with different scales or units are brought to a comparable level. 
- Prevents features with larger scales or more abundant counts from dominating the analysis. 
- Helps in making accurate comparisons and interpretations by standardizing data. 
- Options:\
        1. <u> Total Sum Scaling (TSS) </u>: Each count is divided by the total count for that sample, often multiplied by a constant to transform it into a percentage or proportion. 
        2. <u> Centered Log-Ratio (CLR) </u>: Each feature count is divided by the geometric mean of the counts in the same sample, and then the logarithm is taken. Useful for data where ratios between features are of interest, and it helps deal with the compositional nature of microbiome data.
        3. <u> Cumulative Sum Scaling (CSS) </u>: Does the same basic conversion as TSS but it might include extra adjustments to deal with specific data patterns, giving a potentially more accurate normalization.
        4. <u> Trimmed Mean of M-values(TMM) </u>: TMM normalizes data so you can accurately compare gene or feature counts across samples that may have
different total counts or distributions.\
For each feature (like a gene), TMM computes the log-fold change (M-value) between each sample and  a metadata sample.\
It then removes extreme values (outliers) that could skew the results. This trimming helps focus on more  typical values and reduces the impact of any unusual data points.\
Weighted mean of the remaining M-values is calculated to determine the overall adjustment factor for each sample.\
Finally, this adjustment factor is used to normalize the counts in each sample, making them more comparable.\
                                              
8. **transform** [ Default: "LOG" ] [Options: "LOG", "LOGIT", "AST", "NONE" ] 
- The transform to apply to the datasets.
- This is done to make the data more suitable for the linear models used in MaAslin2, helping to improve the accuracy and reliability of the results.
- Options: \
        1. <u>LOG</u> : The log transformation applies the natural logarithm (log base e) to the data. Used When your data has a wide range of values or
is heavily skewed, such as microbiome abundance data where some taxa are much more abundant than others.
        2. <u>LOGIT</u>:   The logit transformation is used for data that represent proportions or probabilities, where the values lie between 0 and 1. It is defined as logit(x) = log(x / (1 - x)). Used when dealing with data that represents proportions, such as relative abundances that are expressed as fractions or percentages.\
The logit transformation is only applicable to data within the open interval (0, 1), so values exactly at 0 or 1 need to be adjusted                     (e.g., adding a small constant like 0.001).
        3.  <u>Arcsine Square Root Transformation (AST)</u> : It is a statistical transformation used primarily on proportion or percentage data.\
The transformation starts by taking the square root of the proportion value. This step reduces the impact of extreme values.\
Next, it applies the arcsine function (the inverse of the sine function) to the square root result. The arcsine function helps to normalize the distribution further.\
Used when you are working with proportion data, such as relative abundances in microbiome studies, where the values are bounded between 0 and 1. 
                                 
9. **analysis_method** [ Default: "LM" ] [ Options: "LM", "CPLM", "ZICP", "NEGBIN", "ZINB" ]
- The analysis method to apply.\
- Options: \
        1. <u>Linear Model (LM)</u>: Determines how changes in metadata are associated with changes in the taxonomy data. 
        2. <u>Compositional Proportional Linear Model (CPLM)</u>: used for analyzing compositional data, where the taxa abundances are proportions or percentages that sum to 1.
        3. <u>Zero-Inflated Count Model (ZCIP)</u>:used when there are many zero counts in the microbiome data. It handles datasets where a large number of taxa are absent in may samples. 
        4. <u>Negative Binomial Model (NEGBIN)</u>: used for count data where there is overdispersion (variance exceeds the mean).
        5. <u>Zero-Inflated Negative Binomial Model (ZIND)</u> : combines features of both zero-inflation and negative binomial models, useful for count data with both excess zeros and overdispersion.   
                                                          
10. **correction or adjustment methods** [ Default: "BH" ] : 
- When performing numerous statistical tests simultaneously, like testing the association of many microbial taxa with various metadata variables, the risk of finding false positives increases. 
- Correction methods help control this risk to ensure that the results are reliable and that significant findings are not due to random chance.
- This is done by computing the q-value,which is a measure of how many false positives are expected among the significant results. 
- Options:\
      1. <u>Benjamini & Hochberg(BH)(aka false discovery rate(fdr))</u>: A common method used for FDR correction. It ranks the p-values from smallest to largest and adjusts them base on their rank and the total number of tests.
      2. <u>Benjamini & Yekutieli(BY)</u> :  Similar to Benjamini-Hochberg but includes a correction factor that accounts for the correlation between tests.
      3. <u> Bonferroni correction</u> :   Divides the significance threshold (alpha level) by the number of tests performed and then compare each p-value to this adjusted significance level to determine if it is statistically significant. 
      4. <u>Hochberg </u>:  It is similar to the Bonferroni correction but is often more powerful, meaning it has better statistical power to detect true effects while controlling for false positives.
      5. <u>Hommel</u>: controls the Family-Wise Error Rate (FWER) by adjusting p-values in a step-down fashion, starting from the smallest p-value and progressively increasing the threshold. It is more powerful than Bonferroni while maintaining strict error control.
      6. <u>Holm</u>:  controls the Family-Wise Error Rate (FWER) by sequentially adjusting p-values from smallest to largest, comparing ch step. This stepwise approach is less conservative than the Bonferroni correction, offering greater statistical power.\
**FWER** is the probability of finding at least one false positive among all the tests performed, assuming all null hypotheses are true.
**FWER Control** is used to minimize the risk of incorrectly claiming significant results when there are none, thus maintaining the overall reliability of the results.\
For more information on correction methods , [click here](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust).

11.  **standardize** : Apply z-score so continuous metadata are on the same scale [ Default: TRUE ]
12.  **plot_heatmap** : Generate a heatmap for the significant associations [ Default: TRUE ]
13. **heatmap_first_n** : In heatmap, plot top N features with significant associations [ Default: 50 ]
14.  **plot_scatter** : Generate scatter plots for the significant associations [ Default: TRUE ]
15.   **cores** : The number of R processes to run in parallel [ Default: 1 ]


# Reading Output Files
The tool generate the following five major files:
- **Data output files**
    1. `residuals.rds`
            This file contains a data frame with residuals for each feature.
    2. `significant_results.tsv`
            Provides the most important output from MaAsLin2 which is the list of significant associations.
    3.  `all_results.tsv`
            Same format as significant_results.tsv, but include all association results (instead of just the significant ones).

- **Visualization output files**
    4.  `heatmap.pdf`
        This file contains a heatmap of the significant associations.
        ![heatmap](../../images/heatmap_maaslin2.png "heatmap of significatn associations")

    5.  `plots :`
        A plot is generated for each significant association.
        Scatter plots are used for continuous metadata.
        Box plots are for categorical data.
        Data points plotted are after normalization, filtering, and transform.

   > <question-title></question-title>
    >
    > 1. Open the heatmap.pdf output file and observe the heatmap. What do you notice?
    > 2. What do you observe in the significant.tsv file?
    >
    > > <solution-title></solution-title>
    > > 1. You will observe the association of the microbes with the fixed effects variables as `age` and `diagnosis`.
    > > - Observe how setting the reference value as `CD` for the categorical variable `diagnosis` in MaAsLin2 implies that this reference level will be used as the baseline for comparison against other levels of the variable, i.e, `nonIBD` and `UC`.
    > > - The effects of other levels will be interpreted relative to this reference level, helping to understand their impact on microbial features.
    > > - The **colors** of the heatmap represent the magnitude and direction of associations between microbial features and metadata variables.
    >       > - **Color Intensity** : The intensity of the color indicates the strength of the association. Darker or more vivid colors usually represent stronger associations.
    >       > - **Color Hue** : The hue (e.g., red, blue) typically indicates the direction of the association. For instance, red represents positive associations (where an increase in the metadata variable is associated with an increase in the microbial feature) and another color blue represents negative associations (where an increase in the metadata variable is associated with a decrease in the microbial feature).\
    >       >
    >       > For example, if you look for `Bifidobacterium longum` in the heatmap, you'll notice that its occurrence in the human gut is least affected by the individual's age and shows a neutral effect in relation to their diagnosis of UC (Ulcerative Colitis) and non-IBD (non-Inflammatory Bowel Disease).
    > > 2. The significant.tsv file shows statistically significant associations between microbial features and metadata variables that meet a specified threshold (in our case, the default `Maximum significance = 0.25`). It includes effect sizes, p-values, and adjusted p-values (q-values) to indicate the strength, direction, and reliability of each association. This file helps identify meaningful relationships in the microbiome data.
    > {: .solution }
    {: .question}


# Studies involving MaAslin2 tool for analysis

- [**Integrating Dietary Data into Microbiome Studies: A Step Forward for Nutri-Metaomics**](https://d1wqtxts1xzle7.cloudfront.net/86457398/pdf-libre.pdf?1653486881=&response-content-disposition=inline%3B+filename%3DIntegrating_Dietary_Data_into_Microbiome.pdf&Expires=1724359088&Signature=gX3tWDORon-KCwFoaPZjJSGjlE6zE5QsQLrEPsB5exvs75mlu5Tk0P9T4lMmXO4Yb-8oVApN9SpM3zLLvchssL99Ps4I5wPri-YN-zwen8tcotQa10KYClxmaELe5VeR3qa-d3WIgu3leoM6rlkjk32eO9sjK3uo6enF~MnxB5yKnfvj2onoou~CrbxA~f712ik~c-E6Q3g6~yhtIqawFElMRhMZvUVMiTRnyeA4U8qI8tRpoT05Ng-plQDkcWOV33pUB8jiqM4I1Qkfltz4TBMfd2APn5X3UtSXvZtU3OVv6eGWWfbU6W1TZ2NU2VCnfg5EBt1iI6yNYEOo84iW4g__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA) : \
The study explores the enhancement of microbiome research through the incorporation of dietary data. The research emphasizes that integrating detailed dietary information with microbiome analyses provides a more comprehensive understanding of how diet influences gut microbiota composition and function. By applying advanced techniques in nutri-metaomics, the study aims to link specific dietary patterns with microbial changes, revealing insights into the interactions between diet, the microbiome, and health outcomes. This approach improves the ability to identify diet-related biomarkers and tailor personalized nutrition interventions based on microbial profiles.\
MaAsLin2 was used to assess how specific dietary patterns influence the abundance and diversity of gut microbiota by integrating detailed dietary data with microbiome profiles. \
MaAsLin2 was set up with the following parameters: \
1. **normalization** : TMM
2. **transform**:  LOG
3. **correction** : BH
4. **analysis_method** : LM
5. **max_significance** : 0.25 (default significance threshold)
6. **min_abundance** :  0.0001
7. **min_prevalence** : 0.1
8. **fixed effects**: Age, gender, and other characteristics of the participants as well as dietary data were added as fixed effects.
9. **random effects** : as participant samples from two timepoints were included, the participant identification number was added as a random effect.\
All models were adjusted for gender. \
Results with a false-discovery rate (FDR) lower than 0.25 were considered significant.

- [**Longitudinal profiling of the intestinal microbiome in children with cystic fibrosis treated with elexacaftor-tezacaftor-ivacaftor**](https://journals.asm.org/doi/full/10.1128/mbio.01935-23) :\
The study investigated how the intestinal microbiome of children with cystic fibrosis (CF) evolves over time while being treated with the elexacaftor-tezacaftor-ivacaftor (ETI) combination therapy. MaAsLin2 was employed to detect shifts in microbial abundance associated with this treatment. By performing longitudinal microbiome profiling, the research tracked changes in gut microbial diversity and composition in response to ETI treatment. The findings highlighted significant shifts in the microbiome, which could impact gut health and inform future CF treatment approaches. 
MaAsLin2 was set up with the following parameters:

1. **fixed effects**: treatment(ELX/TEZ/IVA), age, and recent antibiotic exposure
2. **random effects**: Subject ID was specified as a random effect due to multiple samples from the same subject
3. **min_prevalence**: The minimum prevalence threshold was set to 0.1, indicating that features must be present in at least 10% of the samples to be included.
4. **transform**: LOG transformed
5. **Analysis method** : The general linear “LM” model was used.
6. **Correction method** : The Benjamini-Hochberg procedure was used to correct P values
7. **Normalization method**:A Centered Log-Ratio (CLR) normalization approach was used instead of default normalization methods



- [**The infant gut resistome is associated with E. coli and early-life exposures** ](https://link.springer.com/article/10.1186/s12866-021-02129-x):\
The study investigated how the infant gut resistome—the collection of antibiotic resistance genes (ARGs) in the gut microbiome—associates with E. coli and early-life exposures using MaAsLin2. The analysis, which utilized additive boosting of generalized linear models for feature reduction, revealed significant associations between ARGs and E. coli presence, as well as early-life factors such as antibiotic use and other exposures. Key parameters included CLR normalization of compositional abundance data, no standardization of continuous variables, and a strict significance threshold (q-value < 0.01) using Benjamini-Hochberg correction. This approach highlighted how early exposures influence the resistome and its relationship with E. coli in the infant gut.


# Conclusion

In essence, uncovering associations between microbial features and metadata variables through tools like MaAsLin2 not only deepens our understanding of microbiome dynamics but also holds promise for clinical applications, personalized health strategies, and advancing the field of microbiome research.

Hurray! You have successfully completed the tutorial. Now try to run the tool with some other parameter values and have fun :)
