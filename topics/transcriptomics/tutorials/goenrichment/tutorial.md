---
layout: tutorial_hands_on

title: "GO Enrichment Analysis"
zenodo_link: "https://zenodo.org/record/1255038#.Wx4qTBwh3CI"
questions:
  - "How can I functionally interpret a list of genes of interest that I obtained from my experiment?"
objectives:
  - "How to perform a GO Enrichment Analysis"
  - "How to interpret and simplify the results"
time_estimation: "1h"
key_points:
  - "The goenrichment tool can be used to perform GO Enrichment analysis"
  - "One needs to be careful when chosing the background population"
  - "There are several methods to simplify the output of the GO Enrichment analysis"
contributors:
  - igcbioinformatics
  - mblue9
---

# Introduction
{:.no_toc}

When we have a large list of genes of interest, such as a list of differentially expressed genes obtained from an RNA-Seq experiment, how do we extract biological meaning from it?

One way to do so is to perform functional enrichment analysis. This method consists of applying statistical tests to verify if genes of interest are more often associated to certain biological functions than what would be expected in a random set of genes. In this tutorial you will learn about enrichment analysis and how to perform it.

**What is the Gene Ontology?** <br/>
The [Gene Ontology](http://www.geneontology.org/) (GO) is a structured, controlled vocabulary for the classification of gene function at the molecular and cellular level. It is divided in three separate sub-ontologies or GO types: biological process (e.g., signal transduction), molecular function (e.g., ATPase activity) and cellular component (e.g., ribosome). These sub-ontologies are structured as directed acyclic graphs (a hierarchy with multi-parenting) of GO terms.

![GO Example](../../images/goenrichment_GOexample.png)
>
**Figure 1** QuickGO - http://www.ebi.ac.uk/QuickGO

The GO Ontology, like other ontologies, are usually coded in the [OBO](http://www.geneontology.org/faq/what-obo-file-format) or the [OWL](http://www.geneontology.org/faq/what-owl-file) formats. It can be downloaded from the [Gene Ontology website](http://geneontology.org/page/download-ontology) or from the [OBO foundry](http://www.obofoundry.org/). You can also find in Galaxy tools that allow you to manipulate and extract information from OBO files, but this is outside the scope of this tutorial.

> ### {% icon comment %} Comment:
> Take note of when and where you obtained your ontology file, as these are constantly being updated.
>
{: .comment}


**What are GO annotations?** <br/>
Genes are associated to GO terms via GO annotations. Each gene can have multiple annotations, even of the same GO type. An important notion to take into account when using GO is that, according to the **true path rule**, a gene annotated to a term is also implicitly annotated to each ancestor of that term in the GO graph. GO annotations have evidence codes that encode the type of evidence supporting them: only a small minority of genes have experimentally verified  annotations; the large majority have annotations inferred electronically based on sequence homology or known patterns.

GO annotations can be obtained from the [Gene Ontology website](http://geneontology.org/page/download-go-annotations), or from species-specific databases. One useful resource to obtain GO annotations is [Ensembl biomart](http://www.ensembl.org/biomart/martview). Again, take note to when and from where you obtained your annotations. For example, if you obtained your data from Ensembl, record the release you used.

> ### Overview
>
> In this tutorial, we will deal with:
>
> 1. Functional Enrichment Analysis
> 2. Methods to simplify the results
> 3. Interpretation of the results
> {:toc}
>
{: .agenda}


# Functional Enrichment Analysis

To perform functional enrichment analysis, we need to have:
- A set of genes of interest (e.g., differentially expressed genes): **study set**
- A set with all the genes to consider in the analysis: **population set** (which must contain the study set)
- **GO annotations**, associating the genes in the population set to GO terms
- The **GO ontology**, with the description of GO terms and their relationships

For each GO term, we need to count the frequency (**k**) of genes in the study set (**n**) that are associated to the term, and the frequency (**K**) of genes in the population set (**N**) that are associated to the same term. Then we test how likely would it be to obtain at least **k** genes associated to the term if **n** genes would be randomly sampled from the population, given the frequency **K** and size **N** of the population.

The appropriate statistical test is the one-tailed variant of Fisher’s exact test, also known as the hypergeometric test for over-representation. When the one-tailed version is applied, this test will compute the probability of observing at least the sample frequency, given the population frequency.  The [hypergeometric distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution) measures precisely the probability of **k** successes in **n** draws, without replacement, from a finite population of size **N** that contains exactly **K** successful objects:
>
![Go Enrichment Formula](../../images/goenrichment_formula.png)

For this first exercise we will use data from [Trapnell et al. 2014](https://www.ncbi.nlm.nih.gov/pubmed/22383036 "Trapnell et al. data"). In this work, the authors created an artificial dataset of gene expression in *Drosophila melanogaster*, where 300 random genes were set (insilico) to be differentially expressed between two conditions.

> ### {% icon hands_on %} Hands-on:
> The data for this tutorial is available at [Zenodo](https://zenodo.org/record/1255038#.Wx4qTBwh3CI) to download. For convenience and reproducibility of results, we already added the GO ontology and annotations in the Zenodo repository.
> 1. **Create a new history**
>
> 2. **Upload to the Galaxy** the following files:
>    - go.obo
>    - drosophila_gene_association.fb
>    - trapnellPopulation.tab
>
>    > ### {% icon tip %} Tip: Upload data to Galaxy [<span style="color:red">[1]</span>](https://galaxyproject.github.io/training-material/topics/introduction/tutorials/galaxy-intro-peaks2genes/tutorial.html)
>    > * **Click** on the upload button in the upper left of the interface.
>    >
>    > ![Galaxy Upload](../../images/goenrichment_galaxy_upload.png)
>    >
>    > * Press **Choose local file** and search for your file.
>    > * Press **Start** and wait for the upload to finish.
>    {: .tip}
>
> 3. **Rename** the *go.obo* file to `GO` and *drosophila_gene_association.fb* file to `GO annotations Drosophila melanogaster`.
>
>    > After you upload the files, and if you press the {% icon galaxy-eye %} (eye) icon of `trapnellPopulation.tab` it should look something like this:
>    > ![Trapnell File](../../images/goenrichment_trapnellFile.png)
> **Figure 2** Trapnell file
>
>
>    > ### {% icon comment %} Comments
>    > The study set represents the differentially expressed genes. These were chosen as having an adjusted p-value for the differential expression test (last column) smaller than a given threshold. In this case, we want to select the genes with an adjusted p-value < 0.05.
>    {: .comment}
>
> 4. **Filter data on any column using simple expressions** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `trapnellPopulation.tab`
>    - {% icon param-text %} *"With following condition"*: `c7 < 0.05`
> ![Galaxy Filter](../../images/goenrichment_galaxyFilter.png)
>
> 5. This generates one file. **Rename** it to `trapnellStudy`.
>
>    > ### {% icon comment %} Comments
>    > Both files have the same type of information, the difference between them being the number of genes, as the genes in the study sample are a subset of the population.
>    {: .comment}
>
> 6. **GOEnrichment** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Gene Ontology File"*: `GO`
>    - {% icon param-file %} *"Gene Product Annotation File"*: `Go Annotations Drosophila melanogaster`
>    - {% icon param-file %} *"Study set File"*: `trapnellStudy`
>    - {% icon param-file %} *"Population set File"*: `trapnellPopulation.tab`
>    - Use the default options for the rest.
> ![GO Enrichment Trapnell](../../images/goenrichment_galaxyTrapnell.png)
>
>
>    > ### {% icon question %} Questions
>    >
>    > What were the results from running GOEnrichment?
>    > <details>
>    >
>    > <summary>Click to view answers</summary>
>    > This will generate 6 files with the respective default names: `goenrichment on trapnellStudy MF Table', `goenrichment on trapnellStudy BP Table', `goenrichment on trapnellStudy MF Table', `goenrichment on trapnellStudy MF Graph', `goenrichment on trapnellStudy BP Graph' and `goenrichment on trapnellStudy CC Graph'. The three Table files list the results of the statistical test for each GO Term, ordered by p-value, and the three Graph files are image files displaying a graph view of the enriched GO terms.
>    > </details>
>    {: .question}
>
>    > ### {% icon comment %} Comments
>    > For each GO term we obtain a p-value corresponding to a single, independent test. Since we are making multiple similar tests, the probability of at least one of them being a false positive increases. Therefore we need to make a correction for multiple testing.
>    {: .comment}
>
>    > ### {% icon question %} Questions
>    >
>    > How many significant terms do we get?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > When we ask how many significant terms, we want to see GO terms that have a p-value < 0.05. According to the results, for Molecular Function we have 5 GO terms, Biological Process we have 43 GO terms and Component Cellular we have 10 GO terms.
>    > </details>
>    {: .question}
>
> If you press the {% icon galaxy-eye %} (eye) icon of the *Molecular Function* file (`MF Trapnell`) you should see something like this:
>
> ![GO Enrichment Trapnell MF](../../images/goenrichment_mfTrapnell.png)
>
>    > ### {% icon question %} Questions
>    >
>    > Did you expect to see significant terms?
>    >
>    > <details>
>    > <summary>Click to view answers</summary>
>    > The ~300 genes should be random, so we wouldn't expect to see any enriched term. Nonetheless we still have significant terms.
>    > </details>
>    {: .question}
{: .hands_on}

> ### {% icon comment %} Comments
>
> Let's go back a little bit, and reopen the `trapnellPopulation.tab` file. If you go through the file, you'll see genes with 'NA' as an adjusted p-value. This means that there are genes in our background population for which the differential expression test was not even performed (usually genes that were not expressed in any sample). These genes are irrelevant for this functional enrichment analysis. We also need to note that the study set does not include those genes!
{: .comment}

Let's remove the irrelevant genes from the background population (`trapnellPopulation.tab`), to see the differences in results.

> ### {% icon hands_on %} Hands-on:
>
> 1. **Filter data on any column using simple expressions** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `trapnellPopulation.tab`
>    - {% icon param-text %} *"With following condition"*: `c7 != 'NA'`
> ![Galaxy Filter NA](../../images/goenrichment_galaxyFilterNA.png)
>
> 2. This generates one file. **Rename** to `trapnellFilteredPopulation`.
> 3. **GOEnrichment** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Gene Ontology File"*: `GO`
>    - {% icon param-file %} *"Gene Product Annotation File"*: `Go Annotations Drosophila melanogaster`
>    - {% icon param-file %} *"Study set File"*: `trapnellStudy`
>    - {% icon param-file %} *"Population set File"*: `trapnellFilteredPopulation`
>    - Use the default options for the rest.
> ![GO Enrichment Trapnell new Population](../../images/goenrichment_galaxyTrapnellNewPop.png)
>
> 4. **Rename** all 6 output files by appending `FilteredPop` to the name, to distinguish them from the previous outputs.
> Let's check the new graph **goenrichment on trapnellStudy MF Graph FilteredPop**.
>
> ![GO Enrichment Trapnell new Population MF](../../images/goenrichment_mfTrapnellNew.png)
>
>
>    > ### {% icon question %} Question
>    >
>    > 1. How many significant terms do we get now?
>    > 2. Why do you see these differences?
>    > <details>
>    >
>    > <summary>Click to view answers</summary>
>    > <ol type="1">
>    > <li> According to the results, in *Molecular Function* and *Biological Process* we have 0 GO terms and *Component Cellular* just 1 GO term. </li>
>    > <li> The background population genes that we removed are not random, they are usually genes that are expressed in specific conditions, tissues or time points. If they are included in the test, we will obtain false enrichments, as we saw.</li>
>    > </ol>
>    > </details>
>    {: .question}
{: .hands_on}


# Simplification of graphs

Graphical views are essential, but sometimes the graph view can become overwhelming due to the size of the results. To exemplify this issue, we will next perform functional enrichment analysis using a more realistic dataset from a study using the mouse model organism. The original dataset can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30352). In this [study](https://www.nature.com/articles/nature10532), the authors compared the gene expression of several tissues. Here, we will use results from the comparison between heart and brain.

> ### {% icon hands_on %} Hands-on:
>
> For the first exercise we will use as a study set the differential genes (padjusted<0.05).
> 1. **Upload to Galaxy** the mouse_brain_vs_heart.txt, Mus_musculus_annotations_biomart_e92.tab and mouse_brain_vs_heart.difgenes.txt files.
> 2. **Rename** the *mouse_brain_vs_heart.txt* file to `Mouse population`, *Mus_musculus_annotations_biomart_e92.tab* file to `GO annotations Mus musculus` and *mouse_brain_vs_heart.difgenes.txt* file to `Mouse diff`.
> 3. **GOEnrichment** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Gene Ontology File"*: `GO`
>    - {% icon param-file %} *"Gene Product Annotation File"*: `GO annotations Mus musculus`
>    - {% icon param-file %} *"Study set File"*: `Mouse diff`
>    - {% icon param-file %} *"Population set File"*: `Mouse population`
>    - Select **'No'** in the Summarize Output option.
> ![GO Enrichment Mouse Diff](../../images/goenrichment_galaxyMouseDiff.png)
>
> 4. This will generate 6 files with the names: `goenrichment on Mouse diff MF Table`, `goenrichment on Mouse diff BP Table`, `goenrichment on Mouse diff CC Table`, `goenrichment on Mouse diff MF Graph`, `goenrichment on Mouse diff BP Graph` and `goenrichment on Mouse diff CC Graph`.
> 5. Analyze the table and graph from *Biological Process*.
>
>
> ![GO Enrichment Mouse Diff BP](../../images/goenrichment_bpMouseDiff.png)
>
>    > ### {% icon comment %} Comments
>    > As you can see the three graphs are very complex and difficult to analyze.
>    {: .comment}
{: .hands_on}

As you see, the number of enriched GO Terms is very high, with graphs that are too extensive to analyze manually. And this is despite the fact that GOEnrichment ignores singletons and skips dependent tests by default, precisely to avoid enrichment results that are too extensive and not informative.

The Summarize Output option in the GOEnrichment tool addresses this problem by conflating branches/families of enriched GO terms and selecting the most representative term(s) from them (usually 1-2 term per family). The greatly simplifies the results while retaining branch information, and thus ensuring that every enriched family of functions is present in the results. Some specificity is necessarily lost, but the trade-off is that the results become easier and more intuitive to analyze.

> ### {% icon hands_on %} Hands-on:
>
> 1. **GOEnrichment** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Gene Ontology File"*: `GO`
>    - {% icon param-file %} *"Gene Product Annotation File"*: `GO annotations Mus musculus`
>    - {% icon param-file %} *"Study set File"*: `Mouse diff`
>    - {% icon param-file %} *"Population set File"*: `Mouse population`
>    - Use the default options for the rest (notice that by default the Summarize Option is on).
> ![GO Enrichment Mouse Diff Summarize](../../images/goenrichment_galaxyMouseDiffSum.png)
> 2. Analyze again the table and graph from *Biological Process*.
>
>
> ![GO Enrichment Mouse Diff Summarize BP](../../images/goenrichment_bpMouseDiffSum.png)
>
>    > ### {% icon question %} Questions
>    > 1. Are there differences in complexity comparing the graph with and without the summarize output option?
>    > <details>
>    >
>    > <summary>Click to view answers</summary>
>    > 1. Yes, there are differences. As you can see, the activation of the Summarize option reduces the size of the graph because this parameter causes families of GO terms to be conflated. Each major branch in the full results is still present in the summarized results, but now is reduced to 1 or 2 most representative terms, leading to a graph that is much easier to interpret while still containing all the key functional information.
>    > </details>
>    {: .question}
{: .hands_on}

Another approach to reduce the complexity of the results is to use a shallower version of GO, the GO slims. GO slims are transversal cuts of GO that cover all key branches but lack specific terms. Thus, using them leads to much simpler results than using the full GO, but also leads to a substantial loss in specificity, which is greater than that of the Summarize Output option. You can download slimmed versions of GO from the [Gene Ontology website](http://geneontology.org/page/go-subset-guide).

To test the GO slim approach, let us use the mouse dataset again. First, however, we need to use GOSlimmer tool to convert the annotations file from full GO to GO slim (as GO annotations are typically made to terms that are too specific to be in the GO slim, and thus need to be extended by the true path rule).

> ### {% icon hands_on %} Hands-on:
>
> 1. **Upload to the Galaxy** the goslim_generic.obo file.
> 2. **Rename** the *goslim_generic.obo* file to `GO Slim`.
> 3. Run **GOSlimmer** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Full Gene Ontology File"*: `GO`
>    - {% icon param-file %} *"GOSlim File"*: `GO Slim`
>    - {% icon param-file %} *"Gene Product Annotation File"*: `GO annotations Mus musculus`
>
> ![Galaxy GO Slimmer](../../images/goenrichment_galaxySlimer.png)
>
> This will generate one file called `Slim Annotations`.
{: .hands_on}
>
>

Now we will go use the GOEnrichment tool with the new Slim Annotations file and the same study set.

> ### {% icon hands_on %} Hands-on:
>
> 1. **GOEnrichment** {% icon tool %} with the following parameters:
>    - Use the `GO Slim`, `Slim Annotations`, `Mouse diff` and `Mouse population` files.
>    - Use the default options.
> ![GO Enrichment Mouse Slim](../../images/goenrichment_galaxyMouseSlim.png)
>
>    > ### {% icon question %} Questions
>    >
>    > 1. What differences do you observe when comparing the results obtained with the GO Slim to those obtained with the full GO, with the Summarize Output option?
>    >
>    > *Component Cellular* with full GO
>    > ![GO Enrichment Mouse Summarize CC](../../images/goenrichment_ccSum.png)
>    > *Component Cellular* with GO Slim
>    > ![GO Enrichment Mouse Slim CC](../../images/goenrichment_ccSlimSum.png)
>    > <details>
>    >
>    > <summary>Click to view answers</summary>
>    > 1. The differences that you observe are due to the ontology used. When we apply the summarize option with the full GO, the GOEnrichment tool will return a summarized output (as we have seen previously). When we opted for GO Slim, the original annotation was already summarized, resulting in an even more summarized output, but with a consequent loss of specificity.
>    > </details>
>    {: .question}
>
{: .hands_on}


# Interpretation of the results
The interpretation of the results will depend on the biological information that we intend to extract. Enrichment analysis can be used in validation (e.g., of a protocol for extracting membrane proteins), characterization (e.g., of the effects of a stress in a organism) and elucidation (e.g., of the functions impacted by the knock-out of a transcription factor).

There is one important point to keep in mind during the analysis: statistically significant is different from biologically meaningful. That said, it is typically possible to obtain some biological or technical insight about the underlying experiment from statistically enriched terms, even if it isn’t readily apparent.

Terms that are very generic tend to be difficult to interpret, because the meaning they convey is shallow. On the other hand, very specific terms are generally not integrative and thus not useful in interpreting a gene set collectively. The interesting terms are those that are sufficiently specific to transmit substantial biological meaning, yet generic enough to integrate multiple genes.

For the second exercice, we will continue to work with the same study set as before but now we analyze separately genes that are over- and under-expressed, and see the enriched GO terms presents in the brain and heart from the mouse.

> ### {% icon hands_on %} Hands-on:
>
> 1. **Upload to Galaxy** the mouseOverexpressed.txt and the mouseUnderexpressed.txt files.
>
>    > ### {% icon comment %} Comments
>    >
>    > The differentially expressed genes can be identified using the adjusted p-value (also known as FDR). The logFC values indicate whether genes are more expressed (logFC>0) or less expressed (logFC<0) in one condition when comparing with another condition.
>    >
>    {: .comment}
>
> 2. **GOEnrichment** {% icon tool %} with the following parameters for the both study files (*mouseOverexpressed.txt* and the *mouseUnderexpressed.txt*).
>    - Use the `GO`, the `GO annotations Mus musculus` and the `Mouse population` files.
>    - Use the default options.
> ![GO Enrichment Mouse Overexpress](../../images/goenrichment_galaxyMouseOver.png)
> ![GO Enrichment Mouse Underexpress](../../images/goenrichment_galaxyMouseUnder.png)
>
> 3. This will generate 12 files, 6 for each sample file, like in previous cases.
>
>    > ### {% icon question %} Questions
>    >
>    > Analyze both Biological Process tables. According to the study, which tissues are over- and underexpressed?
>    > <details>
>    >
>    > <summary>Click to view answers</summary>
>    > The samples correspond to the expressions that occur in the tissues referring to the brain and heart, so the results in the tables (and also in the graphs) should correspond to the specific functions of each organ. When we analyze the tables of enriched functional terms, we can see that the results from underexpressed genes reveal functions related to the brain. While in the case of the genes overexpressed, we identify functions related to muscle / heart function.
>    > </details>
>    {: .question}
{: .hands_on}


# Conclusion
{:.no_toc}
Functional enrichment is a good way to look for patterns in gene lists, but interpretation of results can become a complicated process. One way to reduce this complexity is to use the GOEnrichment tool. This tool not only performs the GO Enrichment test, showing us enriched GO terms from our sets, but also contains functionality to simplify the results and make them more easily interpretable. Independently of this, we need to be careful when choosing our genes of interest, but also the background set of genes against which we want to compare.

