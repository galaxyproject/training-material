---
layout: tutorial_hands_on

title: "Variant calling: Diploid case"
enable: "false"
zenodo_link: ""
questions:
  - "What is the pipeline of the processes for variant calling and processing it?"
objectives:
  - "Identification of the genetic variation using the variant calling"
  - "Using FreeBayes calls for variants generating"
  - "Quering with GEMINI"
time_estimation: "1H"
key_points:
  - "The modern variant callers attempt to assign a reliability estimate for each genotype call. This is done using Bayes reasoning."
  - "FreeBayes variant caller looks at a haplotype window"
  - "After variants have been annotated and post-processed one can manipulate the data with GEMINI queries"
contributors:
  - bebatut
  - torhou
  - nekrut
---

> Much of Galaxy-related features described in this section have been developed by Björn Grüning (@bgruening) and configured by Dave Bouvier (@davebx).

# Introduction
{:.no_toc}

Variant calling is a complex field that was significantly propelled by advances in DNA sequencing and efforts of large scientific consortia such as the [1000 Genomes](http://www.1000genomes.org). Here we summarize basic ideas central to Genotype and Variant calling. First, let's contrast the two things although they often go together:

* **Variant calling** - identification of positions where the sequenced sample is different from the reference sequence (or [reference genome graph](https://github.com/vgteam/vg));
* **Genotype calling** - identifying individual's genotype at variable sites.

A typical workflow for variation discovery involves the following steps (e.g., see Nielsen et al. [2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3593722/)):

1. Mapping reads against the reference genome;
2. Thresholding BAM datasets by, for example, retaining paired, properly mapped reads;
3. Performing quality score recalibration;
4. Performing realignment;
5. Performing variant calling/genotype assignment;
6. Performing filtering and genotype quality score recalibration;
7. Annotating variants and performing downstream analyses.

However, continuing evolution of variant detection methods has made some of these steps obsolete. For instance, omitting quality score recalibration and re-alignment (steps 3 and 4 above) when using haplotype-aware variant callers such as [FreeBayes](https://github.com/ekg/freebayes) does not have an effect on the resulting calls (see Brad Chapman's methodological comparisons at [bcbio](https://bit.ly/1S9kFJN)

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# How does SNP calling and genotyping work?

Before going forward with an actual genotype calling in Galaxy let's take a look as some basic ideas behind modern variant callers.

### SNP calling and genotyping examples:

Consider a set of sequencing reads derived from a diploid individual:

```
REFERENCE: atcatgacggcaGtagcatat
--------------------------------
READ1:     atcatgacggcaGtagcatat
READ2:         tgacggcaGtagcatat
READ3:     atcatgacggcaAtagca
READ4:            cggcaGtagcatat
READ5:     atcatgacggcaGtagc
```

The capitalized position contains a G->A [transition](https://en.wikipedia.org/wiki/Transition_(genetics)). So, in principle this can be a heterozygous site with two alleles **G** and **A**. A commonly used naïve procedure would define a site as *heterozygous* if there is a non-reference allele with frequency between 20% and 80%. In this case **A** is present in 1/5 or 20% of the cases, so we can say that this is a heterozygous site. Yet it is only represented by a single read and is hardly reliable. Here are some of the possibilities that would explain this *variant*. It can be:

* A true variant;
* A library preparation artifact (e.g., a PCR error);
* A base call error;
* A misalignment (though unlikely in the above example).

The modern variant callers attempt to assign a reliability estimate for each genotype call. This is done using Bayes reasoning.
> ### {% icon tip %} Tip: Further reading on
>
> For a great visual explanation see [blog](https://oscarbonilla.com/2009/05/visualizing-bayes-theorem/) by Oscar Bonilla). Here we present a >SNP-relevant "translation" on this explanation (with inspiration from [Erik Garrison](https://github.com/ekg)).
{: .tip}

Suppose you have **A** samples with a variant in a population. You are performing re-sequencing and observe a variant in **B** of your sequencing reads. We want to estimate the probability of having the real polymorphism **A** given our observations in **B**. The logic is as follows:
* The probability of having polymorphism **A** in the population is *P(A) = |A|/|U|*;
* The probability of seeing a variant given our identification approach (i.e., sequencing) is *P(B) = |B|/|U|*;
* Now, the probability of having a variant and it being observed in our sequencing data is the overlap between **A** and **B** sets *P(AB) = |AB|/|U|*.

| Polymorphisms | Variant Calls | Polymorphisms and Variant Calls |
|---------------|---------------|---------------------------------|
|![P(A) polymorphisms](../../images/pA.png)|![P(B) variant calls](../../images/pB.png)|![P(AB) polymorphisms and variant calls](../../images/pAB.png)|

> ### {% icon question %} Questions
>
> 1. What is the probability of a having a real polymorphism **A** given our observation of variants in reads **B**? In other words what is the probability of A given **B**? Or, as stated in the original [blog](https://oscarbonilla.com/2009/05/visualizing-bayes-theorem/): "given that we are in region **B** what is the probability that we are in the region **AB**?"
> 2. Now, let's ask an opposite question. Given a true polymorphism **A** what are the chances that we do detect it (i.e., find ourselves in **AB**)?
>
> > ### {% icon solution %} Solution
> > 1. $$P(A|B) = \frac{\vert AB \vert}{|B|}$$
> >
> >    Dividing by $$\vert U \vert$$: $$P(A \vert B) = \frac{\frac{\vert AB \vert}{\vert U \vert}}{\frac{\vert B \vert}{\vert U \vert}}$$
> >
> >    Because we know that $$P(AB) = \frac{\vert AB \vert}{\vert U \vert}$$ and $$P(B) = \frac{\vert B \vert}{\vert U \vert}$$, we can rewrite the equation in the previous bullet point as $$P(A \vert B) = \frac{P(AB)}{P(B)}$$
> >
> > 2. It will be $$P(B \vert A) = \frac{P(AB)}{P(A)}$$.
> >
> >    So, because we know that $$P(A \vert B) = \frac{P(AB)}{P(B)}$$ and we just reasoned that $$P(B \vert A) = \frac{P(AB)}{P(A)}$$, we can say that $$P(A \vert B)P(B) = P(B \vert A)P(A)$$ leading us to the Bayes formula $$P(A \vert B) = \frac{P(B \vert A)P(A)}{P(B)}$$.
> >
> >    Translating this into "genomics terms" the probability of having a genotype G given reads R is: $$P(G \vert R) = \frac{P(R \vert G)P(G)}{P(R)}$$. Because in a given calculation of $$P(G \vert R)$$ reads are fixed we can re-write the Bayes formula in the following way $$P(G \vert R) \sim P(R \vert G)P(G)$$ with $$P(R)$$ becoming a constant.
> {: .solution }
{: .question}

:exclamation: This leaves us with the need to estimate two things:

1. $$P(R \vert G)$$, the data likelihood;
2. $$P(G)$$, the prior probability for the variant.

In the simplest case we can estimate these as follows:

#### *P(R|G)*

Suppose $$R_{i}$$ is a base in read $$i$$ corresponding to a genome position with genotype $$G$$. The probability of seeing $$R_{i}$$ given $$G$$, $$P(R_{i} \vert G)$$, is given by the quality score of $$R_{i}$$ (the quality scores are given by base calling software and reported as [phred scores](https://en.wikipedia.org/wiki/Phred_quality_score)). Thus the genotype likelihood $$P(R \vert G)$$ is the product of $$P(R_{i} \vert G)$$ over all $$i$$. In reality however there are many other sources of uncertainty (in addition to base qualities) that are incorporated in the calculation of data likelihoods including NGS technology-related issues, dependency of error rates on substitution type (e.g., transitions versus transversions), sequencing context etc...

#### *P(G)* - a single sample case

One can assign an equal probability to all possible genotypes, or to source this information based on previously obtained knowledge containing in a database, such as [dbSNP](https://www.ncbi.nlm.nih.gov/SNP/). In this case (as exemplified in [Nielsen et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3593722/) we may, for instance, have a site with a **G/T** polymorphism and genotypes **GG**, **TT**, and **GT** having frequencies of 0.45, 0.45, 0.09, respectively. We will use these values as priors.

#### *P(G)* - a multi-sample case

Genotype calling reliability can be significantly improved when analyzing multiple samples jointly. In this case genotype frequencies can be inferred from allele frequencies using Hardy-Weinberg equilibrium ([HWE](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle)). The following example (again from [Nielsen et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3593722/)) illustrates this idea: suppose you are calling genotypes for a single individual using a combination of multiple samples. There are two genotypes, **AT** and **AA**, with equally large genotype likelihoods. If, however, in our collection of multiple samples the frequency of **A** is 1% (*p* = 0.01; *q* = 1 - *p* = 0.99), then from the HWE we have:

| AA (*p<sup>2</sup>*)   |   AT (2*pq*) |   TT (*q<sup>2</sup>*)|
|---------|---------|--------|
| 0.0001 | 0.0198 | 0.9801 |

This makes it highly unlikely that **AA** is a true genotype of this individual.

# Calling with FreeBayes

[FreeBayes](https://github.com/ekg/freebayes) is an open source variant caller that has been battle tested by the 1000 Genomes community and is extensively used today (also see [bcbio](https://bcbio.wordpress.com/)). It has a number of features that simply variant discovery workflows. These include (from FreeBayes github page):

> * **Indel realignment is accomplished internally** using a read-independent method, and issues resulting from discordant alignments are dramatically reducedy through the direct detection of haplotypes;
> * **The need for base quality recalibration is avoided** through the direct detection of haplotypes. Sequencing platform errors tend to cluster (e.g. at the ends of reads), and generate unique, non-repeating haplotypes at a given locus;
> * **Variant quality recalibration is avoided** by incorporating a number of metrics, such as read placement bias and allele balance, directly into the Bayesian model;
> * **Ability to incorporate non-diploid cases** such as pooled datasets or data from polyploid samples.

Freebayes is a *haplotype-based* variant caller. This implies that instead of looking at an individual positions within an alignment of reads to the reference genome, it looks at a haplotype window, length of which is dynamically determined (see section 3.2. in [FreeBayes manuscript](https://arxiv.org/pdf/1207.3907v2.pdf)):

|Haplotype-based calling |
|------------------------|
|![Haplotype-based calling](../../images/freebayes.png)|
|<sub>Looking at a haplotype window makes misalignments tolerable. In this case a low complexity poly(A) stretch is misaligned. As a result looking at individual positions will result in calling multiple spurious varians. In the case of FreeBayes looking at a haplotype identifies two alleles (this is a diploid example) `A(7)` and `A(6)`, while `A(8)` is likely an error. Image by [Erik Garrison](https://github.com/ekg/freebayes)</sub>|

# Let's try it

## The data

In this example we will perform variant calling and annotation using [genome in the bottle data](https://sites.stanford.edu/abms/content/giab-reference-materials-and-data). Specifically, we will use Ashkenazim Father-Mother-Son trio data from the Personal Genome Project:

* HG002- NA24385 - huAA53E0 (son)
* HG003- NA24149 - hu6E4515 (father)
* HG004- NA24143 - hu8E87A9 (mother)

Yet for a quick tutorial these datasets are way too big, so we created a [downsampled dataset](https://doi.org/10.5281/zenodo.60520). This dataset was produced by mapping the trio reads against `hg19` version of the human genome, merging the resulting bam files together (we use readgroups to label individual reads so they can be traced to each of the original individuals), and restricting alignments to a small portion of chromosome 19 containing the [*POLRMT*](http://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=Graphics&list_uids=5442) gene.

> ### {% icon hands_on %} Hands-on: Variant calling
>
> 1. Import [the data](https://zenodo.org/record/60520/files/GIAB-Ashkenazim-Trio.txt) into your history {% icon tool %}
> 2. Specify the used genome for mapping {% icon tool %}
>     1. Click on **Edit attributes** (pencil icon on right panel)
>     2. Select `Human Feb 2009` on **Database/Build**
>     3. Save it
> 3. Import the reference genome {% icon tool %}
>     1. Go on **Data Libraries** in **Shared data** (top panel on Galaxy's interface)
>     2. Click on **Training Data**
>     3. Select `hg19`
>     4. Click on **Import selected datasets into history** (just below the top panel)
>     5. Import it
>     6. Convert it from 2bit to fasta with **twoBitToFa** from **Convert Formats**
>
{: .hands_on}

## Generating and post-processing FreeBayes calls

> ### {% icon hands_on %} Hands-on: Generating FreeBayes calls
>
> 1. Select **FreeBayes** from **Phenotype Association** section of the tool menu (left pane of Galaxy's interface) {% icon tool %}
> 2. Make sure the top part of the interface looks like shown below. Here we selected `GIAB-Ashkenazim-Trio-hg19` as input and set **Using reference genome** to `hg19` and **Choose parameter selection level** to `5. Complete list of all options`:
>
>   ![FreeBayes input and parameters](../../images/FreeBayes_settings.png)
>
> 3. Scrolling down to **Tweak algorithmic features?** click `Yes` and set **Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output** to `Yes`. This would help us evaluating the quality of genotype calls:
>
>   ![FreeBayes additional features](../../images/freebayes_gq.png)
{: .hands_on}

This produces a dataset in [VCF](http://www.1000genomes.org/wiki/Analysis/variant-call-format) format containing 35 putative variants. Before we can continue we need to post-process this dataset by breaking compound variants into multiple independent variants with **VcfAllelicPrimitives** tool found within **VCF Tools** section. This is necessary for ensuring the smooth sailing through downstream analyses:

> ### {% icon hands_on %} Hands-on: Post-processing
>
> 1. Select FreeBayes output as the input for this tool {% icon tool %}
> 2. Make sure **Maintain site and allele-level annotations when decomposing** and **Maintain genotype-level annotations when decomposing** are set to `Yes`
>
>   ![VcfAllelicPrimitives input and parameters](../../images/vcfallelicprimitives.png)
>
> **VCFAllelicPrimitives** generated a VCF files containing 37 records (the input VCF only contained 35). This is because a multiple nucleotide polymorphism (`TAGG|CAGA`) at position 618851 have been converted to two:
>
> | Before | After |
> |--------|---------|
> |`chr19 618851 . TAGG CAGA 81.7546`<br> | `chr19 618851 . T C 81.7546`<br>`chr19 618854 . G A 81.7546`|
>
{: .hands_on}

## Annotating variants with SnpEff

At this point we are ready to begin annotating variants using [**SnpEff**](http://snpeff.sourceforge.net/SnpEff.html). **SnpEff**, a project maintained by [Pablo Cingolani](https://www.linkedin.com/in/pablocingolani) "*...annotates and predicts the effects of variants on genes (such as amino acid changes)...*" and so is critical for functional interpretation of variation data.

### {% icon hands_on %} Annotating variants

1. Download `hg19` database with **SnpEff Download** {% icon tool %}
2. Launch annotation of your variants with **SnpEff** from **Annotation**, using the downloaded database (reference genome from your history)

SnpEff will generate two outputs: (1) an annotated VCF file and (2) an HTML report. The report contains a number of useful metrics such as distribution of variants across gene features

  ![SnpEff output chart](../../images/snpeff_chart.png)

or changes to codons

  ![SnpEff output codons](../../images/snpeff_codons.png)

## Manipulating variation data with GEMINI

Now that we have an annotated VCF file it is time to peek inside our variation data. [Aaron Quinlan](http://quinlanlab.org/), creator of [GEMINI](http://gemini.readthedocs.org/en/latest/index.html), calls it *Detective work*.

### Loading data into GEMINI

The first step is to convert a VCF file we would like to analyze into a GEMINI database. For this we will use **GEMINI Load** tool from **Gemini** section. GEMINI takes as input a VCF file and a [PED](https://www.cog-genomics.org/plink2/formats#ped) file describing the relationship between samples. In the case of our dataset the PED file looks like this (second imported file):

```
#family_id sample_id            paternal_id          maternal_id         sex phenotype ethnicity
family1    HG004_NA24143_mother -9                   -9                   2  1         CEU
family1	   HG003_NA24149_father -9                   -9                   1  1         CEU
family1	   HG002_NA24385_son	HG003_NA24149_father HG004_NA24143_mother 1  2         CEU
```

> ### {% icon hands_on %} Hands-on: Loading data {% icon tool %}
>
> So let's load data into GEMINI by setting VCF and PED inputs
>
>   ![GEMINI load parameters](../../images/gemini_load.png)
>
>
> This creates a sqlite database. To see the content of the database use **GEMINI_db_info**:
>
>   ![GEMINI_db_info parameters](../../images/gemini_db_info.png)
>
{: .hands_on}

This produce a list of all tables and fields in the database.

### Querying GEMINI database

GEMINI database is queried using the versatile SQL language (more on SQL [here](https://swcarpentry.github.io/sql-novice-survey)). In Galaxy's version of GEMINI this is done using **GEMINI_query** tool. Within this tool SQL commands are typed directly into the **The query to be issued to the database** text box. Let's begin getting information from some of the tables we discovered with **GEMINI_db_info** tool above.

> ### {% icon tip %} Tip: Gemini tutorials
>
> The examples below are taken from "[Intro to Gemini](https://s3.amazonaws.com/gemini-tutorials/Intro-To-Gemini.pdf)" tutorial. For extensive documentation see "[Querying GEMINI](https://gemini.readthedocs.org/en/latest/content/querying.html)".
{: .tip}

> ### {% icon hands_on %} Hands-on: Selecting "novel" variants that are not annotated in dbSNP database
>
> Type `SELECT count(*) FROM variants WHERE in_dbsnp == 0` into **The query to be issued to the database**
>
>     ![GEMINI query](../../images/gemini_query1.png)
>  
> As we can see from output there are 21 variants that are not annotated in dbSNP
>
{: .hands_on}

> ### {% icon hands_on %} Find variants in POLRMT gene
>
> The query `SELECT * FROM variants WHERE filter is NULL and gene = 'POLRMT'` will produce output with very large number of columns. To restrict the number of columns to a manageable set let's use this command: `SELECT rs_ids, aaf_esp_ea, impact, clinvar_disease_name, clinvar_sig FROM variants WHERE filter is NULL and gene = 'POLRMT'` (column definitions can be found [here](https://gemini.readthedocs.org/en/latest/content/database_schema.html))

Output shows variants found within the *POLRMT* gene.

### Querying genotypes

**GEMINI** provides access to genotype, sequencing depth, genotype quality, and genotype likelihoods for each individual (`subjectID`):

#### Querying

* `gt_types.subjectID` - three types of genotype types: `HOM_REF`, `HET`, `HOM_ALT`;
* `gt_quals.subjectID` - genotype quality
* `gt_depths.subjectID` - total number of reads in this subject at position
* `gt_ref_depths.subjectID` -  number of reference allele reads in this subject at position
* `gt_alt_depths.subjectID` - number of alternate allele reads in this subject at position

> ### {% icon question %} Questions
>  
> 1. At how many sites does child have a non-reference allele?
> 2. At how many sites both father and son have non reference alleles?
> 3. List genotypes for father and son where they have non-reference alleles.
>
> > ### {% icon solution %} Solution
> > 1. To answer this question we will use two fields of GEMINI_query interface:
> >
> >    - Type `SELECT * from variants` into **The query to be issued to the database**
> >    - Type `gt_types.HG002_NA24385_son <> HOM_REF` into **Restrictions to apply to genotype values**
> >
> >    Here is an [example](../../images/gemini_query2.png). It will generate this output
> >
> > 2. Typing the same expression `SELECT * from variants` into **The query to be issued to the database** and `(gt_types.HG002_NA24385_son <> HOM_REF AND gt_types.HG003_NA24149_father <> HOM_REF)` into **Restrictions to apply to genotype values** will generate this output.
> > 3. Typing `SELECT gts.HG002_NA24385_son, gts.HG003_NA24149_father from variants` into **The query to be issued to the database** and `(gt_types.HG002_NA24385_son <> HOM_REF AND gt_types.HG003_NA24149_father <> HOM_REF)` into **Restrictions to apply to genotype values** will generate this output.
> {: .solution }
{: .question}

### Using wildcards

Wildcards simply writing SQL expressions when searching across multiple terms. The syntax for genotype filter wilcards is `COLUMN).(SAMPLE_WILDCARD).(SAMPLE_WILDCARD_RULE).(RULE_ENFORCEMENT).`. Let's look at few examples.

> ### {% icon question %} Question
>
> At which variants are every sample heterozygous?
>
> > ### {% icon solution %} Solution
> >
> > Type `SELECT chrom, start, end, ref, alt, gene, impact, (gts).(*) FROM variants` into **The query to be issued to the database** and `(gt_types).(*).(==HET).(all)` into **Restrictions to apply to genotype values**. Here we use wildcards for the query (`(gts.*)` = get genotypes for **all** samples) and genotype filtering (`(gt_types).(*).(==HET).(all)`, the [all operator](https://gemini.readthedocs.org/en/latest/content/querying.html#the-all-operator) implies that want results for **all** afftected individuals). It will generate this output.
> >
> {: .solution }
{: .question}

# Going further
{:.no_toc}

This short tutorial should give you an overall idea on how generate variant data in Galaxy and process it with GEMINI. Yet there is much more to learn. Below we list GEMINI tutorials and links to Galaxy libraries with relevant data.

| Tutorial | PDF | Data | Galaxy History |
|:----------|:-----:|:------:|:----------------:|
| Introduction | [ :notebook: ](https://s3.amazonaws.com/gemini-tutorials/Intro-To-Gemini.pdf) | [:open_file_folder:](https://usegalaxy.org/library/list#folders/F0283ca691a41c352) | [:arrow_forward:]( https://usegalaxy.org/u/aun1/h/gemini-introduction) |
| Identifying *de novo* mutations underlying Mendelian disease | [ :notebook: ](https://s3.amazonaws.com/gemini-tutorials/Gemini-DeNovo-Tutorial.pdf) | [:open_file_folder:](https://usegalaxy.org/library/list#folders/F775008f45cbbf010) | [:arrow_forward:]( https://usegalaxy.org/u/aun1/h/gemini-de-novo-mutations) |
| Identifying autosomal recessive variants underlying Mendelian disease | [ :notebook: ](https://s3.amazonaws.com/gemini-tutorials/Gemini-Recessive-Tutorial.pdf) | [:open_file_folder:](https://usegalaxy.org/library/list#folders/F35b262f5ac8aa63a) | [:arrow_forward:](https://usegalaxy.org/u/aun1/h/gemini-autosomal-recessive) |
| Identifying autosomal dominant variants underlying Mendelian disease | [ :notebook: ](https://s3.amazonaws.com/gemini-tutorials/Gemini-Dominant-Tutorial.pdf) | [:open_file_folder:](https://usegalaxy.org/library/list#folders/F1c4722ad56892a31) | [:arrow_forward:]( https://usegalaxy.org/u/aun1/h/gemini-autosomal-dominant ) |

## How to use these tutorials?

* Right click on the :notebook: symbol and open tutorial in a new browser tab;
* Right click on :arrow_forward: symbol and open Galaxy history in another new browser tab;
* When Galaxy history interface opens you will need to click **Import history** link highlighted with a red border in the following figure:

  ![Importing history](../../images/import_history.png)

* If you have a wide screen arrange browsers tabs side by side:

  ![Arranging browsers side by side](../../images/side-by-side.png)

* Proceed with tutorial. For example, to repeat the following command from GEMINI tutorial:

  ![GEMINI command](../../images/gemini_command.png)

* Use Galaxy's **GEMINI_load** tool:

  ![GEMINI load tool queries](../../images/galaxy_command.png)

* and so on....


> ### {% icon tip %} Tip: If things don't work...
>
> ...you need to complain. Use [Galaxy's BioStar Channel](https://usegalaxy.org/biostar/biostar_redirect) to do this.
{: .tip}
