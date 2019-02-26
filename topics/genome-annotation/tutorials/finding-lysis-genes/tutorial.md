---
layout: tutorial_hands_on
topic_name: genome-annotation
tutorial_name: finding-lysis-genes
---

> ### Agenda
>
> 1. A Background on Lysis Genes
>    > * Genes Involved
>    > * Expected Genetic Context
> 2. Finding the Endolysin
>    > * Conserved Domain (InterPro)
>    > * BLAST Hits
>    > * SAR Endolysins
>    > * Types of Endolysins
> 3. Finding Holins and Antiholins
>    > * TMHMM Tracks - Transmembrane Prediction
>    > * BLAST Homology
>    > * Basic Considerations When Looking for Holins
>    > * Annotating Holoin/Antiholin Pairs
> 4. Spanin Genes
>    > * Gram-Negative Hosts
>    > * Spanin Architectures: i-spanin/o-spanin, u-spanin
>    > * What is an SPII signal? What is a lipobox?
> 5. Missing Lysis Genes
> {:toc}
>
{: .agenda}

# A Background on Lysis Genes

> ### {% icon question %} What are the different bacterial cell envelope structures?
>    > ### {% icon solution %}
>    > ![](../../images/finding-lysis-genes-screenshots/1_bacterial_cell_envelopes.png)
> {: .solution}
{: .question}

## Genes Involved

![](../../images/finding-lysis-genes-screenshots/2_holin_triggering.png)

![](../../images/finding-lysis-genes-screenshots/3_holin_endolysin_pathway.png)

![](../../images/finding-lysis-genes-screenshots/4_pinholin_sar_endolysin_pathway.png)


> ### {% icon tip %} Gram-positive & Mycolata Bacteria
>
{: .tip}

> ### {% icon comment %} U-Spanins
>
{: .comment}

![](../../images/finding-lysis-genes-screenshots/5_spanins_trapped.png)

![](../../images/finding-lysis-genes-screenshots/7_endolysin_activity.png)

![](../../images/finding-lysis-genes-screenshots/10_spanins_oligomerize.png)

![](../../images/finding-lysis-genes-screenshots/11_membrane_disruption.png)

## Expected Genetic Context

![](../../images/finding-lysis-genes-screenshots/12_lysis_cassette_not_conserved.png)

> ### {% icon comment %} Distributed Lysis Genes
>
> ![](../../images/finding-lysis-genes-screenshots/14_distributed_lysis_genes.png)
>
{: .comment}

![](../../images/finding-lysis-genes-screenshots/15_holin_antiholin_pairs.png)

# Finding the Endolysin

## Conserved Domain (InterPro)

![](../../images/finding-lysis-genes-screenshots/21_endolysin_interpro.png)

## BLAST Hits

## SAR Endolysins

## Types of Endolysins

> * **Glycosidase** (T4 E)

> * **Transglycosidase (lambda R)

> * **Amidase** (T7 3.5)

> * Endopeptidase (T5 Lys)

# Finding Holins and Antiholins

## TMHMM Tracks - Transmembrane Prediction

![](../../images/finding-lysis-genes-screenshots/16_tmhmm_tracks.png)

> ### {% icon tip %} An Important Note
> If the holin and antiholin are **NOT NEXT TO** the endolysin, then it can only be identified if there is *ONLY ONE* small TMD-containing protein in the entire genome, or via BLAST homology.
{: .tip}

## BLAST Homology

## Basic Considerations When Looking for Holins

> * Holins usually have 1-4 TMDs

![](../../images/finding-lysis-genes-screenshots/6_holin_classes.png)

> * Holins are small (70 - 220 amino acids)

## Annotating Holin/Antiholin Pairs

![](../../images/finding-lysis-genes-screenshots/17_holin_antiholin_dual_starts.png)

# Spanins

## Gram-Negative Hosts

The [CPT](https://cpt.tamu.edu/) has curated a [spanin database](https://cpt.tamu.edu/spanindb/#/phages) containing information on hundreds of annotated spanins.

## Spanin Architectures: i-spanin/o-spanin, u-spanin

> * **i-spanin**
>    > * 1 N-terminal TMD
>    > * Small (100-200 amino acids)
>    > * Functionally linked to *o-spanins*

> * o-spanin
>    > * N-terminal SPII signal
>    > * 1 N-terminal lipobox
>    > * Genetically starts completely/partially embedded or just downstream of the i-spanin

> ### {% icon comment %} An Important Note on Spanin Annotation
> The i-spanin is **never** embedded in the o-spanin.
{: .comment}

![](../../images/finding-lysis-genes-screenshots/13_different_spanin_contexts.png)

> * u-spanin
>    > * 1 N-terminal SPII signal
>    > * 1 C-terminal TMD 

![](../../images/finding-lysis-genes-screenshots/9_lambda_t1_spanin_comparison.png)

## What is an SPII signal? What is a lipobox?

![](../../images/finding-lysis-genes-screenshots/18_spi_signal.png)

![](../../images/finding-lysis-genes-screenshots/19_lipoprotein.png)

![](../../images/finding-lysis-genes-screenshots/20_spii_signal_lipobox)

# Missing Lysis Genes

