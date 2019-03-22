---
layout: tutorial_hands_on
topic_name: additional-analyses
tutorial_name: hhpred-analysis-of-proteins
---

> ### Agenda
>
> Background: Using Protein Structure Prediction Tools to Predict Function
>    > * Protein Structure Prediction
>    > * Tools That Predict Protein Structure
>    > * Caveats
>
> Using & Interpreting HHPred with Caution
>
> Annotation in Apollo Based on HHPred Results
>    > * Directions on Naming the Protein
>    > * What to Include in the Notes
>
{: .agenda}

# Background: Using Protein Structure Prediction Tools to Predict Function

## Protein Structure Prediction

## Tools That Predict Protein Structure

> ### {% icon tip %} A Relevant Read
> [ A Completely Reimplemented MPI Bioinformatics Toolkit with a New HHpred Server at its Core. J Mol Biol. 2017 Dec 16.](https://www.ncbi.nlm.nih.gov/pubmed/29258817)
{: .comment}

> ### {% icon comment %} Note that...
> Jobs are stored for 3 weeks.
{: .comment}

## Caveats

# Using & Interpreting HHPred with Caution

> * **Input: Getting your protein sequence from Apollo**

> * **Parameters to run the tool**

> * **Looking at the results**

*Probability* takes into consideration secondary structure, and uses the weighted conservation of amino acids at each position for homology determination for prediction.

*E-value* indicates how many chance hits with a better score than this would be expected in a fully unrelated database (I. E.: Smaller is better, and smaller than one is **absolutely** required).

> ### {% icon tip %} Alignment Guide
> * In the alignments, the amino acids are colored based on their physio-chemical properties (E. G.: positive charge = red, negative charge = blue, aliphatic = green).
>
> * Query-related rows start with 'Q' and target-related ones with 'T'.
>
> * While ss_pred indicates secondary structure (SS) states as predicted by PSIPRED, ss_dssp indicates SS states calculated using DSSP from known structures (e = beta-strand, h = alpha-helix).
>
> * Consensus rows show the conservation pattern of residues in the query and target alignments.
{: .tip}

> * **Sorting through results critically**

HHPred on hypothetical proteins

> * 95% probability = in the bag
> * 50% probability = should be checked for reasonableness
> * 30% probability AND in the top three hits = should be checked for reasonableness

> ### {% icon question %} Questions About Interpreting Results?
>    > ### {% icon solution %}
>    > For answers to questions about looking at the results and more ("What is the meaning of the symbols (+.-|~)?" or, "What is the significance of upper vs. lowercase letters in both the consensus and ss_pred lines?"), see this document: [https://github.com/soedinglab/hh-suite/wiki](https://github.com/soedinglab/hh-suite/wiki)
>    {: .solution}
{: .question}

<!-- NEED CLARIFICATION ON LINK, PREVIOUS LINK LED TO ERROR 404 -->

# Annotation in Apollo Based on HHPred Results

## Directions on Naming the Protein

> ### {% icon comment %} A Reminder of Naming Guidelines
>
{: .comment}

## What to Include in Notes

**HHPred** *(at minimum)*:
> * IDs of hits
> * Scores

