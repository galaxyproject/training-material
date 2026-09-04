---
layout: tutorial_hands_on

title: "Transcriptome Complexity Analysis with TSENAT"

tags:
    - rna-seq
    - isoform
    - entropy
    - transcriptomics
    - complexity
    - information-theory

level: Advanced
time_estimation: 6h

questions:
    - What is Tsallis entropy and how does it relate to isoform complexity?
    - How can I quantify transcriptome complexity changes between conditions?
    - What are scale-dependent diversity measures and why do they matter?
    - How do I detect isoform switching at different scales?
    - What's the difference between Shannon entropy and Tsallis entropy?

objectives:
    - Understand the mathematical foundations of Tsallis entropy
    - Explain the role of the entropic index (q) as a sensitivity parameter
    - Describe the concept of scale-dependent isoform complexity
    - Interpret q-curves and identify complexity patterns
    - Understand the advantages of entropy-based analysis over traditional abundance measures
    - Design and interpret scale-adaptive interaction models
    - Identify scale-dependent isoform switching in RNA-seq data

key_points:
    - Tsallis entropy captures isoform complexity across multiple scales via the q parameter
    - Traditional RNA-seq tools miss isoform reorganization without abundance changes
    - Scale-adaptive interaction models detect condition-specific diversity patterns
    - The q-spectrum reveals which biological scales are affected (rare vs abundant isoforms)
    - Entropy-based approaches complement abundance-focused tools like DESeq2

follow_up_training:
    -
        type: "internal"
        topic_name: transcriptomics
        tutorials:
            - ref-based
            - rna-seq-reads-to-counts

contributions:
  authorship:
    - gallardoalba

lang: en

---

This tutorial covers the theory and practical application of **TSENAT (Tsallis Entropy Analysis Toolbox)**, a  R package for quantifying and analyzing isoform-usage complexity in RNA-seq data using information-theoretic principles.


> <agenda-title></agenda-title>
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background: The Isoform Complexity Problem

## Why Transcriptome Complexity Matters

Standard RNA-seq analysis measures whether **transcript abundance** changes between conditions:
- DESeq2 and edgeR ask: "Do genes produce more or less RNA?"
- These tools focus on **total gene abundance**

However, a critical complementary question remains largely unexplored:
- **How does the diversity of isoforms change?**

A gene may show **little change in total abundance** while dramatically **reshuffling its isoform repertoire**—a phenomenon often missed by traditional methods. This **isoform switching** reflects strategic shifts in protein function driven by splicing regulation, and represents genuine biological signal invisible to abundance-focused analyses.

## Motivation: Evidence from the Literature

### Cancer Heterogeneity and Network Entropy

Recent evidence demonstrates that transcriptome complexity constitutes a measurable biological phenomenon:

- **Tarabichi et al. (2013)**: {% cite tarabichi2013systems %} Cancer cells show "increased entropy of signaling and gene interaction networks" as cancer progresses, contributing to cellular heterogeneity and evolvability
- **Nijman (2020)**: {% cite nijman2020perturbation %} Cancer-associated perturbations increase network entropy, driving phenotypic heterogeneity through disruption of gene regulatory networks
- **Cao et al. (2017)**: {% cite cao2017comprehensive %} Single-cell transcriptomics reveals that transcript-level complexity varies systematically across cell types and developmental states

**The biological principle**: Changes in isoform organization—not just abundance—are valid biological signals worth measuring.

## Abundance-Based Tools vs. Entropy-Based Tools

| Aspect | DESeq2/edgeR | DRIMSeq | TSENAT |
|--------|--------------|---------|--------|
| **What it measures** | Total gene abundance | Differential transcript usage | Isoform heterogeneity/complexity |
| **Detection capability** | Gene-level fold changes | Transcript switching | Reorganization of isoform landscape |
| **A gene that doubles total count but reshuffles isoforms** | ✓ Detected as DE | Maybe detected | ✓ Detected as complexity change |
| **A gene that maintains count but remixes isoforms** | ✗ Missed | Maybe detected | ✓ Detected as complexity change |
| **Interpretability** | Fold-change magnitude | Proportional shift | Diversity pattern (rare vs abundant) |

---

# Theoretical Foundations: From Shannon to Tsallis Entropy

## Historical Development

### Claude Shannon (1948): Birth of Information Theory

Claude Shannon's landmark 1948 work established that **information content could be quantified mathematically**. Shannon entropy measures the uncertainty when drawing a single observation from a probability distribution:

$$H(p) = -\sum_{i=1}^{n} p_i \log_2(p_i)$$

Where:
- $p_i$ is the probability of outcome $i$
- Higher entropy = more unpredictable outcome
- In transcriptomics: higher entropy = more evenly distributed isoforms

**Key insight**: If you draw a transcript from a gene's isoform pool, how predictable is the outcome?

### Limitations of Shannon Entropy

Shannon entropy treats all outcomes equally, regardless of frequency:
- **Rare isoforms** and **common isoforms** receive identical weight
- Cannot be tuned to emphasize one scale over another
- Cannot separately examine complexity driven by rare variants vs. dominant isoforms

### Generalized Entropy Families

To address this limitation, mathematicians developed parametric entropy families:
- **Rényi entropy** (1961): Introduced a tunable parameter α {% cite van2014renyi %}
- **Tsallis entropy** (1988): Alternative generalization with parameter q {% cite tsallis2017foundations %}
- **Hill numbers** (2010): Modern ecological reinterpretation unifying all approaches {% cite chao2010phylogenetic %}

#### Remarkable fact

Despite appearing different mathematically, Rényi and Tsallis entropies can be unified through generalized logarithmic and exponential functions—both answer the fundamental question: **What organizational scales matter?**

## Mathematical Definition of Tsallis Entropy

For a discrete probability vector **p** = (p₁, …, pₙ) representing isoform proportions within a gene, Tsallis entropy is defined as:

$$S_q(p) = \frac{1 - \sum_{i=1}^{n} p_i^q}{q - 1}$$

### Properties and Special Cases

#### The q parameter acts as a sensitivity dial

- **Limit as q → 1**: Recovers Shannon entropy (balanced sensitivity)
- **q = 0**: Richness (number of distinct isoforms present)
- **q = 2**: Simpson index (emphasizes dominant isoforms)
- **q < 1** (e.g., 0.5): Emphasizes **rare isoforms**
- **q > 1** (e.g., 2.0): Emphasizes **abundant isoforms**

### Biological Interpretation: Richness and Evenness

**True diversity** decomposes into two independent components:

1. **Richness**: How many distinct isoforms exist
2. **Evenness**: How evenly transcripts are distributed

Two genes can have identical Shannon entropy yet differ dramatically:

| Scenario | Isoforms | Shannon Entropy | Interpretation |
|----------|----------|-----------------|-----------------|
| **Gene A** | 1 dominant + 99 rare | High (many isoforms contribute) | Heterogeneous: complexity across scales |
| **Gene B** | 100 equally abundant | High (same total entropy) | Uniform: isoforms evenly distributed |

By varying q, TSENAT reveals these differences: Gene A shows **high entropy at low q, declining at high q**, while Gene B shows **consistent entropy across all q-values**.

---


# Quality Control Approaches

## M-Estimation for Sample Influence

Before group-level comparisons, assess whether individual samples exert **disproportionate influence** on entropy estimates.

### Leave-One-Out Influence Assessment

1. **Robust M-estimation**: Iteratively re-weighted least squares with Huber loss
2. **Design**: Quantifies how much sample removal affects location estimates
3. **Independence from group assignment**: Sample-level QC at biological level
4. **Interpretation**: Samples with high influence scores may indicate:
   - Technical quality issues
   - Subject-specific outliers
   - Genuine biological heterogeneity

### Flagging Criteria

Samples warrant review if:
- Distance from centroid > 1.5
- Proportion of affected genes > 0.85
- Outlier status flagged in QC analysis

---

# The Q-Spectrum: Visualization Framework

## What is a Q-Curve?

A **q-curve** plots Tsallis entropy across the entropic index spectrum from q=0 (rare-emphasizing) to q=∞ (abundant-emphasizing). 

### Interpretation by Shape

Three fundamental patterns emerge:

1. **Flat curve (Balanced)**
   - Entropy constant across q-values
   - All isoforms shift proportionally
   - Suggests systematic rebalancing
   - Example: developmental transition

2. **Declining curve (Rare-driven)**
   - D(q=0.5) >> D(q=2.0)
   - Low-abundance isoforms are condition-specific
   - Treatment preferentially expresses rare transcripts
   - Example: exploratory isoform expression

3. **Rising curve (Abundant-driven)**
   - D(q=0.5) << D(q=2.0)
   - High-abundance isoforms shift between conditions
   - Treatment remodels dominant transcript landscape
   - Example: functional specialization

### Multi-Scale Analysis

Different biological processes prioritize different organizational scales:

- **Isoform switching** (discrete change): Detected at specific q-values
- **Wholesale reorganization** (widespread change): Detected at all q-values
- **Cryptic isoform activation**: Revealed at low q (rare scale)
- **Core architecture change**: Revealed at high q (abundant scale)

---

# Scale-Adaptive Interaction Models (SAIT)

## The Core Question

For each gene, we have a **q-curve**: entropy trajectory across entropic scales. We want to test: **Does this curve differ between groups?**

Statistically, this is a **q × condition interaction**: the effect of condition on entropy depends on which q-value you examine.

---

# Transcript Switching Analysis

## The Challenge: From Genes to Transcripts

SAIT tells us **which genes** have scale-dependent complexity changes. But:
- Which **specific transcripts** drive those changes?
- Do the same transcripts switch across all scales?
- Or do different regulatory mechanisms dominate at rare vs. abundant scales?

## Jackknife Resampling Approach

TSENAT uses **leave-one-out jackknife resampling** to identify robust transcript switches:

### Delta Influence Metric

For each transcript:
- **Positive delta influence**: More important in control condition
- **Negative delta influence**: More important in treatment condition
- **Magnitude**: Robustness across bootstrap iterations
- **Visualization**: Color intensity reflects consistency

### Switching Classification

1. **Consistent switching**: Same transcripts show strong delta influence across all q-values
   - Robust, scale-independent splicing shift
   - Suggests coordinated regulatory mechanism

2. **Mixed/scale-dependent switching**: Different transcripts matter at different q-values
   - Regulatory complexity varies by scale
   - Different mechanisms for rare vs. abundant isoforms


---

# Tsallis Divergence: Effect Size Quantification

## Definition

For two probability distributions P (control) and Q (treatment) representing isoform proportions, **Tsallis divergence** measures information-theoretic distance:

$$D_q(P||Q) = \frac{1 - \sum_i p_i^q \cdot q_i^{1-q}}{q - 1}$$

Where:
- D > 0 always (distance measure)
- D = 0 only when distributions identical
- Asymmetry: $D_q(P||Q) \ne D_q(Q||P)$ (information-theoretic property)

## Biological Interpretation

**Effect size threshold**: D > 0.1 indicates meaningful information-theoretic separation between conditions

### Pattern Classification

#### Rare-driven complexity change
- D(q=0.5) >> D(q=2.0)
- Divergence concentrated at low q
- Interpretation: Condition differences driven by rare transcripts

#### Abundant-driven complexity change
- D(q=0.5) << D(q=2.0)
- Divergence concentrated at high q
- Interpretation: Condition differences in dominant isoforms

#### Balanced complexity change
- D(q=0.5) ≈ D(q=2.0)
- Divergence consistent across scales
- Interpretation: All isoforms reorganize proportionally

---

# Method Validation: Non-Parametric Concordance

## Why Two Statistical Methods?

Statistical testing in transcriptomics faces a fundamental trade-off:
- **Parametric methods (GAM)**: Maximum power when assumptions hold, but sensitive to violations
- **Non-parametric methods (SRH)**: Robust to any distribution, but reduced power

### Scheirer-Ray-Hare (SRH) Rank-Based Test

#### Advantages for entropy analysis

1. **Robustness to distribution violations**: Operates on ranks (uniformly distributed), valid for ANY continuous distribution
   - Non-normality: No effect on rank ordering
   - Skewness: Rank transformation corrects automatically
   - Bimodality: Inherently handled

2. **Paired design for q-ordered structure**: 
   - Treats q-measurements as within-subject repeated measures
   - Captures sequential nature of q-spectrum
   - AR(1) autocorrelation doesn't affect rank ordering (exchangeability property)

3. **Interaction testing**:
   - Two-way ANOVA on ranks
   - Tests q × condition effects directly
   - F-statistic follows F-distribution under null

### Concordance Interpretation

High agreement between GAM and SRH:
- ✓ Robust discoveries across statistical frameworks
- ✓ Results likely generalizable

GAM-only discoveries:
- Potential artifacts if GAM assumptions violated
- Alternative: Trust only high-confidence genes significant in both methods

SRH-only discoveries:
- Rare; indicates distributional violations affecting GAM
- May represent genuine signals from highly skewed distributions
