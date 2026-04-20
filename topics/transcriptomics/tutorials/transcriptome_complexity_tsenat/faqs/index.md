# Frequently Asked Questions about TSENAT and Transcriptome Complexity Analysis

## Conceptual Questions

### What's the difference between Tsallis entropy and Shannon entropy?

**Shannon entropy** (q=1 in Tsallis framework) treats all isoforms equally. **Tsallis entropy** adds a parameter `q` that acts as a sensitivity dial:
- When q < 1: emphasizes **rare isoforms**
- When q = 1: recovers Shannon entropy (balanced)
- When q > 1: emphasizes **abundant isoforms**

This allows you to examine complexity at different biological scales in a single analysis.

### Why is transcriptome complexity important if DESeq2 already detects changes?

DESeq2 and edgeR measure **total gene abundance** changes. They miss situations where:
- A gene maintains constant total expression but **reorganizes its isoforms**
- Rare isoforms become more common (or vice versa) without affecting overall counts
- The **diversity** of the isoform landscape changes

These are valid biological signals (isoform switching via splicing regulation) that complement abundance-based methods.

## Technical Questions

### What sample size do I need?

TSENAT assumes **≥6-8 samples per group** for paired/longitudinal designs:
- Smaller samples are underpowered to detect entropy shifts
- Larger samples increase power to detect subtle complexity changes
- Biological replicates are essential (not technical replicates)

### How do I choose the q-spectrum (range and granularity)?

**Standard configuration**: q = 0 to 2.0 with Δq = 0.05 (41 q-values)
- Range: 0-2 covers extreme rare-emphasis to abundant-emphasis
- Granularity: Δq=0.05 provides good resolution without computational burden
- For visualization: Can downsample to q = (0, 0.5, 1, 1.5, 2) for quick exploration

### Does read depth affect entropy estimates?

Yes. TSENAT uses **automatic pseudocount selection** to handle this:
- Pseudocounts scale with library size
- Ensures entropy estimates are comparable across samples with different sequencing depths
- Breaks ties in rank-based validation methods
- Standard practice in RNA-seq analysis

### What's the difference between "consistent" vs. "scale-dependent" transcript switching?

- **Consistent switching**: Same transcripts show high delta-influence values across all q-values
  - Robust, scale-independent splicing shift
  - Suggests coordinated regulatory mechanism
  
- **Scale-dependent switching**: Different transcripts matter at different q-values
  - Regulatory complexity varies by scale
  - Different mechanisms for rare vs. abundant isoforms

## Methodological Questions

### What does "q × condition interaction" mean?

It means the effect of **condition on entropy depends on which q-value** you examine:
- At q=0.5, one condition might have higher entropy
- At q=2, the pattern might reverse
- This reveals **scale-dependent** biological signals

### What's the advantage of Tsallis divergence over fold-change?

**Fold-change**:
- Simple 1:1 ratio between conditions
- Misses complexity patterns
- Unbounded (can be arbitrarily large)

**Tsallis divergence**:
- Information-theoretic distance between distributions
- Captures pattern type (rare-driven, balanced, abundant-driven)
- Symmetric [0] to asymmetric, bounded interpretability
- Directly quantifies isoform complexity differences

## Practical Workflow Questions

### Where do I get transcript-level counts for this analysis?

TSENAT expects **SALMON or Kallisto output** including:
- NumReads: raw fragment counts per transcript
- TPM: length- and library-normalized estimates
- EffectiveLength: read-length corrected transcript length

These are standard outputs from pseudoalignment-based quantification.


## Interpretation Questions

### A gene shows significant q×condition interaction but small effect size—should I report it?

Guidelines for meaningful results:
- **Effect size (Tsallis divergence) D > 0.1**: ~10% information divergence (notable)
- **D > 0.2**: Substantial divergence (recommended threshold)
- **Small D with p < 0.05**: May reflect statistical noise; prioritize high-effect genes

Report both for completeness:
- High-effect genes (D > 0.2): Strong biological signal
- Borderline genes (D ~ 0.1, p < 0.05): Mention as candidates for validation


## Troubleshooting

### My results show no significant genes. What could be wrong?

1. **Sample size too small**: Need ≥6-8 per group (check design)
2. **Entropy differences too subtle**: Try visualizing q-curves manually
3. **Filtering too stringent**: Relaxing filters may reveal additional patterns
4. **Method choice**: Try alternative statistical method (GAM vs. LMM vs. SRH)
5. **Data quality**: Run M-estimation QC; check for outlier samples

### Some q-values give NA or NaN results. Why?

Possible causes:
1. **Zero counts**: Pseudocount selection should handle; check filtering
2. **Numerical instability**: Very small divergence values can underflow
3. **Method-specific**: SRH (rank test) may have issues with many tied values

Solution: Ensure pseudocount regularization is enabled; try alternative statistical method.

### Q-curves show erratic patterns. Is this normal?

**Expected**: Small fluctuations in entropy across q-values
**Concerning**: Large spikes or reversals
- May indicate normalization issues
- Check: Are housekeeping genes (GAPDH, ACTB) showing balanced patterns?
- If not: Revisit pseudocount and filtering parameters

### Can I use TSENAT with very long reads (PacBio, Nanopore)?

Conceptually yes, but:
- TSENAT expects pseudoalignment-based quantification (SALMON, Kallisto)
- Long reads typically use alignment-based tools (Minimap2, etc.)
- Workflow would need adaptation to parse alignment-based counts
- Overall approach (entropy framework) remains valid
