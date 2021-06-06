# Consensus genome characterization
### Call consensus from the aligned BAM file with `ivar consensus` from tool bar.
> 1. Input is the trimmed BAM reads.
> 2. Enter **20** for Minimum quality score threshold
> 3. Enter **0.7** for Minimum frequency threshold
> 4. Enter **20** for Minimum depth
> 5. Select **No** for Exclude regions with smaller depth than the minimum threshold
> 6. Select **Yes** for Use N instead of - for regions with less than minimum coverage
### Rename the FASTA header with `Text transformation with sed`
> Paste this `/^>/s/Consensus_(.)threshold./\1/` sed code in the SED Program window.
