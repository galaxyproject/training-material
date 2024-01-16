---
title: What does it mean to normalize the LFQ intensities?
area: analysis
layout: faq
box_type: question
contributors: [foellmelanie]
---


Median normalization typically refers to subtracting the median of all intensities within one sample from all of the intensities (e.g. Intensity of Protein A - Median of all intensities from Sample 1) , to account for measurement variations. Before normalization log2 transformation is required since many statistical tests demand that the data is actually normal distributed. (Non log intensities show very high values but have a minimum (limit of quantification) leading to a somehow right skewed distribution, after log-transformation the intensity distribution is more like a gaussian distribution. Beside the median (or median-polish) normalization there is also other e.g. the quantile  normalization.

