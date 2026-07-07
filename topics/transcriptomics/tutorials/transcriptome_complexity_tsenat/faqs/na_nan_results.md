---
title: Some q-values give NA or NaN results. Why?
area: troubleshooting
layout: faq
box_type: tip
contributors: [gallardoalba]
---

Possible causes:
1. **Zero counts**: Pseudocount selection should handle; check filtering
2. **Numerical instability**: Very small divergence values can underflow
3. **Method-specific**: SRH (rank test) may have issues with many tied values

Solution: Ensure pseudocount regularization is enabled; try alternative statistical method.
