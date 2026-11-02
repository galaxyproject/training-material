---
title: What's the difference between Tsallis entropy and Shannon entropy?
area: conceptual
layout: faq
box_type: tip
contributors: [gallardoalba]
---

**Shannon entropy** (q=1 in Tsallis framework) treats all isoforms equally. **Tsallis entropy** adds a parameter `q` that acts as a sensitivity dial:
- When q < 1: emphasizes **rare isoforms**
- When q = 1: recovers Shannon entropy (balanced)
- When q > 1: emphasizes **abundant isoforms**

This allows you to examine complexity at different biological scales in a single analysis.
