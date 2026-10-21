---
title: Can I use TSENAT with very long reads (PacBio, Nanopore)?
area: troubleshooting
layout: faq
box_type: tip
contributors: [gallardoalba]
---

Conceptually yes, but:
- TSENAT expects pseudoalignment-based quantification (SALMON, Kallisto)
- Long reads typically use alignment-based tools (Minimap2, etc.)
- Workflow would need adaptation to parse alignment-based counts
- Overall approach (entropy framework) remains valid
