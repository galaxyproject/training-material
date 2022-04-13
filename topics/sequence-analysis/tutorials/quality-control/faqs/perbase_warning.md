---
title: When I get a warning for base per sequence content, what should I do?
box_type: tip
layout: faq
contributors: [heylf]
---

So far it does not mean that your data is bad. Your protocol or your data might have a bias that you normally expect. Check first the following things:
- Adapter content (maybe some adapters are still in your data)
- Kmer content/Over represented sequences (this would indicate a contamination or a protocol/sequence bias)
- Per base quality plot. If the overall quality is not good, then probably the sequencing was poorly performed.
- Read about your protocol, e.g., ChIP-Seq and ATAC-Seq typically have a nucleotide bias. For example this article about [ATAC-Seq](https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-3-3).


