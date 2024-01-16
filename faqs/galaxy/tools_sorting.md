---
title: Sorting Tools
area: tools
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

Sometimes input errors are caused because of non-sorted inputs. Try using these:
- **Picard SortSam:** Sort SAM/BAM by coordinate or queryname.
- **Samtools Sort:** Alternate for SAM/BAM, best when used for coordinate sorting only.
- **SortBED order the intervals:** Best choice for BED/Interval.
- **Sort data in ascending or descending order:** Alternate choice for Tabular/BED/Interval/GTF.
- **VCFsort:** Best choice for VFC.
- **Tool Form Options for Sorting:** Some tools have an option to sort inputs during job execution. Whenever possible, sort inputs before using tools, especially if jobs fail for not having enough memory resources.