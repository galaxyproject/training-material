---
layout: tutorial_hands_on

title: Identification of RNA regions protected by non-ribosomal protein complexes
zenodo_link: 'https://figshare.com/s/4cd0d0cd4c81705fdf92'
questions:
- How to identify RNA fragments from non-ribosomal protein-RNA complexes?
objectives:
- Learning the method to identify RNA fragments from non-ribosomal protein-RNA complexes
time_estimation: '1H'
key_points:
- RNA fragments from non-ribosomal protein-RNA complexes may have potential functions in the cell, but we need to remove them when processing Ribo-Seq data analysis.
contributors: 
- ldyang14
- IceApink
---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

Sequencing reads from Ribosome profiling are not all from actively translated regions of RNA. Because it is inevitable that non-ribosomal protein-RNA complexes still remained during the process of extracting ribosome-protected fragments. However, it is not that these reads are useless and need to be discard. On the contrary, these reads may come from some potential regulatory regions such as lncRNAs, miRNAs. Hence, we introduce how to identify RNA regions protected by non-ribosomal protein complexes (the region enclosed by blue circle in figure below) in this tutorial.

![Fragments from non-ribosomal protein-RNA complexes](../../images/foot-non-ribosomes/foot_from_non-ribosomes.png "Fragments from non-ribosomal protein-RNA complexes (cited from {% cite ingolia2019ribosome %})")

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}



# Import data

> ### {% icon hands_on %} Hands-on: Upload data
>
> 1. Create a new history and give it a name
>
>    {% include snippets/create_new_history.md %}
>    {% include snippets/rename_history.md %}
>
> 2. Upload the required files
>
>    {% include snippets/import_via_link.md %}
>
>    - the links of these files are shown below
>
>      ```
>      https://ndownloader.figshare.com/files/20028590?private_link=4cd0d0cd4c81705fdf92
>      https://ndownloader.figshare.com/files/20028587?private_link=4cd0d0cd4c81705fdf92
>      ```
>    {% include snippets/import_from_data_library.md %}
{: .hands_on}

# Detecting non-ribosomal reads 

The difference between non-ribomsome protein-RNA complexes and ribosome complexes is shown below. Reads from ribosome-RNA complexes can be observed obviously 3-nt periodicity while reads from non-ribosomal protein-RNA complexes are highly localized distribution. 

![Diff non-ribosome and ribosome](../../images/foot-non-ribosomes/diff-of-non-ribosome-foot.png "The difference between non-ribosome protein-RNA complexes and ribosome complexes (cited from {% cite ji2018rfoot %} )")

> ### {% icon hands_on %} Hands-on: Detect reads from non-ribosomes
>
> - Run **Rfoot** {% icon tool %} with following parameters:
>   - {% icon param-collection %} *"input read mapping file in SAM format"*: `sub_RPF_WT_1.sorted.q20.sam`
>   - {% icon param-file %} *"transcript annotation file in genePred format"*: `gencode.v32.annotation.genePred.txt`
>
{: .hands_on}

The results from Rfoot's analysis revealed potential non-ribosomal protein-RNA regions. Some parts of results are shown in below. Each row represents a region from non-ribosomal protein-RNA and columns show some basic information of this region, such as  the start and end site, length, and so on, according to which we can discover potential gene regulatory regions. Moreover, we can explore the regulatory mechanism of translation of these regions and then provide evidence for relavant researches.

| transcriptID       | chrom | strand | start    | end      | length | read.num | max.pos  | max.num | positions                                                    |
| ------------------ | ----- | ------ | -------- | -------- | ------ | -------- | -------- | ------- | ------------------------------------------------------------ |
| ENST00000202773.13 | chr12 | -      | 1.12E+08 | 1.12E+08 | 9      | 88       | 1.12E+08 | 21      | 112408645:1\|112408646:10\|112408647:3\|112408648:7\|112408649:21\|112408650:17\|112408651:17\|112408652:8\|112408653:4\| |
| ENST00000215375.7  | chr19 | +      | 1241812  | 1241814  | 3      | 21       | 1241813  | 14      | 1241812:1\|1241813:14\|1241814:6\|                           |
| ENST00000217133.1  | chr20 | +      | 59024173 | 59024180 | 8      | 15       | 59024174 | 10      | 59024173:2\|59024174:10\|59024175:1\|59024179:1\|59024180:1\| |
| ENST00000219821.9  | chr16 | +      | 19498729 | 19498738 | 10     | 781      | 19498734 | 384     | 19498729:1\|19498730:5\|19498731:11\|19498732:10\|19498733:23\|19498734:384\|19498735:52\|19498736:238\|19498738:57\| |
| ENST00000221265.8  | chr19 | -      | 39391023 | 39391026 | 4      | 13       | 39391025 | 10      | 39391023:1\|39391024:1\|39391025:10\|39391026:1\|            |

# Conclusion

{:.no_toc}

We get a table contained detailed information of RNA fragments from non-ribosomal protein-RNA complexes through `Rfoot`. Then we can further explore functions of them if we interested some specific genomic regions. Besides, we also removed some obstacles for the subsequent analysis.