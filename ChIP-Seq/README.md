ChIP-seq data analysis
======================

ChIP-sequencing is a method used to analyze protein interactions with DNA.
Here, we propose some materials to learn how to analyze ChIP-seq data.

# Slides

Several deck of slides are available for this topic:

- [General introduction about ChIP-seq data analysis](http://bgruening.github.io/training-material/ChIP-Seq/slides)

# Tutorials

A tutorial with hands-on is available for this topic:

- [Identification of the binding sites of the T-cell acute lymphocytic leukemia protein 1 (TAL1) with ChIP-sequencing](tutorials/TAL1_binding_site_identification.md.md)

## Input datasets

The input datasets for the tutorial are available on [Zenodo](https://doi.org/10.5281/zenodo.197100
).

## Galaxy instance

For these tutorials, you can use the [dedicated Docker image](docker/README.md):

```
docker run -d -p 8080:80 bgruening/galaxy-chip-seq-training
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080 ](http://localhost:8080).

# References

## Papers

**Stephen G. Landt et al:** [ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia](http://genome.cshlp.org/content/22/9/1813.long)

> A very useful "encyclopedic" paper with many details about the tools the (mod)ENCODE consortia use and also contains a long section about antibody validation etc..

**Gabriel E Zentner and Steven Henikoff:** [Surveying the epigenomic landscape, one base at a time](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-250)

> Overview of popular sequencing techniques with very nice descriptions of DNase-seq, MNase-seq, FAIRE-seq.

**Benjamin L Kidder et al:** [ChIP-Seq: technical considerations for obtaining high-quality data](http://www.nature.com/ni/journal/v12/n10/abs/ni.2117.html)

> Nice, readable introduction into all aspects of ChIP-seq experiments (from antibodies to cell numbers to replicates to data analysis)

**Marion Leleu et al:** [Processing and analyzing ChIP-seq data](http://bfg.oxfordjournals.org/content/9/5-6/466)

> Fairly detailed review of key concepts of ChIP-seq data processing (less detailed on analysis)

**Peter J. Park:** [ChIP-seq: Advantages and challenges of a maturing technology](http://www.nature.com/nrg/journal/v10/n10/full/nrg2641.html)

**Peter V Kharchenko et al:** [Design and analysis of ChIP-seq experiments for DNA-binding proteins](http://www.nature.com/nbt/journal/v26/n12/full/nbt.1508.html)

**Edison T Liu et al:** [Q&A: ChIP-seq technologies and the study of gene regulation](http://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-56)

> Short overview of several (typical) issues of ChIP-seq analysis

**Thomas S. Carroll et al:**  [Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data](http://journal.frontiersin.org/article/10.3389/fgene.2014.00075/full)

**Shirley Pepke et al:** [Computation for ChIP-seq and RNA-seq studies](http://www.nature.com/nmeth/journal/v6/n11s/full/nmeth.1371.html)

> First comparison of peak callers, focuses on the explanation of basic principles of ChIP-seq data processing and general workflows of peak calling algorithms

**Elizabeth G. Wilbanks & Marc T. Facciotti** [Evaluation of Algorithm Performance in ChIP-Seq Peak Detection](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011471)

> Another comparison of peak callers - focuses more on the evaluation of the peak callers performances than Shirley Pepke et al.

**Mariann Micsinai et al:** [Picking ChIP-seq peak detectors for analyzing chromatin modification experiments](http://nar.oxfordjournals.org/content/40/9/e70.full)

> How to choose the best peak caller for your data set - their finding: default parameters, surprisingly, yield the most reproducible results regardless of the data set type

**Jianxing Fen et al:** [Identifying ChIP-seq enrichment using MACS](http://www.nature.com/nprot/journal/v7/n9/abs/nprot.2012.101.html)

> How to use MACS

**Yong Zhang et al:** [Model-based Analysis of ChIP-Seq (MACS)](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137)

> Original publication of MACS

**Modan K Das & Ho-Kwok Dai:** [A survey of DNA motif finding algorithms](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-S7-S21)

> Review of motif analysis tools

**Philip Machanick and Timothy L. Bailey:** [MEME-ChIP: motif analysis of large DNA datasets](http://bioinformatics.oxfordjournals.org/content/27/12/1696.short)

> MEME-ChIP-paper

**Timothy L. Bailey and Philip Machanick:** [Inferring direct DNA binding from ChIP-seq](http://nar.oxfordjournals.org/content/40/17/e128)

> Centrimo: position-specific motif analysis, especially useful for ChIP-seq data

**Morgane Thomas-Chollier et al:** [Transcription factor binding predictions using TRAP for the analysis of ChIP-seq data and regulatory SNPs](http://www.nature.com/nprot/journal/v6/n12/abs/nprot.2011.409.html),

> How to use TRAP

**Helge G. Roider et al:** [Predicting transcription factor affinities to DNA from a biophysical model.](http://bioinformatics.oxfordjournals.org/content/23/2/134.short)

> Theoretical background of TRAP

# Contributors

This material is maintained by:

- Mallory Freeberg
- Mo Heydarian

For any question related to this topic and the content, you can contact them.

The following individuals have contributed to this training material:

- Mallory Freeberg
- Mo Heydarian
- Friederike Dündar
- Bérénice Batut
