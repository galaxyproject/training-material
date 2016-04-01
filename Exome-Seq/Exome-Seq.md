# GALAXY WORKSHOP on Exome-seq DATA ANALYSIS

Exome sequencing means that all protein-coding genes in a genome are sequenced. In Humans there are ~180,000 exons that makes up 1% of the human genome which contain ~30 million base pairs. Mutations in the exome have usually a higher impact and severe consequences, than in the remaining 99% of the genome. For the identication of genetic variations exome sequencing is cheaper than whole-genome sequencing. With a high coverage rate of 100+ DP, 98% of all exons are covered.

(Also have a look at this tutorial: https://github.com/nekrut/galaxy/wiki/Diploid-variant-calling)


Things exome sequencing can't do:

- Not all genes are in your exon data, especially those buried in stretches of repeats out towards the chromosome tips, aren’t part of exome sequencing chips.
- Mutations in the handful of genes that reside in mitochondria, rather than in the nucleus.
- “Structural variants,” such as translocations and inversions, that move or flip DNA but don’t alter the base sequence.
- Triplet repeat disorders, such as Huntington’s disease and fragile X syndrome can't be detected.
- Other copy number variants will remain beneath the radar, for they too don’t change the sequence, but can increase disease risk.
- Genes in introns. A mutation that jettisons a base in an intron can have dire consequences: inserting intron sequences into the protein, or obliterating the careful stitching together of exons, dropping gene sections. For example, a mutation in the apoE4 gene, associated with Alzheimer’s disease risk, puts part of an intron into the protein.
- “Uniparental disomy.” Two mutations from one parent, rather than one from each, appear the same in an exome screen: the kid has two mutations. But whether mutations come from only mom, only dad, or one from each has different consequences for risk to future siblings. In fact, a case of UPD reported in 1988 led to discovery of the cystic fibrosis gene.
- Control sequences. Much of the human genome tells the exome what to do, like a gigantic instruction manual for a tiny but vital device. For example, mutations in microRNAs cause cancer by silencing various genes, but the DNA that encodes about half of the 1,000 or so microRNAs is intronic – and therefore not on exome chips.
- Epigenetic changes. Environmental factors can place shielding methyl groups directly onto DNA, blocking expression of certain genes. Starvation during the “Dutch Hunger Winter” of 1945, for example, is associated with schizophrenia in those who were fetuses at the time, due to methylation of certain genes. Exome sequencing picks up DNA sequences – not gene expression.
- Gene-gene (epistatic) interactions. One gene affecting the expression of another can explain why siblings with the same single-gene disease suffer to a different extent. For example, a child with severe spinal muscular atrophy, in which an abnormal protein shortens axons of motor neurons, may have a brother who also inherits SMA but has a milder case thanks to a variant of a second gene that extends axons. Computational tools will need to sort out networks of interacting genes revealed in exome sequencing.

## Freiburger Exome-seq pipeline

![This is our Exome-seq pipeline in Galaxy](https://github.com/bgruening/presentations/raw/master/shared/resources/img/genVAST.png)

#Hands on!

1. This time we will start with a few BAM files. You will get one BAM file with all aligned reads from one family. The mother/father are healthy but the child has a yet unknown disease. Import all 3 BAM's into a new history:

[Father](https://github.com/bgruening/training_data/raw/master/Father)  
[Mother](https://github.com/bgruening/training_data/raw/master/Mother)  
[Child = Patient](https://github.com/bgruening/training_data/raw/master/Patient)

2. To call our variants we will use FreeBayes. Please search for the latest version of FreeBayes and have a look at the advanced settings, to get an impression how powerful it can be. Choose basic setting and call all variants from the father.
3. Congratulations, you have created you first VCF file, one of most complicated file formats in bioinformatics. In such a file your called variants are stored. One variant per line (+header lines).
4. Post-Process your data with: *VcfAllelicPrimitives*; Maintain site/allele-level and genotype-level annotations when decomposing.
5. Filter your VCF file with *SnpSift Filter* (version 4.0.0). You only want to have SNPs with a Quality >= 30 and a Coverage >= 10. Have a look at the examples to help you construct the correct expression.

Annotate your variants: 

6. Get the [dbSNP.vcf](https://github.com/bgruening/training_data/raw/master/dbSNP_138.hg19.vcf)
7. At first assign know variants ID's form dbSNP to your variants. Use *SnpSift Annotate* (version 4.0.0) for this task. SnpEff (version 3.4) will annotate your variants with some functional information. Have a look at your INFO column again. You will get some gene names for your variants, but also a predicted impact and if your variant is located inside of a known gene.

At this stage you have your first annotated variants and in theory everything you need for your further studies included in your VCF file. 

8. Please extract your history to a workflow and apply this workflow to the other BAM files.
9. You should now have 3 annotated variant files, from mother, father and the patient. It might be a good idea to rename them accordingly.
10. Combine all 3 files into one with the tool *VCFcombine*
11. Create a pedigree file (PED) with the content given in the table at the bottom.
12. Use the tool *GEMINI load* (version 0.18.1) to create a database out of your combined VCF file and the PED file. This can take some time, especially if we had a bigger vcf file.
13. Import the new 
[GEMINI test database 0.18.1](https://github.com/bgruening/training_data/raw/master/GEMINI%20test%20database.tar.gz)
14. Either way you have know a database with all your variants, with pedigree relations, additional annotations and most importantly its fast to search. Have a look at all different GEMINI tools and run as many tools as possible on your test GEMINI databases. 

Try to get a feeling of what is possible with a variant database in GEMINI.

- GEMINI query is the most versatile of all the GEMINI tools. You can use it to ask 'interesting' questions in simple SQL (see the GEMINI handbook on its usage). Switch on the --header parameter.
- "select chrom, start, end from variants" will show you some information on all variants that were found in any of the three samples
- "select chrom, start, end, (gts).(*) from variants" will also show you the genotype of each sample also with the help wildcards
- "select chrom, start, end, gene, impact, (gts).(*) from variants v where v.impact_severity='HIGH'" will show you some more information and filter out only those variants that have a high impact
- "select chrom, vcf_if, start, end, ref, alt, gene, impact, (gts).(*) from variants v where v.impact_severity='HIGH'" also shows you the reference allele and the alternative allele and the RSID for the snp if it exists
- Pick one of your found variants and look at this position in IGV with all 3 BAM files loaded.





