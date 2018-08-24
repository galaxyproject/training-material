---
layout: tutorial_hands_on

title: "RAD-Seq to construct genetic maps"
zenodo_link: "https://doi.org/10.5281/zenodo.1219888"
questions:
  - "How to analyze RAD sequencing data for a genetic map study?"
objectives:
  - "SNP calling from RAD sequencing data"
  - "Find and correct haplotypes"
  - "Create input files for genetic map building software"
time_estimation: "12H"
key_points:
contributors:
  - yvanlebras
---

# Introduction
{:.no_toc}

This tutorial is based on the analysis described in [publication](http://www.genetics.org/content/188/4/799). 
Further information about the pipeline is available from [the official STACKS website](http://catchenlab.life.illinois.edu/stacks).
The authors developed a genetic map in the spotted gar and presented data from a single linkage group.
The gar genetic map is an F1 pseudotest cross between two parents and 94 of their F1 progeny. They took the markers that
appeared in one of the linkage groups and worked backwards to provide the raw reads from all of the stacks contributing to that linkage group.

This tutorial re-analyzes these data through to genotype determination. These data do not require demultiplexing and do not need processing though `Process Radtags tool`.

> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Pretreatments

## Data upload

The original data is available at [STACKS website](http://catchenlab.life.illinois.edu/stacks/) and the subset used here is findable on [Zenodo](https://zenodo.org/record/1219888#.WtZlK5c6-00).

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this RAD-seq exercise. 
> 2. Import Fasta files from parents and 20 progeny.
>
>    > ### {% icon comment %} Comments
>    > If you are using the [GenOuest Galaxy instance](https://galaxy.genouest.org), you can load the dataset using 'Shared Data' -> 'Data Libraries' -> '1 Galaxy teaching folder' -> 'EnginesOn' -> 'RADseq' -> 'Genetic map'
>    {: .comment}
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the following links into the text field
>    >  ```
>    > https://zenodo.org/record/1219888/files/female
>    > https://zenodo.org/record/1219888/files/male
>    > https://zenodo.org/record/1219888/files/progeny_1
>    > https://zenodo.org/record/1219888/files/progeny_2
>    > https://zenodo.org/record/1219888/files/progeny_3
>    > https://zenodo.org/record/1219888/files/progeny_4
>    > https://zenodo.org/record/1219888/files/progeny_5
>    > https://zenodo.org/record/1219888/files/progeny_6
>    > https://zenodo.org/record/1219888/files/progeny_7
>    > https://zenodo.org/record/1219888/files/progeny_8
>    > https://zenodo.org/record/1219888/files/progeny_9
>    > https://zenodo.org/record/1219888/files/progeny_10
>    > https://zenodo.org/record/1219888/files/progeny_11
>    > https://zenodo.org/record/1219888/files/progeny_12
>    > https://zenodo.org/record/1219888/files/progeny_13
>    > https://zenodo.org/record/1219888/files/progeny_14
>    > https://zenodo.org/record/1219888/files/progeny_15
>    > https://zenodo.org/record/1219888/files/progeny_16
>    > https://zenodo.org/record/1219888/files/progeny_17
>    > https://zenodo.org/record/1219888/files/progeny_18
>    > https://zenodo.org/record/1219888/files/progeny_19
>    > https://zenodo.org/record/1219888/files/progeny_20
>    > ```
>    > * Press **Start**  
>    {: .tip}
>
>    As default, Galaxy takes the link as name. It does not link the dataset to a database or a reference genome.
>
{: .hands_on}

# Building loci using STACKS

Run `Stacks: De novo map` Galaxy tool. This program will run `ustacks`, `cstacks`, and `sstacks` on each individual, accounting for the alignments of each read.

> ### {% icon comment %} Comment
>
> Information on `denovo_map.pl` and its parameters can be found online: https://creskolab.uoregon.edu/stacks/comp/denovo_map.php.
{: .comment}


> ### {% icon hands_on %} Hands-On: Stacks: De novo map
> **Stacks: De novo map** {% icon tool %}: Run Stacks selecting the Genetic map usage. Specify each parent as a sample in the appropriate box, then each of the 20 progeny and specify a CP Cross type, 3 for the Minimum number of identical raw reads required to create a stack, 3 for minimum number of identical raw reads required to create a stack in 'progeny' individuals, 3 for the number of mismatches allowed between loci when building the catalog and activate the option "remove, or break up, highly repetitive RAD-Tags in the ustacks program".
>
>    ![De novo map input](../../images/RAD2_Genetic_Map/denovomap_in.png)
>
>    Once Stacks has completed running, you will see 5 new data collections and 8 datasets.
>
>    ![The output of de novo map](../../images/RAD2_Genetic_Map/denovomap_out.png)
>
>    Investigate the output files: `result.log` and `catalog.*` (snps, alleles and tags).
>
>    Looking at the first file, `denovo_map.log`, you can see the command line used and the start as end execution time.
>
>    ![De novo map log file](../../images/RAD2_Genetic_Map/denovo_map_log_top.png)
>
>    Then are the different STACKS steps:
>
>    `ustacks`
>
>    ![De novo map:ustacks log](../../images/RAD2_Genetic_Map/denovo_map_log_ustacks.png)
>
>    `cstacks`
>
>    ![De novo map: cstacks](../../images/RAD2_Genetic_Map/denovo_map_log_cstacks.png)
>
>
>    > ### {% icon question %} Question
>    >
>    > 1. Can you identify the meanning of the number 425?
>    > 2. Looking at the catalog.tags file, identify specific and shared loci from each individual. Count the number of catalog loci coming from the first individual, from the second, and find on both parents.
>    >
>    > > ### {% icon solution %} Solution
>    > > 1. Here, the catalog is made with 459 tags, 425 coming from the "reference individual", a female. Some of these 425 can be shared with the other parent.
>    > > 2. 3500
>    > {: .solution }
>    {: .question}
>    `sstacks`
>
>    ![De novo map: sstacks log](../../images/RAD2_Genetic_Map/denovo_map_log_sstacks.png)
>
>
>    Lastly, `genotypes` is executed. It searches for markers identified on the parents and the associate progenies' haplotypes. If the first parent have a GA (ex: aatggtgtGgtccctcgtAc) and AC (ex: aatggtgtAgtccctcgtCc) haplotypes, and the second parent only a GA (ex: aatggtgtGgtccctcgtAc) haplotype, STACKS declares an ab/aa marker for this locus. Genotypes program then associate GA to a and AC to b and then scan progeny to determine which haplotype is found on each of them.
>
>    ![De novo map: genotypes log](../../images/RAD2_Genetic_Map/denovo_map_log_genotypes1.png)
>
>
>    ![De novo map: genotypes log](../../images/RAD2_Genetic_Map/denovo_map_log_genotypes2.png)
>
>
>    ![De novo map: genotypes log](../../images/RAD2_Genetic_Map/denovo_map_log_end.png)
>
>    Finally, 447 loci, markers, are kept to generate the `batch_1.genotypess_1.tsv` file. 459 loci are stored on the observed haplotype file `batch_1.haplotypes_1.tsv`.
>
{: .hands_on}

### Matches files

Here are `sample1.snps` (left) and `sample2.snps` (right)

![De novo map matches files](../../images/RAD2_Genetic_Map/denovo_map_matches1.PNG)

Catalog_ID (= catalog Stacks_ID) is composed by the `Stack_ID` from the "reference" individual, here sample1, but number is different from sample2 `Stack_ID`. Thus, in the `catalog.alleles.tsv`, the `Stack_ID` number 3 corresponds to the `Stack_ID` number 16 from sample2!

You can inspect the matches files (you maybe have to change the tsv datatype to a tabular one to correctly display the datasets).

![Male and female matches files](../../images/RAD2_Genetic_Map/denovo_map_matches2.png)

Consider catalog SNPs 27 & 28, on the 302 catalog locus:

![De novo map matches considering catalog SNPs](../../images/RAD2_Genetic_Map/denovo_map_matches_snps.png)

nd the corresponding catalog haplotypes, 3 on the 4 possible (AA, AT, GT but no GA):

![De novo map matches considering catalog haplotypes](../../images/RAD2_Genetic_Map/denovo_map_matches_alleles_haplotypes.png)

heterozygosity is observed on each parent (one ab, the other ac) and there are 19 genotypes for the 22 individuals.

![De novo map macthes: markers](../../images/RAD2_Genetic_Map/denovo_map_matches_markers.png)

We can then see that Stack_ID 330 for female corresponds to the 39 for male:


![De novo map matches: male and female](../../images/RAD2_Genetic_Map/denovo_map_matches_alleles_male_female.png)

# Genotypes determination

> ### {% icon hands_on %} Hands-on: Stacks: Genotypes
> **Stacks: genotypes** {% icon tool %}: Re-Run the last step of `Stacks: De novo map` pipeline specifying more options as:
>    > 1. The genetic map type (ie F1, F2 (left figure, F1xF1), Double Haploid, Back Cross (F1xF0), Cross Pollination (right figure, F1 or F2 but resulting from the cross of pure homozygous parents))
>    >
>    >    ![The genetic map type F2](../../images/RAD2_Genetic_Map/Genetic_map_F2.png)    ![The genetic map CrossPollination](../../images/RAD2_Genetic_Map/Genetic_map_CrossPollination.png)
>    >
>    >
>    > 2. Genotyping options output file type for input in genetic mapper tools (ie JoinMap, R/qtl, ...). Observe that the R/qtl format for an F2 cross type can be an input for MapMaker or Carthagene.
>    >
>    > 3. Thresholds concerning a minimal number of progeny and/or minimum stacks depth to consider a locus
>    >
>    > 4. Make Automated Corrections to the Data. This option allows the user to have the program automatically correct some types of errors. This setting can correct errors with the homozygous tags verification in the progeny by confirming the presence or absence of the SNP. If SNP detection model can't identify a site as heterygous or homozygous, that site is temporarily tagged as homozygous to facilitate the search, by sstacks, in concordance with the loci catalog. If a second allele is detected on the catalog (ie, in parents) and is found on a progeny with a weak frequency (<10% of a stack reads number), the genotypes program can correct the genotype. Additionally, it will delete a homozygous genotype on a particular individual if locus genotype is supported by less than 5 reads. Corrected genotypes are marked uppercase.
>
>    Here is an example of a locus originally marked as homozygous before automatic correction because an allele is supported by less than 5 reads. After correction, this locus is marked as heterozygous.
>
>    ![Automatic correction of genotypes](../../images/RAD2_Genetic_Map/genotypes_automatic_correction.png)
>
>    You can re-run **Stacks: genotypes** {% icon tool %}: modifying the number of genotyped progeny to consider a marker and thus be more or less stringent. Compare results.
>
{: .hands_on}

### Genotypes.tsv files

One line by locus, one column by individual (aa, ab, AB if automatic correction applied, bb, bc, ...) with observed genotype for each locus:

![Genotypes.tsv file overview](../../images/RAD2_Genetic_Map/genotypes_tsv.png)

### Genotypes.txt files

One line by individual, and for each individual, for each catalog locus, genotype:

![Genotypes.txt file overview](../../images/RAD2_Genetic_Map/genotypes_txt.png)

### Haplotypes.tsv files

One line by locus, one column by individual (aa, ab, AB if automatic correction applied, bb, bc, ...) with observed genotype for each locus:

![Haplotypes.tsv file overview](../../images/RAD2_Genetic_Map/haplotypes_tsv.png)

> ### {% icon question %} Question
>
> 1. The use of the deleverage algorithm allows to not consider loci obtained from merging more than 3 stacks. Why 3 if biologically, you are waiting something related to 2 for diploid organisms?
> 2. Re-execute **Stacks: De novo map** pipeline modifying the p-value treshold for the SNP model. What is the difference regarding to unverified haplotypes ?
>
> > ### {% icon solution %} Solution
> > 1. This value of 3 is important to use if we don't want to blacklist loci for whom 99.9% of individuals have one and/or the alt allele and 0.01% have a third one, resulting of a sequencing error.
> > 2. We see a moficiation of the number of unverified haplotypes
> {: .solution }
{: .question}


# Conclusion
{:.no_toc}

In this tutorial, we have analyzed real RAD sequencing data to extract useful information, such as genotypes and haplotypes to generate input files for downstream genetic map creation. This approach can be summarized with the following scheme:

![The genetic map tutorial workflow](../../images/genetic_map_workflow.PNG)
