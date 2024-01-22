---
layout: tutorial_hands_on

title: "Introduction to sequencing with python 1"
questions:
  - What are the origins of Sanger sequencing
  - How did sequencing machines evolve?
  - Haw can we simuate Sanger sequecning with Python?
objectives:
  - Have a basic understanding of history of sequencing
  - Understand Python basics
time_estimation: "1h"
key_points:
  - Sanger sequedncing is sequning by synthesis  
  - Python is powerful
contributions:
  authorship:
  - nekrut

subtopic: gnmx
draft: true
---

![](https://i.imgur.com/1PCleoW.png)

# The problem

The difficulty with sequencing nucleic acids is nicely summarized by [Hutchinson:2007](http://dx.doi.org/10.1093/nar/gkm688):

1. The chemical properties of different DNA molecules were so similar that it appeared difficult to separate them.
2. The chain length of naturally occurring DNA molecules was much greater than for proteins and made complete sequencing seems unapproachable.
3. The 20 amino acid residues found in proteins have widely varying properties that had proven useful in the separation of peptides. The existence of only four bases in DNA therefore seemed to make sequencing a more difficult problem for DNA than for protein.
4. No base-specific DNAases were known. Protein sequencing had depended upon proteases that cleave adjacent to certain amino acids.

It is therefore not surprising that protein-sequencing was developed before DNA sequencing by [Sanger and Tuppy:1951](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1197535/). 

tRNA was the first complete nucleic acid sequenced (see pioneering work of [Robert Holley and colleagues](http://www.jstor.org/stable/1715055) and also [Holley's Nobel Lecture](https://www.nobelprize.org/nobel_prizes/medicine/laureates/1968/holley-lecture.pdf)). Conceptually, Holley's approach was similar to Sanger's protein sequencing: break molecule into small pieces with RNases, determine sequences of small fragments, use overlaps between fragments to reconstruct (assemble) the final nucleotide sequence. 

The work on finding approaches to sequencing DNA molecules began in late 60s and early 70s. One of the earliest contributions has been made by Ray Wu (Cornell) and Dave Kaiser (Stanford), who used _E. coli_ DNA polymerase to [incorporate radioactively labelled nucleotides into protruding ends of bacteriphage lambda](http://www.sciencedirect.com/science/article/pii/S0022283668800129?via%3Dihub). It took several more years for the development of more "high throughput" technologies by Sanger and Maxam/Gilbert. The Sanger technique has ultimately won over Maxam/Gilbert's protocol due to its relative simplicity (once dideoxynucleotides has become commercially available) and the fact that it required smaller amount of starting material as the polymerase was used to generate fragments necessary for sequence determination. 

# [Sanger/Coulson](http://www.sciencedirect.com/science/article/pii/0022283675902132?via%3Dihub) plus/minus method


>[Fred Sanger](https://www.nature.com/articles/505027a) is one of only [four people](https://en.wikipedia.org/wiki/Category:Nobel_laureates_with_multiple_Nobel_awards), who received two Nobel Prizes in their original form (for scientific, not societal, breakthroughs).
{: .comment}

This methods builds on idea of Wu and Kaiser (for *minus* part) and on special property of DNA polymerase isolated from phage T4 (for *plus* part). The schematics of the method is given in the following figure:

-----

![](https://i.imgur.com/COhizBe.png)

**Plus/minus method**. From [Sanger & Coulson: 1975](http://www.sciencedirect.com/science/article/pii/0022283675902132?via%3Dihub)

----

In this method a primer and DNA polymerase is used to synthesize DNA in the presence of P<sup>32</sup>-labeled nucleotides (only one of four is labeled). This generates P<sup>32</sup>-labeled copies of DNA being sequenced. These are then purified and (without denaturing) separated into two groups: *minus* and *plus*. Each group is further divided into four equal parts. 

In the case of *minus* polymerase and a mix of nucleotides minus one are added to each of the four aliquotes: ACG (-T), ACT (-G), CGT (-A), AGT (-C). As a result in each case DNA strand is extended up to a missing nucleotide. 


In the case of plus only one nucleotide is added to each of the four aliquotes (+A, +C, +G, and +T) and T4 DNA polymerase is used. T4 DNA polymerase acts as an exonuclease that would degrade DNA from 3'-end up to a nucleotide that is supplied in the reaction. 

The products of these are loaded into a denaturing polyacrylamide gel as a eight tracks (four for minus and four for plus):

----

![](https://i.imgur.com/EEAKSA8.png)

**Plus/minus method gel radiograph**. From [Sanger & Coulson: 1975](http://www.sciencedirect.com/science/article/pii/0022283675902132?via%3Dihub)

-----

# Maxam/Gilbert chemical cleavage method

In this method DNA is terminally labeled with P<sup>32</sup>, separated into four equal aliquotes.  Two of these are treated with [Dimethyl sulfate (DMSO)](https://en.wikipedia.org/wiki/Dimethyl_sulfate) and remaining two are treated with [hydrazine](https://en.wikipedia.org/wiki/Hydrazine). 

DMSO methylates G and A residues. Treatment of DMSO-incubated DNA with alkali at high temperature will break DNA chains at G and A with Gs being preferentially broken, while treatment of DMSO-incubated DNA with acid will preferentially break DNA at As. Likewise treating hydrazine-incubated DNA with [piperidine](https://en.wikipedia.org/wiki/Piperidine) breaks DNA at C and T, while DNA treated with hydrazine in the presence of NaCl preferentially brakes at Cs. The four reactions are then loaded on a gel generating the following picture:

------

![](https://i.imgur.com/wETzOTp.jpg)

**Radiograph of Maxam/Gilbert gel**. From [Maxam & Gilbert: 1977](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC392330/pdf/pnas00024-0174.pdf)

-----

# Sanger dideoxy method

The original Sanger +/- method was not popular and had a number of technical limitations. In a new approach Sanger took advantage of inhibitors that stop the extension of a DNA strand at particular nucleotides. These inhibitors are dideoxy analogs of normal nucleotide triphosphates:

-----

![](https://i.imgur.com/VUa996S.png)

**Sanger ddNTP gel**. From [Sanger:1977](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC431765/pdf/pnas00043-0271.pdf)

----

# Original approaches were laborious

In the original [Sanger paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC431765/pdf/pnas00043-0271.pdf) the authors sequenced bacteriophage phiX174 by using its own restriction fragments as primers. This was an ideal set up to show the proof of principle for the new method. This is because phiX174 DNA is homogeneous and can be isolated in large quantities. Now suppose that you would like to sequence a larger genome (say _E. coli_). Remember that the original version of Sanger method can only sequence fragments up to 200 nucleotides at a time. So to sequence the entire _E. coli_ genome (which by-the-way was not sequenced until [1997](http://science.sciencemag.org/content/277/5331/1453)) you would need to split the genome into multiple pieces and sequence each of them individually. This is hard, because to produce a readable Sanger sequencing gel each sequence must be amplified to a suitable amount (around 1 nanogram) and be homogeneous (you cannot mix multiple DNA fragments in a single reaction as it will be impossible to interpret the gel). Molecular cloning enabled by the availability of commercially available restriction enzymes and cloning vectors simplified this process. Until the onset of next generation sequencing in 2005 the process for sequencing looked something like this:

* (**1**) - Generate a collection of fragments you want to sequence. It can be a collection of fragments from a genome that was mechanically sheared or just a single fragment generated by PCR.
* (**2**) - These fragment(s) are then cloned into a plasmid vector (we will talk about other types of vectors such as BACs later in the course).
* (**3**) - Vectors are transformed into bacterial cells and positive colonies (containing vectors with insert) are picked from an agar plate.
* (**4**) - Each colony now represents a unique piece of DNA. 
* (**5**) - An individual colony is used to seed a bacterial culture that is grown overnight.
* (**6**) - Plasmid DNA is isolated from this culture and now can be used for sequencing because it is (1) homogeneous and (2) we now a sufficient amount.
* (**7**) - It is sequenced using universal primers. For example the image below shows a map for pGEM-3Z plasmid (a pUC18 derivative). Its multiple cloning site is enlarged and sites for **T7** and **SP6** sequencing primers are shown. These are the **pads** I'm referring to in the lecture. These provide universal sites that can be used to sequence any insert in between. 

-----

![](https://i.imgur.com/bkfrwVJ.png)

**pGEM-3Z**. Figure from Promega, Inc.

-----

Until the invention of NGS the above protocol was followed with some degree of automation. But you can see that it was quite laborious if the large number of fragments needed to be sequenced. This is because each of them needed to be subcloned and handled separately. This is in part why Human Genome Project, a subject of our next lecture, took so much time to complete. 

# Evolution of sequencing machines

The simplest possible sequencing machine is a [gel rig with polyacrylamide gel](https://en.wikipedia.org/wiki/Polyacrylamide_gel_electrophoresis). Sanger used it is his protocol obtaining the following results:

-----

![](https://i.imgur.com/7mRCGEi.png)

Figure from [Sanger et al. 1977](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC431765/pdf/pnas00043-0271.pdf).

----

Here for sequencing each fragment four separate reactions are performed (with ddA, ddT, ggC, and ddG) and four lanes on the gel are used. [One simplification of this process](http://www.jstor.org/stable/pdf/2879539.pdf) that came in the 90s was to use fluorescently labeled dideoxy nucleotides. This is easier because everything can be performed in a single tube and uses a single lane on a gel:

-----

![](https://i.imgur.com/0eC2lLT.png)

Figure from Applied Biosystems [support site](https://www3.appliedbiosystems.com/cms/groups/mcb_support/documents/generaldocuments/cms_041003.pdf).

------

However, there is still substantial labor involved in pouring the gels, loading them, running machines, and cleaning everything post-run. A significant improvement was offered by the development of capillary electrophoresis allowing automation of liquid handling and sample loading. Although several manufacturers have been developing and selling such machines a _de facto_ standard in this area was (and still is) the Applied Biosystems (ABI) Genetics and DNA Anlayzer systems. The highest throughput ABI system, 3730_xl_, had 96 capillaries and could automatically process 384 samples. 

## NGS!

384 samples may sound like a lot, but it is nothing if we are sequencing an entire genome. The beauty of NGS is that these technologies are not bound by sample handling logistics. They still require preparation of libraries, but once a library is made (which can be automated) it is processed more or less automatically to generate multiple copies of each fragment (in the case of 454, Illumina, and Ion Torrent) and loaded onto the machine, where millions of individual fragments are sequenced simultaneously. The following videos and slides explains how these technologies work.

## Historic NGS technologies: 454

454 Technology is a massively parallel modification of [pyrosequencing](http://genome.cshlp.org/content/11/1/3) technology. Incorporation of nucleotides are registered by a [CCD](https://en.wikipedia.org/wiki/Charge-coupled_device) camera as a flash of light generated from the interaction between ATP and Luciferin. The massive scale of 454 process is enabled by generation of a population of beads carrying multiple copies of the same DNA fragment. The beads are distributed across a microtiter plate where each well of the plate holding just one bead. Thus every unique coordinate (a well) on the plate generates flashes when a nucleotide incorporation event takes plate. This is "monochrome" technique: flash = nucleotide is incorporated; lack of flash = no incorporation. Thus to distinguish between A, C, G, and T individual nucleotides are "washed" across the microtiter plate at discrete times: if **A** is being washed across the plate and a flash of light is emitted, this implies that A is present in the fragment being sequenced. 

454 can generated fragments up 1,000 bases in length. Its biggest weakness is inability to precisely determine the length of [homopolymer runs](https://www.broadinstitute.org/crd/wiki/index.php/Homopolymer). Thus the main type if sequencing error generated by 454 are insertions and deletions (indels).


### Reading

* 2001 | [Overview of pyrosequencing methodology - Ronaghi](http://genome.cshlp.org/content/11/1/3)
* 2005 | [Description of 454 process - Margulies et al.](http://www.nature.com/nature/journal/v437/n7057/pdf/nature03959.pdf)
* 2007 | [History of pyrosequencing - Pål Nyrén](http://link.springer.com/protocol/10.1385/1-59745-377-3:1)
* 2007 | [Errors in 454 data - Huse et al. ](http://genomebiology.com/content/pdf/gb-2007-8-7-r143.pdf)
* 2010 | [Properties of 454 data - Balzer et al.](http://bioinformatics.oxfordjournals.org/content/26/18/i420.full.pdf+html)

## A few classical papers

In a series of now classical papersp Philip Green and co-workers has developed a quantitative framework for the analysis of data generated by automated DNA sequencers. 

-----

[![](https://i.imgur.com/gvaAfpa.png)
](http://genome.cshlp.org/content/8/3/175.full)

[![](https://i.imgur.com/fUJbwXH.png)
](http://genome.cshlp.org/content/8/3/186.full)
**Two papers** were published back to back in *Genome Research*.

-----

## Myers - Green debate

![](https://i.imgur.com/F7WriuG.png)

In particular they developed a standard metric for describing the reliability of base calls:

>An important technical aspect of our work is the use of log-transformed error probabilities rather than untransformed ones, which facilitates working with error rates in the range of most importance (very close to 0). Specifically, we define the quality value q assigned to a base-call to be
>$q = -10\times log_{10}(p)$
> where p is the estimated error probability for that base-call. Thus a base-call having a probability of 1/1000 of being incorrect is assigned a quality value of 30. Note that high quality values correspond to low error probabilities, and conversely.

We will be using the concept of "*quality score*" or "*phred-scaled quality score*" repeatedly in this course. 

