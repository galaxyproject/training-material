---
layout: tutorial_hands_on

title: Genome annotation (Funannotate), comparison (OrthoFinder) and visualization
  (JBrowse2 & Circos) [Galaxy Training Material]
zenodo_link: https://zenodo.org/records/17612216
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- contributor1
- contributor2

---


# Introduction

<!-- This is a comment. -->

General introduction about the topic and then an introduction of the
tutorial (the questions and the objectives). It is nice also to have a
scheme to sum up the pipeline used during the tutorial. The idea is to
give to trainees insight into the content of the tutorial and the (theoretical
and technical) key concepts they will learn.

You may want to cite some publications; this can be done by adding citations to the
bibliography file (`tutorial.bib` file next to your `tutorial.md` file). These citations
must be in bibtex format. If you have the DOI for the paper you wish to cite, you can
get the corresponding bibtex entry using [doi2bib.org](https://doi2bib.org).

With the example you will find in the `tutorial.bib` file, you can add a citation to
this article here in your tutorial like this:
{% raw %} `{% cite Batut2018 %}`{% endraw %}.
This will be rendered like this: {% cite Batut2018 %}, and links to a
[bibliography section](#bibliography) which will automatically be created at the end of the
tutorial.


**Please follow our
[tutorial to learn how to fill the Markdown]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html)**

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Title for your first section

Give some background about what the trainees will be doing in the section.
Remember that many people reading your materials will likely be novices,
so make sure to explain all the relevant concepts.

## Title for a subsection
Section and subsection titles will be displayed in the tutorial index on the left side of
the page, so try to make them informative and concise!

# Hands-on Sections
Below are a series of hand-on boxes, one for each tool in your workflow file.
Often you may wish to combine several boxes into one or make other adjustments such
as breaking the tutorial into sections, we encourage you to make such changes as you
see fit, this is just a starting point :)

Anywhere you find the word "***TODO***", there is something that needs to be changed
depending on the specifics of your tutorial.

have fun!

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

# Title of the section usually corresponding to a big step in the analysis

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

![Alternative text](../../images/image_name "Legend of the image")

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> <details-title> More details about the theory </details-title>
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **gfastats**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [gfastats](toolshed.g2.bx.psu.edu/repos/bgruening/gfastats/gfastats/1.3.11+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `output` (Input dataset)
>    - *"Specify target sequences"*: `Disabled`
>    - *"Tool mode"*: `Summary statistics generation`
>        - *"Report mode"*: `Genome assembly statistics (--nstar-report)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **GC Skew**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [GC Skew](toolshed.g2.bx.psu.edu/repos/iuc/circos/circos_gc_skew/0.69.8+galaxy9) %} with the following parameters:
>    - *"Source for reference genome"*: `Use a genome from history`
>        - {% icon param-file %} *"Select a reference genome"*: `output` (Input dataset)
>    - *"Window size"*: `5000`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **RNA STAR**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [RNA STAR](toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.11b+galaxy0) %} with the following parameters:
>    - *"Single-end or paired-end reads"*: `Paired-end (as individual datasets)`
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, forward reads"*: `output` (Input dataset)
>        - {% icon param-file %} *"RNA-Seq FASTQ/FASTA file, reverse reads"*: `output` (Input dataset)
>    - *"Custom or built-in reference genome"*: `Use reference genome from history and create temporary index`
>        - {% icon param-file %} *"Select a reference genome"*: `output` (Input dataset)
>        - *"Length of the SA pre-indexing string"*: `11`
>        - *"Build index with or without known splice junctions annotation"*: `build index without gene-model`
>            - *"Per gene/transcript output"*: `No per gene or transcript output as no GTF was provided`
>        - *"Diploid mode"*: `No`
>    - *"Use 2-pass mapping for more sensitive novel splice junction discovery"*: `No`
>    - In *"Output filter criteria"*:
>        - *"Would you like to set additional output filters?"*: `No`
>    - In *"Algorithmic settings"*:
>        - *"Configure seed, alignment and limits options"*: `Use Defaults`
>    - *"Compute coverage"*: `No coverage`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **OrthoFinder**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [OrthoFinder](toolshed.g2.bx.psu.edu/repos/iuc/orthofinder_onlygroups/orthofinder_onlygroups/2.5.5+galaxy1) %} with the following parameters:
>    - *"Orthofinder starting point"*: `From fasta files`
>        - {% icon param-collection %} *"Select input fasta files"*: `output` (Input dataset collection)
>        - *"Sequence search program"*: `Diamond (faster)`
>    - *"Orthofinder run mode"*: `Full run (including gene trees)`
>        - *"Method for gene tree inference"*: `Dendroblast (faster)`
>        - *"Outputs"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Circos: bigWig to Scatter**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Circos: bigWig to Scatter](toolshed.g2.bx.psu.edu/repos/iuc/circos/circos_wiggle_to_scatter/0.69.8+galaxy9) %} with the following parameters:
>    - {% icon param-file %} *"Data file"*: `output` (output of **GC Skew** {% icon tool %})
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Funannotate predict annotation**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Funannotate predict annotation](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_predict/funannotate_predict/1.8.17+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Assembly to annotate"*: `output` (Input dataset)
>    - In *"Organism"*:
>        - *"Name of the species to annotate"*: `Mucor mucedo`
>        - *"Strain name"*: `muc1`
>    - In *"Evidences"*:
>        - {% icon param-file %} *"RNA-seq mapped to genome to train Augustus/GeneMark-ET"*: `mapped_reads` (output of **RNA STAR** {% icon tool %})
>        - *"Select protein evidences"*: `Custom protein sequences`
>            - {% icon param-file %} *"Proteins to map to genome"*: `output` (Input dataset)
>    - In *"EVM settings (advanced)"*:
>        - *"Split contigs into partitions for EVM processing?"*: `Yes`
>    - *"Which outputs should be generated"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **eggNOG Mapper**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [eggNOG Mapper](toolshed.g2.bx.psu.edu/repos/galaxyp/eggnog_mapper/eggnog_mapper/2.1.8+galaxy4) %} with the following parameters:
>    - *"Basis for annotation"*: `Seed orthologs computed with Diamond (diamond)`
>        - {% icon param-file %} *"Fasta sequences to annotate"*: `fasta_proteins` (output of **Funannotate predict annotation** {% icon tool %})
>        - *"Type of sequences"*: `proteins`
>        - *"Scoring matrix and gap costs"*: `BLOSUM62`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **InterProScan**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [InterProScan](toolshed.g2.bx.psu.edu/repos/bgruening/interproscan/interproscan/5.59-91.0+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Protein FASTA File"*: `fasta_proteins` (output of **Funannotate predict annotation** {% icon tool %})
>    - *"Use applications with restricted license, only for non-commercial use?"*: `Yes`
>    - *"Output format"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Funannotate functional**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Funannotate functional](toolshed.g2.bx.psu.edu/repos/iuc/funannotate_annotate/funannotate_annotate/1.8.17+galaxy5) %} with the following parameters:
>    - *"Input format"*: `GenBank (from 'Funannotate predict annotation' tool)`
>        - {% icon param-file %} *"Genome annotation in genbank format"*: `annot_gbk` (output of **Funannotate predict annotation** {% icon tool %})
>    - {% icon param-file %} *"Eggnog-mapper annotations file"*: `annotations` (output of **eggNOG Mapper** {% icon tool %})
>    - {% icon param-file %} *"InterProScan5 XML file"*: `outfile_xml` (output of **InterProScan** {% icon tool %})
>    - *"Strain name"*: `muc1`
>    - *"locus_tag from NCBI to rename GFF gene models with"*: `MMUCEDO_`
>    - *"Which outputs should be generated"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Busco**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Busco](toolshed.g2.bx.psu.edu/repos/iuc/busco/busco/5.8.0+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Sequences to analyse"*: `fa_proteins` (output of **Funannotate functional** {% icon tool %})
>    - *"Cached database with lineage"*: `Busco v5 Linage Datasets (all+2024-03-21)`
>    - *"Mode"*: `Annotated gene sets (protein)`
>    - *"Auto-detect or select lineage?"*: `Select lineage`
>        - *"Lineage"*: ``
>    - *"Which outputs should be generated"*: ``
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **JBrowse2**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [JBrowse2](toolshed.g2.bx.psu.edu/repos/fubar/jbrowse2/jbrowse2/3.6.5+galaxy1) %} with the following parameters:
>    - *"Action"*: `New JBrowse Instance`
>    - In *"Genome Assembly"*:
>        - {% icon param-repeat %} *"Insert Genome Assembly"*
>            - *"Reference genome to display"*: `Use a genome from history`
>                - {% icon param-file %} *"Select the reference genome"*: `output` (Input dataset)
>            - In *"Track Category"*:
>                - {% icon param-repeat %} *"Insert Track Category"*
>                    - *"Track Category Label"*: `Annotation`
>                    - In *"Track"*:
>                        - {% icon param-repeat %} *"Insert Track"*
>                            - *"Track Type"*: `GFF/GFF3/BED Features`
>                                - *"GFF/GFF3/BED Track Data"*: `Use data from history`
>                                - In *"Advanced: styling options"*:
>                                    - *"Display style"*: `LinearBasicDisplay`
>                - {% icon param-repeat %} *"Insert Track Category"*
>                    - *"Track Category Label"*: `RNA-seq`
>                    - In *"Track"*:
>                        - {% icon param-repeat %} *"Insert Track"*
>                            - *"Track Type"*: `BAM Pileups`
>                                - *"BAM Track Data"*: `Use data from history`
>                                - In *"Advanced: styling options"*:
>                                    - *"Display style"*: `Alignments (Pileup + SNP Coverage)`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Circos: Interval to Tiles**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Circos: Interval to Tiles](toolshed.g2.bx.psu.edu/repos/iuc/circos/circos_interval_to_tile/0.69.8+galaxy9) %} with the following parameters:
>    - *"Data Format"*: `GFF3`
>        - {% icon param-file %} *"GFF3 File"*: `gff3` (output of **Funannotate functional** {% icon tool %})
>        - *"GFF3 Attribute to pull value from"*: `ID`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}

## Sub-step with **Circos**

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Circos](toolshed.g2.bx.psu.edu/repos/iuc/circos/circos/0.69.8+galaxy9) %} with the following parameters:
>    - In *"Karyotype"*:
>        - *"Reference Genome Source"*: ` FASTA File from History (can be slow, generate a length file to improve execution time.)`
>            - {% icon param-file %} *"Source FASTA Sequence"*: `output` (Input dataset)
>    - In *"General"*:
>        - *"Plot Background"*: `Solid Color`
>    - In *"2D Data Tracks"*:
>        - In *"2D Data Plot"*:
>            - {% icon param-repeat %} *"Insert 2D Data Plot"*
>                - *"Plot Type"*: `Tiles`
>                    - In *"Plot Format Specific Options"*:
>                        - *"Fill Color"*: `#678de4`
>                        - *"Stroke Thickness"*: `0`
>                        - *"Overflow Behavior"*: `Hide: overflow tiles are not drawn`
>                - *"Minimum / maximum options"*: `Plot all values`
>                - In *"Rules"*:
>                    - In *"Rule"*:
>                        - {% icon param-repeat %} *"Insert Rule"*
>                            - In *"Conditions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Conditions to Apply"*
>                                    - *"Condition"*: `Based on qualifier value (when available)`
>                                        - *"Qualifier name"*: `strand`
>                                        - *"Condition"*: `Less than (numeric)`
>                                        - *"Qualifier value to compare against"*: `0`
>                            - In *"Actions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Actions to Apply"*
>                                    - *"Action"*: `Change Fill Color for all points`
>                                        - *"Fill Color"*: `#f79308`
>            - {% icon param-repeat %} *"Insert 2D Data Plot"*
>                - *"Outside Radius"*: `0.79`
>                - *"Inside Radius"*: `0.7`
>                - *"Plot Type"*: `Histogram`
>                    - {% icon param-file %} *"Histogram Data Source"*: `output` (output of **Circos: bigWig to Scatter** {% icon tool %})
>                    - In *"Plot Format Specific Options"*:
>                        - *"Fill Color"*: `#000000`
>                        - *"Stroke Thickness"*: `0`
>                - *"Minimum / maximum options"*: `Plot all values`
>                - In *"Rules"*:
>                    - In *"Rule"*:
>                        - {% icon param-repeat %} *"Insert Rule"*
>                            - In *"Conditions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Conditions to Apply"*
>                                    - *"Condition"*: `Apply to Every Point`
>                            - In *"Actions to Apply"*:
>                                - {% icon param-repeat %} *"Insert Actions to Apply"*
>                                    - *"Action"*: `Change Fill Color based on Value`
>                                        - *"Fill Color"*: `Sequential: Blue - Purple`
>    - In *"Limits"*:
>        - *"Maximum number of ideograms to draw"*: `1500`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > <comment-title> short description </comment-title>
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> <question-title></question-title>
>
> 1. Question1?
> 2. Question2?
>
> > <solution-title></solution-title>
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# Conclusion

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.