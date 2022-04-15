---
title: Adjusting identifiers or input source for any mixed sourced data
area: datasets
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


- The inputs are a match for sequence content but simply adding "chr" will not make all chromosomes identifiers synch up between the inputs. How to fix or replace the inputs so that a match is possible.

- The underlying genome sequence content may or may not be identical.* Read method descriptions carefully to learn if that method is right for your usage case.

**Option 1**

**Sequence content is a match but adding "chr" is not enough to obtain an exact identifier match. You want to try to fix the identifiers anyway!!**

1. If the data is in a tabular format (BED, Interval, GTF -- with any headers removed first), and a suitable identifier mapping file can be obtained or created, the tool *Replace column by values which are defined in a convert file* can be used. Note that this will NOT work with BAM, VCF, Wiggle or other structured formats, as these are not tabular formatted data.
1. Manipuations with tools can often be used to split up a dataset, perform text substitutions and additions, concatinate datasets, and most other common operations one could do with command-line shell tools.
1. The dataset could also be downloaded locally to your computer and manipulated there using command-line tools or the text editor of choice.

**Option 2**

**Sequence content is a match but adding "chr" is not enough to obtain an exact identifier match. You DO NOT want to try to fix the identifiers or it is overly complicated or it is simply not possible to fix the data without an external reference mapping file (not always available).**

1. Obtain a reference annotation dataset that is a match for the reference genome used
2. Sometimes the source is the same for both
3. Sometimes the source is the same, but the content of the reference annotation is not ideal for the tools used
    - Example: The tool Cuffdiff makes use of specific attributes in the reference annotation (p_id, tss_id, gene_name). If these attributes are not present in the GTF dataset, the resuls will not be fully annotated and some calculations will be skipped
    - Use the iGenomes version of the reference annotation, as described below
    - Using Cuffdiff and the Gene ID is not present? Check your GTF file - the attribute gene_name is probably missing
4. Sometimes the source can be **iGenomes**, which does contain the extra specific attributes needed for RNA-seq and certain other operationsar
    - Example: The tool htseq-count is used and the attributes gene_id and transcript_id need to be distinct values (also true for the tool Cuffdiff for the best results)
    - Two sources: https://support.illumina.com/sequencing/sequencing_software/igenome.html and http://cole-trapnell-lab.github.io/cufflinks/igenome_table/index.html
    - Download the .tar file locally, uncompress it, then upload just the genes.gtf dataset to Galaxy
        - Note: the compression format .tar is not accepted by the Upload tool
        - If a .tar dataset is attempted to be uploaded, the load may fail or just the first file in the archive is uploaded (and it will not be the genes.gtf file)
5. [Genecode Genes](https://www.gencodegenes.org) is also an annotation source for some genome builds.
6. In summary:
    - For **Gencode**, copy the link to the GTF and paste it into the *Upload* tool. Hg38 data is here https://www.gencodegenes.org/. After it is loaded, remove the headers (lines that start with a "#") with the *Select* tool using the options "NOT Matching" with the regular expression `^#`. Once the formatting is fixed, change the datatype to be `gft` under Edit Attributes (pencil icon). The data will be given the datatype `gff` by default, which works fine with some tools and but not with others. Avoid the `gff3` version of this particular data (contains duplicated IDs and several RNA-seq tools do not work with annotation in that format anyway).
    - For **iGenomes**, the archive corresponding to the target genome/build needs to be locally downloaded, the tar archive unpacked, and then just the `genes.gtf` data uploaded to Galaxy (browse the local file, or use FTP). Find all available genome/builds here: https://support.illumina.com/sequencing/sequencing_software/igenome.html


**Option 3**

**Sequence content is NOT a match or you want to try using a different reference genome instead of a different reference annotation source (reverse of *Method 6* above).**

1. Map against the same reference genome that the reference annotation is based on
2. Where and if this reference genome is available will depend on the genome build
3. In most cases, the source will be the same for both
4. If loading your own genome, make sure it is formatted correctly as a Custom Genome
5. Promote the Custom Genome to a Custom Build and assign the genome/build metadata attribute to datasets
6. [Custom Reference Genome help](/learn/custom-genomes/)
7. Be aware that if the genome is large, this option may result in a memory failure. Try *Method 2* or consider moving to a local or cloud Galaxy where you can control the resources
