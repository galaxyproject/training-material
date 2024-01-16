---
title: How to format fastq data for tools that require .fastqsanger format?
area: datatypes
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---


- Most tools that accept FASTQ data expect it to be in a specific FASTQ version: `.fastqsanger`. The `.fastqsanger` datatype must be assigned to each FASTQ dataset.

In order to do that:

- Watch the **[FASTQ Prep Illumina](http://vimeo.com/galaxyproject/fastqprep)** video for a complete walk-through.
- Run **FastQC** first to assess the type.
    - Run **FASTQ Groomer** if the data needs to have the quality scores rescaled.
    - If you are certain that the quality scores are already scaled to Sanger Phred+33 (the result of an Illumina 1.8+ pipeline), the datatype `.fastqsanger` can be directly assigned. Click on the pencil icon {% icon galaxy-pencil %} to reach the **Edit Attributes** form. In the center panel, click on the "Datatype" tab, enter the datatype `.fastqsanger`, and save.
- Run **FastQC** again on the entire dataset *if any changes were made to the quality scores* for QA.

**Other tips**

- If you are not sure what type of FASTQ data you have (maybe it is not Illumina?), see the help directly on the **FASTQ Groomer** tool for information about types.
    - For _Illumina_, first run **FastQC** on a sample of your data ([how to read the full report](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)). The output report will note the quality score type interpreted by the tool. If not `.fastqsanger`, run **FASTQ Groomer** on the entire dataset. If `.fastqsanger`, just assign the datatype.
    - For _SOLiD_, run **NGS: Fastq manipulation → AB-SOLID DATA → Convert**, to create a `.fastqcssanger` dataset. If you have uploaded a color space fastq sequence with quality scores already scaled to Sanger Phred+33 (`.fastqcssanger`), first confirm by running **FastQC** on a sample of the data. Then if you want to double-encode the color space into psuedo-nucleotide space (required by certain tools), see the instructions on the tool form **Fastq Manipulation** for the conversion.
    - If your data is _FASTA_, but you want to use tools that require FASTQ input, then using the tool **NGS: QC and manipulation → Combine FASTA and QUAL**. This tool will create "placeholder" quality scores that fit your data. On the output, click on the pencil icon {% icon galaxy-pencil %} to reach the **Edit Attributes** form. In the center panel, click on the "Datatype" tab, enter the datatype `.fastqsanger`, and save.
