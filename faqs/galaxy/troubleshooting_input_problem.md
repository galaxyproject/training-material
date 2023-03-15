---
redirect_from: [/faqs/galaxy/analysis_job_failure_input_problem]
title: Understanding input error messages
area: troubleshooting
box_type: tip
layout: faq
contributors: [jennaj, garimavs]
---

Input problems are very common across any analysis that makes use of programmed tools.

- Causes:
    - No quality assurance or content/formatting checks were run on the first datasets of an analysis workflow.
    - Incomplete dataset Upload.
    - Incorrect or unassigned [datatype](https://training.galaxyproject.org/training-material/faqs/galaxy/#datasets) or [database](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_change_dbkey.html).
    - Tool-specific formatting requirements for inputs were not met.
    - Parameters set on a tool form are a mismatch for the input data content or format.
    - Inputs were in an error state (red) or were putatively successful (green) but are empty.
    - Inputs do not meet the datatype specification.
    - Inputs do not contain the exact content that a tool is expecting or that was input in the form.
    - Annotation files are a mismatch for the selected or assigned reference genome build.
    - **Special case:** Some of the data were generated outside of Galaxy, but later a built-in indexed genome build was assigned in Galaxy for use with downstream tools. This scenario can work, but only if those two reference genomes are an exact match.

- Solutions:
    - Review our [Troubleshooting Tips](https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_troubleshooting.html) for what and where to check.
    - Review the [GTN](https://training.galaxyproject.org/) for related tutorials on tools/analysis plus FAQs.
    - Review [Galaxy Help](https://help.galaxyproject.org/) for prior discussion with extended solutions.
    - Review [datatype FAQs](https://training.galaxyproject.org/training-material/faqs/galaxy/datatypes_understanding_datatypes.html).
    - Review the tool form.
        - Input selection areas include usage help.
        - The help section at the bottom of a tool form often has examples. Does your own data match the format/content?
        - See the links to publications and related resources.
    - Review the inputs.
        - All inputs must be in a success state (green) and actually contain content.
        - Did you [directly assign the datatype](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_change_datatype.html) or [convert the datatype](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_convert_datatype.html)? What results when the [datatype is detected](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_detect_datatype.html) by Galaxy? If these differ, there is likely a content problem.
        - For most analysis, allowing Galaxy to detect the datatype during Upload is best and adjusting a datatype later should rarely be needed. If a datatype is modified, the change has a specific purpose/reason.
        - Does your data have headers? Is that in specification for the datatype? Does the tool form have an option to specify if the input has headers or not? Do you need to remove headers first for the correct datatype to be detected? Example [GTF](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_working_with_reference_annotation.html).
        - Large inputs? Consider modifying your inputs to be smaller. Examples: [FASTQ](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_working_with_fastq.html) and [FASTA](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_working_with_fasta.html).
    - Run quality checks on your data.
        - Search [GTN](https://training.galaxyproject.org/) tutorials with the keyword “qa-qc” for examples.
        - Search [Galaxy Help](https://help.galaxyproject.org/) with the keywords “qa-qc” and your datatype(s) for more help.
    - Reference annotation tips.
        - In most cases, [GTF is preferred over GFF3](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_working_with_reference_annotation.html).
        - Search [Galaxy Help](https://help.galaxyproject.org/) with the keywords “gtf” and “gff3” for more help.
    - Input mismatch tips.
        - Do the chromosome/sequence identifiers exactly match between all inputs? Search [Galaxy Help](https://help.galaxyproject.org/) for more help about how to correct build/version identifier mismatches between inputs.
        - "Chr1" and "chr1" and "1" do not mean the same thing to a tool.
    - Custom genome transcriptome exome tips. See [FASTA](https://training.galaxyproject.org/training-material/faqs/galaxy/datasets_working_with_fasta.html).
