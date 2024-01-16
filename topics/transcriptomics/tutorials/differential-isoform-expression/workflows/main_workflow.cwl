#!/usr/bin/env cwl-runner

cwlVersion: v1.2.0-dev2
class: Workflow
doc: |-
  Abstract CWL Automatically generated from the Galaxy workflow file: Genome-wide alternative splicing analysis: human

inputs:
  Active sites dataset:
    doc: Active sites dataset.
    type: File
    format: data
  CPAT header:
    doc: Galaxy compatible CPAT header.
    type: File
    format: data
  Control IDs:
    doc: IDs of control sequences, corresponding to filenames.
    type: File
    format: data
  Genome annotation:
    doc: Reference genome annotation in GTF format.
    type: File
    format: data
  Pfam-A HMM Stockholm file:
    doc: Stockholm file.
    type: File
    format: data
  Pfam-A HMM library:
    doc: HMM library.
    type: File
    format: data
  RNA-seq data collection:
    doc: Collection of paired-end RNA-seq datasets.
    type: File
    format: data
  Reference genome:
    doc: Genome reference in FASTA format.
    type: File
    format: data

outputs: {}

steps:
  10_fastp:
    in:
      single_paired|paired_input: RNA-seq data collection
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_fastp_fastp_0_23_2+galaxy0
      class: Operation
      inputs:
        single_paired|paired_input:
          type: File
          format: Any
      outputs:
        output_paired_coll:
          doc: input
          type: File
        report_json:
          doc: json
          type: File
    out:
    - output_paired_coll
    - report_json
  11_Flatten collection:
    in:
      input: RNA-seq data collection
    run:
      id: __FLATTEN__
      class: Operation
      inputs:
        input:
          type: File
          format: Any
      outputs:
        output:
          doc: input
          type: File
    out:
    - output
  12_Search in textfiles:
    in:
      infile: Genome annotation
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_grep_tool_1_1_1
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        output:
          doc: input
          type: File
    out:
    - output
  13_Search in textfiles:
    in:
      infile: Genome annotation
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_grep_tool_1_1_1
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        output:
          doc: input
          type: File
    out:
    - output
  14_Convert GTF to BED12:
    in:
      gtf_file: Genome annotation
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_gtftobed12_gtftobed12_357
      class: Operation
      inputs:
        gtf_file:
          type: File
          format: Any
      outputs:
        bed_file:
          doc: bed12
          type: File
    out:
    - bed_file
  15_FastQC:
    in:
      input_file: 11_Flatten collection/output
    run:
      id: toolshed_g2_bx_psu_edu_repos_devteam_fastqc_fastqc_0_73+galaxy0
      class: Operation
      inputs:
        input_file:
          type: File
          format: Any
      outputs:
        html_file:
          doc: html
          type: File
        text_file:
          doc: txt
          type: File
    out:
    - html_file
    - text_file
  16_gffread:
    in:
      input: 12_Search in textfiles/output
      reference_genome|genome_fasta: Reference genome
    run:
      id: toolshed_g2_bx_psu_edu_repos_devteam_gffread_gffread_2_2_1_3+galaxy0
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        reference_genome|genome_fasta:
          type: File
          format: Any
      outputs:
        output_exons:
          doc: fasta
          type: File
    out:
    - output_exons
  17_gffread:
    in:
      input: 13_Search in textfiles/output
      reference_genome|genome_fasta: Reference genome
    run:
      id: toolshed_g2_bx_psu_edu_repos_devteam_gffread_gffread_2_2_1_3+galaxy0
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        reference_genome|genome_fasta:
          type: File
          format: Any
      outputs:
        output_exons:
          doc: fasta
          type: File
    out:
    - output_exons
  18_MultiQC:
    in:
      results_0|software_cond|output_0|input: 15_FastQC/text_file
      results_1|software_cond|input: 10_fastp/report_json
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_multiqc_multiqc_1_11+galaxy1
      class: Operation
      inputs:
        results_0|software_cond|output_0|input:
          type: File
          format: Any
        results_1|software_cond|input:
          type: File
          format: Any
      outputs:
        html_report:
          doc: html
          type: File
        stats:
          doc: input
          type: File
    out:
    - stats
    - html_report
  19_Search in textfiles:
    in:
      infile: 15_FastQC/text_file
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_grep_tool_1_1_1
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        output:
          doc: input
          type: File
    out:
    - output
  20_Concatenate datasets:
    in:
      inputs: 19_Search in textfiles/output
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_cat_0_1_1
      class: Operation
      inputs:
        inputs:
          type: File
          format: Any
      outputs:
        out_file1:
          doc: input
          type: File
    out:
    - out_file1
  21_Search in textfiles:
    in:
      infile: 20_Concatenate datasets/out_file1
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_grep_tool_1_1_1
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        output:
          doc: input
          type: File
    out:
    - output
  22_Cut:
    in:
      input: 21_Search in textfiles/output
    run:
      id: Cut1
      class: Operation
      inputs:
        input:
          type: File
          format: Any
      outputs:
        out_file1:
          doc: tabular
          type: File
    out:
    - out_file1
  23_Text reformatting:
    in:
      infile: 22_Cut/out_file1
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        outfile:
          doc: input
          type: File
    out:
    - outfile
  24_Parse parameter value:
    in:
      input1: 23_Text reformatting/outfile
    run:
      id: param_value_from_file
      class: Operation
      inputs:
        input1:
          type: File
          format: Any
      outputs:
        integer_param:
          doc: expression.json
          type: File
    out:
    - integer_param
  25_RNA STAR:
    in:
      refGenomeSource|GTFconditional|sjdbGTFfile: Genome annotation
      refGenomeSource|GTFconditional|sjdbOverhang: 24_Parse parameter value/integer_param
      refGenomeSource|genomeFastaFiles: Reference genome
      singlePaired|input: 10_fastp/output_paired_coll
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_rgrnastar_rna_star_2_7_10b+galaxy3
      class: Operation
      inputs:
        refGenomeSource|GTFconditional|sjdbGTFfile:
          type: File
          format: Any
        refGenomeSource|GTFconditional|sjdbOverhang:
          type: File
          format: Any
        refGenomeSource|genomeFastaFiles:
          type: File
          format: Any
        singlePaired|input:
          type: File
          format: Any
      outputs:
        mapped_reads:
          doc: bam
          type: File
        output_log:
          doc: txt
          type: File
        splice_junctions:
          doc: interval
          type: File
    out:
    - output_log
    - splice_junctions
    - mapped_reads
  26_Concatenate datasets:
    in:
      inputs: 25_RNA STAR/splice_junctions
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_cat_0_1_1
      class: Operation
      inputs:
        inputs:
          type: File
          format: Any
      outputs:
        out_file1:
          doc: input
          type: File
    out:
    - out_file1
  27_Text reformatting:
    in:
      infile: 26_Concatenate datasets/out_file1
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        outfile:
          doc: input
          type: File
    out:
    - outfile
  28_Cut:
    in:
      input: 27_Text reformatting/outfile
    run:
      id: Cut1
      class: Operation
      inputs:
        input:
          type: File
          format: Any
      outputs:
        out_file1:
          doc: tabular
          type: File
    out:
    - out_file1
  29_Sort:
    in:
      input: 28_Cut/out_file1
    run:
      id: sort1
      class: Operation
      inputs:
        input:
          type: File
          format: Any
      outputs:
        out_file1:
          doc: input
          type: File
    out:
    - out_file1
  30_Unique:
    in:
      infile: 29_Sort/out_file1
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_sorted_uniq_1_1_0
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        outfile:
          doc: input
          type: File
    out:
    - outfile
  31_RNA STAR:
    in:
      refGenomeSource|GTFconditional|sjdbGTFfile: Genome annotation
      refGenomeSource|GTFconditional|sjdbOverhang: 24_Parse parameter value/integer_param
      refGenomeSource|genomeFastaFiles: Reference genome
      singlePaired|input: 10_fastp/output_paired_coll
      twopass|sj_precalculated: 30_Unique/outfile
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_rgrnastar_rna_star_2_7_10b+galaxy3
      class: Operation
      inputs:
        refGenomeSource|GTFconditional|sjdbGTFfile:
          type: File
          format: Any
        refGenomeSource|GTFconditional|sjdbOverhang:
          type: File
          format: Any
        refGenomeSource|genomeFastaFiles:
          type: File
          format: Any
        singlePaired|input:
          type: File
          format: Any
        twopass|sj_precalculated:
          type: File
          format: Any
      outputs:
        mapped_reads:
          doc: bam
          type: File
        output_log:
          doc: txt
          type: File
        signal_unique_str1:
          doc: bedgraph
          type: File
        signal_uniquemultiple_str1:
          doc: bedgraph
          type: File
        splice_junctions:
          doc: interval
          type: File
    out:
    - output_log
    - splice_junctions
    - mapped_reads
    - signal_unique_str1
    - signal_uniquemultiple_str1
  32_StringTie:
    in:
      guide|guide_source|ref_hist: Genome annotation
      input_options|input_bam: 31_RNA STAR/mapped_reads
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_stringtie_stringtie_2_2_1+galaxy1
      class: Operation
      inputs:
        guide|guide_source|ref_hist:
          type: File
          format: Any
        input_options|input_bam:
          type: File
          format: Any
      outputs:
        output_gtf:
          doc: gtf
          type: File
    out:
    - output_gtf
  33_Junction Annotation:
    in:
      input: 31_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    run:
      id: |-
        toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_junction_annotation_5_0_1+galaxy2
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        refgene:
          type: File
          format: Any
      outputs:
        outputjpdf:
          doc: pdf
          type: File
        outputpdf:
          doc: pdf
          type: File
        outputxls:
          doc: tabular
          type: File
        stats:
          doc: txt
          type: File
    out:
    - outputpdf
    - outputjpdf
    - outputxls
    - stats
  34_Junction Saturation:
    in:
      input: 31_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    run:
      id: |-
        toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_junction_saturation_5_0_1+galaxy2
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        refgene:
          type: File
          format: Any
      outputs:
        outputpdf:
          doc: pdf
          type: File
        outputr:
          doc: txt
          type: File
    out:
    - outputpdf
    - outputr
  35_Gene Body Coverage (BAM):
    in:
      batch_mode|input: 31_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    run:
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_geneBody_coverage_5_0_1+galaxy2
      class: Operation
      inputs:
        batch_mode|input:
          type: File
          format: Any
        refgene:
          type: File
          format: Any
      outputs:
        outputcurvespdf:
          doc: pdf
          type: File
        outputtxt:
          doc: txt
          type: File
    out:
    - outputcurvespdf
    - outputtxt
  36_Inner Distance:
    in:
      input: 31_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    run:
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_inner_distance_5_0_1+galaxy2
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        refgene:
          type: File
          format: Any
      outputs:
        outputfreqtxt:
          doc: txt
          type: File
        outputpdf:
          doc: pdf
          type: File
        outputtxt:
          doc: txt
          type: File
    out:
    - outputpdf
    - outputtxt
    - outputfreqtxt
  37_Infer Experiment:
    in:
      input: 31_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    run:
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_infer_experiment_5_0_1+galaxy2
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        refgene:
          type: File
          format: Any
      outputs:
        output:
          doc: txt
          type: File
    out:
    - output
  38_Read Distribution:
    in:
      input: 31_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    run:
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_read_distribution_5_0_1+galaxy2
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        refgene:
          type: File
          format: Any
      outputs:
        output:
          doc: txt
          type: File
    out:
    - output
  39_gffread:
    in:
      input: 32_StringTie/output_gtf
      reference_genome|genome_fasta: Reference genome
    run:
      id: toolshed_g2_bx_psu_edu_repos_devteam_gffread_gffread_2_2_1_3+galaxy0
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        reference_genome|genome_fasta:
          type: File
          format: Any
      outputs:
        output_exons:
          doc: fasta
          type: File
    out:
    - output_exons
  40_StringTie merge:
    in:
      guide_gff: Genome annotation
      input_gtf: 32_StringTie/output_gtf
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_stringtie_stringtie_merge_2_2_1+galaxy1
      class: Operation
      inputs:
        guide_gff:
          type: File
          format: Any
        input_gtf:
          type: File
          format: Any
      outputs:
        out_gtf:
          doc: gtf
          type: File
    out:
    - out_gtf
  41_MultiQC:
    in:
      results_0|software_cond|output_0|type|input: 31_RNA STAR/output_log
      results_1|software_cond|output_0|type|input: 37_Infer Experiment/output
      results_1|software_cond|output_1|type|input: 38_Read Distribution/output
      results_1|software_cond|output_2|type|input: 34_Junction Saturation/outputr
      results_1|software_cond|output_3|type|input: 33_Junction Annotation/stats
      results_1|software_cond|output_4|type|input: 35_Gene Body Coverage (BAM)/outputtxt
      results_1|software_cond|output_5|type|input: 36_Inner Distance/outputfreqtxt
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_multiqc_multiqc_1_11+galaxy1
      class: Operation
      inputs:
        results_0|software_cond|output_0|type|input:
          type: File
          format: Any
        results_1|software_cond|output_0|type|input:
          type: File
          format: Any
        results_1|software_cond|output_1|type|input:
          type: File
          format: Any
        results_1|software_cond|output_2|type|input:
          type: File
          format: Any
        results_1|software_cond|output_3|type|input:
          type: File
          format: Any
        results_1|software_cond|output_4|type|input:
          type: File
          format: Any
        results_1|software_cond|output_5|type|input:
          type: File
          format: Any
      outputs:
        html_report:
          doc: html
          type: File
        stats:
          doc: input
          type: File
    out:
    - stats
    - html_report
  42_rnaQUAST:
    in:
      gene_coordinates|gtf: Genome annotation
      reference: Reference genome
      transcripts: 39_gffread/output_exons
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_rnaquast_rna_quast_2_2_3+galaxy0
      class: Operation
      inputs:
        gene_coordinates|gtf:
          type: File
          format: Any
        reference:
          type: File
          format: Any
        transcripts:
          type: File
          format: Any
      outputs:
        fasta_files:
          doc: input
          type: File
        short_report_pdf:
          doc: pdf
          type: File
        stats:
          doc: txt
          type: File
    out:
    - fasta_files
    - stats
    - short_report_pdf
  43_gffread:
    in:
      input: 40_StringTie merge/out_gtf
      reference_genome|genome_fasta: Reference genome
    run:
      id: toolshed_g2_bx_psu_edu_repos_devteam_gffread_gffread_2_2_1_3+galaxy0
      class: Operation
      inputs:
        input:
          type: File
          format: Any
        reference_genome|genome_fasta:
          type: File
          format: Any
      outputs:
        output_exons:
          doc: fasta
          type: File
    out:
    - output_exons
  44_StringTie:
    in:
      guide|guide_source|ref_hist: 40_StringTie merge/out_gtf
      input_options|input_bam: 31_RNA STAR/mapped_reads
    run:
      id: toolshed_g2_bx_psu_edu_repos_iuc_stringtie_stringtie_2_2_1+galaxy1
      class: Operation
      inputs:
        guide|guide_source|ref_hist:
          type: File
          format: Any
        input_options|input_bam:
          type: File
          format: Any
      outputs:
        exon_expression:
          doc: tabular
          type: File
        exon_transcript_mapping:
          doc: tabular
          type: File
        intron_expression:
          doc: tabular
          type: File
        intron_transcript_mapping:
          doc: tabular
          type: File
        output_gtf:
          doc: gtf
          type: File
        transcript_expression:
          doc: tabular
          type: File
    out:
    - output_gtf
    - exon_expression
    - intron_expression
    - transcript_expression
    - exon_transcript_mapping
    - intron_transcript_mapping
  45_Filter collection:
    in:
      how|filter_source: Control IDs
      input: 44_StringTie/transcript_expression
    run:
      id: __FILTER_FROM_FILE__
      class: Operation
      inputs:
        how|filter_source:
          type: File
          format: Any
        input:
          type: File
          format: Any
      outputs:
        output_discarded:
          doc: input
          type: File
        output_filtered:
          doc: input
          type: File
    out:
    - output_filtered
    - output_discarded
  46_IsoformSwitchAnalyzeR:
    in:
      functionMode|genomeAnnotation: Genome annotation
      functionMode|tool_source|averageSize: 24_Parse parameter value/integer_param
      functionMode|tool_source|first_factor|factorLevel: 5_Input parameter/output
      functionMode|tool_source|first_factor|trans_counts: 45_Filter collection/output_filtered
      functionMode|tool_source|novoisoforms|stringtieAnnotation: 40_StringTie merge/out_gtf
      functionMode|tool_source|second_factor|factorLevel: 4_Input parameter/output
      functionMode|tool_source|second_factor|trans_counts: 45_Filter collection/output_discarded
      functionMode|transcriptome: 43_gffread/output_exons
    run:
      id: |-
        toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      class: Operation
      inputs:
        functionMode|genomeAnnotation:
          type: File
          format: Any
        functionMode|tool_source|averageSize:
          type: File
          format: Any
        functionMode|tool_source|first_factor|factorLevel:
          type: File
          format: Any
        functionMode|tool_source|first_factor|trans_counts:
          type: File
          format: Any
        functionMode|tool_source|novoisoforms|stringtieAnnotation:
          type: File
          format: Any
        functionMode|tool_source|second_factor|factorLevel:
          type: File
          format: Any
        functionMode|tool_source|second_factor|trans_counts:
          type: File
          format: Any
        functionMode|transcriptome:
          type: File
          format: Any
      outputs:
        switchList:
          doc: rdata
          type: File
    out:
    - switchList
  47_IsoformSwitchAnalyzeR:
    in:
      functionMode|robject: 46_IsoformSwitchAnalyzeR/switchList
    run:
      id: |-
        toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      class: Operation
      inputs:
        functionMode|robject:
          type: File
          format: Any
      outputs:
        isoformAA:
          doc: fasta
          type: File
        isoformNT:
          doc: fasta
          type: File
        switchList:
          doc: rdata
          type: File
        switchSummary:
          doc: tabular
          type: File
    out:
    - switchList
    - isoformAA
    - isoformNT
    - switchSummary
  48_CPAT:
    in:
      c: 17_gffread/output_exons
      gene: 47_IsoformSwitchAnalyzeR/isoformNT
      n: 16_gffread/output_exons
      r: Reference genome
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_cpat_cpat_3_0_4+galaxy0
      class: Operation
      inputs:
        c:
          type: File
          format: Any
        gene:
          type: File
          format: Any
        n:
          type: File
          format: Any
        r:
          type: File
          format: Any
      outputs:
        no_orf_seqs:
          doc: txt
          type: File
        orf_seqs:
          doc: fasta
          type: File
        orf_seqs_prob:
          doc: tsv
          type: File
        orf_seqs_prob_best:
          doc: tsv
          type: File
    out:
    - orf_seqs
    - orf_seqs_prob
    - orf_seqs_prob_best
    - no_orf_seqs
  49_PfamScan:
    in:
      active_sites|active_file: Active sites dataset
      fasta: 47_IsoformSwitchAnalyzeR/isoformAA
      pfam_data: Pfam-A HMM Stockholm file
      pfam_library: Pfam-A HMM library
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_pfamscan_pfamscan_1_6+galaxy0
      class: Operation
      inputs:
        active_sites|active_file:
          type: File
          format: Any
        fasta:
          type: File
          format: Any
        pfam_data:
          type: File
          format: Any
        pfam_library:
          type: File
          format: Any
      outputs:
        output:
          doc: tabular
          type: File
    out:
    - output
  4_Input parameter:
    in: {}
    run:
      id:
      class: Operation
      inputs: {}
      outputs: {}
    out: []
  50_Remove beginning:
    in:
      input: 48_CPAT/orf_seqs_prob_best
    run:
      id: Remove_beginning1
      class: Operation
      inputs:
        input:
          type: File
          format: Any
      outputs:
        out_file1:
          doc: input
          type: File
    out:
    - out_file1
  51_Text reformatting:
    in:
      infile: 50_Remove beginning/out_file1
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      class: Operation
      inputs:
        infile:
          type: File
          format: Any
      outputs:
        outfile:
          doc: input
          type: File
    out:
    - outfile
  52_Concatenate datasets:
    in:
      inputs: CPAT header
      queries_0|inputs2: 51_Text reformatting/outfile
    run:
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_cat_0_1_1
      class: Operation
      inputs:
        inputs:
          type: File
          format: Any
        queries_0|inputs2:
          type: File
          format: Any
      outputs:
        out_file1:
          doc: input
          type: File
    out:
    - out_file1
  53_IsoformSwitchAnalyzeR:
    in:
      functionMode|coding_potential|analyzeCPAT: 52_Concatenate datasets/out_file1
      functionMode|protein_domains|analyzePFAM: 49_PfamScan/output
      functionMode|robject: 47_IsoformSwitchAnalyzeR/switchList
    run:
      id: |-
        toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      class: Operation
      inputs:
        functionMode|coding_potential|analyzeCPAT:
          type: File
          format: Any
        functionMode|protein_domains|analyzePFAM:
          type: File
          format: Any
        functionMode|robject:
          type: File
          format: Any
      outputs:
        single_gene:
          doc: pdf
          type: File
        switchList:
          doc: rdata
          type: File
    out:
    - switchList
    - single_gene
  54_IsoformSwitchAnalyzeR:
    in:
      functionMode|coding_potential|analyzeCPAT: 52_Concatenate datasets/out_file1
      functionMode|protein_domains|analyzePFAM: 49_PfamScan/output
      functionMode|robject: 47_IsoformSwitchAnalyzeR/switchList
    run:
      id: |-
        toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      class: Operation
      inputs:
        functionMode|coding_potential|analyzeCPAT:
          type: File
          format: Any
        functionMode|protein_domains|analyzePFAM:
          type: File
          format: Any
        functionMode|robject:
          type: File
          format: Any
      outputs:
        consequencesEnrichment:
          doc: tabular
          type: File
        consequencesSummary:
          doc: tabular
          type: File
        consequences_fulldata:
          doc: tabular
          type: File
        genes_consequences:
          doc: input
          type: File
        genes_wo_consequences:
          doc: input
          type: File
        isoformFeatures:
          doc: tabular
          type: File
        mostSwitching:
          doc: tabular
          type: File
        plots_summary:
          doc: input
          type: File
        splicingEnrichment:
          doc: tabular
          type: File
        splicingSummary:
          doc: tabular
          type: File
        splicing_fulldata:
          doc: tabular
          type: File
        switchList:
          doc: rdata
          type: File
    out:
    - plots_summary
    - genes_consequences
    - genes_wo_consequences
    - switchList
    - mostSwitching
    - consequencesSummary
    - consequencesEnrichment
    - splicingSummary
    - splicingEnrichment
    - splicing_fulldata
    - consequences_fulldata
    - isoformFeatures
  5_Input parameter:
    in: {}
    run:
      id:
      class: Operation
      inputs: {}
      outputs: {}
    out: []
