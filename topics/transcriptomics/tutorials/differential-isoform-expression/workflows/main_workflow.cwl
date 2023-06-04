class: Workflow
cwlVersion: v1.2.0-dev2
doc: 'Abstract CWL Automatically generated from the Galaxy workflow file: Genome-wide
  alternative splicing analysis'
inputs:
  Active sites dataset:
    format: data
    type: File
  CPAT header:
    format: data
    type: File
  Control IDs:
    format: data
    type: File
  Genome annotation:
    format: data
    type: File
  Pfam-A HMM Stockholm file:
    format: data
    type: File
  Pfam-A HMM library:
    format: data
    type: File
  Protein coding transcripts:
    format: data
    type: File
  RNA-seq data collection:
    format: data
    type: File
  Reference genome:
    format: data
    type: File
  lncRNA transcripts:
    format: data
    type: File
outputs: {}
steps:
  12_Flatten collection:
    in:
      input: RNA-seq data collection
    out:
    - output
    run:
      class: Operation
      id: __FLATTEN__
      inputs:
        input:
          format: Any
          type: File
      outputs:
        output:
          doc: input
          type: File
  13_fastp:
    in:
      single_paired|paired_input: RNA-seq data collection
    out:
    - output_paired_coll
    - report_json
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_fastp_fastp_0_23_2+galaxy0
      inputs:
        single_paired|paired_input:
          format: Any
          type: File
      outputs:
        output_paired_coll:
          doc: input
          type: File
        report_json:
          doc: json
          type: File
  14_Convert GTF to BED12:
    in:
      gtf_file: Genome annotation
    out:
    - bed_file
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_gtftobed12_gtftobed12_357
      inputs:
        gtf_file:
          format: Any
          type: File
      outputs:
        bed_file:
          doc: bed12
          type: File
  15_FastQC:
    in:
      input_file: 12_Flatten collection/output
    out:
    - html_file
    - text_file
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_devteam_fastqc_fastqc_0_73+galaxy0
      inputs:
        input_file:
          format: Any
          type: File
      outputs:
        html_file:
          doc: html
          type: File
        text_file:
          doc: txt
          type: File
  16_Gene BED To Exon_Intron_Codon BED:
    in:
      input1: 14_Convert GTF to BED12/bed_file
    out:
    - out_file1
    run:
      class: Operation
      id: gene2exon1
      inputs:
        input1:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: bed
          type: File
  17_MultiQC:
    in:
      results_0|software_cond|output_0|input: 15_FastQC/text_file
      results_1|software_cond|input: 13_fastp/report_json
    out:
    - stats
    - html_report
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_multiqc_multiqc_1_11+galaxy1
      inputs:
        results_0|software_cond|output_0|input:
          format: Any
          type: File
        results_1|software_cond|input:
          format: Any
          type: File
      outputs:
        html_report:
          doc: html
          type: File
        stats:
          doc: input
          type: File
  18_Search in textfiles:
    in:
      infile: 15_FastQC/text_file
    out:
    - output
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_grep_tool_1_1_1
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        output:
          doc: input
          type: File
  19_Text reformatting:
    in:
      infile: 16_Gene BED To Exon_Intron_Codon BED/out_file1
    out:
    - outfile
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        outfile:
          doc: input
          type: File
  20_Concatenate datasets:
    in:
      inputs: 18_Search in textfiles/output
    out:
    - out_file1
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_cat_0_1_1
      inputs:
        inputs:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: input
          type: File
  21_Text reformatting:
    in:
      infile: 19_Text reformatting/outfile
    out:
    - outfile
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        outfile:
          doc: input
          type: File
  22_Text reformatting:
    in:
      infile: 19_Text reformatting/outfile
    out:
    - outfile
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        outfile:
          doc: input
          type: File
  23_Cut:
    in:
      input: 20_Concatenate datasets/out_file1
    out:
    - out_file1
    run:
      class: Operation
      id: Cut1
      inputs:
        input:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: tabular
          type: File
  24_Parse parameter value:
    in:
      input1: 21_Text reformatting/outfile
    out:
    - integer_param
    run:
      class: Operation
      id: param_value_from_file
      inputs:
        input1:
          format: Any
          type: File
      outputs:
        integer_param:
          doc: expression.json
          type: File
  25_Parse parameter value:
    in:
      input1: 22_Text reformatting/outfile
    out:
    - integer_param
    run:
      class: Operation
      id: param_value_from_file
      inputs:
        input1:
          format: Any
          type: File
      outputs:
        integer_param:
          doc: expression.json
          type: File
  26_Text reformatting:
    in:
      infile: 23_Cut/out_file1
    out:
    - outfile
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        outfile:
          doc: input
          type: File
  27_Parse parameter value:
    in:
      input1: 26_Text reformatting/outfile
    out:
    - integer_param
    run:
      class: Operation
      id: param_value_from_file
      inputs:
        input1:
          format: Any
          type: File
      outputs:
        integer_param:
          doc: expression.json
          type: File
  28_RNA STAR:
    in:
      algo|params|align|alignIntronMax: 24_Parse parameter value/integer_param
      algo|params|align|alignIntronMin: 25_Parse parameter value/integer_param
      refGenomeSource|GTFconditional|sjdbGTFfile: Genome annotation
      refGenomeSource|GTFconditional|sjdbOverhang: 27_Parse parameter value/integer_param
      refGenomeSource|genomeFastaFiles: Reference genome
      singlePaired|input: 13_fastp/output_paired_coll
    out:
    - output_log
    - splice_junctions
    - mapped_reads
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_rgrnastar_rna_star_2_7_10b+galaxy3
      inputs:
        algo|params|align|alignIntronMax:
          format: Any
          type: File
        algo|params|align|alignIntronMin:
          format: Any
          type: File
        refGenomeSource|GTFconditional|sjdbGTFfile:
          format: Any
          type: File
        refGenomeSource|GTFconditional|sjdbOverhang:
          format: Any
          type: File
        refGenomeSource|genomeFastaFiles:
          format: Any
          type: File
        singlePaired|input:
          format: Any
          type: File
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
  29_Concatenate datasets:
    in:
      inputs: 28_RNA STAR/splice_junctions
    out:
    - out_file1
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_cat_0_1_1
      inputs:
        inputs:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: input
          type: File
  30_Text reformatting:
    in:
      infile: 29_Concatenate datasets/out_file1
    out:
    - outfile
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        outfile:
          doc: input
          type: File
  31_Cut:
    in:
      input: 30_Text reformatting/outfile
    out:
    - out_file1
    run:
      class: Operation
      id: Cut1
      inputs:
        input:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: tabular
          type: File
  32_Sort:
    in:
      input: 31_Cut/out_file1
    out:
    - out_file1
    run:
      class: Operation
      id: sort1
      inputs:
        input:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: input
          type: File
  33_Unique:
    in:
      infile: 32_Sort/out_file1
    out:
    - outfile
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_sorted_uniq_1_1_0
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        outfile:
          doc: input
          type: File
  34_RNA STAR:
    in:
      algo|params|align|alignIntronMax: 24_Parse parameter value/integer_param
      algo|params|align|alignIntronMin: 25_Parse parameter value/integer_param
      refGenomeSource|GTFconditional|sjdbGTFfile: Genome annotation
      refGenomeSource|GTFconditional|sjdbOverhang: 27_Parse parameter value/integer_param
      refGenomeSource|genomeFastaFiles: Reference genome
      singlePaired|input: 13_fastp/output_paired_coll
      twopass|sj_precalculated: 33_Unique/outfile
    out:
    - output_log
    - splice_junctions
    - mapped_reads
    - signal_unique_str1
    - signal_uniquemultiple_str1
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_rgrnastar_rna_star_2_7_10b+galaxy3
      inputs:
        algo|params|align|alignIntronMax:
          format: Any
          type: File
        algo|params|align|alignIntronMin:
          format: Any
          type: File
        refGenomeSource|GTFconditional|sjdbGTFfile:
          format: Any
          type: File
        refGenomeSource|GTFconditional|sjdbOverhang:
          format: Any
          type: File
        refGenomeSource|genomeFastaFiles:
          format: Any
          type: File
        singlePaired|input:
          format: Any
          type: File
        twopass|sj_precalculated:
          format: Any
          type: File
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
  35_Read Distribution:
    in:
      input: 34_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    out:
    - output
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_read_distribution_5_0_1+galaxy2
      inputs:
        input:
          format: Any
          type: File
        refgene:
          format: Any
          type: File
      outputs:
        output:
          doc: txt
          type: File
  36_Junction Saturation:
    in:
      input: 34_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    out:
    - outputpdf
    - outputr
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_junction_saturation_5_0_1+galaxy2
      inputs:
        input:
          format: Any
          type: File
        refgene:
          format: Any
          type: File
      outputs:
        outputpdf:
          doc: pdf
          type: File
        outputr:
          doc: txt
          type: File
  37_Junction Annotation:
    in:
      input: 34_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    out:
    - outputpdf
    - outputjpdf
    - outputxls
    - stats
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_junction_annotation_5_0_1+galaxy2
      inputs:
        input:
          format: Any
          type: File
        refgene:
          format: Any
          type: File
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
  38_Inner Distance:
    in:
      input: 34_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    out:
    - outputpdf
    - outputtxt
    - outputfreqtxt
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_inner_distance_5_0_1+galaxy2
      inputs:
        input:
          format: Any
          type: File
        refgene:
          format: Any
          type: File
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
  39_Infer Experiment:
    in:
      input: 34_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    out:
    - output
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_infer_experiment_5_0_1+galaxy2
      inputs:
        input:
          format: Any
          type: File
        refgene:
          format: Any
          type: File
      outputs:
        output:
          doc: txt
          type: File
  40_Gene Body Coverage (BAM):
    in:
      batch_mode|input: 34_RNA STAR/mapped_reads
      refgene: 14_Convert GTF to BED12/bed_file
    out:
    - outputcurvespdf
    - outputtxt
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_nilesh_rseqc_rseqc_geneBody_coverage_5_0_1+galaxy2
      inputs:
        batch_mode|input:
          format: Any
          type: File
        refgene:
          format: Any
          type: File
      outputs:
        outputcurvespdf:
          doc: pdf
          type: File
        outputtxt:
          doc: txt
          type: File
  41_StringTie:
    in:
      guide|guide_source|ref_hist: Genome annotation
      input_options|input_bam: 34_RNA STAR/mapped_reads
    out:
    - output_gtf
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_stringtie_stringtie_2_2_1+galaxy1
      inputs:
        guide|guide_source|ref_hist:
          format: Any
          type: File
        input_options|input_bam:
          format: Any
          type: File
      outputs:
        output_gtf:
          doc: gtf
          type: File
  42_MultiQC:
    in:
      results_0|software_cond|output_0|type|input: 34_RNA STAR/output_log
      results_1|software_cond|output_0|type|input: 39_Infer Experiment/output
      results_1|software_cond|output_1|type|input: 35_Read Distribution/output
      results_1|software_cond|output_2|type|input: 36_Junction Saturation/outputr
      results_1|software_cond|output_3|type|input: 37_Junction Annotation/stats
      results_1|software_cond|output_4|type|input: 40_Gene Body Coverage (BAM)/outputtxt
      results_1|software_cond|output_5|type|input: 38_Inner Distance/outputfreqtxt
    out:
    - stats
    - html_report
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_multiqc_multiqc_1_11+galaxy1
      inputs:
        results_0|software_cond|output_0|type|input:
          format: Any
          type: File
        results_1|software_cond|output_0|type|input:
          format: Any
          type: File
        results_1|software_cond|output_1|type|input:
          format: Any
          type: File
        results_1|software_cond|output_2|type|input:
          format: Any
          type: File
        results_1|software_cond|output_3|type|input:
          format: Any
          type: File
        results_1|software_cond|output_4|type|input:
          format: Any
          type: File
        results_1|software_cond|output_5|type|input:
          format: Any
          type: File
      outputs:
        html_report:
          doc: html
          type: File
        stats:
          doc: input
          type: File
  43_StringTie merge:
    in:
      guide_gff: Genome annotation
      input_gtf: 41_StringTie/output_gtf
    out:
    - out_gtf
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_stringtie_stringtie_merge_2_2_1+galaxy1
      inputs:
        guide_gff:
          format: Any
          type: File
        input_gtf:
          format: Any
          type: File
      outputs:
        out_gtf:
          doc: gtf
          type: File
  44_StringTie:
    in:
      guide|guide_source|ref_hist: 43_StringTie merge/out_gtf
      input_options|input_bam: 34_RNA STAR/mapped_reads
    out:
    - output_gtf
    - exon_expression
    - intron_expression
    - transcript_expression
    - exon_transcript_mapping
    - intron_transcript_mapping
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_stringtie_stringtie_2_2_1+galaxy1
      inputs:
        guide|guide_source|ref_hist:
          format: Any
          type: File
        input_options|input_bam:
          format: Any
          type: File
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
  45_gffread:
    in:
      input: 43_StringTie merge/out_gtf
      reference_genome|genome_fasta: Reference genome
    out:
    - output_exons
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_devteam_gffread_gffread_2_2_1_3+galaxy0
      inputs:
        input:
          format: Any
          type: File
        reference_genome|genome_fasta:
          format: Any
          type: File
      outputs:
        output_exons:
          doc: fasta
          type: File
  46_Filter collection:
    in:
      how|filter_source: Control IDs
      input: 44_StringTie/transcript_expression
    out:
    - output_filtered
    - output_discarded
    run:
      class: Operation
      id: __FILTER_FROM_FILE__
      inputs:
        how|filter_source:
          format: Any
          type: File
        input:
          format: Any
          type: File
      outputs:
        output_discarded:
          doc: input
          type: File
        output_filtered:
          doc: input
          type: File
  47_IsoformSwitchAnalyzeR:
    in:
      functionMode|genomeAnnotation: Genome annotation
      functionMode|tool_source|averageSize: 27_Parse parameter value/integer_param
      functionMode|tool_source|first_factor|factorLevel: 5_Input parameter/output
      functionMode|tool_source|first_factor|trans_counts: 46_Filter collection/output_filtered
      functionMode|tool_source|novoisoforms|stringtieAnnotation: 43_StringTie merge/out_gtf
      functionMode|tool_source|second_factor|factorLevel: 4_Input parameter/output
      functionMode|tool_source|second_factor|trans_counts: 46_Filter collection/output_discarded
      functionMode|transcriptome: 45_gffread/output_exons
    out:
    - switchList
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      inputs:
        functionMode|genomeAnnotation:
          format: Any
          type: File
        functionMode|tool_source|averageSize:
          format: Any
          type: File
        functionMode|tool_source|first_factor|factorLevel:
          format: Any
          type: File
        functionMode|tool_source|first_factor|trans_counts:
          format: Any
          type: File
        functionMode|tool_source|novoisoforms|stringtieAnnotation:
          format: Any
          type: File
        functionMode|tool_source|second_factor|factorLevel:
          format: Any
          type: File
        functionMode|tool_source|second_factor|trans_counts:
          format: Any
          type: File
        functionMode|transcriptome:
          format: Any
          type: File
      outputs:
        switchList:
          doc: rdata
          type: File
  48_IsoformSwitchAnalyzeR:
    in:
      functionMode|robject: 47_IsoformSwitchAnalyzeR/switchList
    out:
    - switchList
    - isoformAA
    - isoformNT
    - switchSummary
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      inputs:
        functionMode|robject:
          format: Any
          type: File
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
  49_CPAT:
    in:
      c: Protein coding transcripts
      gene: 48_IsoformSwitchAnalyzeR/isoformNT
      n: lncRNA transcripts
      r: Reference genome
    out:
    - orf_seqs
    - orf_seqs_prob
    - orf_seqs_prob_best
    - no_orf_seqs
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_cpat_cpat_3_0_4+galaxy0
      inputs:
        c:
          format: Any
          type: File
        gene:
          format: Any
          type: File
        n:
          format: Any
          type: File
        r:
          format: Any
          type: File
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
  4_Input parameter:
    in: {}
    out: []
    run:
      class: Operation
      id: null
      inputs: {}
      outputs: {}
  50_PfamScan:
    in:
      active_sites|active_file: Active sites dataset
      fasta: 48_IsoformSwitchAnalyzeR/isoformAA
      pfam_data: Pfam-A HMM Stockholm file
      pfam_library: Pfam-A HMM library
    out:
    - output
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_pfamscan_pfamscan_1_6+galaxy0
      inputs:
        active_sites|active_file:
          format: Any
          type: File
        fasta:
          format: Any
          type: File
        pfam_data:
          format: Any
          type: File
        pfam_library:
          format: Any
          type: File
      outputs:
        output:
          doc: tabular
          type: File
  51_Remove beginning:
    in:
      input: 49_CPAT/orf_seqs_prob_best
    out:
    - out_file1
    run:
      class: Operation
      id: Remove_beginning1
      inputs:
        input:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: input
          type: File
  52_Text reformatting:
    in:
      infile: 51_Remove beginning/out_file1
    out:
    - outfile
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_awk_tool_1_1_2
      inputs:
        infile:
          format: Any
          type: File
      outputs:
        outfile:
          doc: input
          type: File
  53_Concatenate datasets:
    in:
      inputs: CPAT header
      queries_0|inputs2: 52_Text reformatting/outfile
    out:
    - out_file1
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_bgruening_text_processing_tp_cat_0_1_1
      inputs:
        inputs:
          format: Any
          type: File
        queries_0|inputs2:
          format: Any
          type: File
      outputs:
        out_file1:
          doc: input
          type: File
  54_IsoformSwitchAnalyzeR:
    in:
      functionMode|coding_potential|analyzeCPAT: 53_Concatenate datasets/out_file1
      functionMode|protein_domains|analyzePFAM: 50_PfamScan/output
      functionMode|robject: 48_IsoformSwitchAnalyzeR/switchList
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
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      inputs:
        functionMode|coding_potential|analyzeCPAT:
          format: Any
          type: File
        functionMode|protein_domains|analyzePFAM:
          format: Any
          type: File
        functionMode|robject:
          format: Any
          type: File
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
  55_IsoformSwitchAnalyzeR:
    in:
      functionMode|coding_potential|analyzeCPAT: 53_Concatenate datasets/out_file1
      functionMode|protein_domains|analyzePFAM: 50_PfamScan/output
      functionMode|robject: 48_IsoformSwitchAnalyzeR/switchList
    out:
    - switchList
    - single_gene
    run:
      class: Operation
      id: toolshed_g2_bx_psu_edu_repos_iuc_isoformswitchanalyzer_isoformswitchanalyzer_1_20_0+galaxy5
      inputs:
        functionMode|coding_potential|analyzeCPAT:
          format: Any
          type: File
        functionMode|protein_domains|analyzePFAM:
          format: Any
          type: File
        functionMode|robject:
          format: Any
          type: File
      outputs:
        single_gene:
          doc: pdf
          type: File
        switchList:
          doc: rdata
          type: File
  5_Input parameter:
    in: {}
    out: []
    run:
      class: Operation
      id: null
      inputs: {}
      outputs: {}

