layout: learning-pathway
tags: [intermediate, immunopeptidomics, cancer, proteogenomics, label-free]
type: use

title: Prediction of potential neoantigens
description: |
  This learning path introduces a comprehensive immunopeptidogenomics workflow for neoantigen discovery using label-free mass spectrometry data. The modules guide you through fusion and variant database generation, peptide identification with FragPipe, peptide validation using PepQuery2, and immunogenicity assessment through HLA binding predictions and IEDB screening.

cover-image: shared/images/proteomics.png
cover-image-alt: image of a 3D protein folding structure

editorial_board:
- subinamehta

pathway:
  - section: "Neoantigen 1: Fusion-Database Generation"
    description: |
      Learn how to generate a personalized fusion peptide database using RNA-seq data. This step sets the foundation for identifying tumor-specific fusion peptides in downstream analyses.
    tutorials:
      - name: neoantigen-1-fusion-db-generation
        topic: proteomics

  - section: "Neoantigen 2: Non-Reference Database Generation"
    description: |
      Construct a non-reference proteogenomic database incorporating somatic mutations, indels, and other genomic alterations from VCF data.
    tutorials:
      - name: neoantigen-2-non-reference-database-generation
        topic: proteomics

  - section: "Neoantigen 3: Database Merge and FragPipe Discovery"
    description: |
      Merge the fusion and non-reference databases, then use FragPipe for mass spectrometry-based discovery of putative neopeptides.
    tutorials:
      - name: neoantigen-3-fragpipe-discovery
        topic: proteomics

  - section: "Neoantigen 4: PepQuery2 Verification"
    description: |
      Perform targeted verification of neoantigen candidates using PepQuery2 for peptide-spectrum match validation.
    tutorials:
      - name: neoantigen-4-peptide-verification
        topic: proteomics

  - section: "Neoantigen 5: Variant Annotation"
    description: |
      Annotate validated neopeptides with their corresponding genomic variants and protein context.
    tutorials:
      - name: neoantigen-5-variant-annotation
        topic: proteomics

  - section: "Neoantigen 6: Predicting HLA Binding"
    description: |
      Predict MHC binding affinity of validated neopeptides using tools such as NetMHCpan or similar.
    tutorials:
      - name: neoantigen-6-predicting-hla-binding
        topic: proteomics

  - section: "Neoantigen 7: IEDB Binding of PepQuery Validated Neopeptides"
    description: |
      Assess immunogenic potential of neopeptides by checking their binding predictions against immune epitope databases such as IEDB.
    tutorials:
      - name: neoantigen-7-hla-binding-novel-peptides
        topic: proteomics

---

New to immunopeptidogenomics or neoantigen prediction? Follow this learning path to discover how label-free mass spectrometry and proteogenomic integration can accelerate identification of clinically relevant neoantigens and help in personalized immunotherapy.
