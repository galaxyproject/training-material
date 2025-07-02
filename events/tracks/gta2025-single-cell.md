---
layout: event-track

title: Single Cell
description: Learn all about Single Cell on Galaxy. Start with the tutorial at your own pace. If you need support contact us via the Slack Channel [gta_single-cell](https://gtnsmrgsbord.slack.com/channels/{{page.slack_channel}}).

slack_channel: gta_single-cell

contributions:
    organisers:
        - pavanvidem
        - nomadscientist
    instructors:
        - annasyme
        - Nilchia
        - dianichj
        - hrhotz
        - igormakunin
        - jennaj
        - MarisaJL
        - pavanvidem
        - nomadscientist


program:
  - section: "Introduction"
    description: |
      Introduction lecture to single-cell data analysis and data formats.
    tutorials:
      - name: scrna-intro
        topic: single-cell
      - name: scrna-data-formats
        topic: single-cell

  - section: "Raw sequencing data to count matrix"
    description: |
      Generate a cell-by-gene matrix from droplet-based single-cell RNA sequencing data.
    tutorials:
      - name: scrna-preprocessing-tenx
        topic: single-cell
      - name: scrna-case_alevin
        topic: single-cell

  - section: "Your first analysis based on Scanpy toolkit"
    description: |
      This section includes standard single-cell data analysis, including preprocessing, clustering and the identification of cell types via known marker genes; and pseudobulk analysis to help in understanding cell-type-specific gene expression changes. The whole analysis uses Anndata objects, Scanpy toolkit and decoupler.
    tutorials:
      - name: scrna-scanpy-pbmc3k
        topic: single-cell
      - name: pseudobulk-analysis
        topic: single-cell

  - section: "Your first analysis based on Seurat tools"
    description: |
      This section includes single-cell data analysis including preprocessing, clustering and the identification of cell types via known marker genes. The whole analysis uses Seurat objects and Seurat R package.
    tutorials:
      - name: scrna-seurat-pbmc3k
        topic: single-cell

---
