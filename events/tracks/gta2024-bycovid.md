---
layout: event-track

title: BeYond-COVID
description: |
  Learn about viral pathogen data analysis with Galaxy. Content here has been developed over the course of the COVID-19 pandemic and to a large part under the umbrella of the BY-COVID project as explained in a recent [blog post](https://galaxyproject.org/news/2024-09-16-by-covid-eol/).

  Work through the material at your own pace. If you need support contact us via the Slack Channel [#gta_BY-COVID](https://gtnsmrgsbord.slack.com/archives/C07NGS6BSLA).

slack_channel: gta_BY-COVID

contributions:
    organisers:
        - wm75
    instructors:
        - annasyme
        - bgruening
        - GarethPrice-Aus 
        - igormakunin
        - mateeullah905
        - wm75



program:
  - section: "Viral genomes sequencing data analysis " 
    description: |
      Here you find three hands-on tutorials demonstrating sequencing data analysis with Galaxy for three, rather different viral pathogens in increasing order of analysis complexity. Feel free to work through all three of them or just pick the one covering the virus you are most interested in, but consider watching first the introductory video, which explains some of the particularities of each virus and the corresponding analysis challenges.
    tutorials:
      - type: custom
        name: "[Sequencing spectrum viral genomes](https://gallantries.github.io/video-library/videos/virology/sequencing-spectrum-viral-genomes)"
        description: |
          [Lecture Video](https://gallantries.github.io/video-library/videos/virology/sequencing-spectrum-viral-genomes)
      - name: sars-cov-2
        topic: variant-analysis
      - name: aiv-analysis
        topic: variant-analysis
      - name: pox-tiled-amplicon
        topic: variant-analysis

  - section: "SARS-CoV-2 genome surveillance system"
    description: |
      Learn about how to scale up your Galaxy-based viral sequencing data analysis to batches of samples and how to automate the processing of such batches for highest sample throughput. Watch the introductory video for an overview of the contents of this section and how the pieces are connected. 
    tutorials:
      - type: custom
        name: "An automated SARS-CoV-2 genome surveillance system built around Galaxy"
        description: |
          [Lecture Video](https://gallantries.github.io/video-library/videos/one-health/galaxy-pathogen-surveillance); [Showcase page mentioned in the video](https://www.infectious-diseases-toolkit.org/showcase/covid19-galaxy)
      - name: sars-cov-2-variant-discovery
        topic: variant-analysis
      - name: workflow-automation
        topic: galaxy-interface
      - type: custom
        name: "[The usegalaxy.* SARS-CoV-2 Bot in Action](https://gallantries.github.io/video-library/videos/sars-cov2/usegalaxy-star-bot/)"
        description: |
          [Demo Video](https://gallantries.github.io/video-library/videos/sars-cov2/usegalaxy-star-bot/)
      - type: custom
        name: "[Upload to ENA](https://gallantries.github.io/video-library/videos/sars-cov2/upload-ena)"
        description: |
          [Demo Video](https://gallantries.github.io/video-library/videos/sars-cov2/upload-ena)
---
