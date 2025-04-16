---
layout: learning-pathway
tags: [advanced, single-cell]
cover-image: assets/images/wab-lp-deconvolution.png
cover-image-alt: "XY plot of actual versus inferred proportions, with coloured dots representing clusters and largely falling in a 1-1 slope"
type: use

editorial_board:
- nomadscientist
- hexhowells
- mtekman

funding:
  - elixir-uk-dash

title: "Bioinformatics Projects: Using deconvolution to get new insights from old bulk RNA-seq data"
description: |
  Are you an educator looking for project ideas for students to practice independent enquiry and research skills? Are you a student looking for a project idea? Look no more - here, you will find a learning pathway of tutorials that can guide you through the skills to find old data and transform it into new results!

  To be clear, we will only provide the methods - you will need to come up with your own research question by exploring the literature and available public datasets, apply these analyses, and interpret the results. Your research question will take the form of, **"How does `variable X` impact the cell type proportions in `issue/sample/organism Y`?"**

  Note: You will need to be familiar with the Galaxy interface and single-cell RNA-seq analysis in general to follow this Learning Pathway. You can do so by completing the [Introduction to single-cell analysis learning pathway]({% link learning-pathways/intro_single_cell.md %}). It would be a bonus to also complete the [Beyond single cell learning pathway]({% link learning-pathways/beyond_single_cell.md %}) to reinforce that knowledge.

  For support throughout these tutorials, join our Galaxy [single cell chat group on Matrix](https://matrix.to/#/#Galaxy-Training-Network_galaxy-single-cell:gitter.im) to ask questions!



pathway:
  - section: "Module 1: What's deconvolution?"
    description: |
      First, you will learn about the concept of deconvolution. This will help you formulate your question and identify datasets next.
    tutorials:
      - name: bulk-music
        topic: single-cell

  - section: "Module 2: Picking & importing a dataset"
    description: |
      Next, you will need to pick a bulk RNA-seq dataset, along with a corresponding single-cell dataset as a reference to perform deconvolution. You will need to then transform these datasets into ESet objects. We have set up these tutorials to work with datasets from the European Bioinformatics Institute, because these are carefully curated and work with our workflows. You can try others, but you may experience challenges.
    tutorials:
      - name: bulk-music-3-preparebulk
        topic: single-cell
      - name: EBI-retrieval
        topic: single-cell
      - name: bulk-music-2-preparescref
        topic: single-cell

  - section: "Module 3: Does my reference work well?"
    description: |
              Next, you will benchmark your reference dataset to see how accurate it is at inferring proportions. If it does not work well, you may need to pick a different dataset and try again!
    tutorials:
      - name: bulk-deconvolution-evaluate
        topic: single-cell

  - section: "Module 4: Analysing your data!"
    description: |
              At long last, you've done all the hard work of learning about deconvolution, picking your datasets, reformatting them for use, and making sure your reference is of a high quality. You can now finally infer cell proportions from your bulk RNA-seq samples, and compare them across a variable of interest!
    tutorials:
     - name: bulk-music-4-compare
       topic: single-cell

  - section: "The End!"
    description: |
      And now you're done! We hope that you generated interesting results that you are able to write up in a fantastic project. We would love to hear from you if you have! Contact us via our Galaxy [single cell chat group on Matrix](https://matrix.to/#/#Galaxy-Training-Network_galaxy-single-cell:gitter.im). Alternatively, if you prefer Slack, join the [GTN's Slack workspace](https://gxy.io/gtn-slack) and message our #single-cell-users channel.

      You will find more features, tips and tricks in our general [Galaxy Single-cell Training page](/training-material/topics/single-cell/index.html).
---

Need a short bioinformatics project idea? Follow this learning path to create new insights from old data!
