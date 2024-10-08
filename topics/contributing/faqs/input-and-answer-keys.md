---
title: Input Histories & Answer Keys
area: gtn
box_type: tip
layout: faq
contributors: [hexylena, nomadscientist]
---

Tutorials sometimes require significant amounts of data or data prepared in a very specific manner which often is shown to cause errors for learners that significantly affect downstream results. Input histories are an answer to that:

```yaml
input_histories:
  - label: "UseGalaxy.eu"
    history: https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/cs1pre-processing-with-alevin---input-1
    date: "2021-09-01"
```

Additionally once the learner has gotten started, tutorials sometimes feature tools which produce stochastic outputs, or have very long-running steps. In these cases, the tutorial authors may provide answer histories to help learners verify that they are on the right track, or to enable them to catch up if they fall behind or something goes wrong.

```yaml
answer_histories:
    - label: "UseGalaxy.eu"
    history: https://humancellatlas.usegalaxy.eu/u/j.jakiela/h/generating-a-single-cell-matrix-using-alevin-3
    - label: "Older Alevin version"
    history: https://humancellatlas.usegalaxy.eu/u/wendi.bacon.training/h/cs1pre-processing-with-alevin---answer-key
    date: 2024-01-01
```
