---
title: Datasets versus collections
description: Explanation of why collections are needed and what they are
area: collections, histories
box_type: tip
layout: faq
contributors: [nekrut]
---

**Datasets versus collections**

In Galaxy's history datasets can be present as individual entries or they can be combined into *Collections*. Why do we need collections? Collections combine multiple individual
datasets into a single entity which is easy to manage. Galaxy tools can use collections directly as inputs. Collection can be **simple** or **nested**.

**Simple collections**

Imagine that you've uploaded a hundred FASTQ files corresponding to a hundred samples. These will appear as a hundred individual datasets in your history making it very long.
But the chances are that when you analyze these data you will do the same thing on each dataset.

To simplify this process you can combine all hundred datasets into a single entity called a *dataset collection* (or simply a *collection* or a *list*). It will appear as a single box in your history making it much easier to understand. Galaxy tools are designed to take collections as inputs. So, for example, if you want to map each of these datasets against a reference genome using, say, {% tool Minimap2 %}, you will need to provide `minmap2` with just one input, the collection, and it will automatically start 100 jobs behind the scenes and will combine all outputs into a single collection containing BAM files.

![A simple collection is a container containing individual datasets]({% link shared/images/simple_collection.svg %})

There is a number of situations when simple collections are not sufficient to reflect the complexity of the data. To deal with this situation Galaxy allows for **nested** collections.

**Nested collections**

Probably the most common example of this is pared end data when each sample is represented by two files: one containing forward reads and another containing reverse reads. In Galaxy you can create **nested** collection that reflects the hierarchy of the data. In the case of paired data Galaxy supports **paired** collections.

![A paired collection is a container containing individual datasets and preserving their hierarchy]({% link shared/images/paired_collection.svg %})

<!-- Original editable image for simple collections = https://docs.google.com/drawings/d/1A-tRerNLzC4FJfShUFT327wMvSSX4Y8AAdaxD_Fwaa0/edit?usp=sharing -->
<!-- Original editable image for paired collection = https://docs.google.com/drawings/d/1Bbx4UmIYdDAqK3KSm6LtLQ8zwEXxDglHmVbb0MK-mbQ/edit?usp=sharing -->
