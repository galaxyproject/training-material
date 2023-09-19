---
layout: tutorial_hands_on
title: FAIR and RDM from a selfish perspective
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
  GTN: Galaxy Training Network
  RDM: Research Data Management
level: Introductory
zenodo_link: ''
questions:
- Why should other researchers adhere to FAIR principles?
- How can I profit from good RDM?
objectives:
- Learn about the benefits of good data management
time_estimation: "30M"
key_points:
- RDM can lead save time and money.
tags:
- fair
- gtn
- rdm
contributions:
  authorship:
    - korneelhens
    - nomadscientist
subtopic: fair-data

requirements:
  - type: "internal"
    topic_name: fair
    tutorials:
      - fair-intro
---


# Introduction

**You** are obviously too important to bother with this research data management (RDM) nonsense. Just copy the data management plan from a previous grant application and stop wasting time.
However, there might be some value in persuading your colleagues to do their RDM properly and to adhere to FAIR standards. In this tutorial, we explore the purely selfish reasons for encouraging good data management practices based on some real-life examples.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}




## Excercise 1: Find the experiment

Imagine the following situation: An undergrad student wants to do a lab-based project with you. This is great because your post-doc has recently left the lab and there are several promising leads from their work that can be followed up. You propose the following topic to the student: 

**Effect of mutations in the gene CHES-I-like on feeding behaviour in _Drosophila melanogaster_**

The only problem is, you can't remember what your post-doc has done so far. This is no problem, as you have followed your data management plan to the letter: _"lab books will be preserved for 10 years in the PI’s office"._ Let's see how effective this strategy actually is.

<question-title></question-title>
The PDF below contains a scan of the lab books in this example. Somewhere, there is an experiment to show the effect of the knock-down of a gene called CHES-I-like on the feeding behaviour of fruit flies. Try to find this information in the lab book.

> [5 years worth of lab books](Lab-book_excercise1.pdf)

> > <solution-title></solution-title>
> > The information can be found on page XXX of the pdf.
> {: .solution}
>
{: .question}

Too difficult? Maybe try something a bit more easy

<question-title></question-title>

The information can be found on page XXX of the lab book. Can you tell me how this experiment was designed?

> [5 years worth of lab books](Lab-book_excercise1.pdf)

> > <solution-title></solution-title>
> > No you can't. Neither can I. Obviously our data management plan did not do us any favours here.
> {: .solution}
>
{: .question} 


## Excercise 2: Reproduce the experiment

Imagine the following situation: You have got great data. The editors of Nature, Science and Cell are all knocking on your door to publish in their journals. The only thing standing between you and a Nobel prize is one picture of an immunostaining that is not publication quality and the immunostaining needs to be redone. 

![Wonky immunostaining](flybrain.png "Immunostaining")

<question-title></question-title>
This is an immunostaining of a fruitfly brain. Can you list 15 more pieces of information that you would need to be able to reproduce this image?

> > <solution-title></solution-title>
> > 1. Which fly strain
> > 2. Which developmental stage
> > 3. Which genotype
> > 4. Time of day
> > 5. Feeding conditions
> > 6. Sex
> > 7. Mating condition
> > 8. Fixation time and conditions
> > 9. Primary antibody
> > 10. Secondary antibody
> > 11. Dilutions of antibodies
> > 12. Staining protocol
> > 13. Which microscope
> > 14. Which laser lines
> > 15. What exposure time
> > 16. Filter settings
> > 17. Magnification
> > 18. Pinhole size
> > 19. Number of z slices
> > 20. ...
> > 
> {: .solution}
>
{: .question}

Too difficult? Maybe try something a bit more easy



## Excercise 4: Get the data


## Conclusion

The Galaxy Training Network is an example of a robust, effective Community of Practice.
For more information please look at this great article {% cite hiltemann2023galaxy %}, the corresponding FAIR guidelines {% cite fair-training-materials %} and follow [short introduction to FAIR data stewardship](http://fellowship.elixiruknode.org/).
