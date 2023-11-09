---
layout: tutorial_hands_on
title: FAIR in a nutshell
abbreviations:
  FAIR: Findable, Accessible, Interoperable, Reusable
zenodo_link: ''
questions:
- What do the FAIR principles stand for [Wilkinson et al. 2016](https://www.nature.com/articles/sdata201618)?
- How to make data FAIR? 
objectives:
- Learn the FAIR principles
- Recognise the relationship between FAIR and Open data
time_estimation: "10M"
key_points:
- FAIR data are data which meet principles of findability, accessibility, interoperability, and reusability (FAIR).
- FAIR data are as open as possible, and as closed as necessary.
- The main objective of FAIR data is to increase data reuse by researchers. 
tags:
- fair
- open
- data stewardship
priority: 1
contributions:
  authorship:
    - kkamieniecka
    - poterlowicz-lab
  editing:
    - hexylena
  funding:
      - ELIXIR-UK-DaSH
subtopic: fair-data

---


The FAIR (Findable, Accessible, Interoperable, Reusable) principles emphasize machine-actionability. The main objective of FAIR is to increase data reuse by researchers. The core concepts of the FAIR principles are based on good scientific practice and intuitively grounded. 

This tutorial is a short introduction to the FAIR principles and their origin. You can find out more at the [FAIR Pointers](https://elixir-uk-dash.github.io/FAIR-Pointers/ep1/index.html) online course.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# FAIR and its origins

The FAIR Guiding Principles aid in designing data publishing platforms for easier manual and automated deposition, exploration, sharing, and reuse {% cite wilkinson %}. FAIR stands for specific improvements in data management and archival practices as it outlines clear, high-level, domain-independent principles that can be used to create a variety of scholarly outputs:


| The FAIR Guiding Principles |                                                                                                                                                                                                                                                                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| To be Findable:             | F1. (meta)data are assigned a globally unique and persistent identifier<br>F2. data are described with rich metadata (defined by R1 below)<br>F3. metadata clearly and explicitly include the identifier of the data it describes <br>F4. (meta)data are registered or indexed in a searchable resource                                   |
| To be Accessible:           | A1. (meta)data are retrievable by their identifier using a standardized communications protocol <br>A1.1 the protocol is open, free, and universally implementable<br>A1.2 the protocol allows for an authentication and authorization procedure, where necessary <br>A2. metadata are accessible, even when the data are no longer available |
| To be Interoperable:        | I1. (meta)data use a formal, accessible, shared, and broadly applicable language for knowledge representation. <br>I2. (meta)data use vocabularies that follow FAIR principles<br>I3. (meta)data include qualified references to other (meta)data                                                                                         |
| To be Reusable:             | R1. meta(data) are richly described with a plurality of accurate and relevant attributes <br>R1.1. (meta)data are released with a clear and accessible data usage license<br>R1.2. (meta)data are associated with detailed provenance<br>R1.3. (meta)data meet domain-relevant community standards                                        |

Table 1: The FAIR guiding principles as described in Wilkinson, M., Dumontier, M., Aalbersberg, I. et al. The FAIR Guiding Principles for scientific data management and stewardship {% cite wilkinson %}.

A report from the European Commission Expert Group on [FAIR data](https://zenodo.org/record/1285272#.ZGc58exByha) describes the origins of FAIR and its development in 2014-2015 by a [FORCE11](https://force11.org/groups/) Working Group. The following exercise dips into this report and asks you to investigate some of FAIR’s history and foundation.

## Open data and FAIR
The level of accessibility and usability criteria distinguish open from FAIR data. Open data is accessible without limitations, whereas FAIR data specifies certain requirements for access and use.

![text reading fair does not equal open](../../images/fair_open.png)

Open data can be modified, and distributed for any reason. Although extensively used and accessible, FAIR data additionally includes the following usability standards that go beyond permission alone: 

- In order to be found and cited, FAIR data must be identified and deposited into online public records. 

- It is necessary to make FAIR data available so that it may be accessed, read, and processed. 

- FAIR data needs to be captured and presented in a form that can be used.

The published FAIR Guiding Principles: Wilkinson, M., Dumontier, M., Aalbersberg, I. et al. The FAIR Guiding Principles for scientific data management and stewardship {% cite wilkinson %}.

## FAIRification and FAIRness of data
Making your data FAIR compatible involves adopting the 15 guiding principles listed in Table 1, which is the process known as "FAIRification." The degree to which you adhere to these criteria determines how FAIR your data are. In other words, FAIRness refers to how FAIR your data is to a certain level and indicates a tacit method of evaluating its compliance.

Documentation and frameworks for data FAIRification. Each of the 15 FAIR principles is put into context with real data examples: [GO FAIR](https://www.go-fair.org/fair-principles/), [FAIR Cookbook](https://faircookbook.elixir-europe.org/content/home.html).

# Conclusion
The FAIR Principles place a strong focus on encouraging individual data reuse while also supporting the capacity of machines to automatically discover and utilise the data. This short introduction aims to develop and disseminate guidance and processes needed to make and keep data FAIR. 
