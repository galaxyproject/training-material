---
title: How does the GTN implement the "Ten simple rules for collaborative lesson development"
area: contributors
layout: faq
box_type: tip
contributors: [hexylena]
---


The GTN framework is inherently collaborative and community-driven, and comprises a growing number of contributors with expertise in a wide range of scientific and technical domains. Given this highly collaborative nature of a community with very different skill sets, the GTN framework has evolved over the years to facilitate the contribution and maintenance of the tutorials. We aim to adhere to best-practice guidelines for collaborative lesson development described in {% cite Devenyi_2018 %}. The structure of the tutorials and repository has been made modular with unified syntax and use of snippets enabling easy access for authors to add common tips and tricks new users might need to know. This system allows for easy updating of all tutorials, if there is a change in tools or interface. More generally, we continually strive to lower contribution barriers for content creators by providing a framework that is easy to use for training developers regardless of their level of knowledge of the underlying technical framework.

Implementation of the "Ten simple rules for collaborative lesson development" ({% cite Devenyi_2018 %}) in the training material:

Rules                                            | Implementation in the GTN framework
----                                             | ---
Clarify audience                                 | Tutorial metadata includes level indicators (introductory, intermediate, advanced) and a list of prerequisite tutorials as recommended prior knowledge. This information is rendered at the top of each tutorial.
Make lessons modular                             | Development of small tutorials linked together via learning paths
Teach best practice lesson development           | We maintain the topic [*Contributing to the Galaxy Training Material*]({% link topics/contributing/index.md %}) including numerous tutorials describing how to create new content. Furthermore, quarterly online collaboration fest (CoFests) are organized, where contributors can get direct support. Development of a Train the Trainer program and a mentoring program for instructors, in which lesson development is taught
Encourage and empower contributors               | Involve them in reviews. Mentor them. Encourage them to become maintainers.
Build community around lessons                   | Quarterly online collaboration fest (CoFests) and Community calls. Chat on our Gitter/Matrix channel.
Publish periodically and recognize contributions | Author listed on tutorials. Hall of fame listing all contributors. Full tutorial citation at the end of the tutorial. Tweet about new or updated tutorials. List of new or updated tutorials in Galaxy Community newsletter. Soon: publication of tutorials via article
Evaluate lessons at several scales               | Tutorial change (Pull Request) review. Embedded feedback form in tutorials for trainee feedback. Instructor feedback. Automatic workflow testing
Reduce, re-use, recycle                          | Sharing content between tutorials, specially using snippets. Development of small modular tutorials linked by learning paths
Link to other resources                          | Links to original paper, documentation, external tutorials and other material
You can't please everyone                        | but we can try (several different Galaxy introduction tutorials for different audience). Aim to clearly state what the tutorial does and does not cover, at the start.
