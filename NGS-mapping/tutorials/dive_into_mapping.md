Dive into mapping
=================

:grey_question: ***Questions***

- *What is sequence alignment and how does it work?*
- *Why do we need it?*
- *What are typical troubles?*

:dart: ***Objectives***

- *First objective of this tutorial (It is a single sentence describing what a learner will be able to do once they have sat through the lesson. The objectives must be technical, but also theoretical, objectives. You can check [SWC lessons](http://swcarpentry.github.io/instructor-training/19-lessons/) to help you writing learning objectives.)*
- *Second objective*
- *Third objective*
- *...*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *Quality control of NGS data*
- *Third requirement*
- *...*

:hourglass: ***Time estimation*** *1d/3h/6h*

# Introduction

An important role in NGS data analysis is to search for the correct locations of a read in the genome. The basic concept to do this is to 
align or to map two sequences of DNA in the computer. The idea is to compute how similar two given sequences are. The more similar they are, the liklier it is
that these two sequences are having something in common. This can be a relationship in a phylogentic tree or a potential position in the genome. 

# The basic idea
Given two sequences:
AACCGCCTT
AGGGGCCTT

To find out how similar two sequences are a similarity or a distance score can be computed.

## Similarity
A match is given if at the same position the same nucleotide in both sequences is given, e.g.
at position 1 both have an "A". The more matches two sequences are having, the more similar they are.

## Distance
To compute the distance the mismatches needs to be counted. A mismatch is given if at the same position in both sequences two different nucliotides are given
e.g. at position 2 "A" in the first sequence and "G" in the second. The less mismatches are given, the better.


## Gaps
The basic idea is lacking of two important real world szenarios:
1. Sequences of different length

2. Deletion and or insertions in a sequence.

Lets have a look at the following two sequnces:
AAT
AAAT
How to handle the different sizes? Should it count as a mismatch? Or would it be more realistic to assume that at position 3 of the second sequences an additional "A"
was included respectivly a deletion at position 3 of sequence one? To overcome this issue, gaps are introduced. 

A gap is noted as "-". It is inserted in one sequence to show that at this position a deletion or insertion in one of the sequences happend:
AA-T
AAAT

A gap counts now as a mismatch.
# Parameters, parameters, parameters

# Quality

# State of the art algorithms
The algorithms that are used for real-world data are all based on the idea of precomputation of an index. 

## Hashing

## BWT


H

Introduction about this part

## Subpart 1

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

## Subpart 2

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Part 2

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

## Subpart 2

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

Some blabla

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Conclusion

Conclusion about the technical key points. And then relation between the technics and the biological question to end with a global view.

:grey_exclamation: ***Key Points***

- *Simple sentence to sum up the first key point of the tutorial (Take home message)*
- *Second key point*
- *Third key point*
- *...*

# :clap: Thank you
