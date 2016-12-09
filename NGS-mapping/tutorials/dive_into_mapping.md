---
layout: tutorial_hands_on
topic_name: NGS-mapping
tutorial_name: dive_into_mapping
---

Dive into mapping
=================

:grey_question: ***Questions***

- *What is sequence alignment and how does it work?*
- *Why do we need it?*
- *What needs to be considered*

:dart: ***Objectives***

- *Understanding of the basic principles of sequence alignment*
- *How mapping of NGS-data works*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *Quality control of NGS data*

:hourglass: ***Time estimation*** *1h*

# Introduction

An important role in NGS data analysis is to search for the correct locations of a read in the genome. The basic concept to do this is to 
align or to map two sequences of DNA in the computer. The idea is to compute how similar two given sequences are. The more similar they are, the liklier it is
that these two sequences are having something in common. This can be a relationship in a phylogentic tree or a potential position in the genome. 

# The basic idea
Given two sequences:

```
A A C C G C C T T

A G G G G C C T T
```


To find out how similar two sequences are a similarity or a distance score can be computed.

## Similarity
A match is given if at the same position the same nucleotide in both sequences is given, e.g.
at position 1 both have an "A". The more matches two sequences are having, the more similar they are. In the example below the two
sequences are having a similarity score of 6.

```
A A C C G C C T T
|       | | | | |
A G G G G C C T T

1 2 3 4 5 6 7 8 9
```

## Distance
To compute the distance the mismatches needs to be counted. A mismatch is given if at the same position in both sequences two different nucliotides are given
e.g. at position 2 "A" in the first sequence and "G" in the second. The less mismatches are given, the better. In the example below the two
sequences are having a similarity score of 3.

```
A A C C G C C T T
  X X X 
A G G G G C C T T

1 2 3 4 5 6 7 8 9
```

The two scores can give different results if they are used to compare the quality of alignments:

```
A A C G 
|     | 
A G G G 

1 2 3 4
```

```
A A C C G
|       |
A G G G G 

1 2 3 4 5 
```

The first example is having a similarity and distance score of 2, but the second example is having a similarity score of 2 and a distance
score of 3. Which alignment is now better? In terms of the similarity both are equal, if the distance score is considered, the first alignment is better.

## Gaps
The basic idea is lacking of two important real world szenarios:

1. Sequences of different length

2. Deletion and or insertions in a sequence.

Lets have a look at the following two sequnces:

```
A A T

A A A T

1 2 3 4
```
How to handle the different sizes? Should it count as a mismatch? Or would it be more realistic to assume that at position 3 of the second sequences an additional "A"
was included respectivly a deletion at position 3 of sequence one? To overcome this issue, gaps are introduced. 

A gap is noted as "-". It is inserted in one sequence to show that at this position a deletion or insertion in one of the sequences happend:

```
A A - T
| |   |
A A A T

1 2 3 4
```
A gap counts now as a mismatch.


# NGS-data mapping
The use case for mappers in NGS-mapping is to find the position of a short sequence in the genome. This is difficult to a achieve.
Consider the example below:


```
>Genome:    A A G G G C C T C T  G  G  G C T C G G G G C C T C T A T A G C G C G C
>Sequence:  G G C C T C T C G
>Sequence:    G G C C T C T C G
>Sequence:      G G C C T C T C  G
>Sequence:        G G C C T C T  C  G
            1 2 3 4 5 6 7 8 9 10 11 12
```

1. At which position the sequence is aligned best?
2. How to handle the case multiple positions are the best?
3. What about the environment of a match?
4. Should gaps be allowed and if yes, how long they should be?

Even without a computer scientist background it is easy to see that a lot of computing resources (Runtime, memory) are required to solve this issue. 
And for a biologist even more important is the quality of the mappings.
approches handle the case different, they precompute not every match to the genome everytime from scratch. 

## State of the art algorithms 
The algorithms that are used for real-world data are all based on the idea of precomputation of an index. 

### Hashing

### BWT

# Conclusion

1. The more you know about your data the better the mapping will be.
2. Mapping is not a trivial task and it can happen that the second or tenth best alignment is the alignment you are searching for.
