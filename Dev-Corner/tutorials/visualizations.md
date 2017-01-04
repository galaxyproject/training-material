Visualizations in Galaxy
============================================

:grey_question: ***Questions***

- *Major question that would be addressed in this tutorial (mostly general biological questions)*
- *Second question*
- *Third question*
- *...*

:dart: ***Objectives***

- *First objective of this tutorial (It is a single sentence describing what a learner will be able to do once they have sat through the lesson. The objectives must be technical, but also theoretical, objectives. You can check [SWC lessons](http://swcarpentry.github.io/instructor-training/19-lessons/) to help you writing learning objectives.)*
- *Second objective*
- *Third objective*
- *...*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *Javascript Knowledge*
- *...*

:hourglass: ***Time estimation*** *1d/3h/6h*

# Introduction

Visualizations may be very helpful in understanding data better.
These can be simple barplots but also projections of high dimensional data or even genomes.
What is characteristical about visualizations is that they often require endless tweaking and change in settings like zoom and colors.
Most often visualizations are ideally interactive, where changing the settings is effectively the browsing through the data.
For this reason it may be inconvenient to make use of a static galaxy tool because it misses interactive features.
Therefore galaxy offers the option to create *visualizations plugins*, file format specific javascripts that intergrate with the history menu.

In this tutorial we shall go through how this works and create a simple tool ourselves.
The tool will create a visualization of the number of aligned reads per chromosome in a BAM file, and we will discuss possible optimizations.

Because visualizations are written in HTML5 and JavaScript it is essential to have a good understand if you want to make visualizations ready for production.
However, for this tutorial we will keep it very basic.


Documentation about Galaxy visualizations can be found here. It can be used as base:

- https://wiki.galaxyproject.org/VisualizationsRegistry
- https://wiki.galaxyproject.org/VisualizationsRegistry/Cookbook
- https://wiki.galaxyproject.org/VisualizationsRegistry/Code
- https://wiki.galaxyproject.org/DataProviders
- https://wiki.galaxyproject.org/DataProviders/Cookbook
- https://wiki.galaxyproject.org/Develop/Visualizations

# Part 1

The visualization we want to write is a tool that shows the number of aligned reads per chromosome, of a BAM file.
The first thing we need to do is to come up with a name.
Let's call it **.

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
