Dive into quality control
=========================

:grey_question: ***Questions***

- *How to control quality of NGS data?*
- *What are the quality parameters to check for each dataset?*
- *How to improve the quality of a sequence dataset?*

:dart: ***Objectives***

- *Manipulate FastQ files*
- *Control quality from a FastQ file*
- *Use FastQC tool*
- *Understand FastQC output*
- *Use tools for quality correction*

:heavy_check_mark: ***Requirements***

- *Galaxy introduction*
- *Second requirement*
- *Third requirement*
- *...*

:hourglass: ***Time estimation*** *1d/3h/6h*

# Introduction

General introduction about the topic and then an introduction of the tutorial (the questions and the objectives). It is nice also to have a scheme to sum up the pipeline used during the tutorial. The idea is to give to trainees insight into the content of the tutorial and the (theoretical and technical) key concepts they will learn.

# Get your sequence dataset

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

# Check the quality of your sequences

Chacking sequences quality, and raw data quality in general, is one of the most important part of an analysis. Taking time to evaluate deeply each part of your data can give a lot of benefit on the next steps and avoids loosing time re-executing several times the entire workflow.

## Run FastQC

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

## Check the FastQC output

FastQC outputs are an html report and a raw data file. If the raw data file allows the user to start from original data to create its own grpahs and stats and the html report shows directly a graphical web page, the two outputs are divided in sections :

:pencil2: ***Basics statistics***

We find here simple statistics on the analysed file. We can thus find here:
-the ASCII encoding type, relevant when using FastQ groomer next to FastQC to convert your fastQ files to the fastqsanger format.
-number of sequences
-sequence length
-the overall GC content

:pencil2: ***Per base sequence quality ***

A red line indicates the mean value of the quality score when the yellow box allows the user to determine inter-quartiles values (25% - 75%). box plot top and bottom show respectively the 90% and 10% treshold. The 3 parts of the graph (green, orange, red), represent respectively very good, quite good and not so good ;) quality. It's classical to see the overall score decreasing with le sequence length as sequencing qulity is decreasing with time.

A warning indicates that one of the first quartile is under 10 or if a median value is under 25. 
An error indicates that one of the first quartile is under 5 or if a median value is under 20.

:pencil2: ***Per sequence quality score ***

Allows to detect if a subset of sequences show weak quality. This can for example happen when sequencing step have difficulties to capture sequences on images because on a corner of the field of view. Thus, this can indicate a systematic problem on a flowcell part.

A warning indicates an error rate of 0.2% when the peak is under 27. 
An error indicates an error rate of 1% when the peak is under 20.

:grey_exclamation: ***Key Points***

- *Parameters to check during quality control*
    1. *Parameter*
    2. *Parameter*
    3. *...*

# Improve the quality of your sequences

Short introduction about this subpart.

:pencil2: ***Hands on!***

1. First step
2. Second step
3. Third step

## Run TrimGalore or PRINSEQ?

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

## Check the improvement in quality score

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

- *Quality control steps*
    1. *Check the quality of your sequences*
    2. *Improve the quality of your sequences*
    3. *Check the quality improvement*

# :clap: Thank you
