---
layout: tutorial_hands_on

title: How to turn ancient calendar sources into absolute dates
zenodo_link: ''
questions:
- How can a sequence of preserved lunar observations be turned into candidate absolute dates for an ancient Egyptian ruler?
- Why is the *last* visibility of the lunar crescent the relevant event for the ancient Egyptian calendar?
- How can a chain of generic Galaxy text-manipulation tools implement a chronology pipeline?
objectives:
- Restrict a four-millennia table of computed lunar visibilities to the historically plausible window of one ruler.
- Convert a lunar feast date reported in the Egyptian civil calendar into a search term for the last-visibility table using the −2-day rule.
- Combine the searches for several preserved lunar dates and read the chronologically sorted output as a list of viable dating options.
time_estimation: 1H
key_points:
- A reported first lunar day (`psḏntyw`) is sought in the last-visibility table two days earlier (the −2-day rule).
- ΔT (column 3) decreases monotonically with time, so sorting on it orders the candidate dates chronologically.
- Each preserved lunar date yields several candidate absolute dates; intersecting several of them pins down a single chronology (here: year 1 of Senwosret III ≈ 1872 BCE).
contributors:
- Rita Gautschy
- Johannes Nussbaum

---


# Introduction

This tutorial shows how a sequence of preserved ancient lunar observations can be turned into
candidate **absolute calendar dates** for an ancient ruler. The worked example follows the
chronology of the Egyptian Middle Kingdom (12th Dynasty) and the reign of pharaoh
**Senwosret III**, based on the lunar dates from the temple archive of el-Lahun (Illahun)
{% cite Gautschy2011 %}.

## The astronomical background

In many cultures, the first or last visibility of the Moon close to New Moon announced the
beginning of a new lunar month. Ancient Egypt was unusual: there, the new lunar month began with
the **invisibility of the old crescent** — the day called *psḏntyw*. For that reason the relevant
table for Egypt is the one of **last visibilities** of the waning crescent before New Moon
(`Thebenletzte` = *Theben, letzte Sichtbarkeit*, "Thebes, last visibility").

The input data are precomputed last visibilities of the Moon at Thebes, one row per lunar month,
covering the four millennia from 2000 BCE to 2000 CE. Each row gives the date in the Julian
calendar, the time, ΔT (the difference between dynamical and universal time), and — crucially — the
corresponding date in the **Egyptian civil calendar**. The Egyptian civil year was strictly
schematic: three seasons (*akhet*, *peret*, *shemou*) of four 30-day months each, plus five
epagomenal days, giving 365 days every year with no leap day. In this table a civil date is written
as `<month 1–4> <season> <day>`, e.g. `1 peret    5` for "month I of peret, day 5".

The table can be retrieved as a text file (`Thebenletzte.txt`) from Rita Gautschy's collection of
lunar data {% cite Gautschy-mond %} and converted into the tab-separated `Thebenletzte.tsv` used here
with [OpenRefine](https://openrefine.org/).

## From a preserved date to an absolute date

A preserved lunar feast date (a *psḏntyw*, the start of a lunar month) is recorded in the Egyptian
civil calendar, but its **absolute** (Julian) date is unknown — that is exactly what we want to
recover. Because the Egyptian day began at dawn and the lunar month began with the invisibility of
the old crescent, the **last visibility falls two days earlier** than the reported first lunar day:

> If a first lunar day is reported on I peret 7, one has to seek in the table of the last
> visibilities for the date I peret 5.
{: .quote cite="{% cite_url Gautschy-mond %}"}

So a reported feast date is converted into a **search term = reported date − 2 days**, which is then
searched in the last-visibility table. Every matching row is a **candidate absolute date**. A single
preserved date yields several candidates (because the civil calendar drifts through the seasons over
the centuries); combining the searches for **several** preserved dates and keeping only chronologies
consistent with the known regnal years pins the sequence down to one absolute chronology.

## The pipeline

1. Because the original file spans four millennia, it is first cut down to the historically
   reasonable window for the ruler: lines are removed from the beginning
   (**Remove beginning of a file**) and the end (**Select first lines from a dataset**).
2. The window is then searched for each preserved lunar date (**Search in textfiles (grep)**).
3. The individual search results are concatenated
   (**Concatenate multiple datasets or collections**).
4. Duplicate rows are removed (**Unique occurrences of each record**).
5. Finally the rows are sorted chronologically (**Sort data in ascending or descending order**).

The result is a text file listing every remaining viable dating option for the preserved lunar
dates, given simultaneously in the Julian and the Egyptian civil calendar.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Get data

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import `Thebenletzte.tsv` (the precomputed last visibilities of the Moon at Thebes) from
>    [Zenodo]({{ page.zenodo_link }}) or from the shared data library (`GTN - Material` ->
>    `{{ page.topic_name }}` -> `{{ page.title }}`):
>
>    ```
>    Thebenletzte.tsv
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the dataset to `Thebenletzte.tsv` if needed
> 4. Check that the datatype is set to `tabular`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

> <comment-title> The structure of the input file </comment-title>
>
> Each row is one computed last visibility of the lunar crescent at Thebes. The columns are:
> `date | best time | deltat[s] | sigma | Egyptian date | TSun | TMoon | q | code | new moon date | new moon time`.
> Two columns matter for this tutorial:
> - **column 5 — Egyptian date**, written as `<month 1–4> <season> <day>` (e.g. `1 peret    5`),
>   which we search for; and
> - **column 3 — ΔT** (seconds), which decreases monotonically with time and is therefore used to
>   sort the results chronologically.
>
> The file covers 2000 BCE to 2000 CE, so the year in column 1 is negative for BCE dates
> (e.g. `-1872` ≈ 1872 BCE).
{: .comment}

# Restrict the table to the ruler's window

The input table spans four millennia. Searching all of it would return matches from every epoch,
most of them historically meaningless. We therefore first cut the table down to the window in which
the chosen ruler could plausibly have reigned.

For this example we focus on **Senwosret III** of the 12th Dynasty. The el-Lahun lunar dates are
usually placed in the 19th–18th centuries BCE; Gautschy's best-fitting solution puts **year 1 of
Senwosret III around 1872 BCE** {% cite Gautschy2011 %}. To be safe we keep a generous window of
roughly **1900 BCE to 1748 BCE**, which comfortably brackets Senwosret III and his successor
Amenemhat III.

The two cuts below are tuned to this file: in `Thebenletzte.tsv`, line 1240 is the first row of
1900 BCE, and the next 1900 rows reach down to about 1748 BCE.

## Remove the lines before the window

> <hands-on-title> Cut off the early millennia </hands-on-title>
>
> 1. {% tool [Remove beginning](Remove beginning1) %} with the following parameters:
>    - *"Remove first"*: `1240`
>    - {% icon param-file %} *"from"*: `output` (Input dataset)
>
>    This drops everything before about 1900 BCE.
>
{: .hands_on}

> <question-title></question-title>
>
> 1. Why do we remove lines by *count* (1240) rather than by date?
>
> > <solution-title></solution-title>
> >
> > 1. The tools used here (**Remove beginning** / **Select first**) are generic text tools that
> >    operate on line numbers, not on calendar values. Because the file has exactly one row per
> >    lunar month in chronological order, a line count is a reliable proxy for a date: line 1240 is
> >    the first row of 1900 BCE.
> >
> {: .solution}
>
{: .question}

## Keep only the window

> <hands-on-title> Cut off everything after the window </hands-on-title>
>
> 1. {% tool [Select first](Show beginning1) %} with the following parameters:
>    - *"Select first"*: `1900`
>    - {% icon param-file %} *"from"*: `out_file1` (output of **Remove beginning** {% icon tool %})
>
>    Keeping the first 1900 of the remaining rows reaches down to about 1748 BCE, leaving a
>    ~150-year window around the reign of Senwosret III.
>
{: .hands_on}

# Search for the preserved lunar dates

We now search the ~150-year window for each preserved lunar date. Recall the **−2-day rule**: a
*psḏntyw* reported on a given civil date is sought in the last-visibility table **two days earlier**.

Each search term is written as a small regular expression, because the **Search in textfiles** tool
in this workflow runs `grep` with the `-i` (case-insensitive) and `-P` (Perl regex) options:

- `\s+` matches the variable alignment spaces between the season name and the day in column 5;
- `\b` after the day stops a search for day `5` from also matching `15`, `25` or `50`.

So a reported `psḏntyw` on **I peret 7** becomes the search term `1 peret\s+5\b` (I peret 7 − 2 days
= I peret 5). We keep **0** preceding and **0** trailing lines, because every matching row is itself
one candidate date — context lines would only pollute the result.

> <comment-title> Where the search terms come from </comment-title>
>
> The first term below is the worked example documented by {% cite Gautschy2011 %}. The other two
> are illustrative *psḏntyw* dates chosen to demonstrate the method; they are verified to occur in
> `Thebenletzte.tsv` but are **not** transcribed from a specific papyrus. When adapting this tutorial
> to real research, replace them with the exact attested lunar dates of the ruler you are studying.
{: .comment}

## First lunar date — `I peret 5`

> <hands-on-title> Search for the first preserved date </hands-on-title>
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `out_file1` (output of **Select first** {% icon tool %})
>    - *"Regular Expression"*: `1 peret\s+5\b`
>    - *"Show lines preceding the matched line"*: `0`
>    - *"Show lines trailing the matched line"*: `0`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. How many candidate dates does this search return, and what are they?
> 2. Why does one civil date match several different years?
>
> > <solution-title></solution-title>
> >
> > 1. Six rows, in the years −1896, −1835, −1821, −1810, −1796 and −1771 (i.e. 1896–1771 BCE).
> > 2. The Egyptian civil year has a fixed length of 365 days and ignores the ~365.25-day solar
> >    year, so it slowly drifts through the seasons. A given civil date therefore coincides with the
> >    last visibility of the crescent only in particular years, recurring at intervals as the
> >    calendar drifts. The regnal year recorded with the papyrus, together with the other lunar
> >    dates, selects which of these candidates is the historical one.
> >
> {: .solution}
>
{: .question}

## Second lunar date — `II shemou 18`

> <hands-on-title> Search for the second preserved date </hands-on-title>
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `out_file1` (output of **Select first** {% icon tool %})
>    - *"Regular Expression"*: `2 shemou\s+18\b`
>    - *"Show lines preceding the matched line"*: `0`
>    - *"Show lines trailing the matched line"*: `0`
>
>    This corresponds to a reported `psḏntyw` on II shemou 20 (− 2 days). It returns four candidate
>    years: −1878, −1853, −1828 and −1803.
>
{: .hands_on}

## Third lunar date — `I shemou 7`

> <hands-on-title> Search for the third preserved date </hands-on-title>
>
> 1. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Select lines from"*: `out_file1` (output of **Select first** {% icon tool %})
>    - *"Regular Expression"*: `1 shemou\s+7\b`
>    - *"Show lines preceding the matched line"*: `0`
>    - *"Show lines trailing the matched line"*: `0`
>
>    This corresponds to a reported `psḏntyw` on I shemou 9 (− 2 days). It returns five candidate
>    years: −1877, −1863, −1852, −1838 and −1752.
>
{: .hands_on}

# Combine, deduplicate and sort the candidates

## Concatenate the three searches

> <hands-on-title> Merge the search results into one table </hands-on-title>
>
> 1. {% tool [Concatenate multiple datasets or collections](cat1) %} with the following parameters:
>    - {% icon param-file %} *"Concatenate Dataset"*: `output` (output of the **second** search, `II shemou 18`)
>    - In *"Dataset"*:
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `output` (output of the **first** search, `I peret 5`)
>        - {% icon param-repeat %} *"Insert Dataset"*
>            - {% icon param-file %} *"Select"*: `output` (output of the **third** search, `I shemou 7`)
>
{: .hands_on}

## Remove duplicate rows

> <hands-on-title> Keep only unique candidate dates </hands-on-title>
>
> 1. {% tool [Unique](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sorted_uniq/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"File to scan for unique values"*: `out_file1` (output of **Concatenate multiple datasets or collections** {% icon tool %})
>    - *"Advanced Options"*: `Hide Advanced Options`
>
>    Identical rows that happen to be returned by more than one search are collapsed into one.
>
{: .hands_on}

## Sort the candidates chronologically

> <hands-on-title> Order the candidate dates by time </hands-on-title>
>
> 1. {% tool [Sort](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sort_header_tool/9.5+galaxy2) %} with the following parameters:
>    - {% icon param-file %} *"Sort Query"*: `outfile` (output of **Unique** {% icon tool %})
>    - In *"Column selections"*:
>        - {% icon param-repeat %} *"Insert Column selections"*
>            - *"on column"*: `c3`
>            - *"in"*: `Descending order`
>
{: .hands_on}

> <question-title></question-title>
>
> 1. The data are sorted on **column 3 (ΔT)** in *descending* order. Why does that order the
>    candidate dates chronologically, from oldest to youngest?
>
> > <solution-title></solution-title>
> >
> > 1. ΔT (the difference between dynamical time and universal time, in seconds) shrinks steadily as
> >    we approach the present — from roughly 46,000 s around 2000 BCE to about 64 s in 2000 CE. It is
> >    therefore a monotonic proxy for absolute time. Sorting on it in *descending* order puts the
> >    largest ΔT (the oldest date) first, giving a clean chronological list without having to parse
> >    the calendar dates themselves.
> >
> {: .solution}
>
{: .question}

# Reading the result

The final dataset is the deduplicated union of the three searches, ordered chronologically — 15
candidate dates spanning 1896–1752 BCE. Each row gives the candidate's Julian date (column 1) and its
Egyptian civil date (column 5):

```
13. 4.-1896   ...   1 peret    5
19. 9.-1878   ...   2 shemou  18
 9. 8.-1877   ...   1 shemou   7
 5. 8.-1863   ...   1 shemou   7
13. 9.-1853   ...   2 shemou  18
 2. 8.-1852   ...   1 shemou   7
30. 7.-1838   ...   1 shemou   7
29. 3.-1835   ...   1 peret    5
 6. 9.-1828   ...   2 shemou  18
26. 3.-1821   ...   1 peret    5
23. 3.-1810   ...   1 peret    5
31. 8.-1803   ...   2 shemou  18
19. 3.-1796   ...   1 peret    5
13. 3.-1771   ...   1 peret    5
 8. 7.-1752   ...   1 shemou   7
```

This list is the workflow's contribution to the chronological problem: every astronomically possible
absolute date for each preserved lunar observation. Choosing the historically correct one is the
historian's step — only one combination of these candidates is mutually consistent with the regnal
years recorded on the papyri and with the known reign length of the ruler. For Senwosret III, that
consistent solution places his year 1 around 1872 BCE {% cite Gautschy2011 %}.

# Conclusion

Starting from a four-millennia table of computed last visibilities of the Moon, we restricted it to a
ruler's plausible window, searched it for several preserved lunar feast dates (each converted with
the −2-day rule), and combined the results into a single chronologically sorted list of candidate
absolute dates. The same pipeline of generic Galaxy text tools can be reused for any sequence of
preserved lunar observations and any ruler, simply by adjusting the window (the two line cuts) and
the search terms.