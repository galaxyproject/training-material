---
layout: tutorial_hands_on
title: Analysing Solar Eclipse Frequency with a Galaxy Workflow
lang: en
level: Intermediate
questions:
- How can we identify groups of solar eclipses above a user-defined observable magnitude
  within a specified time interval?
- How can runtime parameters make the analysis reusable and reproducible?
objectives:
- Import a pre-calculated solar-eclipse dataset into Galaxy.
- Convert a text file to tab-separated tabular data.
- Split the eclipse magnitude field into observable (mag1) and theoretical (mag2)
  magnitudes.
- Filter eclipses using a user-defined minimum observable magnitude.
- Sort eclipse records chronologically.
- Identify groups of N consecutive eclipses occurring within a user-defined maximum
  number of years.
- Remove duplicate eclipse records while retaining one occurrence of each eclipse.
- Run the analysis reproducibly with workflow parameters.
requirements:
  - type: "internal"
    topic_name: digital-humanities
    tutorials:
      - introduction_to_dh
time_estimation: 45m
key_points:
- The observable magnitude mag1 is used for filtering.
- 'The workflow exposes three runtime parameters: minimum magnitude, number of eclipses,
  and maximum interval in years.'
- Galaxy workflow parameters separate the scientific question from hard-coded values.
- AWK steps can perform compact, reproducible transformations and selections on tabular
  astronomical data.
- The final result contains each qualifying eclipse once and is sorted chronologically.
contributions:
  authorship:
  - gautschr
  data:
  - dasch-swiss
  funding:
  - swissuniversities
  infrastructure:
  - GALAXY-CH
  - dasch-swiss
---

# Introduction
Ancient sources contain numerous astronomical observations that can be used to establish absolute dates. Solar eclipses are particularly important in this context. In antiquity, eclipses were often perceived as a threat and bad omen, causing fear and potentially political disturbances. By far not every partial solar eclipse has been noticed by an unprepared observer. Changes in ambient illumination remain relatively inconspicuous until a large part of the solar disk is covered. Even when about three quarters of the Sun's area is obscured, the reduction in daylight may still pass largely unnoticed with the most dramatic changes in illumination occuring only during the very deep partial phases close to totality. Historical visibility must be understood primarily in terms of the conspicuous environmental effects of a deep eclipse—diminished and altered daylight, unusual shadows, and, near totality, a pronounced twilight-like appearance. These facts make eclipses of high magnitude particularly relevant when considering events that could have attracted widespread attention without advance knowledge of the phenomenon. A series of such high magnitude solar eclipses at a certain location within a limited number of years is a rare phenomenon. Such eclipse series have been interpreted to being the cause of political disturbances.

This tutorial will show you how to identify such sequences of notable eclipses step-by-step based on pre-calculated tables which are available for the following locations relevant for classicale antiquity: Alexandria, Amarna, Assur, Athens, Babylon, Jerusalem, Knossos, Mari, Memphis, Rome, and Thebes. The data are deposited in DaSCH's repository DSP ({% cite GautschyDaSCHData %}). The dataset contains calculations of notable solar eclipses for the period between 2500 BCE and 1000 CE. It is a solar eclipse canon with identifications of the eclipses recorded in historical sources. {% cite Gautschy2012 %} discusses two examples of the solar eclipse canon's application. If you want to know more about the basics of eclipses, {% cite GautschyDaSCHDocu %} provides an introduction to eclipse geometry and the underlying calculations.

<agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Analysing Solar Eclipse Frequency with a Galaxy Workflow

Solar eclipses can be relevant to historical and chronological
investigations. For example, one may ask whether several conspicuous
eclipses occurred within a relatively short period of time. In antiquity,
eclipses were often perceived as a bad omen, as signs of the dissatisfaction
of the gods with a ruler. Therefore, an unusual sequence of
remarkable eclipses could have had drastic consequences.

In this tutorial, we use pre-calculated solar-eclipse data and a Galaxy
workflow to answer a configurable question:

> **Which eclipses with an observable magnitude at or above a chosen
> threshold belong to a group of N consecutive eclipses occurring within
> at most X years?**

The published workflow is deliberately parameterised. Instead of fixing a
particular threshold, number of eclipses, or time interval in the
workflow, the person running it supplies these values at runtime.

## The input data

The source data are pre-calculated lists of solar eclipses for locations
relevant to classical antiquity, covering 2500 BCE to 1000 CE available
in DaSCH's DSP repository. Files can be fetched from DSP or first
downloaded there and uploaded from the local computer.

Examples of suitable input files are:

- Amarna: [https://ingest.dasch.swiss/projects/0868/assets/0V3H7UR6naZ-rOIyKnuP8An/original](https://ingest.dasch.swiss/projects/0868/assets/0V3H7UR6naZ-rOIyKnuP8An/original)
- Athens: [https://ingest.dasch.swiss/projects/0868/assets/5DeOeAhy2Mf-HZBdK4jDWHS/original](https://ingest.dasch.swiss/projects/0868/assets/5DeOeAhy2Mf-HZBdK4jDWHS/original)
- Babylon: [https://ingest.dasch.swiss/projects/0868/assets/4hh1eegWr9U-SkJdggNo4fu/original](https://ingest.dasch.swiss/projects/0868/assets/4hh1eegWr9U-SkJdggNo4fu/original)
- Jerusalem: [https://ingest.dasch.swiss/projects/0868/assets/6AaA42isNKi-xuIFgPapETm/original(https://ingest.dasch.swiss/projects/0868/assets/6AaA42isNKi-xuIFgPapETm/original)
- Rome: [https://ingest.dasch.swiss/projects/0868/assets/6UlGbASVDAP-oQoqtM0sJT4/original](https://ingest.dasch.swiss/projects/0868/assets/6UlGbASVDAP-oQoqtM0sJT4/original)
- Thebes: [https://ingest.dasch.swiss/projects/0868/assets/44KbkIeJKAV-c1ZniSserLQ/original](https://ingest.dasch.swiss/projects/0868/assets/44KbkIeJKAV-c1ZniSserLQ/original)

Input files are expected to have the following structure:

`Y  M  D  Type  deltaT  mag  Tmax  Tend`

The original text file contains the following columns:

 | Column | Meaning | 
 | ------ | --------------------------------------------------- |
 | Y | Year |
 | M | Month |
 | D | Day |
 | Type | Eclipse type |
 | deltaT | Delta T |
 | mag | Two magnitude values separated by `/` |
 | Tmax | Two times of maximum eclipse values separated by `/` |
 | Tend | End time |

The `mag` and `Tmax` fields contain two values. This workflow separates the two `mag` values into:

-   **mag1**: the maximum magnitude actually observable at the selected
    location.
-   **mag2**: the theoretical maximum magnitude; the Sun may still or already be below
    the horizon when this value is reached.

For this analysis, **mag1 is the relevant magnitude**.

Relevant columns for the workflow are:
- Column 1: `Y`(year; years are in astronomical notation including a historically non-existent year 0, year -1919 corresponds to 1920 BCE)
- Column 2: `M` (month)
- Column 3: `D` (day)
- Column 6: `mag` (maximum magnitude of the solar eclipse)

# Workflow overview

The workflow analysis consists of nine processing stages:

1.  Upload / import the eclipse data.
2.  Convert the original text file to tab-separated data.
3.  Split `mag` into `mag1` and `mag2`.
4.  Select eclipses whose observable magnitude meets the user-defined
    threshold.
5.  Sort the selected eclipses chronologically.
6.  Find groups of N consecutive eclipses occurring within X years.
7.  Sort the eclipses chronologically.
8.  Remove duplicate eclipse rows while retaining one occurrence.
9.  Sort the final result chronologically.

![Visulisation of the workflow](Visualisierung_Solec_GalaxyWorkflow.png "Visulisation of the workflow")

If you run the published workflow rather than following the tutorial step-by-step on your own, you need to upload a file with the source data and to supply the following three values at runtime:

 | Workflow parameter | Example | Meaning |
 | ------------------------- | ------- | ------------------------------------------------- |
 | Minimum magnitude | 0.8 | Minimum observable magnitude (`mag1`) |
 | Number of eclipses | 3 | Number of consecutive eclipses, N |
 | Maximum interval in years | 10 | Maximum span between the first and Nth eclipse, X |

# Step 1: Import the eclipse data

You can upload data in various ways. Here are some examples:

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the desired `.txt` eclipse file from your local computer, or retrieve it directly from the [DaSCH repository](https://ingest.dasch.swiss/projects/0868/assets/0V3H7UR6naZ-rOIyKnuP8An/original):
>
>    ```
>    https://ingest.dasch.swiss/projects/0868/assets/0V3H7UR6naZ-rOIyKnuP8An/original
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>    
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
{: .hands_on}

Be sure to choose a source file whose pre-calculated lower magnitude limit does not exclude eclipses you
want to investigate. For example, a study threshold of `0.85` requires a source list calculated with a threshold of `0.8` or lower.

Once you click start, your upload should begin. It will first turn orange while in progress, then turn green once it is successfully uploaded. After this step you should have one file in your History.


# Step 2: Convert the text file to tabular data

> <hands-on-title> Convert the file </hands-on-title>
>
> 1. {% tool [Convert](Convert characters1) %} with the following parameters:
>    - {% icon param-file %} *"in Dataset"*: `Solar eclipse data (.txt)` (Input dataset)
>
>    > <comment-title> The original `.txt` file is converted into a tab-separated format suitable for subsequent Galaxy tools by changing the spaces/whitespaces to tabs. The output is a tabular dataset.</comment-title>
>
>
{: .hands_on}

# Step 3: Separate observable and theoretical magnitude

> <hands-on-title> Separate observable and theoretical magnitude </hands-on-title>
>
> The original `mag` and `Tmax` columns contain two values separated by `/`. A **Replace text** step splits the `mag` field and creates the columns
`mag1` and `mag2`. The two `Tmax` values correspond to the two `mag` values, they are split as well, but not used further during the workflow.
> 
> 1. {% tool [Replace Text in entire line](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file1` (output
of **Convert** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `mag`
>            - *"Replace with:"*: `mag1 \t mag2`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `Tmax`
>            - *"Replace with:"*: `Tmax1 \t Tmax2`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `/`
>            - *"Replace with:"*: `\t`
>    - Click *"Run Tool"*.
>
>    > <comment-title>The first two replacements prepare the header for new columns. The third replacement then actually replaces all "/" with a tab separator.</comment-title>
>
{: .hands_on}

The resulting table has the structure:

`Y  M  D  Type  deltaT  mag1  mag2  Tmax1 Tmax2  Tend`

# Step 4: Select eclipses above the magnitude threshold

> <hands-on-title> Select eclipses above the magnitude threshold </hands-on-title>
>
>The workflow uses a **Filter** step to select all eclipses with a magnitude greater or equal to a user-supplied **Magnitude threshold**.
>
> 1. {% tool [Filter](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `out_file1` (output of **Replace Text in entire line** {% icon tool %})
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `mag`
>            - *"Replace with:"*: `mag1 \t mag2`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `Tmax`
>            - *"Replace with:"*: `Tmax1 \t Tmax2`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `/`
>            - *"Replace with:"*: `\t`
>    - Click *"Run Tool"*.
>
>    > <comment-title>The first two replacements prepare the header for new columns. The third replacement then actually replaces all "/" with a tab separator.</comment-title>
>
{: .hands_on}



> ### {% icon hands_on %} Filter by observable magnitude
>
> **Workflow input:** `Magnitude threshold`
>
> Enter a numerical value, for example:
>
> `0.8`
>
> **Tool:** `Text reformatting` (AWK), version `9.5+galaxy3`
>
> ``` awk
> BEGIN {
>     OFS="\t"
>     threshold = VAR1 + 0
> }
>
> NR == 1 {
>     print
>     next
> }
>
> $6 >= threshold {
>     print
> }
> ```
>
> `VAR1` receives the value entered by the user. `$6` is `mag1`, so the
> comparison is performed inside the workflow rather than entered by the
> user.

> ### {% icon comment %} Why `+ 0`?
>
> The Galaxy AWK tool exposes its workflow variables as text values.
> Adding zero explicitly converts the threshold to a numeric value for
> the comparison.

# Step 5: Sort the eclipses chronologically

The remaining eclipses must be in chronological order before consecutive
groups can be identified.

> ### {% icon hands_on %} Sort chronologically
>
> **Tool:** `Sort`, version `1.2.0`
>
> Configure numeric ascending sorts in this order:
>
> 1.  Column 1: Year (`Y`)
> 2.  Column 2: Month (`M`)
> 3.  Column 3: Day (`D`)
>
> Header lines: `1`

This produces a chronologically ordered sequence of eclipses satisfying
the magnitude criterion.

# Step 6: Find groups of eclipses

This is the central analytical step. The user supplies:

-   **Number of eclipses (N)**
-   **Maximum interval in years (X)**

Both are connected to variables in a **Text reformatting with awk**
step.

For each possible sequence of N consecutive eclipses, the workflow
compares the year of the Nth eclipse with the year of the first eclipse:

`year[N] - year[1] <= X`

If the condition is satisfied, all eclipses in that group are returned.

> ### {% icon hands_on %} Select qualifying groups
>
> **Workflow inputs:**
>
> -   Number of eclipses: for example `3`
> -   Maximum interval in years: for example `10`
>
> **Tool:** `Text reformatting` (AWK), version `9.5+galaxy3`
>
> ``` awk
> BEGIN {
>     OFS="\t"
>     n = VAR1
>     years = VAR2
> }
>
> NR == 1 {
>     next
> }
>
> {
>     year[NR-1] = $1
>     line[NR-1] = $0
>     count = NR-1
> }
>
> END {
>     for (i=1; i<=count-n+1; i++) {
>         if (year[i+n-1] - year[i] <= years) {
>             for (j=i; j<=i+n-1; j++) {
>                 print line[j]
>             }
>         }
>     }
> }
> ```


 > <hands-on-title> Split the magnitude column </hands-on-title>
 >
 > 1. {% tool [Text
reformatting](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/9.5+galaxy3)
%} with the following parameters:
 >    - {% icon param-file %} *"File to process"*: `out_file1` (output
of **Convert** {% icon tool %})
 >    - *"AWK Program"*: `BEGIN { OFS=\t }
NR == 1 {
     print $1,$2,$3,$4,$5,mag1,mag2,$7,$8
     next
}
{
     split($6,m,/);
     print $1,$2,$3,$4,$5,m[1],m[2],$7,$8
}`
 >
{: .hands_on}

> <question-title></question-title>
>
> 1. What has to change in the AWK expression if the an import file has no header?
> 2. What has to change in the AWK expression if column 5 would contain the magnitudes and not column 6?
>
> > <solution-title></solution-title>
> >
> > 1. NR == 0
> > 2. split($5,m,/);
> >
> {: .solution}
>
{: .question}


> ### {% icon comment %} What counts as a group?
>
> With `N = 3` and `X = 10`, the workflow examines every three
> consecutive eclipses in the magnitude-filtered chronological list. A
> group qualifies when the year difference between its first and third
> eclipse is no more than ten years.

## Overlapping groups

An eclipse can belong to more than one qualifying group. For example, if
eclipses A-B-C qualify and B-C-D also qualify, B and C are emitted twice
by the group-selection step.

This is intentional: the group-selection step first identifies all
qualifying windows. The next step collapses repeated eclipse rows.

# Step 7: Remove duplicate eclipse rows

A further **Text reformatting with awk** step removes repeated rows but
retains the first occurrence of each eclipse.

> ### {% icon hands_on %} Keep one occurrence of each eclipse
>
> **Tool:** `Text reformatting` (AWK), version `9.5+galaxy3`
>
> ``` awk
> NR == 1 {
>     print
>     next
> }
>
> !seen[$0]++
> ```
>
> The expression `!seen[$0]++` uses the complete row as a key. The first
> occurrence is printed; subsequent identical rows are suppressed.

# Step 8: Sort the final result

Because the duplicate-removal step preserves first occurrence rather
than chronological position, the result is sorted once more.

> ### {% icon hands_on %} Sort the final result
>
> **Tool:** `Sort`, version `1.2.0`
>
> Sort numerically in ascending order by:
>
> 1.  Column 1: Year (`Y`)
> 2.  Column 2: Month (`M`)
> 3.  Column 3: Day (`D`)
>
> Header lines: `0`

The final dataset is therefore a chronological list of eclipses that
satisfy the magnitude threshold and belong to at least one qualifying
group.

# Run the workflow

A useful first test is:

  |  Parameter | Value | 
  |  ------------------------- | ------- | 
  |  Minimum magnitude | 0.8 | 
  |  Number of eclipses | 3 | 
  |  Maximum interval in years | 10 | 

> ### {% icon hands_on %} Run
>
> 1.  Select the eclipse data file as the workflow dataset input.
> 2.  Enter `0.8` for **Magnitude threshold**.
> 3.  Enter `3` for **Number of eclipses**.
> 4.  Enter `10` for **Maximum interval in years**.
> 5.  Run the workflow.
> 6.  Inspect the final sorted dataset.

# Interpreting the result

The workflow returns **all eclipses with `mag1` at or above the
threshold that are members of at least one group of N consecutive
qualifying eclipses occurring within the specified maximum interval**.

The result should not be interpreted as evidence that historical
observers actually saw, recorded, or interpreted every eclipse. The
workflow identifies astronomical configurations in the pre-calculated
dataset; historical interpretation requires appropriate contextual
evidence.

# Explore the parameters

The main advantage of the workflow is that the research question can be
changed without editing the workflow.

> ### {% icon hands_on %} Experiment
>
> Run the workflow several times with different parameters.
>
> For example, compare:
>
>     Minimum magnitude     N   Maximum interval
>   ------------------- ----- ------------------
>                 `0.8`   `3`               `10`
>                 `0.8`   `4`               `20`
>                 `0.9`   `3`               `20`
>
> Consider:
>
> -   How does raising the magnitude threshold affect the number of
>     qualifying eclipses?
> -   How does increasing N change the result?
> -   How sensitive are the results to the chosen time interval?
> -   Are particular periods repeatedly selected under different
>     parameter combinations?

# Why use a Galaxy workflow?

The workflow makes several parts of the analysis explicit:

-   the source dataset,
-   the transformations applied to the data,
-   the definition of observable magnitude used for filtering,
-   the chronological sorting,
-   the mathematical criterion used to define an eclipse group,
-   and the user-selected analytical parameters.

This makes the procedure easier to inspect, rerun, share, and modify
than an analysis in which these operations are performed manually.

# Conclusion

You have built and run a parameterised Galaxy workflow for astronomical
chronology. Starting with a pre-calculated eclipse list, the workflow
converts and restructures the data, selects eclipses by observable
magnitude, identifies temporally concentrated groups, removes
overlapping duplicates, and produces a chronologically sorted result.

Because the three scientific criteria are runtime parameters, the same
workflow can be reused to investigate different hypotheses without
changing its internal structure.
