---
title: NCBI SRA sourced fastq data
area: data upload
box_type: tip
layout: faq
contributors: [jennaj, Melkeb]
---

In these FASTQ data:

- The quality score identifier (+) is sometimes not a match for the sequence identifier (@).
- The forward and reverse reads may be interlaced and need to be separated into distinct datasets.
- Both may be present in a dataset. Correct the first, then the second, as explained below.
- Format problems of any kind can cause tool failures and/or unexpected results.
- Fix the problems before running any other tools (including **FastQC**, **Fastq Groomer**, or other QA tools)

**For inconsistent sequence (@) and quality (+) identifiers**

- Correct the format by running the tool **Replace Text in entire line** with these options:

  - Find pattern: `^\+SRR.+`
  - Replace with: `+`

Note: If the quality score line is named like "+ERR" instead (or other valid options), modify the pattern search to match.

**For interlaced forward and reverse reads**

**Solution 1 (reads named /1 and /2)**

- Use the tool **FASTQ de-interlacer on paired end reads**

**Solution 2 (reads named /1 and /2)**

- Create distinct datasets from an interlaced fastq dataset by running the tool **Manipulate FASTQ reads on various attributes** on the original dataset. It will run twice.

Note: The solution does NOT use the **FASTQ Splitter** tool. The data to be manipulated are interlaced sequences. This is different in format from data that are joined into a single sequence.

- Use the Manipulate FASTQ settings to produce a dataset that contains the `/1` reads**

  *Match Reads*

    - Match Reads by `Name/Identifier`
    - Identifier Match Type `Regular Expression`
    - Match by `.+/2`

  *Manipulate Reads*

    - Manipulate Reads by `Miscellaneous Actions`
    - Miscellaneous Manipulation Type `Remove Read`

- Use these Manipulate FASTQ settings to produce a dataset that contains the `/2` reads**

  - Exact same settings as above except for this change: Match by `.+/1`

**Solution 3 (reads named /1 and /3)**

- Use the same operations as in Solution 2 above, except change the first **Manipulate FASTQ** query term to be:
- Match by `.+/3`

**Solution 4 (reads named without /N)**

- If your data has differently formatted sequence identifiers, the "Match by" expression from *Solution 2* above can be modified to suit your identifiers.

Alternative identifiers such as:

```
@M00946:180:000000000-ANFB2:1:1107:14919:14410 1:N:0:1
@M00946:180:000000000-ANFB2:1:1107:14919:14410 2:N:0:1
```
