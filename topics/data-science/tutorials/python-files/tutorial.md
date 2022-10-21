---
layout: tutorial_hands_on

title: Python - Files & CSV
level: Introductory
requirements: []
follow_up_training: []
questions:
- How can I read from a file?
- How can I parse a CSV file?
- How can I write results out

objectives:
- Read data from a file
- Write new data to a file
- Use `with` to ensure the file is closed properly.
- Use the CSV module to parse comma and tab separated datasets.

time_estimation:  1H30M
key_points:
- "File reading requires a mode: read, write, and append"
- Use the CSV module to parse CSV files.
- do NOT attempt to do it yourself, it will encounter edge cases that the CSV module handles for you
- Use a `with` block to open a file.

subtopic: python-modular
contributors:
  - hexylena
  - dirowa
  - bazante1

priority: 7
notebook:
  language: python
---

Here we'll give a quick tutorial on how to read and write files within Python.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Setup

For this tutorial, we assume you're working in a notebook (Jupyter, CoCalc, etc) so we'll run a quick "setup" step to download some CSV data:

```python
import urllib.request
# Download a copy of Hamlet.
urllib.request.urlretrieve("https://gutenberg.org/cache/epub/1524/pg1524.txt", "hamlet.txt")
# Download some COVID data for Europe
urllib.request.urlretrieve("https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/csv/data.csv", "vaccinations.csv")
# And a fastq file
urllib.request.urlretrieve("https://gist.github.com/hexylena/7d249607f8f763301f06c78a48c3bf6f/raw/a100e278cee1c94035a3a644b16863deee0ba2c0/example.fastq", "example.fastq")
```

And now we're ready to get started learning about files!

## Reading from a file

Reading from and writing to files in Python is very straightforward, we use `open()` to open a file to read from it. Files are accessed through something called a **file [handle](https://en.wikipedia.org/wiki/Handle_(computing))**, you're not accessing the file itself, you're opening a connection to the file, and then you can read lines from through that file handle. When you open a file handle, you must specify it's mode:

Mode | Purpose
--- | ---
`r` | We will **read** from this file handle
`w` | We will **write** to this file handle
`a` | We will **append** to this file handle (we cannot access earlier contents!)

You will mostly use `r` and `w`, `a` is especially useful for writing to program logs where you don't really care what was written before, you just want to add your new logs to the end of the file.

```python
with open('hamlet.txt', 'r') as handle:
    # readlines reads every line of a file into a giant list!
    lines = handle.readlines()
```

Here we introduce a new bit of syntax, the `with` block. Technically `with` begins a "context manager" which allows python to setup some things before the block, run some contents in the block, and automatically handle cleanup of this block. When you open a file, you *must* close it when you're done with it (otherwise bad things can happen!) and `with` prevents most of those issues.

In the above code snippet after the second line, the file (referred to by `handle`) is automatically closed.

> > <code-in-title>Using `with`</code-in-title>
> > ```
> > with open('file.txt', 'r') as handle:
> >     print(handle.readlines())
> > ```
> {: .code-in}
> > <code-out-title>Not using `with`</code-out-title>
> > ```
> > handle = open('file.txt', 'r')
> > print(handle.readlines())
> > handle.close() # Important!
> > ```
> {: .code-out}
{: .code-2col}

Let's see what's in our file. `

```python
print(lines[0])
print(lines[1])
print(lines[2])
```

Notice how it prints out a blank line afterwards! This is due to a `\n`, a newline. A newline just tells the computer "please put content on the next line". We can see it by using the `repr()` function:

```python
print(repr(lines[0]))
```

Every line that's read in ends in a newline currently. This is done because if we wanted to write it back out, we would need to preserve those newlines, or all of the content would be on one giant line. Let's try writing out a file, it's *just* like reading in a file!

```python
with open('hamlet-copy.txt', 'w') as handle:
    # readlines reads every line of a file into a giant list!
    for line in lines:
        handle.write(line)
```

Check this file in your folder, does it look right? Is it identical in size to the original?

### Use Case: Summarisation

Let's use the file's contents for something useful. This file specifically is the play Hamlet, by Shakespeare. The contents are formatted with a speaker indicated with all capital letters, followed by their lines (potentially spread over multiple lines of the file.)

```
HAMLET. Madam, how like you this play?

QUEEN. The lady protests too much, methinks.

HAMLET. O, but she’ll keep her word.
```

So let's count up how many times each character speaks! (Roughly)

```python
with open('hamlet.txt', 'r') as handle:
    # readlines reads every line of a file into a giant list!
    lines = handle.readlines()

speakers = {}
```

Here we've initialised the `lines` variable with the contents of the text, and setup `speakers` as a dictionary that will let us track how many times each character speaks. Next let's define a function to check if a character is speaking on that line. It is if it meets two conditions: the first word is all caps, and it ends with a `.`.

```python
def is_speaker(word):
    return word == word.upper() and word[-1] == '.'
```

We can use that function later to check if a line starts with a speaker

```python
# Loop over every line we read in
for line in lines:
    # Split by default splits on whitespace.
    words = line.split()

    # Are there words on this line? We can't access the first word if we
    # haven't any words.
    if len(words) == 0:
        continue

    # Check if the first word is uppercase, and the last character is a `.`,
    # then it's a character speaking.
    if is_speaker(words[0]):
        # Give this an easier to remember and understand name.
        speaker = words[0]

        # Have we seen this speaker before? If not, we should add them to the
        # speakers dictionary. Hint: Try removing this to see why we do this.
        if speaker not in speakers:
            speakers[speaker] = 0

        # Increment the number of times we've seen them speak.
        speakers[speaker] = speakers[speaker] + 1
```

Ok! We've done a couple things here that all fall into the category of **defensive programming**. As programmers, we often accept input from users or from unknown sources. That input may be wrong, it may have bad data, it may be trying to attack us. So we respond by checking very carefully if things match our expectations, and rejecting the input otherwise. We did a couple things here for that:

- Checking if the input matches our expectations exactly (capitals, `.`)
- Checking if the line is empty, using `continue` to skip it if it was
- Checking if the speaker is already known in the dictionary, adding it otherwise.

Let's see who was the most chatty:

```python
for key, value in speakers.items():
    print(key, value)
```

We've clearly caught a number of values that aren't expected, some section headers (the numeric values, and some rare values we don't expect.)

```python
for key, value in speakers.items():
    if value > 1:
        print(key, value)
```

Hamlet, the titular character, has the vast majority of turns speaking throughout the play:

Character | Turns speaking
---       | ---
HAMLET    | 358
HORATIO   | 107
KING      | 102
POLONIUS  | 86
QUEEN     | 69


## Writing files

Writing a file out is exactly like reading a file, we just use the different file mode `w` to indicate we wish to write to a file:

```python
with open('hello.txt', 'w') as handle:
    # Let's write it out a few times
    handle.write("Hello, world!")
    handle.write("Hello, world!")
    handle.write("Hello, world!")
    handle.write("Hello, world!")
```

Check the file's contents in your folder, does it look like you expect? Remember your `\n`s!

### Use Case: Transformation

A common use case is transforming one file's contents into another file format or file type. Let's do a very simple example of that, taking a FASTQ file and transforming it into a FASTA file. Remember first that a FASTQ file looks like:

```
@M00970:337:000000000-BR5KF:1:1102:17745:1557 1:N:0:CGCAGAAC+ACAGAGTT
GTGCCAGCCGCCGCGGTAGTCCGACGTGGCTGTCTCTTATACACATCTCCGAGCCCACGAGACCGAAGAACATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAGAAGCAAATGACGATTCAAGAAAGAAAAAAACACAGAATACTAACAATAAGTCATAAACATCATCAACATAAAAAAGGAAATACACTTACAACACATATCAATATCTAAAATAAATGATCAGCACACAACATGACGATTACCACACATGTGTACTACAAGTCAACTA
+
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGFGGGGGGAFFGGFGGGGGGGGFGGGGGGGGGGGGGGFGGG+38+35*311*6,,31=******441+++0+0++0+*1*2++2++0*+*2*02*/***1*+++0+0++38++00++++++++++0+0+2++*+*+*+*+*****+0**+0**+***+)*.***1**//*)***)/)*)))*)))*),)0(((-((((-.(4(,,))).,(())))))).)))))))-))-(
```

Line | Contents
-- | --
1 | Identifier, prefixed with `@`
2 | Sequence
3 | `+`
4 | Quality scores

And a fasta file looks like this, where `>` indicates a sequence identifier, and it is followed by one or more lines of ACTGs.

```
>M00970:337:000000000-BR5KF:1:1102:17745:1557 1:N:0:CGCAGAAC+ACAGAGTT
GTGCCAGCCGCCGCGGTAGTCCGACGTGGCTGTCTCTTATACACATCTCCGAGCCCACGAGACCGAAGAACATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAGAAGCAAATGACGATTCAAGAAAGAAAAAAACACAGAATACTAACAATAAGTCATAAACATCATCAACATAAAAAAGGAAATACACTTACAACACATATCAATATCTAAAATAAATGATCAGCACACAACATGACGATTACCACACATGTGTACTACAAGTCAACTA
```

In the setup portion we downloaded a FASTQ file, now let's extract all of the sequences from this file, and write them out as a FASTA file. Why would you want to do this? Sometimes after sequencing a sample (especially metagenomics), you want to blast the sequences to figure out which organisms they belong to. A common way to do that is BLAST which accepts fasta formatted sequences. So we'll write something to convert these formats, removing the `+` and quality score lines.

```python
with open('example.fastq', 'r') as handle:
    data = handle.readlines()

i = 0
with open('example-converted.fasta', 'w') as handle:
    for line in data:
        # Since fastq files have groups of 4 lines, we can use that to extract
        # specific lines of every read:
        if i % 4 == 0:
            handle.write(">" + line[1:])
        if i % 4 == 1:
            handle.write(line)
        i = i + 1
```

That's it! Check if your `fasta` file looks correct.

> <tip-title>`enumerate`</tip-title>
> If you want to know which item number you're on while you're looping over a list, you can use the function `enumerate()`
> Try out this code to see how it works:
>
> ```
> for index, item in enumerate(['a', 'b', 'c']):
>     print(index, item)
> ```
>
> Try using it to clean up the above code.
{: .tip}

## Reading CSV data

If you're reading data from a comma separated value (CSV) or tab separated value (TSV) file, you should use the built in `csv` module to do this. You might ask yourself "Why, csv parsing is easy" and that is a common thought! It would be so simple to do something like

```python
# Please don't do this :)
with open('vaccinations.csv', 'r') as handle:
    first_lines = handle.readlines()[0:10]
    for line in first_lines:
        print(line.split(','))
```

But you would be wrong! This code has a subtle bug that you might not see until someone generates data that specifically affects it, with "quoted" columns. If you have a table like

Patient | Location | Disease Indications
--- | --- | ---
Helena | Den Haag, the Netherlands | Z87.890
Bob | Little Rock, Arkansas, USA | Z72.53
Jane | London, UK | Z86.16

This would probably be exported as a CSV file from Excel that looks like:

```
Patient,Location,Disease Indications
Helena,"Den Haag, the Netherlands",Z87.890
Bob,"Little Rock, Arkansas, USA",Z72.53
Jane,"London, UK",Z86.16
```

Note that some columns are quoted. What do you think will happen with the following code?

```python
csv_data = """
Patient,Location,Disease Indications
Helena,"Den Haag, the Netherlands",Z87.890
Bob,"Little Rock, Arkansas, USA",Z72.53
Jane,"London, UK",Z86.16
""".strip().split('\n')

# Please don't do this :)
for line in csv_data:
    print(line.split(','))
```

Does that look right? Maybe not. Instead we can use the `csv` module to work around this and properly process CSV files:

```python
import csv

csv_data = """
Patient,Location,Disease Indications
Helena,"Den Haag, the Netherlands",Z87.890
Bob,"Little Rock, Arkansas, USA",Z72.53
Jane,"London, UK",Z86.16
""".strip().split('\n')

# Please DO this :)
csv_reader = csv.reader(csv_data, delimiter=",", quotechar='"')
for row in csv_reader:
    print(row)
```

That looks a lot better! Now we've properly handled the quoted columns that contain one or more `,` in the middle of our file. This is actually one of the motivating factors in using the TSV format, the <kbd>tab</kbd> character is much more rare in data than <kbd>,</kbd>. There is less chance for confusion with poorly written software.

Let's read in some statistics about vaccinations:

```python
import csv

vax = []
with open('vaccinations.csv', 'r') as handle:
    csv_reader = csv.reader(handle, delimiter=",", quotechar='"')
    for row in csv_reader:
        # Skip our header row
        if row[0] == 'YearWeekISO':
            continue
        # Otherwise load in the data
        vax.append(row)

print(vax[0:10])
```

Here we have a 2 dimensional array, a list of lists. Each row is an entry in the main list, and each column is an entry in each of those children.

Our columns are:

Column | Value
-- | --
0   | YearWeekISO
1   | ReportingCountry
2   | Denominator
3   | NumberDosesReceived
4   | NumberDosesExported
5   | FirstDose
6   | FirstDoseRefused
7   | SecondDose
8   | DoseAdditional1
9   | UnknownDose
10  | Region
11  | TargetGroup
12  | Vaccine
13  | Population

Let's subset the data to make it a bit easier to work with, maybe we'll just use the Dutch data (please feel free to choose another column though!)

```python
country = 'NL'

subset = []
for row in vax:
    # Here we select for a country
    if row[1] == country:
        subset.append(row)
print(subset[0:10])
```

That should be easier to work with, now we only have one country's data. Let's do some exercises with this data:


> <question-title>Which vaccines were given?</question-title>
>
> Which vaccines were given? Use the `subset` to examine which vaccines were given in the Netherlands. *Tip*: if `x` is a list, `set(x)` will return the unique values in that list.
>
> > <solution-title></solution-title>
> > To figure out which vaccines were given, we can look at column 12:
> >
> > ```python
> > vaccines = []
> > for row in subset:
> >     vaccines.append(row[12])
> > print(set(vaccines))
> > ```
> >
> > We can use the `set` function to convert the list to a set, and show only the unique values.
> {: .solution}
{: .question}

```python
# Try things out here!
```

> <question-title>How many of each vaccine were given?</question-title>
>
> How many of each were given?
>
> *Tip*: use the accumulator pattern.
> *Tip*: Columns 5, 7, and 8 have doses being given out to patients.
>
> > <solution-title></solution-title>
> > ```python
> > doses = {}
> > for row in subset:
> >     brand = row[12]
> >     if brand not in doses:
> >         doses[brand] = 0
> >     doses[brand] = doses[brand] + int(row[5]) + int(row[7]) + int(row[8])
> > print(doses)
> > ```
> {: .solution}
{: .question}

```python
# Try things out here!
```

> <question-title>How many of each vaccine were exported? received?</question-title>
>
> How many of each were exported? received?
>
> *Tip*: you only need to loop once.
> *Tip*: you will need to handle an edge case here. Try and find it out!
>
> > <solution-title></solution-title>
> > ```python
> > export = {}
> > received = {}
> > for row in subset:
> >     brand = row[12]
> >     if brand not in export:
> >         export[brand] = 0
> >     if brand not in received:
> >         received[brand] = 0
> >     if row[4]:
> >         export[brand] = export[brand] + int(row[4])
> >     if row[3]:
> >         received[brand] = received[brand] + int(row[3])
> > print(export)
> > print(received)
> > ```
> {: .solution}
{: .question}

```python
# Try things out here!
```

> <question-title>When was the first dose received?</question-title>
>
> *Tip*: use `break`, and check how many doses were given!
>
> > <solution-title></solution-title>
> > ```python
> > for row in subset:
> >     if row[5] and int(row[5]) > 0:
> >         print(f"On {row[0]}, {row[5]} doses of {row[12]} were given")
> >         break
> > ```
> {: .solution}
{: .question}

```python
# Try things out here!
```

> <question-title>Transform the data for plotting</question-title>
>
> Let's say you want to plot the fraction of the population that has been vaccinated by the various points in time.
>
> - Subset the data further for TargetGroup (column 11) set to 'ALL'
> - Create an accumulator, to count how many doses have been given their FirstDose at each week
> - Use column 13 (population) to calculate the fraction of the population that has been given one of those doses at each week
>
> The output should be a list of percentages ranging from [0.0 to 1.0].
>
> > <solution-title></solution-title>
> > ```python
> > total_doses = 0
> > percent_vaccinated_per_week = []
> > for row in subset:
> >     if row[11] != 'ALL':
> >         continue
> >     total_doses = total_doses + int(row[5])
> >     percent_vaccinated_per_week.append(total_doses / int(row[13]))
> > ```
> {: .solution}
{: .question}

```python
# Try things out here!
percent_vaccinated_per_week = []
# Write code here!




# When you're done, you should have a 'results' variable
# You may need to `pip install matplotlib`
%matplotlib inline
import matplotlib.pyplot as plt
plt.plot(percent_vaccinated_per_week)
plt.xlabel('Week')
plt.ylabel('Percent Vaccinated')
```


> <question-title>Which vaccines were given?</question-title>
>
> Write this out to a new file, with two columns. The index and the Percent Vaccinated. It should be a comma separated file, and should have a header. Save this as a csv file named `weekly-percent-vax.csv`
>
> *Tip*: use a `csvwriter`, it works exactly like a `csvreader`. You can use the `writerow()` function to write out a row.
> *Tip*: Use `enumerate()` to get a list of items with indexes.
>
> > <solution-title></solution-title>
> >
> > ```
> > with open('weekly-percent-vax.csv', 'w') as handle:
> >     writer = csv.writer(csv_data, delimiter=",", quotechar='"')
> >     for row in enumerate(percent_vaccinated_per_week):
> >         writer.writerow(row)
> > ```
> >
> {: .solution}
{: .question}

```python
# Try things out here!
```

Congratulations on getting this far! Hopefully you feel more comfortable working with files.
