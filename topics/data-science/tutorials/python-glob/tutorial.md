---
layout: tutorial_hands_on

title: Python - Globbing
level: Intermediate
requirements: []
follow_up_training: []
questions:
- How can I collect a list of files.

objectives:
- Use glob to collect a list of files
- Learn about the potential pitfalls of glob

time_estimation:  15M
key_points:
- If your data is ordering dependent, sort your globs!
enable: false
subtopic: python-modular
contributors:
  - hexylena
  - dirowa
  - bazante1

priority: 2
notebook:
  language: python
---

Globbing is the term used in computer science when we have a bunch of files and we want to list all of them matching some pattern.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Setup

We'll start by creating some files for use in the rest of this tutorial

```python
import os
import subprocess

dirs = ['a', 'a/b', 'c', 'c/e', 'd', '.']
files = ['a.txt', 'a.csv', 'b.csv', 'b.txt', 'e.glm']

for d in dirs:
    # Create some directories
    os.makedirs(d, exist_ok=True)
    # Create some files
    for f in files:
        subprocess.check_output(['touch', os.path.join(d, f)])
```

Now we should have a pretty full folder!

# Finding Files

We can use the glob module to find files:

```python
import glob
print(glob.glob('*.csv'))
print(glob.glob('*.txt'))
```

Here we use an asterisk (`*`) as a wildcard, it matches any bit of text (but not into folders!) to all matching files. Here we list all matching `csv` or `txt` files. This is great to find files matching a pattern.

We can also use asterisks anywhere in the glob, it doesn't just have to be the filename portion:

```python
print(glob.glob('a*'))
```

Here we even see a third entry: the directory.

# Finding files in directories

Until now we've found only files in a single top level directory, but what if we wanted to find files in subdirectories?

Only need a single directory? Just include that!

```python
print(glob.glob('a/*.csv'))
```

But if you need more levels, or want to look in *all* folders, then you need the double wildcard! With two asterisks `**` we can search recursively through directories for files:

```python
print(glob.glob('**/a.csv'))
```

# Exercise

> ### {% icon question %} Question: Where in the world is the CSV?
>
> 1. How would you find all `.csv` files?
> 2. How would you find all `.txt` files?
> 3. How would you find all files starting with the letter 'e'?
>
> > ### {% icon solution %} Solution
> >
> > 1. `glob.glob('**/*.csv')`
> > 2. `glob.glob('**/*.txt')`
> > 3. `glob.glob('**/e*')`
> {: .solution}
{: .question}

```python
# Try things out here!
```

# Pitfalls

Some analyses (especially simultaions) can be dependent on data input order or data sorting. This was recently seen in {% cite Bhandari_Neupane_2019 %} where the data files used were sorted one way on Windows, and another on Linux, resulting in different results for the same code and the same datasets! Yikes!

If you know your analyses are dependent on file ordering, then you can use `sorted()` to make sure the data is provided in a uniform way every time.

```python
print(sorted(glob.glob('**/a.csv')))
```

If you're not sure if your results will be dependent, you can try sorting anyway. Or better yet, randomising the list of inputs to make sure your code behaves properly in any scenario.
