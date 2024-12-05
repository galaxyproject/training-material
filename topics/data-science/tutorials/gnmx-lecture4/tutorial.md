---
layout: tutorial_hands_on

title: "Introduction to sequencing with Python (part three)"
questions:
  - Understanding Python lists
  - Understanding Python dictionaries
  - Learnig about dynamic programming
objectives:
  - Understanding of lists and dictionaries
  - Learning about dynamic programming
  - Learning about how to translate DNA in Python
time_estimation: "1h"
key_points:
  - Lists and dictionaries are key data structures in Python
contributions:
  authorship:
  - nekrut

priority: 4

subtopic: gnmx
draft: true

notebook:
  language: python
  pyolite: true
---


[![XKCD 2483](https://imgs.xkcd.com/comics/seven.png)](https://xkcd.com/2483)


Preclass prep: Chapters [5](https://greenteapress.com/thinkpython2/html/thinkpython2009.html) and [7](https://greenteapress.com/thinkpython2/html/thinkpython2011.html) from "Think Python"

> This material uses examples from notebooks developed by [Ben Langmead](https://langmead-lab.org/teaching-materials/)
{: .quote cite="https://langmead-lab.org/teaching-materials" }

# Prep

1. Start [JupyterLab](https://mybinder.org/v2/gh/jupyterlab/jupyterlab-demo/try.jupyter.org?urlpath=lab)
2. Within JupyterLab start a new Python3 notebook
3. Open [this page](http://cs1110.cs.cornell.edu/tutor/#mode=edit) in a new browser tab

# Lists: Dynamic programming in sequence alignment

## Dynamic programming matrix as a 2D list

An excellent way to illustrate the utility of lists is to implement a dynamic programming algorithm for sequence alignment. Suppose we have two sequences that deliberately have different lengths:

$$
\texttt{G C T A T A C}$
$$

and 

$$
\texttt{G C G T A T G C}$
$$

Let's represent them as the following matrix where the first character $$\epsilon$$ represents an empty string:

$$
\begin{array}{ c | c | c | c | c | c | c}
& \epsilon & G & C & T & A & T & A & C\\
\hline
\epsilon \\
\hline
G\\
\hline
 C\\
 \hline
G\\
\hline
T\\
\hline
A\\
 \hline
T\\
\hline
G\\
\hline
C
\end{array}
\\
\textbf{Note}: sequence\ \texttt{X}\ is\ vertical\ and\ sequence\ \texttt{Y}\ is\ horizontal.
$$

In this matrix, the cells are addressed as shown below. They filled in using the following logic:

$$
D[i,j] = min\begin{cases} 
          \color{green}{D[i-1,j] + 1} & \\
          \color{blue}{D[i,j-1] + 1} & \\
          \color{red}{D[i-1,j-1] + \delta(x[i-1],y[j-1])}
             \end{cases}
$$

where $$\color{green}{green}$$ is *upper* cell, $$\color{blue}{blue}$$ is *left* cell, and $$\color{red}{red}$$ is *upper-left* cell:

$$
\begin{array}{ c | c | c | c | c | c | c}
& \epsilon & G & C & T & A & T & A & C\\
            \hline
          \epsilon &  &  &  &  &  &  &  & \\
          \hline
          G & \\
          \hline
          C &  & & \color{red}{D[2,2]} & \color{green}{D[2,3]}\\
          \hline
          G &  & & \color{blue}{D[3,2]} & D[3,3]\\
          \hline
          T & \\
          \hline
          A & \\
          \hline
          T & \\
          \hline
          G & \\
          \hline
          C & 
\end{array}
\\
\textbf{Note}: sequence\ \texttt{X}\ is\ vertical\ and\ sequence\ \texttt{Y}\ is\ horizontal.
$$

## Initializing the matrix

Let's initialize the first column and first row of the matrix. Because the distance between a string and an empty string is equal to the length of the
string (e.g., a distance between, say, string $$\texttt{TCG}$$ and an empty string is 3) this resulting matrix will look like this:

$$
\begin{array}{ c | c | c | c | c | c | c}
& \epsilon & G & C & T & A & T & A & C\\
            \hline
          \epsilon & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7\\
          \hline
          G & 1\\
          \hline
          C & 2\\
          \hline
          G & 3\\
          \hline
          T & 4\\
          \hline
          A & 5\\
          \hline
          T & 6\\
          \hline
          G & 7\\
          \hline
          C & 8
\end{array}
\\
\textbf{Note}: sequence\ \texttt{X}\ is\ vertical\ and\ sequence\ \texttt{Y}\ is\ horizontal.
$$

This can be achieved with the following code:

```python
D = np.zeros((len(x)+1, len(y)+1), dtype=int)
D[0, 1:] = range(1, len(y)+1)
D[1:, 0] = range(1, len(x)+1)
```

## Filling the matrix

Now we can fill the entire matrix by using two nested loops: one iterating over $$i$$ coordinate (sequence $$x$$) and the other iterating over $$j$$ coordinate (sequence $$y$$):

```python
for i in range(1, len(x)+1):
    for j in range(1, len(y)+1):
        delt = 1 if x[i-1] != y[j-1] else 0
        D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
```

Let's start with the cell between $$\texttt{G}$$ and $$\texttt{G}$$. Recall that:

$$
D[i,j] = min\begin{cases} 
          \color{green}{D[i-1,j] + 1} & \\
          \color{blue}{D[i,j-1] + 1} & \\
          \color{red}{D[i-1,j-1] + \delta(x[i-1],y[j-1])}
             \end{cases}
$$


where $$\delta(x[i-1],y[j-1]) = 0$$ if $$x[i-1] = y[j-1]$$ and $$1$$ otherwise. And now let's color each of the cells corresponding to each part of the above expression:

$$
\begin{array}{ c | c | c | c | c | c | c}
& \epsilon & G & C & T & A & T & A & C\\
\hline
\epsilon & \color{red}0 & \color{green}1 & 2 & 3 & 4 & 5 & 6 & 7\\
\hline
G & \color{blue}1\\
          \hline
          C & 2\\
          \hline
          G & 3\\
          \hline
          T & 4\\
          \hline
          A & 5\\
          \hline
          T & 6\\
          \hline
          G & 7\\
          \hline
          C & 8
\end{array}
\\
\textbf{Note}: sequence\ \texttt{X}\ is\ vertical\ and\ sequence\ \texttt{Y}\ is\ horizontal.
$$

In this specific case:

$$
D[i,j] = min\begin{cases} 
          \color{green}{D[i-1,j] + 1}\ or\ 0+0=0 & \\
          \color{blue}{D[i,j-1] + 1}\ or\ 1+1=2 & \\
          \color{red}{D[i-1,j-1] + \delta(x[i-1],y[j-1])}\ or\ 1+1=2
             \end{cases}
$$


The minimum of 0, 2, and 2 will be 0, so we are putting zero into that cell:

$$
\begin{array}{ c | c | c | c | c | c | c}
& \epsilon & G & C & T & A & T & A & C\\
\hline
\epsilon & \color{red}0 & \color{green}1 & 2 & 3 & 4 & 5 & 6 & 7\\
\hline
G & \color{blue}1 & \color{red}0\\
          \hline
          C & 2\\
          \hline
          G & 3\\
          \hline
          T & 4\\
          \hline
          A & 5\\
          \hline
          T & 6\\
          \hline
          G & 7\\
          \hline
          C & 8
\end{array}
\\
\textbf{Note}: sequence\ \texttt{X}\ is\ vertical\ and\ sequence\ \texttt{Y}\ is\ horizontal.
$$

Using this logic we can fill the entire matrix like this:

$$
\begin{array}{ c | c | c | c | c | c | c}
& \epsilon & G & C & T & A & T & A & C\\
\hline
\epsilon & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7\\
\hline
                 G & 1 & 0 & 1 & 2 & 3 & 4 & 5 & 6\\
          \hline
                 C & 2 & 1 & 0 & 1 & 2 & 3 & 4 & 5\\
          \hline
                 G & 3 & 2 & 1 & 1 & 2 & 3 & 4 & 5\\
          \hline
                 T & 4 & 3 & 2 & 1 & 2 & 2 & 3 & 4\\
          \hline
                 A & 5 & 4 & 3 & 2 & 1 & 2 & 2 & 3\\
          \hline
                 T & 6 & 5 & 4 & 3 & 2 & 1 & 2 & 3\\
          \hline
                 G & 7 & 6 & 5 & 4 & 3 & 2 & 2 & 3\\
          \hline
                 C & 8 & 7 & 6 & 5 & 4 & 3 & 3 & \color{red}2
\end{array}
\\
\textbf{Note}: sequence\ \texttt{X}\ is\ vertical\ and\ sequence\ \texttt{Y}\ is\ horizontal.
$$


The lower rightmost cell highlighted in red is special. It contains the value for the edit distance between the two strings.
The following Python script implements this idea. You can see that it is essentially instantaneous:

```python
import numpy as np

def edDistDp(x, y):
    """ Calculate edit distance between sequences x and y using
        matrix dynamic programming.  Return matrix and distance. """
    D = np.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D,D[len(x),len(y)]
```
A graphical representation of the matrix between `GCGTATGCACGC` and `GCTATGCCACGC` looks like this:

![A DP matrix with Seaborn](https://i.imgur.com/veMfPFt.png)

This image is generated using [Seaborn](https://seaborn.pydata.org/index.html) package using matrix directly:

```python
sns.heatmap(D,annot=True,cmap="crest")
```


------

# Dictionaries: Translating sequences

Perhaps the best way to demonstrate the utility of dictionaries is by using DNA-to-Protein translation as an example. 

## Using dictionaries to translate DNA

The following dictionary maps codons to corresponding amino acid translations. In this case, codon is the *key* and amino acid is the *value*:


```python
table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
```

Let's generate a random DNA sequence:


```python
import random
seq = "".join([random.choice('atcg') for x in range(100)])
```


```python
seq
```




    'agaccgtagcccaagtgcgtttgaatgtggctacttgggaggatttcattgcggtctgtctccgtacttgttattggtcttctttctgcattatgacgca'



To translate this sequence we write a code that uses a `for` loop that iterates over the DNA sequence in steps of 3, creating a codon at each iteration. If the codon is less than 3 letters long, the loop is broken. The resulting amino acid is then added to the `translation` string:


```python
translation = ""
for i in range(0, len(seq), 3):
    codon = seq[i:i+3].upper()
    if len(codon) < 3: break
    if codon in table:
        translation += table[codon]
    else:
        translation += "X"

print("Translation:", translation)
```

    Translation: RP_PKCV_MWLLGRISLRSVSVLVIGLLSAL_R


Note that the code uses the `upper()` method to ensure the codon is in uppercase since the table dictionary is case-sensitive. Additionally, the code checks if the codon is in the table dictionary and if not, it adds the letter "X" to the translation. This is a common way to represent unknown or stop codons in a protein sequence.


```python
translation
```




    'RP_PKCV_MWLLGRISLRSVSVLVIGLLSAL_R'



Now we define a function that would perform translation so that we can reuse it later:


```python
def translate(seq):
    translation = ''
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        if len(codon) < 3: break
        if codon in table:
            translation += table[codon]
        else:
            translation += "X"
    return(translation)
```


```python
translate(seq)
```




    'RP_PKCV_MWLLGRISLRSVSVLVIGLLSAL_R'

## Expanding to multiple phases (frames)

We can further modify the function by adding a `phase` parameter that would allow translating in any of the three phases:


```python
def translate_phase(seq,phase):
    translation = ''
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    assert phase >= 0 and phase <= 3, "Phase parameter can only be set to 0, 1, or 2! You specified {}".format(phase)
    for i in range(phase, len(seq), 3):
        codon = seq[i:i+3].upper()
        if len(codon) < 3: break
        if codon in table:
            translation += table[codon]
        else:
            translation += "X"
    return(translation)
```


```python
translate_phase(seq,2)
```




    'TVAQVRLNVATWEDFIAVCLRTCYWSSFCIMT'



To translate in all six reading frames (three of the "+" strand and three of the "-" strand) we need to be able to create a reverse complement of the sequence. Let's write a simple function for that.
The cell below implements a function `revcomp` that takes a DNA sequence as input and returns its reverse complement. It works by first reversing the sequence using the slice
notation `seq[::-1]`, which returns the sequence in reverse order. Then, the `translate` method is used with the `str.maketrans` function to replace each occurrence
of 'a', 't', 'g', 'c', 'A', 'T', 'G', and 'C' in the reversed sequence with 't', 'a', 'c', 'g', 'T', 'A', 'C', and 'G', respectively:


```python
import string

def revcomp(seq):
    return seq[::-1].translate(str.maketrans('atgcATCG','tagcTACG'))
```

Now let's use this function to create translation in all six reading frames. The cell below uses a `for` loop that iterates over the range `(0, 3)`, representing the different phases
(or starting positions) of the translation. At each iteration, the `translate_phase` function is called with the DNA sequence and the current phase, and the resulting protein sequence
is appended to the `translations` list along with the `phase` and the `strand` orientation (+ or -):


```python
translations = []

for i in range(0,3):
    translations.append((translate_phase(seq,i),str(i),'+'))
    translations.append((translate_phase(revcomp(seq),i),str(i),'-'))
```


```python
translations
```




    [('RP_PKCV_MWLLGRISLRSVSVLVIGLLSAL_R', '0', '+'),
     ('SLIIDKQQE_ELATDRRIKWWELRRLKASSRSL', '0', '-'),
     ('DRSPSAFECGYLGGFHCGLSPYLLLVFFLHYDA', '1', '+'),
     ('R__STNNRNKN_PQTGESNGGNYGD_KRVPVAC', '1', '-'),
     ('TVAQVRLNVATWEDFIAVCLRTCYWSSFCIMT', '2', '+'),
     ('ADNRQTTGIRTSHRQANQMVGTTEIESEFP_P', '2', '-')]



## Finding coordinates of continuous translations

The translation we've generated above contains stops (e.g., `_` symbols). The actual biologically relevant protein sequences are between stops.
We now need to split translation strings into meaningful peptide sequences and compute their coordinates. Let's begin by splitting a string on `_` and computing the start and end positions of each peptide:


```python
string = "aadsds_dsds_dsds"

split_indices = []
for i,char in enumerate(string):
    if char == "_":
        split_indices.append(i)
print(split_indices)
```

    [6, 11]


The code above generates a list of split indices for a string. The list contains the indices of the characters in the string that match a specified character (in this case, the underscore `_` character).

The `enumerate` function is used to loop over the characters in the string, and at each iteration, the current index and character are stored in the variables `i` and `char`, respectively.
If the current character matches the specified character, the index `i` is appended to the `split_indices` list.

After the loop, the `split_indices` list is printed to the console. For the input string `"aadsds_dsds_dsds"`, the output would be `[6, 11]`, indicating that the dashes are located at indices 6 and 11.

But we actually need coordinates of peptides bound by `_` characters. To get to this let's first modify `split_indices` by adding the beginning and end:


```python
string = "aadsds_dsds_dsds"

split_indices = []
for i,char in enumerate(string):
    if char == "_":
        split_indices.append(i)
        
split_indices.insert(0,-1)
split_indices.append(len(string))
print(split_indices)
```

    [-1, 6, 11, 16]


Now, let's convert these to ranges and also stick the peptide sequence in:


```python
string = "aadsds_dsds_dsds"

split_indices = []
for i, char in enumerate(string):
    if char == "_":
        split_indices.append(i)
        
split_indices.insert(0, -1)
split_indices.append(len(string))

orfs = string.split('_')

parts = []
for i in range(len(split_indices)-1):
    parts.append((orfs[i],split_indices[i]+1, split_indices[i+1]))

print(parts)
```

    [('aadsds', 0, 6), ('dsds', 7, 11), ('dsds', 12, 16)]


Now let's convert this to function:


```python
def extract_coords(translation):
    split_indices = []
    for i, char in enumerate(translation):
        if char == "_":
            split_indices.append(i)
        
    split_indices.insert(0, -1)
    split_indices.append(len(translation))

    parts = []
    for i in range(len(split_indices)-1):
        parts.append((translation.split('_')[i], split_indices[i] + 1, split_indices[i + 1]))

    return(parts)
```


```python
extract_coords(string)
```




    [('aadsds', 0, 6), ('dsds', 7, 11), ('dsds', 12, 16)]



And specify the right parameters to make it truly generic:


```python
def extract_coords_with_annotation(separator,translation,phase,strand):
    split_indices = []
    for i,char in enumerate(translation):
        if char == separator:
            split_indices.append(i)
        
    split_indices.insert(0,-1)
    split_indices.append(len(translation))

    parts = []
    for i in range(len(split_indices)-1):
        parts.append((translation.split(separator)[i], phase, strand, split_indices[i]+1, split_indices[i+1]))

    return(parts)
```


```python
extract_coords_with_annotation('_', string, '0', '+')
```




    [('aadsds', '0', '+', 0, 6),
     ('dsds', '0', '+', 7, 11),
     ('dsds', '0', '+', 12, 16)]



## Analyzing all translations of a given sequence

We begin by parsing the `translations` list we defined above:

```python
all_translations = []
for item in translations:
    all_translations.append(extract_coords_with_annotation('_',item[0],item[1],item[2]))
```


```python
all_translations
```




    [[('RP', '0', '+', 0, 2),
      ('PKCV', '0', '+', 3, 7),
      ('MWLLGRISLRSVSVLVIGLLSAL', '0', '+', 8, 31),
      ('R', '0', '+', 32, 33)],
     [('SLIIDKQQE', '0', '-', 0, 9),
      ('ELATDRRIKWWELRRLKASSRSL', '0', '-', 10, 33)],
     [('DRSPSAFECGYLGGFHCGLSPYLLLVFFLHYDA', '1', '+', 0, 33)],
     [('R', '1', '-', 0, 1),
      ('', '1', '-', 2, 2),
      ('STNNRNKN', '1', '-', 3, 11),
      ('PQTGESNGGNYGD', '1', '-', 12, 25),
      ('KRVPVAC', '1', '-', 26, 33)],
     [('TVAQVRLNVATWEDFIAVCLRTCYWSSFCIMT', '2', '+', 0, 32)],
     [('ADNRQTTGIRTSHRQANQMVGTTEIESEFP', '2', '-', 0, 30),
      ('P', '2', '-', 31, 32)]]



Now the problem with this list is that it is nested. However, we need to make it flat:


```python
flat_list = []
for sublist in all_translations:
    for item in sublist:
        flat_list.append(item)
```


```python
flat_list
```




    [('RP', '0', '+', 0, 2),
     ('PKCV', '0', '+', 3, 7),
     ('MWLLGRISLRSVSVLVIGLLSAL', '0', '+', 8, 31),
     ('R', '0', '+', 32, 33),
     ('SLIIDKQQE', '0', '-', 0, 9),
     ('ELATDRRIKWWELRRLKASSRSL', '0', '-', 10, 33),
     ('DRSPSAFECGYLGGFHCGLSPYLLLVFFLHYDA', '1', '+', 0, 33),
     ('R', '1', '-', 0, 1),
     ('', '1', '-', 2, 2),
     ('STNNRNKN', '1', '-', 3, 11),
     ('PQTGESNGGNYGD', '1', '-', 12, 25),
     ('KRVPVAC', '1', '-', 26, 33),
     ('TVAQVRLNVATWEDFIAVCLRTCYWSSFCIMT', '2', '+', 0, 32),
     ('ADNRQTTGIRTSHRQANQMVGTTEIESEFP', '2', '-', 0, 30),
     ('P', '2', '-', 31, 32)]



Now we can load the list into `Pandas` and plot away:


```python
import pandas as pd
df = pd.DataFrame(flat_list,columns=['aa','frame','phase','start','end'])
```


```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>aa</th>
      <th>frame</th>
      <th>phase</th>
      <th>start</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>RP</td>
      <td>0</td>
      <td>+</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PKCV</td>
      <td>0</td>
      <td>+</td>
      <td>3</td>
      <td>7</td>
    </tr>
    <tr>
      <th>2</th>
      <td>MWLLGRISLRSVSVLVIGLLSAL</td>
      <td>0</td>
      <td>+</td>
      <td>8</td>
      <td>31</td>
    </tr>
    <tr>
      <th>3</th>
      <td>R</td>
      <td>0</td>
      <td>+</td>
      <td>32</td>
      <td>33</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SLIIDKQQE</td>
      <td>0</td>
      <td>-</td>
      <td>0</td>
      <td>9</td>
    </tr>
    <tr>
      <th>5</th>
      <td>ELATDRRIKWWELRRLKASSRSL</td>
      <td>0</td>
      <td>-</td>
      <td>10</td>
      <td>33</td>
    </tr>
    <tr>
      <th>6</th>
      <td>DRSPSAFECGYLGGFHCGLSPYLLLVFFLHYDA</td>
      <td>1</td>
      <td>+</td>
      <td>0</td>
      <td>33</td>
    </tr>
    <tr>
      <th>7</th>
      <td>R</td>
      <td>1</td>
      <td>-</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td></td>
      <td>1</td>
      <td>-</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>9</th>
      <td>STNNRNKN</td>
      <td>1</td>
      <td>-</td>
      <td>3</td>
      <td>11</td>
    </tr>
    <tr>
      <th>10</th>
      <td>PQTGESNGGNYGD</td>
      <td>1</td>
      <td>-</td>
      <td>12</td>
      <td>25</td>
    </tr>
    <tr>
      <th>11</th>
      <td>KRVPVAC</td>
      <td>1</td>
      <td>-</td>
      <td>26</td>
      <td>33</td>
    </tr>
    <tr>
      <th>12</th>
      <td>TVAQVRLNVATWEDFIAVCLRTCYWSSFCIMT</td>
      <td>2</td>
      <td>+</td>
      <td>0</td>
      <td>32</td>
    </tr>
    <tr>
      <th>13</th>
      <td>ADNRQTTGIRTSHRQANQMVGTTEIESEFP</td>
      <td>2</td>
      <td>-</td>
      <td>0</td>
      <td>30</td>
    </tr>
    <tr>
      <th>14</th>
      <td>P</td>
      <td>2</td>
      <td>-</td>
      <td>31</td>
      <td>32</td>
    </tr>
  </tbody>
</table>
</div>

Now let's plot it:


```python
import altair as alt

plus = alt.Chart(df[df['phase']=='+']).mark_rect().encode(
    x = alt.X('start:Q'),
    x2 = alt.X2('end:Q'),
    y = alt.Y('frame:N'),
    color='frame',
    tooltip='aa:N'
             ).properties(
    width=600,
    height=100)

minus = alt.Chart(df[df['phase']=='-']).mark_rect().encode(
    x = alt.X('start:Q',sort=alt.EncodingSortField('start:Q', order='descending')),
    x2 = alt.X2('end:Q'),
    y = alt.Y('frame:N'),
    color='frame',
    tooltip='aa:N'
             ).properties(
    width=600,
    height=100)

plus & minus
```

![ORFS](./images/orfs.svg "Translation on plus (top) and minus (bottom) strands.")
