---
layout: tutorial_hands_on

title: "Introduction to sequencing with Python (part four)"
questions:
  - How to manipulate files in Python
  - How to read and write FASTA
  - How to read and write FASTQ
  - How to read and write SAM
objectives:
  - Understand manipulation of files in Python
time_estimation: "1h"
key_points:
  - Python can be used to read and write files
contributions:
  authorship:
  - nekrut

priority: 5

subtopic: gnmx
draft: true

notebook:
  language: python
  pyolite: true
---

# Reading and writing files in Python

First let's download a file we will be using to your notebooks:

```
!wget https://raw.githubusercontent.com/nekrut/BMMB554/master/2023/data/l9/mt_cds.fa
```

In Python, you can handle files using the built-in `open` function. The `open` function creates a file object, which you can use to read, write, or modify the file.

Here's an example of how to open a file for reading:


```python
f = open("mt_cds.fa", "r")
```

In this example, the open function takes two arguments: the name of the file, and the mode in which you want to open the file. The `r` mode indicates that you want to open the file for reading.

After you've opened the file, you can read its contents using the read method:

```python
contents = f.read()
print(contents)
```

You can also read the file line by line using the readline method:

```python
line = f.readline()
print(line)
```

When you're done reading the file, you should close it using the close method:

```python
f.close()
```

You can also use the `with` statement to automatically close the file when you're done:

```python
with open("mt_cds.fa", "r") as f:
    contents = f.read()
    print(contents)
```

You can also write to files using `write` method (note the `"w"` mode):

```python
f = open("sample.txt", "w")
f.write("This is a new line.")
f.close()
```

If you open an existing file in write mode, its contents will be overwritten. If you want to append to an existing file instead, you can use the `"a"` mode:

```python
f = open("sample.txt", "a")
f.write("This is another line.")
f.close()
```

When you're writing to a file, it's important to make sure you close the file when you're done. If you don't close the file, any changes you make may not be saved.

In addition to reading and writing text files, you can also use Python to handle binary files, such as images or audio files.

Let's download an image:

```
!wget https://imgs.xkcd.com/comics/file_extensions.png
```

Here's an example of how to read an image file:

```python
with open("file_extensions.png", "rb") as f:
    contents = f.read()
```

Note that when working with binary files, you must use the `"rb"` mode for reading and the `"wb"` mode for writing.

There are many more features and methods related to file handling in Python, but the basics covered here should be enough to get you started.

## [Fasta](https://en.wikipedia.org/wiki/FASTA_format)

Fasta is a file format that is commonly used to store biological sequences, such as DNA or protein sequences. In Python, you can read a Fasta file by opening the file, reading the lines one by one, and processing the data as needed.

Here's an example of how you might read a Fasta file in Python:

```python
sequences = {}
with open("mt_cds.fa", "r") as file:
  header = ""
  sequence = ""
  for line in file:
    line=line.rstrip()
    if line.startswith('>'):
      if header != "":
        sequences[header] = sequence
        sequence = ""
      header = line[1:]
    else:
      sequence += line
  if header != "":
    sequences[header] = sequence
```
The code above uses a `with` statement to open the file and read the lines one by one. If a line starts with a `">"`, it is assumed to be a header, and the current sequence is stored in the dictionary using the current header as the key. If the line does not start with a `">"`, it is assumed to be part of the current sequence.


## [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)

Let's download a sample fastq file:

```python
!wget https://raw.githubusercontent.com/nekrut/BMMB554/master/2023/data/l9/reads.fq
```

Fastq is a file format that is commonly used to store high-throughput sequencing data. It consists of a series of records, each of which includes a header, a sequence, a quality score header, and a quality score string. In Python, you can read a Fastq file by opening the file, reading the lines four at a time, and processing the data as needed.

Here's an example of how you might read a Fastq file in Python:

```python
def read_fastq(file_path):
    records = []
    with open(file_path, "r") as f:
        while True:
            header = f.readline().strip()
            if header == "":
                break
            sequence = f.readline().strip()
            quality_header = f.readline().strip()
            quality = f.readline().strip()
            records.append((header, sequence, quality))
    return records
```
In this example, the `read_fastq` function takes a file path as an argument, and returns a list of records, where each record is a tuple of four strings: the header, the sequence, the quality score header, and the quality score string. The function uses a `while` loop to read the lines four at a time until the end of the file is reached.

You can use this function to read a Fastq file like this:

```python
records = read_fastq("reads.fq")
for header, sequence, quality_header, quality in records:
    print(header)
    print(sequence)
    print(quality_header)
    print(quality)
```
This will print the headers, sequences, quality score headers, and quality scores in the Fastq file. You can modify the read_fastq function to process the data in any way you need.


## [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) 

Let's download an example SAM file:

```
!wget https://raw.githubusercontent.com/nekrut/BMMB554/master/2023/data/l9/sam_example.sam
```

SAM (Sequence Alignment/Map) is a file format that is used to store the results of DNA sequencing alignments. In Python, you can read a SAM file by opening the file, reading the lines one by one, and processing the data as needed.

Here's an example of how you might read a SAM file in Python:

```python
def read_sam(file_path):
    records = []
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            records.append(fields)
    return records
```
In this example, the `read_sam` function takes a file path as an argument, and returns a list of records, where each record is a list of fields. The function uses a `with` statement to open the file and read the lines one by one. If a line starts with an `"@"`, it is assumed to be a header and is ignored. If the line does not start with an `"@"`, it is assumed to be a record, and the fields are extracted by splitting the line on tabs.

You can use this function to read a SAM file like this:

```python
records = read_sam("sam_example.sam")
for fields in records:
    print(fields)
```

This will print the fields in the SAM file. You can modify the read_sam function to process the data in any way you need. For example, you might want to extract specific fields, such as the reference name, the start position, and the cigar string.