---
layout: tutorial_hands_on

title: Python - Subprocess
level: Intermediate
requirements: []
follow_up_training: []
questions:
- How can I run another program?

objectives:
- Run a command in a subprocess.
- Learn about `check_call` and `check_output` and when to use each of these.
- Read it's output.

time_estimation:  45M
key_points:
- "**DO NOT USE `os.system`**"
- "**DO NOT USE shell=True**"
- 👍Use `subprocess.check_call()` if you don't care about the output, just that it succeeds.
- 👍Use `subprocess.check_output()` if you want the output
- Use `.decode('utf-8')` to read the output of `check_output()`
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

Sometimes you need to run other tools in Python, like maybe you want to
Here we'll give a quick tutorial on how to read and write files within Python.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Subprocesses

Programs can run other programs, and in Python we do this via the `subprocess` module. It lets you run any other command on the system, just like you could at the terminal.

The first step is importing the module

```python
import subprocess
```

You'll primarily use two functions:

```python
help(subprocess.check_call)
```

Which executes a command and checks if it was successful (or it raises an exception), and

```python
help(subprocess.check_output)
```


# Check Call: Downloading Files

Which executes a command returns the output of that command. This is really useful if you're running a subprocess that writes something to stdout, like a report you need to parse. We'll learn how to use these by running two gene callers, augustus and glimmer. You can install both from Conda if you do not have them already.

```
conda create -n subprocess augustus glimmer3
```

Additionally you'll need two files, you *generally should not do this*, but you can use a subprocess to download the file! We'll use `subprocess.check_call` for this which simply executes the program, and continues on. If there is an error in the execution, it will raise an exception and stop execution.

```python
url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/836/945/GCF_000836945.1_ViralProj14044/GCF_000836945.1_ViralProj14044_genomic.fna.gz"
genome = 'Escherichia virus T4.fna.gz'
subprocess.check_call(['wget', url, '-O', genome])
subprocess.check_call(['gzip', '-d', genome])
```

> ### {% icon tip %} Tip: What do these commands look like on the CLI?
> ```
> wget https://ftp.ncbi.nlm.nih.gov/.... -O "Escherichia virus T4.fna.gz"
> gzip -d "Escherichia virus T4.fna.gz"
> ```
{: .tip}

The above segment
- sets a url variable
- sets an output filename, `Escherichia virus T4.fna.gz`
- runs `check_call` with a single argument: a list
    - `wget` a tool we use to download files
    - the URL
    - `-O` indicating the next argument will be the 'output name'
    - what we want the output filename to be called
- runs `check_call` with a single argument: a list
    - `gzip` a tool to decompress files
    -`-d` indicating we want to decompress
    - and the filename.

This list is especially important. When you run commands on the command line, normally you just type in a really bit of text by yourself. It's one big string, and you're responsible for making sure quotation marks appear in the right place. For instance, if you have spaces in your filenames, you have to quote the filename. Python requires you specify a list of arguments, and then handles the quoting for you! Which, honestly, is easier and safer.

> > ### {% icon code-in %} Terminal
> > Here we manually quote the argument
> > ```
> > glimmer3 "bow genome.txt"
> > ```
> {: .code-in}
> > ### {% icon code-out %} Python
> > Here python handles that for us
> > ```
> > subprocess.check_call(['glimmer3', 'bow genome.txt'])
> > ```
> {: .code-out}
{: .code-2col}

> ### {% icon tip %} Tip: Exploitation!
> This is one of the major reasons we don't use `os.system` or older Python interfaces for running commands.
> If you're processing files, and a user supplies a file with a space, if your program isn't expecting that space in that filename, then it could do something dangerous!
> Like exploit your system!
>
> So, **always** use `subprocess` if you run to commands, never any other module, *despite what you see on the internet!*
{: .tip}


There are more functions in the module, but the vast majority of the time, those are sufficient.

# Check Output: Gene Calling with Augustus

```python
gff3 = subprocess.check_output([
    'augustus', # Our command
    '--species=E_coli_K12', # the first argument, the species, we're using a phage so we call genes  based on it's host organism
    'Escherichia virus T4.fna', # The path to our genome file, without the .gz because we decompressed it.
    '--gff3=on' # We would like gff3 formatted output (it's easier to parse!)
])

gff3 = gff3.decode('utf-8')
gff3 = gff3.split('\n')
```

> ### {% icon tip %} Tip: What does this commands look like on the CLI?
> ```
> augustus --species=E_coli_K12 'Escherichia virus T4.fna' --gff3=on
> ```
{: .tip}

If you're using `subprocess.check_output()` python doesn't return plain text `str` to you, instead it returns a `bytes` object. We can decode that into text with `.decode('utf-8')`, a phrase you should memorise as going next to `check_output()`, for 99% of use cases.

Let's look at the results!

```python
print(gff3[0:20])
```

It's a lot of comment lines, starting with `#`. Let's remove those

```python
cleaned_gff3 = []
for line in gff3:
    if line.startswith('#'):
        continue
    cleaned_gff3.append(line)

print(cleaned_gff3[0:20])
```

And now you've got a set of gff3 formatted gene calls! You can use all of your loop processing skills to slice and dice this data into something great!

# Aside: `stdin`, `stderr`, `stdout`

All unix processes have three default file handles that are available to them:

- `stdin`, where data is passed to the program via a pipe. E.g. `generate-data | my-program`, there the program would read the output of `generate-data` from the pipe.
- `stdout`, the default place where things are written. E.g. if you `print()` in python, it goes to `stdout`. People often redirect `stdout` to a file, like `my-program > output.txt` to save the output.
- `stderr`, *generally* if your program produces output on `stdout`, you might still want to log messages (errors, % done, etc.) If you write to `stdout`, it might get mixed in with the user's outputs, so we write to `stderr`, which also gets printed to the screen, and looks identical as any print statement, but it's coming from a separate pipe.

# Pipes

One of the more complicated cases, however, is when you need pipes.

```python
url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/721/125/GCF_001721125.1_ASM172112v1/GCF_001721125.1_ASM172112v1_cds_from_genomic.fna.gz"
cds = 'E. Coli CDSs.fna.gz'
subprocess.check_call(['wget', url, '-O', cds])
subprocess.check_call(['gzip', '-d', cds])
```

With subprocesses, you can control the stdin, and stdout of the process by using file handles.

> > ### {% icon code-in %} Terminal
> > Here we pipe a file to a process named `build-icm` which takes one argument, the output name. It reads sequences from stdin.
> > ```
> > cat seq.fa | build-icm test.icm
> > # OR
> > build-icm test.icm < seq.fa
> > ```
> {: .code-in}
> > ### {% icon code-out %} Python
> > Here we need to do a bit more.
> > 1. Open a file handle
> > 2. Pass that file handle to `check_call` or `check_output`. This determines where stdin comes from.
> > ```
> > with open('seq.fa', 'r') as handle:
> >     subprocess.check_call(['build-icm', 'test.icm'], stdin=handle)
> > ```
> {: .code-out}
{: .code-2col}

We'll do that now:

```python
with open('E. Coli CDSs.fna', 'r') as handle:
    subprocess.check_call(['build-icm', 'test.icm'], stdin=handle)
```

> ### {% icon tip %} Tip: What does this commands look like on the CLI?
> ```
> build-icm test.icm < 'E. Coli CDSs.fna'
> ```
{: .tip}

Here we build a model, based on the sequences of *E. Coli* K-12, that Glimmer3 can use.

```python
output = subprocess.check_output([
    'glimmer3', # Our program
    'Escherichia virus T4.fna', # The input genome
    'test.icm', # The model we just built
    't4-genes'  # The base name for output files. It'll produce t4-genes.detail and t4-genes.predict.
]).decode('utf-8') # And of course we decode as utf-8

print(output)
```

> ### {% icon tip %} Tip: What does this commands look like on the CLI?
> ```
> glimmer3 'Escherichia virus T4.fna' test.icm t4-genes
> ```
{: .tip}

*What happened here?* The output of the program was written to `stderr`, not `stdout`, so Python may print that out to your screen, but `output` will be empty. To solve this common problem we can re-run the program and collect both `stdout` and `stderr`.

```python
output = subprocess.check_output([
    'glimmer3', # Our program
    'Escherichia virus T4.fna', # The input genome
    'test.icm', # The model we just built
    't4-genes'  # The base name for output files. It'll produce t4-genes.detail and t4-genes.predict.
], stderr=subprocess.STDOUT).decode('utf-8') # And of course we decode as utf-8

print(output)
```

Here we've re-directed the `stderr` to `stdout` and mixed both of them together. This isn't always what we want, but here the program produces no output, and we can do that safely, and now we can parse it or do any other computations we need with it! Our Glimmer3 gene calls are in `t4-genes.detail` and `t4-genes.predict` if we want to open and process those as well.
