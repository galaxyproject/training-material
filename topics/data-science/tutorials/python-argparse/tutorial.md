---
layout: tutorial_hands_on

title: Python - Argparse
level: Intermediate
requirements: []
follow_up_training: []
questions:
- How do I make a proper command line script
- How do I use argparse?
- What problems does it solve?

objectives:
- Learn how sys.argv works
- Write a simple command line program that sums some numbers
- Use argparse to make it nicer.

time_estimation: 30M
key_points:
- "If you are writing a command line script, no matter how small, use argparse."
- "`--help` is even written for us, without us writing any special code to handle that case"
- "It handles a lot of cases and input validation for you"
- "It produces a nice `--help` text that can help you if you've forgotten what your tool does"
- "It's nice for users of your scripts! They don't have to read the code to know how it behaves if you document it well."

subtopic: python-modular
contributions:
  authorship:
  - hexylena
  editing:
  - bazante1
  testing:
  - dirowa
  funding:
  - avans-atgm

priority: 10
---

[`argparse`](https://docs.python.org/3/library/argparse.html) is an argument parsing library for Python that's part of the stdlib. It lets you make command line tools significantly nicer to work with.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

Unlike previous modules, this lesson won't use a Jupyter/CoCalc notebook, and that's because we'll be parsing command lines! You'll need to open a code editor on your platform of choice (`nano`, `vim`, `emacs`, VSCode are all options) and use the following blocks of code to construct your command line tool.

## `sys.argv`

In the coding world, whenever you run a Python script on the command line, it has a special variable available to it named `argv`. This is a list of all of the arguments used when you run a command line program.

> <hands-on-title>Print out argv</hands-on-title>
>
> 1. Create / open the file `run.py` in your text editor of choice
> 2. There we'll create a simple Python script that:
>   1. imports `sys`, the system module needed to access argv.
>   2. Prints out `sys.argv`
>
> > <solution-title></solution-title>
> > ```python
> > import sys
> >
> > print(sys.argv)
> > ```
> {: .solution}
>
> 3. Run this with different command line arguments:
>
>    ```
>    python run.py
>    python run.py 1 2 3 4
>    python run.py --help
>    ```
>
{: .hands_on}

> <question-title></question-title>
> What did you notice about the output? There are two main points.
> > <solution-title></solution-title>
> > 1. The name of the script (`run.py`) is included as the first value every time.
> > 2. All of the arguments are passed as strings, no numbers.
> {: .solution}
{: .question}

## Simple tasks

Let's sum up all of the numbers passed on the command line. We'll do this by hand, and then we'll replace it with `argparse` to see how much effort that saves us.

> <hands-on-title>Hands-on</hands-on-title>
> Update your script to sum up every number passed to it on the command line.
>
> It should handle:
> - 1 or more numbers
> - nothing (and maybe print out a message?)
> - invalid values (print out an error message that the value couldn't be processed.)
>
> Hints:
> - Skip the program name
> - Use `try` and `except` to try converting the string to a number.
>
> > <question-title></question-title>
> >
> > How does your updated script look?
> >
> > > <solution-title></solution-title>
> > >
> > > ```python
> > > import sys
> > >
> > > result = 0
> > >
> > > if len(sys.argv) == 1:
> > >     print("no arguments were supplied")
> > > else:
> > >     for arg in sys.argv[1:]:
> > >         try:
> > >             result += float(arg)
> > >         except:
> > >             print(f"Could not parse {arg}")
> > >
> > >     print(result)
> > > ```
> > {: .solution}
> {: .question}
{: .hands_on}

## Argparse

Argparse saves us a lot of work, because it can handle a number of things for us!

- Ensures that the correct number of arguments are provided (and provide a nice error message otherwise)
- Ensure that the correct types of arguments are provided (no strings for a number field)
- Provide a help message describing your program

Argparse is used as follows. First we need to import it

```python
import argparse
```

And then we can define a 'parser' which will parse our command line. Additionally we can provide a description field which tells people what our tool does:

```python
parser = argparse.ArgumentParser(description='Process some integers.')
```

And finally we can define some arguments that are available. Just like we have arguments to functions, we have arguments to command lines. These come in two flavours:

- required (without a `--`)
- optional "flags" (prefixed with `--`)

Here we have an argument named 'integers', which validates that all input values are of the type `int`. `nargs` is the number of arguments, `+` means '1 or more'. And we have some help text as well:

```python
parser.add_argument('integer', type=int, help='an integer parameter')
parser.add_argument('many_integers', type=int, nargs='+', help='an integer parameter')
```

We can also define an optional flag, here it's called `--sum`. We use `store_true` which will set it as true if the flag is used , otherwise false.

```python
parser.add_argument('--sum', action='store_true', help='Should we sum up the integers?')
```

Finally we parse the arguments, which reads `sys.argv` and processes it according to the above rules. The output is stored in `args`.

```python
args = parser.parse_args()
```

We have two main variables we can use now:

```
args.integer # A single integer
args.many_integers # A list of ints
args.sum # A boolean, True or False.
```

## Using argparse

Let's go back to our script, and replace `sys` with argparse.

> <hands-on-title>Replacing argv.</hands-on-title>
>
> 1. Given the following script, replace the use of `argv` with argparse.
>
>    ```python
>    import sys
>
>    result = 0
>
>    if len(sys.argv) == 1:
>        print("no arguments were supplied")
>    else:
>        for arg in sys.argv[1:]:
>            try:
>                result += float(arg)
>            except:
>                print(f"Could not parse {arg}")
>
>        print(result)
>    ```
>    {: .hands_on}
>
>    You should have one argument: numbers (type=float)
>
>    And print out the sum of those numbers.
>
>    > <question-title></question-title>
>    > How does your final script look?
>    > > <solution-title></solution-title>
>    > > ```python
>    > > import argparse
>    > >
>    > > parser = argparse.ArgumentParser(description='Sum some numbers')
>    > > parser.add_argument('integers', type=float, nargs='+',
>    > >                     help='a number to sum up.')
>    > > args = parser.parse_args()
>    > >
>    > > print(sum(args.integers))
>    > > ```
>    > {: .solution}
>    {: .question}
>
> 2. Try running the script with various values
>
>    ```bash
>    python run.py
>    python run.py 1 3 5
>    python run.py 2 4 O
>    python run.py --help
>    ```
{: .hands_on}

Wow that's a lot simpler! We have to learn how `argparse` is invoked but it handles a lot of cases for us:

- No arguments provided
- Responding to `--help`
- Raising an error for invalid values

`--help` is even written for us, without us writing any special code to handle that case! This is why you need to use `argparse`:
- It handles a lot of cases and input validation for you
- It produces a nice `--help` text that can help you if you've forgotten what your tool does
- It's nice for users of your scripts! They don't have to read the code to know how it behaves if you document it well.

There is a lot of documentation in the [`argparse`](https://docs.python.org/3/library/argparse.html) module for all sorts of use cases!


## Why Argparse?

Using argparse can be a big change to your tool but there are some benefits to using it!

1. Standardised interface to your tool that's familiar to everyone who uses command line tools
2. Automatic Help page
3. Automatic Galaxy Tools?

### Generating Automatic Galaxy Tools (Optional)

With the `argparse2tool` project, and eventually `pyGalGen` which will be merged into `planemo`, you can generate Galaxy tools automatically from `argparse` based Python scripts.

> <hands-on-title>Generate a Galaxy tool wrapper from your script</hands-on-title>
> 1. Write out the python script to a file named `main.py`
>
>    ```python
>    import argparse
>
>    parser = argparse.ArgumentParser(description='Sum some numbers')
>    parser.add_argument('integers', type=float, nargs='+',
>                        help='a number to sum up.')
>    args = parser.parse_args()
>
>    print(sum(args.integers))
>    ```
>
> 2. Create a virtual environment, just in case: ``
>
>    ```console
>    python -m venv .venv
>    . .venv/bin/activate
>    ```
>
> 3. Install `argparse2tool` via pip:
>
>    ```console
>    pip install argparse2tool
>    ```
>
> 4. Generate the tool interface:
>
>    > <code-in-title>Command</code-in-title>
>    > ```
>    > PYTHONPATH=$(argparse2tool) python main.py --generate_galaxy_xml
>    > ```
>    {: .code-in }
>
>    > <code-out-title>Galaxy XML</code-out-title>
>    > ```xml
>    > <tool name="main.py" id="main.py" version="1.0">
>    >   <description>Sum some numbers</description>
>    >   <stdio>
>    >     <exit_code range="1:" level="fatal"/>
>    >   </stdio>
>    >   <version_command><![CDATA[python main.py --version]]></version_command>
>    >   <command><![CDATA[python main.py
>    > #set repeat_var_1 = '" "'.join([ str($var.integers) for $var in $repeat_1 ])
>    > "$repeat_var_1"
>    >
>    > > $default]]></command>
>    >   <inputs>
>    >     <repeat title="repeat_title" min="1" name="repeat_1">
>    >       <param label="a number to sum up." value="0" type="float" name="integers"/>
>    >     </repeat>
>    >   </inputs>
>    >   <outputs>
>    >     <data name="default" format="txt" hidden="false"/>
>    >   </outputs>
>    >   <help><![CDATA[TODO: Write help]]></help>
>    > </tool>
>    >
>    > ```
>    {: .code-out }
>
{: .hands_on}
