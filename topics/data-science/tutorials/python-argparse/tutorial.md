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
contributors:
  - hexylena
  - dirowa
  - bazante1

priority: 10
---

[`argparse`](https://docs.python.org/3/library/argparse.html) is an argument parsing library for Python that's part of the stdlib. It lets you make command line tools significantly nicer to work with.

> ### Agenda
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

> ### {% icon hands_on %} Hands-on: Print out argv
>
> 1. Create / open the file `run.py` in your text editor of choice
> 2. There we'll create a simple Python script that:
>   1. imports `sys`, the system module needed to access argv.
>   2. Prints out `sys.argv`
>
> > ### {% icon solution %} Solution
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

> ### {% icon question %} Question
> What did you notice about the output? There are two main points.
> > ### {% icon solution %} Solution
> > 1. The name of the script (`run.py`) is included as the first value every time.
> > 2. All of the arguments are passed as strings, no numbers.
> {: .solution}
{: .question}

## Simple tasks

Let's sum up all of the numbers passed on the command line. We'll do this by hand, and then we'll replace it with `argparse` to see how much effort that saves us.

> ### {% icon hands_on %} Hands-on
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
> > ### {% icon question %} Question
> > How does your updated script look?
> > > ### {% icon solution %} Solution
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
parser.add_argument('integers', type=int, nargs='+',
                    help='an integer for the accumulator')
```

We can also define an optional flag, here it's called `--sum`. Here it goes to a destination named 'accumulate', the name we'll use to access the value of this argument. It has an action of 'store_const' which just tracks if the flag was supplied or not.

The `const` attribute is set to `sum`, which is actually the function `sum()`, this is what the value will be if we run the command with `--sum`. Otherwise it will `default` to the function `max()`. We again have some help text to tell us how it behaves

```python
parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='sum the integers (default: find the max)')
```

Finally we parse the arguments, which reads `sys.argv` and processes it according to the above rules. The output is stored in `args`.

```python
args = parser.parse_args()
```

We have two main variables we can use now:

```
args.integers # A list of integers.
args.accumulate # Actually a function!
```

## Using argparse

Let's go back to our script, and replace `sys` with argparse.

> ### {% icon hands_on %} Hands-on: Replacing argv.
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
> > ### {% icon question %} Question
> > How does your final script look?
> > > ### {% icon solution %} Solution
> > > ```python
> > > import argparse
> > >
> > > parser = argparse.ArgumentParser(description='Sum some numbers')
> > > parser.add_argument('integers', type=float, nargs='+',
> > >                     help='a number to sum up.')
> > > args = parser.parse_args()
> > >
> > > print(sum(args.integers))
> > > ```
> > {: .solution}
> {: .question}
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

