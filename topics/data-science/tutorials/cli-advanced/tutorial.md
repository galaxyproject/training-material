---
layout: tutorial_hands_on

title: Advanced CLI in Galaxy
level: Intermediate
zenodo_link: ""
requirements:
- type: "internal"
  topic_name: data-science
  tutorials:
      - cli-basics
follow_up_training:
- type: "internal"
  topic_name: data-science
  tutorials:
      - cli-bashcrawl

questions:
- "How can I combine existing commands to do new things?"
- "How can I perform the same actions on many different files?"
- "How can I find files?"
- "How can I find things in files?"
objectives:
- "Redirect a command's output to a file."
- "Process a file instead of keyboard input using redirection."
- "Construct command pipelines with two or more stages."
- "Explain what usually happens if a program or pipeline isn't given any input to process."
- "Explain Unix's 'small pieces, loosely joined' philosophy."
- "Write a loop that applies one or more commands separately to each file in a set of files."
- "Trace the values taken on by a loop variable during execution of the loop."
- "Explain the difference between a variable's name and its value."
- "Explain why spaces and some punctuation characters shouldn't be used in file names."
- "Demonstrate how to see what commands have recently been executed."
- "Re-run recently executed commands without retyping them."
- "Use `grep` to select lines from text files that match simple patterns."
- "Use `find` to find files and directories whose names match simple patterns."
- "Use the output of one command as the command-line argument(s) to another command."
- "Explain what is meant by 'text' and 'binary' files, and why many common tools don't handle the latter well."
time_estimation: 2H
key_points:
- "`wc` counts lines, words, and characters in its inputs."
- "`cat` displays the contents of its inputs."
- "`sort` sorts its inputs."
- "`head` displays the first 10 lines of its input."
- "`tail` displays the last 10 lines of its input."
- "`command > [file]` redirects a command's output to a file (overwriting any existing content)."
- "`command >> [file]` appends a command's output to a file."
- "`[first] | [second]` is a pipeline: the output of the first command is used as the input to the second."
- "The best way to use the shell is to use pipes to combine simple single-purpose programs (filters)."
- "A `for` loop repeats commands once for every thing in a list."
- "Every `for` loop needs a variable to refer to the thing it is currently operating on."
- "Use `$name` to expand a variable (i.e., get its value). `${name}` can also be used."
- "Do not use spaces, quotes, or wildcard characters such as '*' or '?' in filenames, as it complicates variable expansion."
- "Give files consistent names that are easy to match with wildcard patterns to make it easy to select them for looping."
- "Use the up-arrow key to scroll up through previous commands to edit and repeat them."
- "Use <kbd>Ctrl</kbd>+<kbd>R</kbd> to search through the previously entered commands."
- "Use `history` to display recent commands, and `![number]` to repeat a command by number."
- "`find` finds files with specific properties that match patterns."
- "`grep` selects lines in files that match patterns."
- "`--help` is an option supported by many bash commands, and programs that can be run from within Bash, to display more information on how to use these commands or programs."
- "`man [command]` displays the manual page for a given command."
- "`$([command])` inserts a command's output in place."
notebook:
  language: bash
subtopic: bash
contributors:
  - carpentries
  - hexylena
  - bazante1
  - erasmusplus
  - avans-atgm
tags:
- bash
---

This tutorial will walk you through the basics of how to use the Unix command line.

> <comment-title></comment-title>
>
> This tutorial is **significantly** based on [the Carpentries](https://carpentries.org) ["The Unix Shell"](https://swcarpentry.github.io/shell-novice/) lesson, which is licensed CC-BY 4.0. Adaptations have been made to make this work better in a GTN/Galaxy environment.
>
{: .comment}


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Pipes and Filtering

Now that we know a few basic commands,
we can finally look at the shell's most powerful feature:
the ease with which it lets us combine existing programs in new ways.
We'll start with the directory called `shell-lesson-data/molecules`
that contains six files describing some simple organic molecules.
The `.pdb` extension indicates that these files are in Protein Data Bank format,
a simple text format that specifies the type and position of each atom in the molecule.

```bash
cd ~/Desktop/shell-lesson-data/
ls molecules
```


Let's go into that directory with `cd` and run an example command `wc cubane.pdb`:

```bash
cd molecules
wc cubane
```

`wc` is the 'word count' command:
it counts the number of lines, words, and characters in files (from left to right, in that order).

If we run the command `wc *.pdb`, the `*` in `*.pdb` matches zero or more characters,
so the shell turns `*.pdb` into a list of all `.pdb` files in the current directory:

```bash
wc *.pdb
```

Note that `wc *.pdb` also shows the total number of all lines in the last line of the output.

If we run `wc -l` instead of just `wc`,
the output shows only the number of lines per file:

```bash
wc -l *.pdb
```

The `-m` and `-w` options can also be used with the `wc` command, to show
only the number of characters or the number of words in the files.

> <tip-title>Why Isn't It Doing Anything?</tip-title>
>
> What happens if a command is supposed to process a file, but we
> don't give it a filename? For example, what if we type:
>
> ```bash
> $ wc -l
> ```
>
> but don't type `*.pdb` (or anything else) after the command?
> Since it doesn't have any filenames, `wc` assumes it is supposed to
> process input given at the command prompt, so it just sits there and waits for us to give
> it some data interactively. From the outside, though, all we see is it
> sitting there: the command doesn't appear to do anything.
>
> If you make this kind of mistake, you can escape out of this state by holding down
> the control key (<kbd>Ctrl</kbd>) and typing the letter <kbd>C</kbd> once and
> letting go of the <kbd>Ctrl</kbd> key.
> <kbd>Ctrl</kbd>+<kbd>C</kbd>
{: .tip}


## Capturing output from commands

Which of these files contains the fewest lines?
It's an easy question to answer when there are only six files,
but what if there were 6000?
Our first step toward a solution is to run the command:

```bash
wc -l *.pdb > lengths.txt
```

The greater than symbol, `>`, tells the shell to **redirect** the command's output
to a file instead of printing it to the screen. (This is why there is no screen output:
everything that `wc` would have printed has gone into the
file `lengths.txt` instead.)  The shell will create
the file if it doesn't exist. If the file exists, it will be
silently overwritten, which may lead to data loss and thus requires
some caution.

> <tip-title>No <code>&gt;</code> on an AZERTY keyboard?</tip-title>
> You can rewrite this using the tee command which writes out a file, while also showing the output to `stdout`.
> ```bash
> wc -l *.pdb | tee lengths.txt
> ```
> Or you can use copy and paste to copy the `>` character from the materials.
{: .tip}

`ls lengths.txt` confirms that the file exists:

```bash
ls lengths.txt
```

We can now send the content of `lengths.txt` to the screen using `cat lengths.txt`.
The `cat` command gets its name from 'concatenate' i.e. join together,
and it prints the contents of files one after another.
There's only one file in this case,
so `cat` just shows us what it contains:

```bash
cat lengths.txt
```

> <tip-title>Output Page by Page</tip-title>
>
> We'll continue to use `cat` in this lesson, for convenience and consistency,
> but it has the disadvantage that it always dumps the whole file onto your screen.
> More useful in practice is the command `less`,
> which you use with `less lengths.txt`.
> This displays a screenful of the file, and then stops.
> You can go forward one screenful by pressing the spacebar,
> or back one by pressing `b`.  Press `q` to quit.
{: .tip}


## Filtering output

Next we'll use the `sort` command to sort the contents of the `lengths.txt` file.
But first we'll use an exercise to learn a little about the sort command:

> <question-title>What Does `sort -n` Do?</question-title>
>
> The file `shell-lesson-data/numbers.txt`
> contains the following lines:
>
> ```
> 10
> 2
> 19
> 22
> 6
> ```
>
> If we run `sort` on this file, the output is:
>
> ```
> 10
> 19
> 2
> 22
> 6
> ```
>
> If we run `sort -n` on the same file, we get this instead:
>
> ```
> 2
> 6
> 10
> 19
> 22
> ```
>
> Explain why `-n` has this effect.
>
> > <solution-title></solution-title>
> > The `-n` option specifies a numerical rather than an alphanumerical sort.
> {: .solution}
{: .question}

We will also use the `-n` option to specify that the sort is
numerical instead of alphanumerical.
This does *not* change the file;
instead, it sends the sorted result to the screen:

```bash
sort -n lengths.txt
```


We can put the sorted list of lines in another temporary file called `sorted-lengths.txt`
by putting `> sorted-lengths.txt` after the command,
just as we used `> lengths.txt` to put the output of `wc` into `lengths.txt`.
Once we've done that,
we can run another command called `head` to get the first few lines in `sorted-lengths.txt`:

```bash
sort -n lengths.txt > sorted-lengths.txt
```

Using `-n 1` with `head` tells it that
we only want the first line of the file;
`-n 20` would get the first 20,
and so on.
Since `sorted-lengths.txt` contains the lengths of our files ordered from least to greatest,
the output of `head` must be the file with the fewest lines.

> <tip-title>Redirecting to the same file</tip-title>
>
> It's a very bad idea to try redirecting
> the output of a command that operates on a file
> to the same file. For example:
>
> ```
> $ sort -n lengths.txt > lengths.txt
> ```
>
> Doing something like this may give you
> incorrect results and/or delete
> the contents of `lengths.txt`.
{: .tip}

> <question-title>What Does `>>` Mean?</question-title>
>
> We have seen the use of `>`, but there is a similar operator `>>`
> which works slightly differently.
> We'll learn about the differences between these two operators by printing some strings.
> We can use the `echo` command to print strings e.g.
>
> > <code-in-title>Bash</code-in-title>
> > ```
> > $ echo The echo command prints text
> > ```
> {: .code-in}
>
> > <code-out-title></code-out-title>
> > ```
> > The echo command prints text
> > ```
> {: .code-out}
>
> Now test the commands below to reveal the difference between the two operators:
>
> > <code-in-title>Bash</code-in-title>
> > ```
> > $ echo hello > testfile01.txt
> > ```
> {: .code-in}
>
> and:
>
> > <code-in-title>Bash</code-in-title>
> > ```
> > $ echo hello >> testfile02.txt
> > ```
> {: .code-in}
>
> **Hint**: Try executing each command twice in a row and then examining the output files.
>
> > ## Solution
> > In the first example with `>`, the string 'hello' is written to `testfile01.txt`,
> > but the file gets overwritten each time we run the command.
> >
> > We see from the second example that the `>>` operator also writes 'hello' to a file
> > (in this case`testfile02.txt`),
> > but appends the string to the file if it already exists
> > (i.e. when we run it for the second time).
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Appending Data</question-title>
>
> We have already met the `head` command, which prints lines from the start of a file.
> `tail` is similar, but prints lines from the end of a file instead.
>
> Consider the file `shell-lesson-data/data/animals.txt`.
> After these commands, select the answer that
> corresponds to the file `animals-subset.txt`:
>
> ```
> $ head -n 3 animals.txt > animals-subset.txt
> $ tail -n 2 animals.txt >> animals-subset.txt
> ```
>
> 1. The first three lines of `animals.txt`
> 2. The last two lines of `animals.txt`
> 3. The first three lines and the last two lines of `animals.txt`
> 4. The second and third lines of `animals.txt`
>
> > <solution-title></solution-title>
> > Option 3 is correct.
> > For option 1 to be correct we would only run the `head` command.
> > For option 2 to be correct we would only run the `tail` command.
> > For option 4 to be correct we would have to pipe the output of `head` into `tail -n 2`
> > by doing `head -n 3 animals.txt | tail -n 2 > animals-subset.txt`
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```


## Passing output to another command
In our example of finding the file with the fewest lines,
we are using two intermediate files `lengths.txt` and `sorted-lengths.txt` to store output.
This is a confusing way to work because
even once you understand what `wc`, `sort`, and `head` do,
those intermediate files make it hard to follow what's going on.
We can make it easier to understand by running `sort` and `head` together:

```bash
sort -n lengths.txt | head -n 1
```

The vertical bar, `|`, between the two commands is called a **pipe**.
It tells the shell that we want to use
the output of the command on the left
as the input to the command on the right.

This has removed the need for the `sorted-lengths.txt` file.

## Combining multiple commands
Nothing prevents us from chaining pipes consecutively.
We can for example send the output of `wc` directly to `sort`,
and then the resulting output to `head`.
This removes the need for any intermediate files.

We'll start by using a pipe to send the output of `wc` to `sort`:

```bash
wc -l *.pdb | sort -n
```

We can then send that output through another pipe, to `head`, so that the full pipeline becomes:

```bash
wc -l *.pdb | sort -n | head -n 1
```

This is exactly like a mathematician nesting functions like *log(3x)*
and saying 'the log of three times *x*'.
In our case,
the calculation is 'head of sort of line count of `*.pdb`'.


The redirection and pipes used in the last few commands are illustrated below:

![Redirects and Pipes of different commands](../../images/carpentries-cli/redirects-and-pipes.svg)

`wc -l *.pdb` will direct the output to the shell. `wc -l *.pdb > lengths` will
direct output to the file lengths. `wc -l *.pdb | sort -n | head -n 1` will
build a pipeline where the output of the wc command is the input to the sort
command, the output of the sort command is the input to the head command and
the output of the head command is directed to the shell

> <question-title>Piping Commands Together</question-title>
>
> In our current directory, we want to find the 3 files which have the least number of
> lines. Which command listed below would work?
>
> 1. `wc -l * > sort -n > head -n 3`
> 2. `wc -l * | sort -n | head -n 1-3`
> 3. `wc -l * | head -n 3 | sort -n`
> 4. `wc -l * | sort -n | head -n 3`
>
> > <solution-title></solution-title>
> > Option 4 is the solution.
> > The pipe character `|` is used to connect the output from one command to
> > the input of another.
> > `>` is used to redirect standard output to a file.
> > Try it in the `shell-lesson-data/molecules` directory!
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```


## Tools designed to work together
This idea of linking programs together is why Unix has been so successful.
Instead of creating enormous programs that try to do many different things,
Unix programmers focus on creating lots of simple tools that each do one job well,
and that work well with each other.
This programming model is called 'pipes and filters'.
We've already seen pipes;
a **filter** is a program like `wc` or `sort`
that transforms a stream of input into a stream of output.
Almost all of the standard Unix tools can work this way:
unless told to do otherwise,
they read from standard input,
do something with what they've read,
and write to standard output.

The key is that any program that reads lines of text from standard input
and writes lines of text to standard output
can be combined with every other program that behaves this way as well.
You can *and should* write your programs this way
so that you and other people can put those programs into pipes to multiply their power.


> <question-title>Pipe Reading Comprehension</question-title>
>
> A file called `animals.txt` (in the `shell-lesson-data/data` folder) contains the following data:
>
> ```
> 2012-11-05,deer
> 2012-11-05,rabbit
> 2012-11-05,raccoon
> 2012-11-06,rabbit
> 2012-11-06,deer
> 2012-11-06,fox
> 2012-11-07,rabbit
> 2012-11-07,bear
> ```
>
> What text passes through each of the pipes and the final redirect in the pipeline below?
>
> ```
> $ cat animals.txt | head -n 5 | tail -n 3 | sort -r > final.txt
> ```
>
> Hint: build the pipeline up one command at a time to test your understanding
>
> > <solution-title></solution-title>
> > The `head` command extracts the first 5 lines from `animals.txt`.
> > Then, the last 3 lines are extracted from the previous 5 by using the `tail` command.
> > With the `sort -r` command those 3 lines are sorted in reverse order and finally,
> > the output is redirected to a file `final.txt`.
> > The content of this file can be checked by executing `cat final.txt`.
> > The file should contain the following lines:
> > ```
> > 2012-11-06,rabbit
> > 2012-11-06,deer
> > 2012-11-05,raccoon
> > ```
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Pipe Construction</question-title>
>
> For the file `animals.txt` from the previous exercise, consider the following command:
>
> ```
> $ cut -d , -f 2 animals.txt
> ```
>
> The `cut` command is used to remove or 'cut out' certain sections of each line in the file,
> and `cut` expects the lines to be separated into columns by a <kbd>Tab</kbd> character.
> A character used in this way is a called a **delimiter**.
> In the example above we use the `-d` option to specify the comma as our delimiter character.
> We have also used the `-f` option to specify that we want to extract the second field (column).
> This gives the following output:
>
> ```
> deer
> rabbit
> raccoon
> rabbit
> deer
> fox
> rabbit
> bear
> ```
>
> The `uniq` command filters out adjacent matching lines in a file.
> How could you extend this pipeline (using `uniq` and another command) to find
> out what animals the file contains (without any duplicates in their
> names)?
>
> > <solution-title></solution-title>
> > ```
> > $ cut -d , -f 2 animals.txt | sort | uniq
> > ```
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Which Pipe?</question-title>
>
> The file `animals.txt` contains 8 lines of data formatted as follows:
>
> ```
> 2012-11-05,deer
> 2012-11-05,rabbit
> 2012-11-05,raccoon
> 2012-11-06,rabbit
> ...
> ```
>
> The `uniq` command has a `-c` option which gives a count of the
> number of times a line occurs in its input.  Assuming your current
> directory is `shell-lesson-data/data/`, what command would you use to produce
> a table that shows the total count of each type of animal in the file?
>
> 1.  `sort animals.txt | uniq -c`
> 2.  `sort -t, -k2,2 animals.txt | uniq -c`
> 3.  `cut -d, -f 2 animals.txt | uniq -c`
> 4.  `cut -d, -f 2 animals.txt | sort | uniq -c`
> 5.  `cut -d, -f 2 animals.txt | sort | uniq -c | wc -l`
>
> > <solution-title></solution-title>
> > Option 4. is the correct answer.
> > If you have difficulty understanding why, try running the commands, or sub-sections of
> > the pipelines (make sure you are in the `shell-lesson-data/data` directory).
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

## Nelle's Pipeline: Checking Files

Nelle has run her samples through the assay machines
and created 17 files in the `north-pacific-gyre/2012-07-03` directory described earlier.
As a quick check, starting from her home directory, Nelle types:

```bash
cd ~/Desktop/shell-lesson-data/north-pacific-gyre/2012-07-03
wc -l *.txt
```

The output is 18 lines that look like this:

```
300 NENE01729A.txt
300 NENE01729B.txt
300 NENE01736A.txt
300 NENE01751A.txt
300 NENE01751B.txt
300 NENE01812A.txt
... ...
```

Now she types this:

```bash
wc -l *.txt | sort -n | head -n 5
```

Whoops: one of the files is 60 lines shorter than the others.
When she goes back and checks it,
she sees that she did that assay at 8:00 on a Monday morning --- someone
was probably in using the machine on the weekend,
and she forgot to reset it.
Before re-running that sample,
she checks to see if any files have too much data:

```bash
wc -l *.txt | sort -n | tail -n 5
```

Those numbers look good --- but what's that 'Z' doing there in the third-to-last line?
All of her samples should be marked 'A' or 'B';
by convention,
her lab uses 'Z' to indicate samples with missing information.
To find others like it, she does this:

```bash
ls *Z.txt
```

Sure enough,
when she checks the log on her laptop,
there's no depth recorded for either of those samples.
Since it's too late to get the information any other way,
she must exclude those two files from her analysis.
She could delete them using `rm`,
but there are actually some analyses she might do later where depth doesn't matter,
so instead, she'll have to be careful later on to select files using the wildcard expressions
`NENE*A.txt NENE*B.txt`.


> <question-title>Removing Unneeded Files</question-title>
>
> Suppose you want to delete your processed data files, and only keep
> your raw files and processing script to save storage.
> The raw files end in `.dat` and the processed files end in `.txt`.
> Which of the following would remove all the processed data files,
> and *only* the processed data files?
>
> 1. `rm ?.txt`
> 2. `rm *.txt`
> 3. `rm * .txt`
> 4. `rm *.*`
>
> > <solution-title></solution-title>
> > 1. This would remove `.txt` files with one-character names
> > 2. This is correct answer
> > 3. The shell would expand `*` to match everything in the current directory,
> > so the command would try to remove all matched files and an additional
> > file called `.txt`
> > 4. The shell would expand `*.*` to match all files with any extension,
> > so this command would delete all files
> {: .solution}
{: .question}


# Loops

**Loops** are a programming construct which allow us to repeat a command or set of commands
for each item in a list.
As such they are key to productivity improvements through automation.
Similar to wildcards and tab completion, using loops also reduces the
amount of typing required (and hence reduces the number of typing mistakes).

Suppose we have several hundred genome data files named `basilisk.dat`, `minotaur.dat`, and
`unicorn.dat`.
For this example, we'll use the `creatures` directory which only has three example files,
but the principles can be applied to many many more files at once. First, go
into the creatures directory.

```bash
# Change directories here!

```

The structure of these files is the same: the common name, classification, and updated date are
presented on the first three lines, with DNA sequences on the following lines.
Let's look at the files:

```bash
head -n 5 basilisk.dat minotaur.dat unicorn.dat
```

We would like to print out the classification for each species, which is given on the second
line of each file.
For each file, we would need to execute the command `head -n 2` and pipe this to `tail -n 1`.
We’ll use a loop to solve this problem, but first let’s look at the general form of a loop:

```
for thing in list_of_things
do
    operation_using $thing    # Indentation within the loop is not required, but aids legibility
done
```

and we can apply this to our example like this:

```bash
for filename in basilisk.dat minotaur.dat unicorn.dat
do
    head -n 2 $filename | tail -n 1
done
```


> <tip-title>Follow the Prompt</tip-title>
>
> The shell prompt changes from `$` to `>` and back again as we were
> typing in our loop. The second prompt, `>`, is different to remind
> us that we haven't finished typing a complete command yet. A semicolon, `;`,
> can be used to separate two commands written on a single line.
{: .tip}

When the shell sees the keyword `for`,
it knows to repeat a command (or group of commands) once for each item in a list.
Each time the loop runs (called an iteration), an item in the list is assigned in sequence to
the **variable**, and the commands inside the loop are executed, before moving on to
the next item in the list.
Inside the loop,
we call for the variable's value by putting `$` in front of it.
The `$` tells the shell interpreter to treat
the variable as a variable name and substitute its value in its place,
rather than treat it as text or an external command.

In this example, the list is three filenames: `basilisk.dat`, `minotaur.dat`, and `unicorn.dat`.
Each time the loop iterates, it will assign a file name to the variable `filename`
and run the `head` command.
The first time through the loop,
`$filename` is `basilisk.dat`.
The interpreter runs the command `head` on `basilisk.dat`
and pipes the first two lines to the `tail` command,
which then prints the second line of `basilisk.dat`.
For the second iteration, `$filename` becomes
`minotaur.dat`. This time, the shell runs `head` on `minotaur.dat`
and pipes the first two lines to the `tail` command,
which then prints the second line of `minotaur.dat`.
For the third iteration, `$filename` becomes
`unicorn.dat`, so the shell runs the `head` command on that file,
and `tail` on the output of that.
Since the list was only three items, the shell exits the `for` loop.

> <tip-title>Same Symbols, Different Meanings</tip-title>
>
> Here we see `>` being used as a shell prompt, whereas `>` is also
> used to redirect output.
> Similarly, `$` is used as a shell prompt, but, as we saw earlier,
> it is also used to ask the shell to get the value of a variable.
>
> If the *shell* prints `>` or `$` then it expects you to type something,
> and the symbol is a prompt.
>
> If *you* type `>` or `$` yourself, it is an instruction from you that
> the shell should redirect output or get the value of a variable.
{: .tip}

When using variables it is also
possible to put the names into curly braces to clearly delimit the variable
name: `$filename` is equivalent to `${filename}`, but is different from
`${file}name`. You may find this notation in other people's programs.

We have called the variable in this loop `filename`
in order to make its purpose clearer to human readers.
The shell itself doesn't care what the variable is called;
if we wrote this loop as:

```bash
for x in basilisk.dat minotaur.dat unicorn.dat
do
    head -n 2 $x | tail -n 1
done
```

or:

```bash
for temperature in basilisk.dat minotaur.dat unicorn.dat
do
    head -n 2 $temperature | tail -n 1
done
```

it would work exactly the same way.

**Don't do this.**

Programs are only useful if people can understand them,
so meaningless names (like `x`) or misleading names (like `temperature`)
increase the odds that the program won't do what its readers think it does.

> <question-title>Variables in Loops</question-title>
>
> This exercise refers to the `shell-lesson-data/molecules` directory.
> `ls` gives the following output:
>
> ```
> cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> ```
>
> What is the output of the following code?
>
> ```
> for datafile in *.pdb
> do
>     ls *.pdb
> done
> ```
>
> Now, what is the output of the following code?
>
> ```
> for datafile in *.pdb
> do
>    ls $datafile
> done
> ```
>
> Why do these two loops give different outputs?
>
> > <solution-title></solution-title>
> > The first code block gives the same output on each iteration through
> > the loop.
> > Bash expands the wildcard `*.pdb` within the loop body (as well as
> > before the loop starts) to match all files ending in `.pdb`
> > and then lists them using `ls`.
> > The expanded loop would look like this:
> > ```
> > $ for datafile in cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > > do
> > >     ls cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > > done
> > ```
> >
> > ```
> > cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > cubane.pdb  ethane.pdb  methane.pdb  octane.pdb  pentane.pdb  propane.pdb
> > ```
> >
> > The second code block lists a different file on each loop iteration.
> > The value of the `datafile` variable is evaluated using `$datafile`,
> > and then listed using `ls`.
> >
> > ```
> > cubane.pdb
> > ethane.pdb
> > methane.pdb
> > octane.pdb
> > pentane.pdb
> > propane.pdb
> > ```
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Limiting Sets of Files</question-title>
>
> What would be the output of running the following loop in thei
> `shell-lesson-data/molecules` directory?
>
> ```
> for filename in c*
> do
>     ls $filename
> done
> ```
>
> 1.  No files are listed.
> 2.  All files are listed.
> 3.  Only `cubane.pdb`, `octane.pdb` and `pentane.pdb` are listed.
> 4.  Only `cubane.pdb` is listed.
>
> > <solution-title></solution-title>
> > 4 is the correct answer. `*` matches zero or more characters, so any file name starting with
> > the letter c, followed by zero or more other characters will be matched.
> {: .solution}
>
> How would the output differ from using this command instead?
>
> ```
> for filename in *c*
> do
>     ls $filename
> done
> ```
>
> 1.  The same files would be listed.
> 2.  All the files are listed this time.
> 3.  No files are listed this time.
> 4.  The files `cubane.pdb` and `octane.pdb` will be listed.
> 5.  Only the file `octane.pdb` will be listed.
>
> > <solution-title></solution-title>
> > 4 is the correct answer. `*` matches zero or more characters, so a file name with zero or more
> > characters before a letter c and zero or more characters after the letter c will be matched.
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Saving to a File in a Loop - Part One</question-title>
>
> In the `shell-lesson-data/molecules` directory, what is the effect of this loop?
>
> ```
> for alkanes in *.pdb
> do
>     echo $alkanes
>     cat $alkanes > alkanes.pdb
> done
> ```
>
> 1.  Prints `cubane.pdb`, `ethane.pdb`, `methane.pdb`, `octane.pdb`, `pentane.pdb` and
>    `propane.pdb`, and the text from `propane.pdb` will be saved to a file called `alkanes.pdb`.
> 2.  Prints `cubane.pdb`, `ethane.pdb`, and `methane.pdb`, and the text from all three files
>     would be concatenated and saved to a file called `alkanes.pdb`.
> 3.  Prints `cubane.pdb`, `ethane.pdb`, `methane.pdb`, `octane.pdb`, and `pentane.pdb`,
>     and the text from `propane.pdb` will be saved to a file called `alkanes.pdb`.
> 4.  None of the above.
>
> > <solution-title></solution-title>
> > 1 is correct. The text from each file in turn gets written to the `alkanes.pdb` file.
> > However, the file gets overwritten on each loop iteration, so the final content of `alkanes.pdb`
> > is the text from the `propane.pdb` file.
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Saving to a File in a Loop - Part Two</question-title>
>
> Also in the `shell-lesson-data/molecules` directory,
> what would be the output of the following loop?
>
> ```
> for datafile in *.pdb
> do
>     cat $datafile >> all.pdb
> done
> ```
>
> 1.  All of the text from `cubane.pdb`, `ethane.pdb`, `methane.pdb`, `octane.pdb`, and
>     `pentane.pdb` would be concatenated and saved to a file called `all.pdb`.
> 2.  The text from `ethane.pdb` will be saved to a file called `all.pdb`.
> 3.  All of the text from `cubane.pdb`, `ethane.pdb`, `methane.pdb`, `octane.pdb`, `pentane.pdb`
>     and `propane.pdb` would be concatenated and saved to a file called `all.pdb`.
> 4.  All of the text from `cubane.pdb`, `ethane.pdb`, `methane.pdb`, `octane.pdb`, `pentane.pdb`
>     and `propane.pdb` would be printed to the screen and saved to a file called `all.pdb`.
>
> > <solution-title></solution-title>
> > 3 is the correct answer. `>>` appends to a file, rather than overwriting it with the redirected
> > output from a command.
> > Given the output from the `cat` command has been redirected, nothing is printed to the screen.
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

Let's continue with our example in the `shell-lesson-data/creatures` directory.
Here's a slightly more complicated loop:

```bash
for filename in *.dat
do
    echo $filename
    head -n 100 $filename | tail -n 20
done
```

The shell starts by expanding `*.dat` to create the list of files it will process.
The **loop body**
then executes two commands for each of those files.
The first command, `echo`, prints its command-line arguments to standard output.
For example:

```bash
echo hello there
```

prints:

```
hello there
```

In this case,
since the shell expands `$filename` to be the name of a file,
`echo $filename` prints the name of the file.
Note that we can't write this as:

```bash
for filename in *.dat
do
    $filename
    head -n 100 $filename | tail -n 20
done
```

because then the first time through the loop,
when `$filename` expanded to `basilisk.dat`, the shell would try to run `basilisk.dat` as a program.
Finally,
the `head` and `tail` combination selects lines 81-100
from whatever file is being processed
(assuming the file has at least 100 lines).

> <tip-title>Spaces in Names</tip-title>
>
> Spaces are used to separate the elements of the list
> that we are going to loop over. If one of those elements
> contains a space character, we need to surround it with
> quotes, and do the same thing to our loop variable.
> Suppose our data files are named:
>
> ```
> red dragon.dat
> purple unicorn.dat
> ```
>
> To loop over these files, we would need to add double quotes like so:
>
> ```
> $ for filename in "red dragon.dat" "purple unicorn.dat"
> > do
> >     head -n 100 "$filename" | tail -n 20
> > done
> ```
>
> It is simpler to avoid using spaces (or other special characters) in filenames.
>
> The files above don't exist, so if we run the above code, the `head` command will be unable
> to find them, however the error message returned will show the name of the files it is
> expecting:
>
> ```
> head: cannot open ‘red dragon.dat’ for reading: No such file or directory
> head: cannot open ‘purple unicorn.dat’ for reading: No such file or directory
> ```
>
> Try removing the quotes around `$filename` in the loop above to see the effect of the quote
> marks on spaces. Note that we get a result from the loop command for unicorn.dat
> when we run this code in the `creatures` directory:
>
> ```
> head: cannot open ‘red’ for reading: No such file or directory
> head: cannot open ‘dragon.dat’ for reading: No such file or directory
> head: cannot open ‘purple’ for reading: No such file or directory
> CGGTACCGAA
> AAGGGTCGCG
> CAAGTGTTCC
> ...
> ```
{: .tip}

We would like to modify each of the files in `shell-lesson-data/creatures`, but also save a version
of the original files, naming the copies `original-basilisk.dat` and `original-unicorn.dat`.
We can't use:

```
cp *.dat original-*.dat
```

because that would expand to:

```bash
cp basilisk.dat minotaur.dat unicorn.dat original-*.dat
```

This wouldn't back up our files, instead we get an error.

This problem arises when `cp` receives more than two inputs. When this happens, it
expects the last input to be a directory where it can copy all the files it was passed.
Since there is no directory named `original-*.dat` in the `creatures` directory we get an
error.

Instead, we can use a loop:
```bash
for filename in *.dat
do
    cp $filename original-$filename
done
```

This loop runs the `cp` command once for each filename.
The first time,
when `$filename` expands to `basilisk.dat`,
the shell executes:

```
cp basilisk.dat original-basilisk.dat
```

The second time, the command is:

```
cp minotaur.dat original-minotaur.dat
```

The third and last time, the command is:

```
cp unicorn.dat original-unicorn.dat
```

Since the `cp` command does not normally produce any output, it's hard to check
that the loop is doing the correct thing.
However, we learned earlier how to print strings using `echo`, and we can modify the loop
to use `echo` to print our commands without actually executing them.
As such we can check what commands *would be* run in the unmodified loop.

The following diagram
shows what happens when the modified loop is executed, and demonstrates how the
judicious use of `echo` is a good debugging technique.

![The for loop 'for filename in *.dat; do echo cp $filename original-$filename; done' will successively assign the names of all '*.dat' files in your current directory to the variable '$filename' and then execute the command. With the files 'basilisk.dat', 'minotaur.dat' and 'unicorn.dat' in the current directory the loop will successively call the echo command three times and print three lines: 'cp basislisk.dat original-basilisk.dat', then 'cp minotaur.dat original-minotaur.dat' and finally 'cp unicorn.dat original-unicorn.dat'](../../images/carpentries-cli/shell_script_for_loop_flow_chart.svg)

## Nelle's Pipeline: Processing Files

Nelle is now ready to process her data files using `goostats.sh` ---
a shell script written by her supervisor.
This calculates some statistics from a protein sample file, and takes two arguments:

1. an input file (containing the raw data)
2. an output file (to store the calculated statistics)

Since she's still learning how to use the shell,
she decides to build up the required commands in stages.
Her first step is to make sure that she can select the right input files --- remember,
these are ones whose names end in 'A' or 'B', rather than 'Z'.
Starting from her home directory, Nelle types:

```bash
cd ~/Desktop/shell-lesson-data/north-pacific-gyre/2012-07-03
for datafile in NENE*A.txt NENE*B.txt
do
    echo $datafile
done
```

Her next step is to decide
what to call the files that the `goostats.sh` analysis program will create.
Prefixing each input file's name with 'stats' seems simple,
so she modifies her loop to do that:

```bash
for datafile in NENE*A.txt NENE*B.txt
do
    echo $datafile stats-$datafile
done
```

She hasn't actually run `goostats.sh` yet,
but now she's sure she can select the right files and generate the right output filenames.

> <tip-title>Top Terminal Tip: Re-running previous commands</tip-title>
> Typing in commands over and over again is becoming tedious,
> though,
> and Nelle is worried about making mistakes,
> so instead of re-entering her loop,
> she presses <kbd>↑</kbd>.
> In response,
> the shell redisplays the whole loop on one line
> (using semi-colons to separate the pieces):
>
> ```
> for datafile in NENE*A.txt NENE*B.txt; do echo $datafile stats-$datafile; done
> ```
{: .tip}

Using the left arrow key,
Nelle backs up and changes the command `echo` to `bash goostats.sh`:

```
for datafile in NENE*A.txt NENE*B.txt; do bash goostats.sh $datafile stats-$datafile; done
```

When she presses <kbd>Enter</kbd>,
the shell runs the modified command.
However, nothing appears to happen --- there is no output.
After a moment, Nelle realizes that since her script doesn't print anything to the screen
any longer, she has no idea whether it is running, much less how quickly.
She kills the running command by typing <kbd>Ctrl</kbd>+<kbd>C</kbd>,
uses <kbd>↑</kbd> to repeat the command,
and edits it to read:

```
for datafile in NENE*A.txt NENE*B.txt; do echo $datafile; bash goostats.sh $datafile stats-$datafile; done
```

> <tip-title>Beginning and End</tip-title>
>
> We can move to the beginning of a line in the shell by typing <kbd>Ctrl</kbd>+<kbd>A</kbd>
> and to the end using <kbd>Ctrl</kbd>+<kbd>E</kbd>.
{: .tip}

When she runs her program now,
it produces one line of output every five seconds or so:

1518 times 5 seconds,
divided by 60,
tells her that her script will take about two hours to run.
As a final check,
she opens another terminal window,
goes into `north-pacific-gyre/2012-07-03`,
and uses `cat stats-NENE01729B.txt`
to examine one of the output files.
It looks good,
so she decides to get some coffee and catch up on her reading.

> <tip-title>Those Who Know History Can Choose to Repeat It</tip-title>
>
> Another way to repeat previous work is to use the `history` command to
> get a list of the last few hundred commands that have been executed, and
> then to use `!123` (where '123' is replaced by the command number) to
> repeat one of those commands. For example, if Nelle types this:
>
> ```
> $ history | tail -n 5
> ```
>
> ```
>   456  ls -l NENE0*.txt
>   457  rm stats-NENE01729B.txt.txt
>   458  bash goostats.sh NENE01729B.txt stats-NENE01729B.txt
>   459  ls -l NENE0*.txt
>   460  history
> ```
>
> then she can re-run `goostats.sh` on `NENE01729B.txt` simply by typing
> `!458`. This number will be different for you, you should check your history before running it!
{: .tip}

> <tip-title>Other History Commands</tip-title>
>
> There are a number of other shortcut commands for getting at the history.
>
> - <kbd>Ctrl</kbd>+<kbd>R</kbd> enters a history search mode 'reverse-i-search' and finds the
> most recent command in your history that matches the text you enter next.
> Press <kbd>Ctrl</kbd>+<kbd>R</kbd> one or more additional times to search for earlier matches.
> You can then use the left and right arrow keys to choose that line and edit
> it then hit <kbd>Return</kbd> to run the command.
> - `!!` retrieves the immediately preceding command
> (you may or may not find this more convenient than using <kbd>↑</kbd>)
> - `!$` retrieves the last word of the last command.
> That's useful more often than you might expect: after
> `bash goostats.sh NENE01729B.txt stats-NENE01729B.txt`, you can type
> `less !$` to look at the file `stats-NENE01729B.txt`, which is
> quicker than doing <kbd>↑</kbd> and editing the command-line.
{: .tip}

> <question-title>Doing a Dry Run</question-title>
>
> A loop is a way to do many things at once --- or to make many mistakes at
> once if it does the wrong thing. One way to check what a loop *would* do
> is to `echo` the commands it would run instead of actually running them.
>
> Suppose we want to preview the commands the following loop will execute
> without actually running those commands:
>
> ```
> cd ~/Desktop/shell-lesson-data/pdb/
> for datafile in *.pdb
> do
>     cat $datafile >> all.pdb
> done
> ```
>
> What is the difference between the two loops below, and which one would we
> want to run?
>
> > <code-in-title>Version 1</code-in-title>
> > ```
> > for datafile in *.pdb
> > do
> >     echo cat $datafile >> all.pdb
> > done
> > ```
> {: .code-in}
>
> > <code-in-title>Version 2</code-in-title>
> > ```
> > for datafile in *.pdb
> > do
> >     echo "cat $datafile >> all.pdb"
> > done
> > ```
> {: .code-in}
>
> > <solution-title>Solution</solution-title>
> > The second version is the one we want to run.
> > This prints to screen everything enclosed in the quote marks, expanding the
> > loop variable name because we have prefixed it with a dollar sign.
> >
> > The first version appends the output from the command `echo cat $datafile`
> > to the file, `all.pdb`. This file will just contain the list;
> > `cat cubane.pdb`, `cat ethane.pdb`, `cat methane.pdb` etc.
> >
> > Try both versions for yourself to see the output! Be sure to open the
> > `all.pdb` file to view its contents.
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Nested Loops</question-title>
>
> Suppose we want to set up a directory structure to organize
> some experiments measuring reaction rate constants with different compounds
> *and* different temperatures.  What would be the
> result of the following code:
>
> ```
> for species in cubane ethane methane
> do
>     for temperature in 25 30 37 40
>     do
>         mkdir $species-$temperature
>     done
> done
> ```
>
> > <solution-title></solution-title>
> > We have a nested loop, i.e. contained within another loop, so for each species
> > in the outer loop, the inner loop (the nested loop) iterates over the list of
> > temperatures, and creates a new directory for each combination.
> >
> > Try running the code for yourself to see which directories are created!
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```


# Finding Things

In the same way that many of us now use 'Google' as a
verb meaning 'to find', Unix programmers often use the
word 'grep'.
'grep' is a contraction of 'global/regular expression/print',
a common sequence of operations in early Unix text editors.
It is also the name of a very useful command-line program.

`grep` finds and prints lines in files that match a pattern.
For our examples,
we will use a file that contains three haiku taken from a
1998 competition in *Salon* magazine. For this set of examples,
we're going to be working in the writing subdirectory:

```bash
cd
cd Desktop/shell-lesson-data/writing
cat haiku.txt
```


> <tip-title>Forever, or Five Years</tip-title>
>
> We haven't linked to the original haiku because
> they don't appear to be on *Salon*'s site any longer.
> As [Jeff Rothenberg said](https://www.clir.org/wp-content/uploads/sites/6/ensuring.pdf),
> 'Digital information lasts forever --- or five years, whichever comes first.'
> Luckily, popular content often [has backups](http://wiki.c2.com/?ComputerErrorHaiku).
{: .tip}

Let's find lines that contain the word 'not':

```bash
grep not haiku.txt
```

Here, `not` is the pattern we're searching for.
The grep command searches through the file, looking for matches to the pattern specified.
To use it type `grep`, then the pattern we're searching for and finally
the name of the file (or files) we're searching in.

The output is the three lines in the file that contain the letters 'not'.

By default, grep searches for a pattern in a case-sensitive way.
In addition, the search pattern we have selected does not have to form a complete word,
as we will see in the next example.

Let's search for the pattern: 'The'.

```bash
grep The haiku.txt
```

This time, two lines that include the letters 'The' are outputted,
one of which contained our search pattern within a larger word, 'Thesis'.

To restrict matches to lines containing the word 'The' on its own,
we can give `grep` with the `-w` option.
This will limit matches to word boundaries.

Later in this lesson, we will also see how we can change the search behavior of grep
with respect to its case sensitivity.

```bash
grep -w The haiku.txt
```

Note that a 'word boundary' includes the start and end of a line, so not
just letters surrounded by spaces.
Sometimes we don't
want to search for a single word, but a phrase. This is also easy to do with
`grep` by putting the phrase in quotes.

```bash
grep -w "is not" haiku.txt
```

We've now seen that you don't have to have quotes around single words,
but it is useful to use quotes when searching for multiple words.
It also helps to make it easier to distinguish between the search term or phrase
and the file being searched.
We will use quotes in the remaining examples.

Another useful option is `-n`, which numbers the lines that match:

```bash
grep -n "it" haiku.txt
```

Here, we can see that lines 5, 9, and 10 contain the letters 'it'.

We can combine options (i.e. flags) as we do with other Unix commands.
For example, let's find the lines that contain the word 'the'.
We can combine the option `-w` to find the lines that contain the word 'the'
and `-n` to number the lines that match:

```bash
grep -n -w "the" haiku.txt
```

Now we want to use the option `-i` to make our search case-insensitive:

```bash
grep -n -w -i "the" haiku.txt
```

Now, we want to use the option `-v` to invert our search, i.e., we want to output
the lines that do not contain the word 'the'.

```bash
grep -n -w -v "the" haiku.txt
```


If we use the `-r` (recursive) option,
`grep` can search for a pattern recursively through a set of files in subdirectories.

Let's search recursively for `Yesterday` in the `shell-lesson-data/writing` directory:

```bash
grep -r Yesterday .
```

`grep` has lots of other options. To find out what they are, we can type:

```bash
grep --help
```

> <question-title>Using `grep`</question-title>
>
> Which command would result in the following output:
>
> ```
> and the presence of absence:
> ```
>
> 1. `grep "of" haiku.txt`
> 2. `grep -E "of" haiku.txt`
> 3. `grep -w "of" haiku.txt`
> 4. `grep -i "of" haiku.txt`
>
> > <solution-title></solution-title>
> > The correct answer is 3, because the `-w` option looks only for whole-word matches.
> > The other options will also match 'of' when part of another word.
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <tip-title>Wildcards</tip-title>
>
> `grep`'s real power doesn't come from its options, though; it comes from
> the fact that patterns can include wildcards. (The technical name for
> these is **regular expressions**, which
> is what the 're' in 'grep' stands for.) Regular expressions are both complex
> and powerful; if you want to do complex searches, please look at the lesson
> on [our website](http://v4.software-carpentry.org/regexp/index.html). As a taster, we can
> find lines that have an 'o' in the second position like this:
>
> ```
> $ grep -E "^.o" haiku.txt
> ```
>
> ```
> You bring fresh toner.
> Today it is not working
> Software is like that.
> ```
>
> We use the `-E` option and put the pattern in quotes to prevent the shell
> from trying to interpret it. (If the pattern contained a `*`, for
> example, the shell would try to expand it before running `grep`.) The
> `^` in the pattern anchors the match to the start of the line. The `.`
> matches a single character (just like `?` in the shell), while the `o`
> matches an actual 'o'.
{: .tip}

```bash
# Explore the possible solutions here!
```

> <question-title>Tracking a Species</question-title>
>
> Leah has several hundred
> data files saved in one directory, each of which is formatted like this:
>
> ```
> 2013-11-05,deer,5
> 2013-11-05,rabbit,22
> 2013-11-05,raccoon,7
> 2013-11-06,rabbit,19
> 2013-11-06,deer,2
> ```
>
> She wants to write a shell script that takes a species as the first command-line argument
> and a directory as the second argument. The script should return one file called `species.txt`
> containing a list of dates and the number of that species seen on each date.
> For example using the data shown above, `rabbit.txt` would contain:
>
> ```
> 2013-11-05,22
> 2013-11-06,19
> ```
>
> Put these commands and pipes in the right order to achieve this:
>
> ```
> cut -d : -f 2
> >
> |
> grep -w $1 -r $2
> |
> $1.txt
> cut -d , -f 1,3
> ```
>
> Hint: use `man grep` to look for how to grep text recursively in a directory
> and `man cut` to select more than one field in a line.
>
> An example of such a file is provided in `shell-lesson-data/data/animal-counts/animals.txt`
>
> > <solution-title></solution-title>
> >
> > ```
> > grep -w $1 -r $2 | cut -d : -f 2 | cut -d , -f 1,3 > $1.txt
> > ```
> >
> > Actually, you can swap the order of the two cut commands and it still works. At the
> > command line, try changing the order of the cut commands, and have a look at the output
> > from each step to see why this is the case.
> >
> > You would call the script above like this:
> >
> > ```
> > $ bash count-species.sh bear .
> > ```
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <question-title>Little Women</question-title>
>
> You and your friend, having just finished reading *Little Women* by
> Louisa May Alcott, are in an argument.  Of the four sisters in the
> book, Jo, Meg, Beth, and Amy, your friend thinks that Jo was the
> most mentioned.  You, however, are certain it was Amy.  Luckily, you
> have a file `LittleWomen.txt` containing the full text of the novel
> (`shell-lesson-data/writing/data/LittleWomen.txt`).
> Using a `for` loop, how would you tabulate the number of times each
> of the four sisters is mentioned?
>
> Hint: one solution might employ
> the commands `grep` and `wc` and a `|`, while another might utilize
> `grep` options.
> There is often more than one way to solve a programming task, so a
> particular solution is usually chosen based on a combination of
> yielding the correct result, elegance, readability, and speed.
>
> > <solution-title></solution-title>
> > ```
> > for sis in Jo Meg Beth Amy
> > do
> > 	echo $sis:
> >	grep -ow $sis LittleWomen.txt | wc -l
> > done
> > ```
> >
> > Alternative, slightly inferior solution:
> > ```
> > for sis in Jo Meg Beth Amy
> > do
> > 	echo $sis:
> >	grep -ocw $sis LittleWomen.txt
> > done
> > ```
> >
> > This solution is inferior because `grep -c` only reports the number of lines matched.
> > The total number of matches reported by this method will be lower if there is more
> > than one match per line.
> >
> > Perceptive observers may have noticed that character names sometimes appear in all-uppercase
> > in chapter titles (e.g. 'MEG GOES TO VANITY FAIR').
> > If you wanted to count these as well, you could add the `-i` option for case-insensitivity
> > (though in this case, it doesn't affect the answer to which sister is mentioned
> > most frequently).
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

While `grep` finds lines in files,
the `find` command finds files themselves.
Again,
it has a lot of options;
to show how the simplest ones work, we'll use the directory tree shown below.

![A file tree under the directory 'writing' contians several sub-directories and files such that 'writing' contains directories 'data', 'thesis', 'tools' and a file 'haiku.txt'; 'writing/data' contains the files 'Little Women.txt', 'one.txt' and 'two.txt'; 'writing/thesis' contains the file 'empty-draft.md'; 'writing/tools' contains the directory 'old' and the files 'format' and 'stats'; and 'writing/tools/old' contains a file 'oldtool'](../../images/carpentries-cli/find-file-tree.svg)

Nelle's `writing` directory contains one file called `haiku.txt` and three subdirectories:
`thesis` (which contains a sadly empty file, `empty-draft.md`);
`data` (which contains three files `LittleWomen.txt`, `one.txt` and `two.txt`);
and a `tools` directory that contains the programs `format` and `stats`,
and a subdirectory called `old`, with a file `oldtool`.

For our first command,
let's run `find .` (remember to run this command from the `shell-lesson-data/writing` folder).

```bash
find .
```

As always,
the `.` on its own means the current working directory,
which is where we want our search to start.
`find`'s output is the names of every file **and** directory
under the current working directory.
This can seem useless at first but `find` has many options
to filter the output and in this lesson we will discover some
of them.

The first option in our list is
`-type d` that means 'things that are directories'.
Sure enough,
`find`'s output is the names of the five directories in our little tree
(including `.`):

```bash
find . -type d
```


Notice that the objects `find` finds are not listed in any particular order.
If we change `-type d` to `-type f`,
we get a listing of all the files instead:

```bash
find . -type f
```

Now let's try matching by name:

```bash
find . -name *.txt
```

We expected it to find all the text files,
but it only prints out `./haiku.txt`.
The problem is that the shell expands wildcard characters like `*` *before* commands run.
Since `*.txt` in the current directory expands to `haiku.txt`,
the command we actually ran was:

```bash
find . -name haiku.txt
```

`find` did what we asked; we just asked for the wrong thing.

To get what we want,
let's do what we did with `grep`:
put `*.txt` in quotes to prevent the shell from expanding the `*` wildcard.
This way,
`find` actually gets the pattern `*.txt`, not the expanded filename `haiku.txt`:

```bash
find . -name "*.txt"
```

> <tip-title>Listing vs. Finding</tip-title>
>
> `ls` and `find` can be made to do similar things given the right options,
> but under normal circumstances,
> `ls` lists everything it can,
> while `find` searches for things with certain properties and shows them.
{: .tip}

As we said earlier,
the command line's power lies in combining tools.
We've seen how to do that with pipes;
let's look at another technique.
As we just saw,
`find . -name "*.txt"` gives us a list of all text files in or below the current directory.
How can we combine that with `wc -l` to count the lines in all those files?

The simplest way is to put the `find` command inside `$()`:

```bash
wc -l $(find . -name "*.txt")
```

When the shell executes this command,
the first thing it does is run whatever is inside the `$()`.
It then replaces the `$()` expression with that command's output.
Since the output of `find` is the four filenames `./data/one.txt`, `./data/LittleWomen.txt`,
`./data/two.txt`, and `./haiku.txt`, the shell constructs the command:

```bash
wc -l ./data/one.txt ./data/LittleWomen.txt ./data/two.txt ./haiku.txt
```

which is what we wanted.
This expansion is exactly what the shell does when it expands wildcards like `*` and `?`,
but lets us use any command we want as our own 'wildcard'.

It's very common to use `find` and `grep` together.
The first finds files that match a pattern;
the second looks for lines inside those files that match another pattern.
Here, for example, we can find PDB files that contain iron atoms
by looking for the string 'FE' in all the `.pdb` files above the current directory:

```bash
grep "FE" $(find .. -name "*.pdb")
```

> <question-title>Matching and Subtracting</question-title>
>
> The `-v` option to `grep` inverts pattern matching, so that only lines
> which do *not* match the pattern are printed. Given that, which of
> the following commands will find all files in `/data` whose names
> end in `s.txt` but whose names also do *not* contain the string `net`?
> (For example, `animals.txt` or `amino-acids.txt` but not `planets.txt`.)
> Once you have thought about your answer, you can test the commands in the `shell-lesson-data`
> directory.
>
> 1.  `find data -name "*s.txt" | grep -v net`
> 2.  `find data -name *s.txt | grep -v net`
> 3.  `grep -v "net" $(find data -name "*s.txt")`
> 4.  None of the above.
>
> > <solution-title></solution-title>
> > The correct answer is 1. Putting the match expression in quotes prevents the shell
> > expanding it, so it gets passed to the `find` command.
> >
> > Option 2 is incorrect because the shell expands `*s.txt` instead of passing the wildcard
> > expression to `find`.
> >
> > Option 3 is incorrect because it searches the contents of the files for lines which
> > do not match 'net', rather than searching the file names.
> {: .solution}
{: .question}

```bash
# Explore the possible solutions here!
```

> <tip-title>Binary Files</tip-title>
>
> We have focused exclusively on finding patterns in text files. What if
> your data is stored as images, in databases, or in some other format?
>
> A handful of tools extend `grep` to handle a few non text formats. But a
> more generalizable approach is to convert the data to text, or
> extract the text-like elements from the data. On the one hand, it makes simple
> things easy to do. On the other hand, complex things are usually impossible. For
> example, it's easy enough to write a program that will extract X and Y
> dimensions from image files for `grep` to play with, but how would you
> write something to find values in a spreadsheet whose cells contained
> formulas?
>
> A last option is to recognize that the shell and text processing have
> their limits, and to use another programming language.
> When the time comes to do this, don't be too hard on the shell: many
> modern programming languages have borrowed a lot of
> ideas from it, and imitation is also the sincerest form of praise.
{: .tip}

The Unix shell is older than most of the people who use it. It has
survived so long because it is one of the most productive programming
environments ever created --- maybe even *the* most productive. Its syntax
may be cryptic, but people who have mastered it can experiment with
different commands interactively, then use what they have learned to
automate their work. Graphical user interfaces may be easier to use at
first, but once learned, the productivity in the shell is unbeatable.
And as Alfred North Whitehead wrote in 1911, 'Civilization advances by
extending the number of important operations which we can perform
without thinking about them.'

> <question-title>`find` Pipeline Reading Comprehension</question-title>
>
> Write a short explanatory comment for the following shell script:
>
> ```
> wc -l $(find . -name "*.dat") | sort -n
> ```
>
> > <solution-title></solution-title>
> > 1. Find all files with a `.dat` extension recursively from the current directory
> > 2. Count the number of lines each of these files contains
> > 3. Sort the output from step 2. numerically
> {: .solution}
{: .question}

# Final Notes

All of the commands you have run up until now were ad-hoc, interactive commands.
