---
layout: tutorial_hands_on

title: Python - Loops
level: Introductory
requirements:
 - type: internal
   topic_name: data-science
   tutorials:
     - python-iterables
     - python-flow
follow_up_training: []
questions:
- "How can I make a program do many things?"

objectives:
- "Explain what for loops are normally used for."
- "Trace the execution of a simple (unnested) loop and correctly state the values of variables in each iteration."
- "Write for loops that use the Accumulator pattern to aggregate values."

time_estimation:  40M
key_points:
- "A *for loop* executes commands once for each value in a collection."
- "A `for` loop is made up of a collection, a loop variable, and a body."
- "The first line of the `for` loop must end with a colon, and the body must be indented."
- "Indentation is always meaningful in Python."
- "Loop variables can be called anything (but it is strongly advised to have a meaningful name to the looping variable)."
- "The body of a loop can contain many statements."
- "Use `range` to iterate over a sequence of numbers."
- "The Accumulator pattern turns many values into one."

enable: false
subtopic: python-modular
contributors:
  - carpentries
  - hexylena
  - dirowa

priority: 1
notebook:
  language: python
---

A *for loop* tells Python to execute some statements once for each value in a list, a character string, or some other collection: "for each thing in this group, do these operations"

> ### {% icon comment %} Comment
> This tutorial is **significantly** based on [the Carpentries](https://carpentries.org) [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation/), [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation/), and [Plotting and Programming in Python](https://swcarpentry.github.io/python-novice-gapminder/), which are licensed CC-BY 4.0.
>
> Adaptations have been made to make this work better in a GTN/Galaxy environment.
{: .comment}

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# For Loops

Which of these would you rather write

> > ### {% icon code-in %} Manually
> > ```
> > print(2)
> > print(3)
> > print(5)
> > print(7)
> > print(11)
> > ```
> {: .code-in}
>
> > ### {% icon code-out %} With Loops
> > ```python
> > for number in [2, 3, 5, 7, 11]:
> >     print(number)
> > ```
> {: .code-out}
{: .code-2col}

It may be less clear here, since you just need to do one operation (`print`) but if you had to do two operations, three, more?

## Structure

## A `for` loop is made up of a collection, a loop variable, and a body.

```python
for number in [2, 3, 5]:
    doubled = number * 2
    print(f"{number} doubled is {doubled}")
```

- `number` - this is the loop variable. It's a new variable, that's assigned to the values from the collection.
- the collection, `[2, 3, 5]` is a `list` of numbers which we can tell from the square brackets used: `[`, `]`
- the loop body, where we double a number and the print out a message. The loop body is what gets executed for every iteration of the loop.

> > ### {% icon code-in %} The loop
> > ```
> > for number in [2, 3, 5]:
> >     doubled = number * 2
> >     print(f"{number} doubled is {doubled}")
> > ```
> {: .code-in}
>
> > ### {% icon code-out %} What's really happening internally
> > ```python
> > # First iteration, number = 2
> > doubled = number * 2
> > print(f"{number} doubled is {doubled}")
> > # Second iteration, number = 3
> > doubled = number * 3
> > print(f"{number} doubled is {doubled}")
> > # Third iteration, number = 5
> > doubled = number * 5
> > print(f"{number} doubled is {doubled}")
> > ```
> {: .code-out}
{: .code-2col}

Writing loops saves us time and makes sure our code is accurate, that we don't accidentally introduce a typo somewhere in the loop body.


## Indentation

The first line of the `for` loop must end with a colon, and the body must be indented.

> ### {% icon tip %} Tip: Blocks in Python
>
> The colon at the end of the first line signals the start of a *block* of statements.
>
> ```
> for x in y:
>     print(x)
> ```
>
> or
>
> ```
> if x > 10:
>     print(x)
> ```
>
> or even further nesting is possible:
>
> ```
> for x in y:
>     if x > 10:
>         print(x)
> ```
>
{: .tip}

The indentation is in fact, quite necessary. Notice how this fails:

```python
for number in [2, 3, 5]:
print(number)
```

And, likewise, this:

```python
patient1 = "z2910"
  patient2 = "y9583"
```

> ### {% icon question %} Question: Correct the errors
>
> ```
> data = [1, 3, 5, 9]
> for i in data:
> if i < 4:
> print(i)
> else:
> print(i * 2)
> print(len(data))
> ```
>
> > ### {% icon solution %} Solution
> >
> >
> > ```python
> > data = [1, 3, 5, 9]
> > for i in data:
> >     if i < 4:
> >         print(i)
> >     else:
> >         print(i * 2)
> > # But what about this line:
> > print(len(data))
> > ```
> >
> > Here this code is actually ambiguous, we don't know how indented `print(len(data))` should be. This very synthetic example lacks context, but there are three places it could be, with three different effects.
> > The first option, no indentation, prints out the list length at the end.
> >
> > ```python
> > [...]
> >     else:
> >         print(i * 2)
> > print(len(data))
> > ```
> >
> > The second, prints out the length once per loop iteration.
> >
> > ```python
> > data = [1, 3, 5, 9]
> > for i in data:
> >     if i < 4:
> >         print(i)
> >     else:
> >         print(i * 2)
> >     print(len(data))
> > ```
> >
> > First the `if` and `else` case resolve themselves, and then `len(data)` is printed every single loop iteration.
> > The last case only prints out the length if `i >= 4`.
> >
> > ```python
> > data = [1, 3, 5, 9]
> > for i in data:
> >     if i < 4:
> >         print(i)
> >     else:
> >         print(i * 2)
> >         print(len(data))
> > ```
> >
> {: .solution}
{: .question}

```python
# Check your answers here!
data = [1, 3, 5, 9]
for i in data:
if i < 4:
print(i)
else:
print(i * 2)
print(len(data))
```

## Variable Naming

Loop variables can be called anything, `i`, `j`, and `k` are very commong defaults due to their long history of use in other programing languages.
As with all variables, loop variables are: Created on demand, and Meaningless; their names can be anything at all.

```python
for kitten in [2, 3, 5]:
    print(kitten)
```

But *meaningless* is bad for variable names, and whenever possible, we should strive to pick useful, accurate variable names that help use remember what's going on:

```
for sequence in sequences:
    print()
for patient in clinic_patients:
    print()
for nucleotide in dna_sequence:
    print()
```

## Range

You can use `range` to iterate over a sequence of numbers. This is a built in function (check `help(range)`!) so it's always available even if you don't `import` anything. The range produced is non-inclusive: `range(N)` is the numbers `0` to `N-1`, so the result will be exactly the length you requested.

```python
for number in range(3):
    print(number)
```

However in python `range` is a special type of iterable: none of the numbers are created until we need them.

```python
print(range(5))
print(range(-3, 8)[0:4])
```

The easiest way to see what numbers are actually in there is to convert it to a `list`:

```python
print(list(range(5)))
print(list(range(-3, 8)))
print(list(range(0, 10, 2)))
```

## Accumulation

In programming you'll often want to accumulate some values: counting things (or "accumulating"). The pattern consists of creating a variable to store your result, running a loop over some data, and in that loop, adding to the variable for your result.

```python
# Sum the first 10 integers.
total = 0
for number in range(1, 11):
    total = total + (number)
print(total)
```

But how did we get that result? We can add some "debugging" lines to the above code to figure out how we got to that result. Try adding the following line in the above loop

```
print(f'Currently {number}, our total is {total}')
```

You can add it before you update `total`, after it, or both! Compare the outputs to understand what's happening on each line.

> ### {% icon question %} Question: Tracing Execution
>
> Create a table showing the numbers of the lines that are executed when this program runs,
> and the values of the variables after each line is executed.
>
> ```
> total = 0
> for char in "tin":
>     total = total + 1
> ```
> > ### {% icon solution %} Solution
> >
> > | Line | Variables            |
> > |------|----------------------|
> > | 1    | total = 0            |
> > | 2    | total = 0 char = 't' |
> > | 3    | total = 1 char = 't' |
> > | 2    | total = 1 char = 'i' |
> > | 3    | total = 2 char = 'i' |
> > | 2    | total = 2 char = 'n' |
> > | 3    | total = 3 char = 'n' |
> {: .solution}
{: .question}

```python
#Test your code here!
```

> ### {% icon question %} Question: Reversing a String
>
> Fill in the blanks in the program below so that it prints "stressed"
> (the reverse of the original character string "desserts").
>
> ```
> original = "stressed"
> result = ____
> for char in original:
>     result = ____
> print(result)
> ```
> > ### {% icon solution %} Solution
> > ```
> > original = "stressed"
> > result = ""
> > for char in original:
> >     result = char + result
> > print(result)
> > ```
> {: .solution}
{: .question}

```python
# Test your code here!
original = "stressed"
result = ____
for char in original:
    result = ____
print(result)
```


> ### {% icon question %} Question: Practice Accumulating
>
> Fill in the blanks in each of the programs below
> to produce the indicated result.
>
> ```
> # Total length of the strings in the list: ["red", "green", "blue"] => 12
> total = 0
> for word in ["red", "green", "blue"]:
>     ____ = ____ + len(word)
> print(total)
> ```
>
> > ### {% icon solution %} Solution
> > ```
> > total = 0
> > for word in ["red", "green", "blue"]:
> >     total = total + len(word)
> > print(total)
> > ```
> {: .solution}
>
> ```
> # List of word lengths: ["red", "green", "blue"] => [3, 5, 4]
> lengths = ____
> for word in ["red", "green", "blue"]:
>     lengths.____(____)
> print(lengths)
> ```
> > ### {% icon solution %} Solution
> > ```
> > lengths = []
> > for word in ["red", "green", "blue"]:
> >     lengths.append(len(word))
> > print(lengths)
> > ```
> {: .solution}
>
> ```
> # Concatenate all words: ["red", "green", "blue"] => "redgreenblue"
> words = ["red", "green", "blue"]
> result = ____
> for ____ in ____:
>     ____
> print(result)
> ```
> > ### {% icon solution %} Solution
> > ```
> > words = ["red", "green", "blue"]
> > result = ""
> > for word in words:
> >     result = result + word
> > print(result)
> > ```
> {: .solution}
>
> __Create an acronym:__ Starting from the list `["red", "green", "blue"]`, create the acronym `"RGB"` using
> a for loop.
>
> __Hint:__ You may need to use a string method to properly format the acronym.
> > ### {% icon solution %} Solution
> > ```
> > acronym = ""
> > for word in ["red", "green", "blue"]:
> >     acronym = acronym + word[0].upper()
> > print(acronym)
> > ```
> {: .solution}
{: .question}

```python
#Test your code here!
```

> ## Cumulative Sum
>
> Reorder and properly indent the lines of code below
> so that they print a list with the cumulative sum of data.
> The result should be `[1, 3, 5, 10]`.
>
> ```
> cumulative.append(total)
> for number in data:
> cumulative = []
> total += number
> total = 0
> print(cumulative)
> data = [1,2,2,5]
> ```
> > ### {% icon solution %} Solution
> > ```
> > total = 0
> > data = [1,2,2,5]
> > cumulative = []
> > for number in data:
> >     total += number
> >     cumulative.append(total)
> > print(cumulative)
> > ```
> {: .solution}
{: .question}

```python
# Test your code here!
```

> ## Identifying Variable Name Errors
>
> 1. Read the code below and try to identify what the errors are
>    **without** running it.
> 2. Run the code and read the error message.
>    What type of `NameError` do you think this is?
>    Is it a string with no quotes, a misspelled variable, or a
>    variable that should have been defined but was not?
> 3. Fix the error.
> 4. Repeat steps 2 and 3, until you have fixed all the errors.
>
> ```
> for number in range(10):
>     # use a if the number is a multiple of 3, otherwise use b
>     if Number % 3 == 0:
>         message = message + a
>     else:
>         message = message + "b"
> print(message)
> ```
>
> > ### {% icon solution %} Solution
> > - Python variable names are case sensitive: `number` and `Number` refer to different variables.
> > - The variable `message` needs to be initialized as an empty string.
> > - We want to add the string `"a"` to `message`, not the undefined variable `a`.
> >
> > ```
> > message = ""
> > for number in range(10):
> >     # use a if the number is a multiple of 3, otherwise use b
> >     if number % 3 == 0:
> >         message = message + "a"
> >     else:
> >         message = message + "b"
> > print(message)
> > ```
> {: .solution}
{: .question}

```python
# Fix me!
for number in range(10):
    # use a if the number is a multiple of 3, otherwise use b
    if Number % 3 == 0:
        message = message + a
    else:
        message = message + "b"
print(message)
```

> ### {% icon question %} Question: Identifying Item Errors
>
> 1. Read the code below and try to identify what the errors are
>    **without** running it.
> 2. Run the code, and read the error message. What type of error is it?
> 3. Fix the error.
>
> ```
> seasons = ['Spring', 'Summer', 'Fall', 'Winter']
> print(f'My favorite season is {seasons[4]}')
> ```
>
> > ### {% icon solution %} Solution
> > This list has 4 elements and the index to access the last element in the list is `3`.
> > ```
> > seasons = ['Spring', 'Summer', 'Fall', 'Winter']
> > print(f'My favorite season is {seasons[3]}')
> > ```
> {: .solution}
{: .question}

```
# Fix me!
seasons = ['Spring', 'Summer', 'Fall', 'Winter']
print(f'My favorite season is {seasons[4]}')
```
