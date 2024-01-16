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

subtopic: python-modular
contributors:
  - carpentries
  - hexylena
  - dirowa
  - bazante1

priority: 6
notebook:
  language: python
  pyolite: true
---

A *for loop* tells Python to execute some statements once for each value in a list, a character string, or some other collection: "for each thing in this group, do these operations"

> <comment-title></comment-title>
> This tutorial is **significantly** based on [the Carpentries](https://carpentries.org) [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation/), [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation/), and [Plotting and Programming in Python](https://swcarpentry.github.io/python-novice-gapminder/), which are licensed CC-BY 4.0.
>
> Adaptations have been made to make this work better in a GTN/Galaxy environment.
{: .comment}

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# For Loops

Which of these would you rather write

> > <code-in-title>Manually</code-in-title>
> > ```
> > print(2)
> > print(3)
> > print(5)
> > print(7)
> > print(11)
> > ```
> {: .code-in}
>
> > <code-out-title>With Loops</code-out-title>
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

- `number` - this is the loop variable. It's a new variable, that's assigned to the values from the collection. It does not need to be defined before the loop.
- the collection, `[2, 3, 5]` is a `list` of numbers which we can tell from the square brackets used: `[`, `]`
- the loop body, where we double a number and the print out a message. The loop body is what gets executed for every iteration of the loop.

> > <code-in-title>The loop</code-in-title>
> > ```python
> > for number in [2, 3, 5]:
> >     doubled = number * 2
> >     print(f"{number} doubled is {doubled}")
> > ```
> {: .code-in}
>
> > <code-out-title>What's really happening internally</code-out-title>
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

## Things You Can Loop Over

You can loop over characters in a string

```python
dna_string = 'ACTGGTCATCG'
for base in dna_string:
    print(base)
```

You can loop over lists:

```python
cast = ['Elphaba', 'Glinda', 'Fiyero', 'Nessarose']
for character in cast:
    print(character)
```

## Indentation

The first line of the `for` loop must end with a colon, and the body must be indented with *four spaces*. Many editors do this automatically for you and even convert <kbd>Tab</kbd>s into 4 spaces.

> <tip-title>Blocks in Python</tip-title>
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
#Fix me!
for number in [2, 3, 5]:
print(number)
```

And, likewise, this:

```python
patient1 = "z2910"
  patient2 = "y9583"
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
for number in range(10):
    print(number)
```

> <tip-title>Iterables can be weird</tip-title>
> In python `range` is a special type of iterable: none of the numbers are created until we need them.
>
> ```python
> print(range(5))
> print(range(-3, 8)[0:4])
> ```
>
> The easiest way to see what numbers are actually in there is to convert it to a `list`:
>
> ```python
> print(list(range(5)))
> print(list(range(-3, 8)))
> print(list(range(0, 10, 2)))
> ```
{: .tip}

## Accumulation

In programming you'll often want to accumulate some values: counting things (or "accumulating"). The pattern consists of creating a variable to store your result, running a loop over some data, and in that loop, adding to the variable for your result.

```python
# Sum the first 10 integers.
total = 0
for number in range(1, 11):
    total = total + (number)
print(f" final: {{ total }}")
```

But how did we get that result? We can add some "debugging" lines to the above code to figure out how we got to that result. Try adding the following line in the above loop

```
print(f'Currently {number}, our total is {total}')
```

You can add it before you update `total`, after it, or both! Compare the outputs to understand what's happening on each line.

> <tip-title>Controlling your loop!</tip-title>
>
> There are multiple ways to efficiently control your loop if you need it.
> these are the inbuilt python functions: continue & break
>
> when python encounters *continue* in your loop it will stop working and goes to the next iteration of the loop.
> ```
> for letter in 'Galaxy':
>     if letter == 'l':
>         continue
>     print(f'The letters are: {letter}')
>```
> with *break* python stops the loop and continues with the next part of the code like nothing happened
> ```
> for letter in 'Galaxy':
>     if letter == 'l':
>         break
>     print(f'The letters are: {letter}')
> print('Done')
> ```
{: .tip}

```python
# Test break and continue here
```

## Exercises

> <question-title>Tracing Execution</question-title>
>
> Create a table showing the numbers of the lines that are executed when this program runs,
> and the values of the variables after each line is executed.
>
> ```
> total = 0
> for char in "tin":
>     total = total + 1
> ```
> > <solution-title></solution-title>
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

> <question-title>Reversing a String</question-title>
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
> > <solution-title></solution-title>
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


> <question-title>Practice Accumulating</question-title>
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
> > <solution-title></solution-title>
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
> > <solution-title></solution-title>
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
> > <solution-title></solution-title>
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
> > <solution-title></solution-title>
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
> > <solution-title></solution-title>
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

> <question-title>A classic programmer test: Fizz Buzz</question-title>
> [FizzBuzz](https://en.wikipedia.org/wiki/Fizz_buzz) is a classic "test" question that is used in some job interviews to remove candidates who really do not understand programming. Your task is this:
>
> Write a for loop that loops over the numbers 1 to 50.
>
> - If the number is divisible by 3, write Fizz instead of the number
> - If the number is divisible by 5, write Buzz instead of the number
> - If the number is divisible by 3 and 5 both, write FizzBuzz instead of the number
> - Otherwise, write the number itself.
>
> > <solution-title></solution-title>
> > ```python
> > for i in range(1, 50):
> >     if i % 3 == 0 and i % 5 == 0:
> >         print("FizzBuzz")
> >     elif i % 3 == 0:
> >         print("Fizz")
> >     elif i % 5 == 0:
> >         print("Buzz")
> >     else:
> >         print(i)
> > ```
> {: .solution}
{: .question}

```python
# Do a FizzBuzz
```

> <question-title>Identifying Item Errors</question-title>
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
> > <solution-title></solution-title>
> > This list has 4 elements and the index to access the last element in the list is `3`.
> > ```
> > seasons = ['Spring', 'Summer', 'Fall', 'Winter']
> > print(f'My favorite season is {seasons[3]}')
> > ```
> {: .solution}
{: .question}

```python
# Fix me!
seasons = ['Spring', 'Summer', 'Fall', 'Winter']
print(f'My favorite season is {seasons[4]}')
```

> <question-title>Correct the errors</question-title>
>
> This code is completely missing indentation, it needs to be fixed. Can you make some guesses at how indented each line should be?
> ```
> data = [1, 3, 5, 9]
> acc = 0
> for i in data:
> if i < 4:
> acc = acc + i * 2
> else:
> acc = acc + i
> print(f'The value at {i} is {acc}')
> print(f'The answer is {acc}')
> ```
>
> > <solution-title></solution-title>
> >
> >
> > ```python
> > data = [1, 3, 5, 9]
> > acc = 0
> > # There is a : character at the end of this line, so you KNOW the next line
> > # must be indented.
> > for i in data:
> >     # Same here, another :
> >     if i < 4:
> >         acc = acc + i * 2
> >     # And again! Another :
> >     else:
> >         acc = acc + i
> > # But what about these lines?
> > print(f'The value at {i} is {acc}')
> > print(f'The answer is {acc}')
> > ```
> >
> > Here this code is actually ambiguous, we don't know how indented the two prints should be. This very synthetic example lacks good context, but there are three places it could be, with three different effects.
> >
> > There are two bits of knowledge we can use, however:
> > - the first print uses `i`, so it must be within the loop
> > - the second print cannot be indented more than the first print (Why? It would require a block like `for ... :` or `if .. :` to indent further.)
> >
> > The first option, no indentation, prints out the value once per loop, that seems good
> >
> > ```python
> > [...]
> >     else:
> >         acc = acc + i
> >     print(f'The value at {i} is {acc}')
> > ```
> >
> > The second, prints out the value only during the else case, not otherwise.
> >
> > ```python
> >     else:
> >         acc = acc + i
> >         print(f'The value at {i} is {acc}')
> > ```
> >
> > So that's probably wrong, and we should take the first option. That leaves two options for the final print, no indentation, or at the same level as our first print statement. We can guess that we probably want to print out the final result of the loop, and that it should not be indented.
> >
> > ```python
> > data = [1, 3, 5, 9]
> > acc = 0
> > for i in data:
> >     if i < 4:
> >         acc = acc + i * 2
> >     else:
> >         acc = acc + i
> >     print(f'The value at {i} is {acc}')
> > print(f'The answer is {acc}')
> > ```
> >
> {: .solution}
{: .question}

```python
# This code accidentally lost it's indentation! Can you fix it?
data = [1, 3, 5, 9]
acc = 0
for i in data:
if i < 4:
acc = acc + i * 2
else:
acc = acc + i
print(f'The value at {i} is {acc}')
print(f'The answer is {acc}')
```

> <question-title>Trimming a FASTQ string</question-title>
> Given a FASTQ string, and a list with quality scores, use `break` to print out just the good bit of DNA and it's quality score.
>
> ```
> # We've got a Read
> read = """
> @SEQ_ID
> GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
> +
> 55CCF>>>>>>CCCCCCC65!''*((((***+))%%%++)(%%%%).1***-+*''))**
> """.strip().split('\n')
>
> def quality_to_percent(q):
>     return 100 * (1 - (10 ** (q / -10)))
>
> sequence = read[1]
> quality_scores = [ord(x) - 33 for x in read[3]]
>
> for i in ... # TODO
> ```
>
> > <solution-title>Hint: Looping over two variables</solution-title>
> > There are two ways to do this, one you might be able to guess, and one that might be new:
> > 1. Loop over a `range()` using `len(sequence)`. Since `len(sequence) == len(quality_scores)`, when we access the Nth position of either, they match up.
> > 2. `zip(sequence, quality_scores)` will loop over both of these lists together. It produces a new list that looks like `[['G', 20], ['A', 20], ['T', 34]]`.
> {: .solution}
>
> > <solution-title></solution-title>
> > The naïve solution is quite easy and readable:
> > ```
> > for i in range(len(sequence)):
> >     if quality_scores[i] < 15:
> >         break
> >     print(f'Base {i} = {sequence[i]} with {quality_to_percent(quality_scores[i])}% accuracy')
> > ```
> > But we can make this a bit prettier using the `zip()` function:
> > ```
> > for base, score in zip(sequence, quality_scores):
> >     if score < 15:
> >         break
> >     print(f'Base = {base} with {quality_to_percent(score)}% accuracy')
> > ```
> > But note that we don't have the position in the list anymore, so we remove it from the print statement.
> {: .solution}
{: .question}

```python
# We've got a Read
read = """
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
55CCF>>>>>>CCCCCCC65!''*((((***+))%%%++)(%%%%).1***-+*''))**
""".strip().split('\n')

def quality_to_percent(q):
    return 100 * (1 - (10 ** (q / -10)))

# Extract the sequence
sequence = read[1]
# And the quality scores, and map those to the correct values.
quality_scores = [ord(x) - 33 for x in read[3]]

# Write something here
# That loops over BOTH the sequence and Quality Scores.
# And prints them out
# If the quality scores are `<15`, then break and quit printing.
for i in ...
```
