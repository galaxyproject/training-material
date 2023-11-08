---
layout: tutorial_hands_on

title: Python - Flow Control
level: Introductory
requirements: []
follow_up_training: []
questions:
- "How can my programs do different things based on data values?"


objectives:
- "Write conditional statements including `if`, `elif`, and `else` branches."
- "Correctly evaluate expressions containing `and` and `or`."

time_estimation:  40M
key_points:
- "Use `if condition` to start a conditional statement, `elif condition` to
   provide additional tests, and `else` to provide a default."
- "The bodies of the branches of conditional statements must be indented."
- "Use `==` to test for equality."
- "`X and Y` is only true if both `X` and `Y` are true."
- "`X or Y` is true if either `X` or `Y`, or both, are true."
- "Zero, the empty string, and the empty list are considered false;
   all other numbers, strings, and lists are considered true."
- "`True` and `False` represent truth values."
- "`not` can be used to invert the condition"

subtopic: python-modular
contributors:
  - carpentries
  - hexylena
  - dirowa
  - bazante1

priority: 5
notebook:
  language: python
  pyolite: true
---

"Flow Control" is how we describe when we change the flow of code's execution, based on some conditions. Here we'll learn how to take different actions depending on what data out program sees, or how to run code only if some condition is true.

> <comment-title></comment-title>
> This tutorial is **significantly** based on [the Carpentries](https://carpentries.org) [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation/) and [Plotting and Programming in Python](https://swcarpentry.github.io/python-novice-gapminder/), which are licensed CC-BY 4.0.
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

## Comparators

In Python we have the following comparators to do compare two values

- `>`: greater than
- `<`: less than
- `==`: equal to
- `!=`: does not equal
- `>=`: greater than or equal to
- `<=`: less than or equal to

They're all "binary" comparators, we can only compare two values at a time.

```python
print(37 < 38)
print(38 < 38)
print(39 < 38)
```

These print out `True` or `False`, these are the two possible values of the [*boolean*](https://en.wikipedia.org/wiki/Boolean_algebra) datatype in Python.

We can use `<=` to check if it's less than or equal to:

```python
print(19 <= 20)
print(20 <= 20)
print(21 <= 20)
```

And we can use `==` for comparing numbers in Python

```
print(11 == 11)
print(11 != 11)
print(22 != 33)
```

And now that we can compare numbers, we can start doing useful things with them!

## Conditionals

We can ask Python to take different actions, depending on a condition, with an `if` statement:

```python
num = 37
if num > 100:
    print('greater')
else:
    print('not greater')
print('done')
```

The second line of this code uses the keyword `if` to tell Python that we want to make a choice.
If the test that follows the `if` statement is true,
the body of the `if`
(i.e., the set of lines indented underneath it) is executed, and "greater" is printed.
If the test is false,
the body of the `else` is executed instead, and "not greater" is printed.
Only one or the other is ever executed before continuing on with program execution to print "done":

![A flowchart diagram of the if-else construct that tests if variable num is greater than 100](../../images/python-flowchart-conditional.png)

Conditional statements don't have to include an `else`. If there isn't one,
Python simply does nothing if the test is false:

```python
num = 53
print('before conditional...')
if num > 100:
    print(f'{num} is greater than 100')
print('...after conditional')
```

> <question-title>If behaviour</question-title>
>
> Try changing the `num` value and see what happens for different values.
>
> What happens if `num` is a:
> 1. 202
> 2. 3.145
> 3. "test"
> 4. 100.000001
>
> > <solution-title></solution-title>
> > 1. Condition is activated!
> > 2. Nothing, but not because it is a float! Because it's less than 100
> > 3. Traceback, a `TypeError`, you cannot compare strings with integers
> > 4. Condition is activated!
> {: .solution}
{: .question}

## Multiple Branches

But what if you want more branches? What if you need to handle more cases? `elif` to the rescue!

We can chain several tests together using `elif`, which is short for "else if".

```
if todays_temperature > 30:
    print("Wear shorts! Remember your sunscreen")
elif todays_temperature > 20:
    print("It's nice weather finally! Gasp!")
elif todays_temperature < 10:
    print("Time to bundle up!")
else:
    print("Dress normally")
```

> <tip-title>If/Elif/Elif/Elif/Else:</tip-title>
> if/elif/else cases follow these rules:
>
> - must start with an `if`
> - can have 0 or more `elif` conditions
> - can have 0 or 1 `else` condition (if no else condition is supplied, it's equivalent to `else: <nothing>`)
{: .tip}

Each of these three sections is a **branch**, the code pauses, and chooses to go down one of the branches based on the conditions.

The following Python code uses `elif` to print the sign of a number.

```python
num = -3

if num > 0:
    print(f'{num} is positive')
elif num == 0:
    print(f'{num} is zero')
else:
    print(f'{num} is negative')
```

**NB**: To test for equality we use a double equals sign `==`
rather than a single equals sign `=` which is used to assign values.

## Combining Tests

We can also combine tests using `and` and `or`.
`and` is only true if both parts are true:

```python
a = 1
b = -1
if (a > 0) and (b <= 0):
    print('both parts are true')
else:
    print('at least one part is false')
```

> <question-title>Predict what happens</question-title>
> Predict the outcomes of the following values of `a` and `b` above. Predicting what you think the code will do is a useful skill to practice
> 1. a = 0; b = -1
> 2. a = 0; b = 10
> 3. a = 4; b = -22
> 4. a = 99; b = 99
>
> > <solution-title></solution-title>
> > 1. at least one part is false
> > 2. at least one part is false
> > 3. both parts are true
> > 4. at least one part is false
> {: .solution}
{: .question}


while `or` is true if at least one part is true:

```python
a = 1
b = -1
if (a < 0) or (b > 0):
    print('at least one test is true')
```


> <tip-title>`True` and `False`</tip-title>
> `True` and `False` are special words in Python called `booleans`,
> which represent truth values. A statement such as `1 < 0` returns
> the value `False`, while `-1 < 0` returns the value `True`.
{: .tip}

`True` and `False` booleans are not the only values in Python that are true and false.
In fact, *any* value can be used in an `if` or `elif`. This is commonly used to
check, for instance, if a string is empty or if some data is provided:

```python
if '':
    print('empty string is true')
if 'word':
    print('word is true')
```

You can also use it to check if a list is empty or full:

```python
if []:
    print('empty list is true')
if [1, 2, 3]:
    print('non-empty list is true')
# The last statement is equivalent to:
if len([1, 2, 3]) > 0:
    print('non-empty list is true')
```

Or you can check if a number is zero, or non-zero:

```python
if 0:
    print('zero is true')
if 1:
    print('one is true')
```

## Inverting Conditions

Sometimes it is useful to check whether some condition is not true.
The Boolean operator `not` can do this explicitly.
After reading and running the code below,
write some `if` statements that use `not` to test the rule
that you formulated in the previous question.
`not` is a `unary` (not `binary`) operator: it only takes a single value

```python
if not '':
    print('empty string is not true')
if not 'word':
    print('word is not true')
if not not True:
    print('not not True is true')
```

## Ranges

Python makes it super easy to check if a number is within a range.

```python
quality_score = 32 # Try out different values!

if quality_score > 40:
    print("Your data is a bit sus")
elif 20 < quality_score <= 40:
    print("Hey that looks ok")
elif 4 < quality_score <= 20:
    print("Oh you did nanopore sequencing")
else:
    print("It shouldn't be *that* bad. Try again.")
```

There are two important points here:

- `20 < x < 40` is equivalent to `20 < x and x < 40`, checking both sides of the condition, to make sure it's greater than one value and smaller than another
- Note that we checked in the second case `20 < x` and then in the third we had to check `x <= 20`. If we had not had a `<=` on one side, what would have happened to 20? It would have gone straight to else!


## Exercises

> <question-title></question-title>
> `if`s, `elif`s and `else`s get evaluated in blocks. Look at the following code and list the lines that are part of a single block.
>
> ```
> 1.  if x:
> 2.      # ..
> 3.  if y:
> 4.      # ..
> 5.  elif z:
> 6.      # ..
> 7.  if q:
> 8.      # ..
> 9.  else:
> 10.     # ..
> 11. elif t:
> 12.     # ..
> 13. else e:
> 14.     # ..
> ```
>
> > <solution-title></solution-title>
> > "Blocks" of if/elif/elses
> > - must start with an `if`
> > - can have 0 or more `elif` conditions
> > - can have 0 or 1 `else` condition (if no else condition is supplied, it's equivalent to `else: <nothing>`)
> >
> > The above blocks are parsed together, you could not insert a `print` anywhere within the blocks, but between the blocks it would work.
> >
> > - 1-2, Just an `if` by itself. There's no elif, or else, so that's the end of that block
> > - 3-6, `if` and `elif` get evaluated, there is no `else`, so that's the end of that block
> > - 7-10, `if` and `else` is fine
> > - 11-14, error! This is missing an `if` case, it will fail with a syntaxerror.
> {: .solution}
{: .question}

```python
# Test code here.
```

> <question-title>How Many Paths?</question-title>
>
> Consider this code:
>
> ```
> if 4 > 5:
>     print('A')
> elif 4 == 5:
>     print('B')
> elif 4 < 5:
>     print('C')
> ```
>
> Which of the following would be printed if you were to run this code?
> Why did you pick this answer?
>
> 1.  A
> 2.  B
> 3.  C
> 4.  B and C
>
> > <solution-title></solution-title>
> > C gets printed because the first two conditions, `4 > 5` and `4 == 5`, are not true,
> > but `4 < 5` is true.
> {: .solution}
{: .question}

```python
# Test code here.
```

> <question-title>Close Enough</question-title>
>
> Write some conditions that print `True` if the variable `a` is within `10` of the variable `b`
> and `False` otherwise.
> Compare your implementation with your partner's:
> do you get the same answer for all possible pairs of numbers?
>
> > <tip-title>abs</tip-title>
> > There is a [built-in function `abs`][abs-function] that returns the absolute value of
> > a number:
> > ```python
> > print(abs(-12))
> > ```
> >
> > ```
> > 12
> > ```
> {: .tip}
>
> > <solution-title>1</solution-title>
> > ```
> > a = 5
> > b = 5.1
> >
> > if abs(a - b) <= 10:
> >     print('True')
> > else:
> >     print('False')
> > ```
> {: .solution}
>
> > <solution-title>2</solution-title>
> > ```
> > print(abs(a - b) <= 10)
> > ```
> >
> > This works because the Booleans `True` and `False`
> > have string representations which can be printed.
> {: .solution}
{: .question}

```python
# Test code here.
```

> <question-title>Pitfalls</question-title>
>
> A *integer* number between 0 and 100 will be provided to this function. Answer these two questions:
>
> - Will it always print something? If not, which value(s) fail?
> - Can you find any numbers the programmer explicitly wanted to handle, that aren't handled as expected?
>
> ```
> num = 42 # Randomly chosen so the code will execute, try changing it around.
> if num > 90:
>     print("great score")
> elif num < 32:
>     print("Very cold")
> elif num >= 86:
>     print("Almost")
> elif num == 86:
>     print("It's exactly this value!")
> elif 32 < num < 58:
>     print("Getting warmer")
> elif 57 < num <= 86:
>     print("Everything else goes here")
> ```
>
> > <solution-title></solution-title>
> > 1. No, it won't. 32 is the only value there that doesn't print anything. You can either do `x < 57` and later `57 <= x` to test the bigger and smaller values, or you can make use `x < 57` and `56 < x`, which have the same results, but **only with integers**. If your code accepted a float, e.g. `56.5`, both of those tests would be true. So `x < 57` and later `57 <= x` is the preferred way to write that.
> > 2. `86` is the most obvious solution to this, the programmer added a check specifically to see if the value was 86, but instead it's caught by the previous case.
> {: .solution}
{: .question}

```python
num = 42 # Randomly chosen so the code will execute, try changing it around.
if num > 90:
    print("great score")
elif num < 32:
    print("Very cold")
elif num >= 86:
    print("Almost")
elif num == 86:
    print("It's exactly this value!")
elif 32 < num < 58:
    print("Getting warmer")
elif 57 < num <= 86:
    print("Everything else goes here")
```

> <tip-title>Why a synthetic example like this?</tip-title>
> Complicated if/elif/else cases are common in code, you need to be able to spot these sort of issues. For example there are large if/else cases in the [Galaxy codebase](https://github.com/galaxyproject/galaxy/blob/9143dd7ca46d150ebfb26febbe187979f682da51/tools/stats/grouping.py#L153-L176), sometimes nested even, and being ale to predict their behaviour is really important to being able to work with the code. [Missing else cases](https://github.com/galaxyproject/galaxy/blob/9143dd7ca46d150ebfb26febbe187979f682da51/tools/stats/grouping.py#L185-L195) are sometimes important, sometimes a bug, sometimes just the code hasn't been implemented yet, which is why we always write good code comments!
{: .tip}
