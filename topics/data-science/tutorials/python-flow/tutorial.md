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

enable: false
subtopic: python-modular
contributors:
  - carpentries
  - hexylena

priority: 1
notebook:
  language: python
---

"Flow Control" is how we describe when we change the flow of code's execution, based on some conditions. Here we'll learn how to take different actions depending on what data out program sees, or how to run code only if some condition is true.

> ### {% icon comment %} Comment
> This tutorial is **significantly** based on [the Carpentries](https://carpentries.org) [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation/) and [Plotting and Programming in Python](https://swcarpentry.github.io/python-novice-gapminder/), which are licensed CC-BY 4.0.
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

```python
not greater
done
```
{: .output}

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
    print(num, 'is greater than 100')
print('...after conditional')
```


> ### {% icon question %} Question: If behaviour
>
> Try changing the `num` value and see what happens for different values.
>
> What happens if `num` is a:
> 1. 202
> 2. 3.145
> 3. "test"
> 4. 100.000001
>
> > ### {% icon solution %} Solution
> > 1. Condition is activated!
> > 2. Nothing, but not because it is a float! Because it's less than 100
> > 3. Traceback, a `TypeError`, you cannot compare strings with integers
> > 4. Condition is activated!
> {: .solution}
{: .question}

We can also chain several tests together using `elif`,
which is short for "else if".
The following Python code uses `elif` to print the sign of a number.

```python
num = -3

if num > 0:
    print(num, 'is positive')
elif num == 0:
    print(num, 'is zero')
else:
    print(num, 'is negative')
```


Note that to test for equality we use a double equals sign `==`
rather than a single equals sign `=` which is used to assign values.

> ### {% icon tip %} Tip: Comparing in Python
>
> Along with the `>` and `==` operators we have already used for comparing values in our
> conditionals, there are a few more options to know about:
>
> - `>`: greater than
> - `<`: less than
> - `==`: equal to
> - `!=`: does not equal
> - `>=`: greater than or equal to
> - `<=`: less than or equal to
{: .tip}

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

> ### {% icon question %} Question: Predict what happens
> Predict the outcomes of the following values of `a` and `b` above. Predicting what you think the code will do is a useful skill to practice
> 1. a = 0; b = -1
> 2. a = 0; b = 10
> 3. a = 4; b = -22
> 4. a = 99; b = 99
>
> > ### {% icon solution %} Solution
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


> ### {% icon tip %} Tip: `True` and `False`
> `True` and `False` are special words in Python called `booleans`,
> which represent truth values. A statement such as `1 < 0` returns
> the value `False`, while `-1 < 0` returns the value `True`.
{: .tip}

`True` and `False` booleans are not the only values in Python that are true and false.
In fact, *any* value can be used in an `if` or `elif`.
After reading and running the code below,
explain what the rule is for which values are considered true and which are considered false.

```python
if '':
    print('empty string is true')
if 'word':
    print('word is true')
if []:
    print('empty list is true')
if [1, 2, 3]:
    print('non-empty list is true')
if 0:
    print('zero is true')
if 1:
    print('one is true')
```

> ## How Many Paths?
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
> > ## Solution
> > C gets printed because the first two conditions, `4 > 5` and `4 == 5`, are not true,
> > but `4 < 5` is true.
> {: .solution}
{: .question}

> ### {% icon tip %} Tip: That's Not Not What I Meant
>
> Sometimes it is useful to check whether some condition is not true.
> The Boolean operator `not` can do this explicitly.
> After reading and running the code below,
> write some `if` statements that use `not` to test the rule
> that you formulated in the previous question.
>
> ```
> if not '':
>     print('empty string is not true')
> if not 'word':
>     print('word is not true')
> if not not True:
>     print('not not True is true')
> ```
{: .tip}

> ### {% icon question %} Question: Close Enough
>
> Write some conditions that print `True` if the variable `a` is within `10` of the variable `b`
> and `False` otherwise.
> Compare your implementation with your partner's:
> do you get the same answer for all possible pairs of numbers?
>
> > ### {% icon tip %} Tip: abs
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
> > ### {% icon solution %} Solution 1
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
> > ### {% icon solution %} Solution 2
> > ```
> > print(abs(a - b) <= 10)
> > ```
> >
> > This works because the Booleans `True` and `False`
> > have string representations which can be printed.
> {: .solution}
{: .question}
