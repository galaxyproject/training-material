---
layout: tutorial_hands_on

title: Python - Math
level: Introductory
requirements: []
follow_up_training: []
questions:
- How do I do math in Python?

objectives:
- Understand the fundamentals of object assignment and math in python and can write simple statements and execute calcualtions in order to be able to summarize the results of calculations and classify valid and invalid statements.
- Translate some known math functions (e.g. euclidean distance, root algorithm) into python to transfer concepts from mathematics lessons directly into Python.

time_estimation:  30M
key_points:
- Converting mathematics equation to python is pretty easy!
- In real life you'll occasionally need to do this, either converting from a formula you read in a paper, or a description of an algorithm, into code that you'll re-use.
- Did you forget how a particular module or function works? Try `help()`

subtopic: python-modular
contributors:
  - hexylena
  - dirowa
  - bazante1

priority: 1
notebook:
  language: python
---

Here we'll learn some of the fundamentals of python and how to do basic maths in Python.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Python Fundamentals

## Variables

Any Python interpreter can be used as a calculator. We can do simple sums

```python
3 + 5 - 4
```

Multiplication

```python
5 * 4
```

Division

```python
5 / 4
```

Just like you've probably learned in maths courses, we can assign values to variables

```python
x = 10
y = 2.5
z = 3
```

From now on, whenever we use `x` or `y` or `z`, Python will substitute the value we assigned to it.
You'll notice when we assign things to a variable, nothing is printed out. If we've assigned a value to a variable we can print it by doing:

```python
print(x)
print(y, z)
```

You can print out multiple variables if you separate them with a comma `,`.


> <tip-title>Variable Names</tip-title>
> In Python, variable names:
>
> - can include letters, digits, and underscores
> - cannot start with a digit
> - are case sensitive
>
> This means that, for example:
>
> - `x` is a valid name
> - `weight0` is a valid variable name, whereas `0weight` is not
> - `weight` and `Weight` are different variables
{: .tip}

```python
x * (y + z)
```

What do you think the following will result in? The best way to find out is trying things for yourself. If it's wrong, can you correct it?

```python
xy
```

Time to check what we've learned!

> <question-title>Exercise 0: Simple Equations</question-title>
> Given the equation:
>
> $$y = x * 921 + 534$$
>
> What is the value of `y` when `x = 452`
>
> You can test out solutions in the box below this question.
>
> > <solution-title></solution-title>
> > ```
> > x = 452
> > y = x * 921 + 534
> > # Remember to print it out!
> > print(y)
> > ```
> >
> > 416826
> {: .solution}
{: .question}

```python
# Test solutions here!
# By the way, lines starting with a # are comment lines!
# You can use that to take notes, and not affect how your code runes
# People use it for documentation: explaining what the function does,
# what a variable means, why they chose this or that algorithm.
```

## Libraries

A library is a collection of files (called modules) that contains functions for use by other programs. It may also contain data values (e.g., numerical constants) and other things. A library's contents are supposed to be related, but there's no way to enforce that. The Python standard library is an extensive suite of modules that comes with Python itself. Many additional libraries are available from [PyPI (the Python Package Index)](https://pypi.org/).

### Libraries and modules

A library is a collection of modules, but the terms are often used interchangeably, especially since many libraries only consist of a single module, so don't worry if you mix them.

### A program must import a library module before using it.

You can use `import` to load a library module into a program's memory, then refer to things from the `module as module_name.thing_name`. Python uses `.` to mean “part of”. For example, using `datetime`, one of the modules in the standard library:

```python
import datetime

# it is currently
datetime.datetime.now()
```


> <tip-title>More complicated importing</tip-title>
>
> First, you can use `help` to learn about the contents of a library module. it works just like help for a function.
>
> ```
> help(datetime)
> ```
>
> You can import specific items from a library module to shorten programs. You can use `from ... import ...` to load only specific items from a library module. Then refer to them directly without library name as prefix.
>
>
> ```
> from datetime import datetime
>
> datetime.now()
> ```
>
> You can create an alias for a library module when importing it to shorten programs. Use `import ... as ...` to give a library a short alias while importing it. Then refer to items in the library using that shortened name.
>
> ```
> from datetime import datetime as dt
>
> dt.now()
> ```
{: .tip}


## Math Module

Let's import the math module:

```python
import math
```

This imports the math module. If you ever don't know what a function does, you can use the `help()` command:

```python
help(math)
```

And you'll see a list of the functions and properties available. Let's try out one of those functions now:

```python
math.sqrt(9)
```

> <solution-title>Tip: Why `math.`?</solution-title>
> When we import a module like `import math`, we need to use that as a prefix. Imagine we had two different modules, `math` and `other_math`, and both have a `sqrt` function. How would Python know which `sqrt` function we wanted? So we use `math.sqrt` to be explicit about which function we need.
{: .solution}

You might also have done powers (e.g. 2 cubed, or $$2^3$$) in the past, too:

```python
math.pow(2, 3)
```

Above we computed 2 cubed, but how did we know what to write? We might have seen `math.pow` in the `help(math)` page above, but how did we know `2, 3`? The same way of course:

```python
help(math.pow)
```

That would tell us that if we want 2 cubed, we need to write `2` and `3`. So let's do some exercises now and practice some of our maths skills.


> <question-title>Exercise 1: Basics</question-title>
> Please convert this function from an equation, into python code:
>
> $$x = 2^8$$
>
> > <solution-title></solution-title>
> > ```
> > x = math.pow(2, 8)
> > print(x)
> > ```
> {: .solution}
{: .question}

```python
# Test solutions here!
x =
print(x)
```


> <question-title>Exercise 2: Averaging two values</question-title>
> See if you can find the average of these two values, using the math operations
>
> - 23484
> - 12345
>
> You can find the average by summing, and dividing by 2:
>
> > <solution-title></solution-title>
> > ```
> > (23484 + 12345) / 2
> > ```
> {: .solution}
{: .question}

```python
# Test solutions here!
```

> <question-title>Exercise 3: Round-trip</question-title>
> Please convert this function from an equation, into python code. Remember, you can assign values to variables, if you want to split this into multiple steps.
>
> $$x = \sqrt{9^2}$$
>
> > <solution-title></solution-title>
> > ```
> > y = math.pow(9, 2)
> > x = math.sqrt(y)
> > ```
> >
> > Or
> >
> > ```
> > x = math.sqrt(math.pow(9, 2))
> > ```
> {: .solution}
{: .question}

```python
# Test solutions here!
```

> <question-title>Exercise 4: Pythagorean Theorem</question-title>
> The formula for a 90° triangle can be expressed as:
>
> $$a^2 + b^2 = c^2$$
>
> Or, expressed in terms of c, $$c = \sqrt{a^2 + b^2}$$
>
> Please convert this to python and find `c` when
>
> ```
> a = 65
> b = 72
> ```
>
> > <solution-title></solution-title>
> > c = 97
> > ```
> > c = math.sqrt(math.pow(a, 2) + math.pow(b, 2))
> > ```
> {: .solution}
{: .question}

```python
# Test solutions here!
```

> <question-title>Exercise 5: Quadratic Roots</question-title>
> Way back in algebra class, you might have been given a quadratic equation, something like:
>
> $$y = 2*x^2 + x - 1$$ and were told to find the roots of this function, using a complicated equation. So challenge time: reproduce this equation in Python:
>
> Given the following variables:
>
> a = 2
> b = 1
> c = -1
>
> Convert the following formulas to Python, and find the answers
>
> $$\dfrac{-b + \sqrt{b^2 - 4ac}}{2a}$$ and $$\dfrac{-b - \sqrt{b^2 - 4ac}}{2a}$$
>
> Since we're checking your ability to write Python, not do math, it should give `-1` and `0.5` to let you check your work. If you got those values, you got it right!
>
> Make sure you save each root as it's own variable, and then print them out.
>
> > <solution-title></solution-title>
> > ```
> > root1 = (-b + math.sqrt(math.pow(b, 2) - 4 * a * c))/(2 * a)
> > root2 = (-b - math.sqrt(math.pow(b, 2) - 4 * a * c))/(2 * a)
> > print(root1)
> > print(root2)
> > ```
> {: .solution}
{: .question}

```python
# Test solutions here!
root1 =
root2 =
print(root1)
print(root2)
```
