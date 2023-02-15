---
layout: tutorial_hands_on

title: Python - Basic Types & Type Conversion
level: Introductory
requirements: []
follow_up_training: []
questions:
- "What kinds of data do programs store?"
- "How can I convert one type to another?"

objectives:
- "Explain key differences between integers and floating point numbers."
- "Explain key differences between numbers and character strings."
- "Use built-in functions to convert between integers, floating point numbers, and strings."

time_estimation:  30M
key_points:
- "Every value has a type."
- "Use the built-in function `type` to find the type of a value."
- "Types control what operations can be done on values."
- "Strings can be added and multiplied."
- "Strings have a length (but numbers don't)."
- "Must convert numbers to strings or vice versa when operating on them."
- "Can mix integers and floats freely in operations."
- "Variables only change value when something is assigned to them."

subtopic: python-modular
contributors:
  - carpentries
  - hexylena
  - dirowa
  - bazante1

priority: 3
notebook:
  language: python
  pyolite: true
---

Python is a typed language, data has a type, and different types of data cannot always be connected immediately and might need some conversion step before they can be used together. For instance if you add a number to a number, what should happen? If you add a number to a message, what do you expect will happen?

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Types

Every value in a program has a specific type.

Name                  | Python Code | Represents
---                   | ---         | ---
Integer               | `int`       | represents positive or negative whole numbers like 3 or -512.
Floating point number | `float`     | represents real numbers like 3.14159 or -2.5.
Character string      | `str`       | text, written with either `'` or `"` quotes (they must match)

### Checking the Type

Use the built-in function `type` to find out what type a value has. This works on values as well as variables. But remember: the *value* has the type --- the *variable* is just a label.

Check the type of values with the [`type()`](https://docs.python.org/3.8/library/functions.html#type) function:

```python
print(type(52))
print(type(3.14159))
```

You can also check the types of variables

```python
fitness = 'average'
print(type(fitness))
```

### Methods

A value's type determines what the program can do to it. Some operations may work

```python
print(5 - 3)
```

And some operations may not work:

```python
print('hello' - 'h')
```

For instance, you can use the `+` and `*` operators on strings.

```python
full_name = 'Ahmed' + ' ' + 'Walsh'
print(full_name)
separator = '=' * 10
print(separator)
```

Some methods only accept specific types, or only work on specific types.

The built-in function `len` returns the length of your data. Which of the following would you expect to work? `len(string)`? `len(int)`?

```python
print(len(full_name))
print(len(52))
```

### Matching Types

Not all types support all operations, adding an integer to a string doesn't make much sense:

```python
print(1 + '2')
```

This does not work because it's ambiguous: should `1 + '2'` be `3` (a number) or `'12'` (a string)? Some types can be converted to other types by using the type name as a function.

```python
print(1 + int('2'))
print(str(1) + '2')
```

### Operation Support

Here is a quick chart showing which operations are allowed for each pair:

Left\Right | `int`  | `float` | `str`
---------- | ---    | -----   | ---
`int`      | `+-*/` | `+-*/`  | `*`
`float`    | `+-*/` | `+-*/`  | ``
`str`      | `*`    | ``      | `+`

As you can see you can do `3 * "test"` and `"test" * 3`, but it doesn't work with floats.

## Can mix integers and floats freely in operations.

Integers and floating-point numbers can be mixed in arithmetic. Python 3 (which we use) automatically converts integers to floats as needed.

```python
print(f'half is {1 / 2.0}')
print(f'three squared is {3.0 ** 2}')
```


## Variables only change value when something is assigned to them.

If we make one cell in a spreadsheet depend on another, and update the latter,
the former updates automatically. However, this does **not** happen in programming languages.

```python
variable_one = 1
variable_two = 5 * variable_one
variable_one = 2
print(f'first is {variable_one} and second is {variable_two}')
```

The computer reads the value of `first` when doing the multiplication, creates
a new value, and assigns it to `second`. After that, `second` does not remember
where it came from. Every computation happens line-by-line.

> <question-title>Fractions</question-title>
>
> What type of value is 3.14159?
> How can you find out?
>
> > <solution-title></solution-title>
> >
> > It is a floating-point number (often abbreviated "float").
> > It is possible to find out by using the built-in function `type()`.
> >
> > ```
> > print(type(3.14159))
> > <class 'float'>
> > ```
> {: .solution}
{: .question}

```python
# Test out solutions here!

```

> <question-title>Automatic Type Conversion</question-title>
>
> What type of value is the result of (3.25 + 4)?
>
> > <solution-title></solution-title>
> >
> > It is a float:
> > integers are automatically converted to floats as necessary.
> >
> > ```
> > result = 3.25 + 4
> > print(f'result is {type(result)}')
> > ```
> >
> > ```
> > 7.25 is <class 'float'>
> > ```
> {: .solution}
{: .question}

```python
# Test out solutions here!

```

> <question-title>Choose a Type</question-title>
>
> What type of value (integer, floating point number, or character string)
> would you use to represent each of the following?  Try to come up with more than one good answer for each problem.  For example, in  # 1, when would counting days with a floating point variable make more sense than using an integer?
>
> 1. Number of days since the start of the year.
> 2. Time elapsed from the start of the year until now in days.
> 3. Serial number of a piece of lab equipment.
> 4. A lab specimen's age
> 5. Current population of a city.
> 6. Average population of a city over time.
>
> > <solution-title></solution-title>
> >
> > The answers to the questions are:
> > 1. Integer, since the number of days would lie between 1 and 365.
> > 2. Floating point, since fractional days are required
> > 3. Character string if serial number contains letters and numbers, otherwise integer if the serial number consists only of numerals
> > 4. This will vary! How do you define a specimen's age? whole days since collection (integer)? date and time (string)?
> > 5. Choose floating point to represent population as large aggregates (eg millions), or integer to represent population in units of individuals.
> > 6. Floating point number, since an average is likely to have a fractional part.
> {: .solution}
{: .question}

> <question-title>Division Types</question-title>
>
> In Python 3, the `//` operator performs integer (whole-number) floor division, the `/` operator performs floating-point
> division, and the `%` (or *modulo*) operator calculates and returns the remainder from integer division:
>
> ```python
> print(f'5 // 3: {5 // 3}')
> print(f'5 / 3: {5 / 3}')
> print(f'5 % 3: {5 % 3}')
> ```
>
> ```
> 5 // 3: 1
> 5 / 3: 1.6666666666666667
> 5 % 3: 2
> ```
>
> If `num_subjects` is the number of subjects taking part in a study,
> and `num_per_survey` is the number that can take part in a single survey,
> write an expression that calculates the number of surveys needed
> to reach everyone once.
>
> > <solution-title></solution-title>
> > We want the minimum number of surveys that reaches everyone once, which is
> > the rounded up value of `num_subjects/ num_per_survey`. This is
> > equivalent to performing a floor division with `//` and adding 1. Before
> > the division we need to subtract 1 from the number of subjects to deal with
> > the case where `num_subjects` is evenly divisible by `num_per_survey`.
> > ```
> > num_subjects = 600
> > num_per_survey = 42
> > num_surveys = (num_subjects - 1) // num_per_survey + 1
> >
> > print(num_subjects, 'subjects,', num_per_survey, 'per survey:', num_surveys)
> > ```
> >
> > ```
> > 600 subjects, 42 per survey: 15
> > ```
> {: .solution}
{: .question}

```python
# Test out solutions here!
```

> <question-title>Strings to Numbers</question-title>
>
> Where reasonable, `float()` will convert a string to a floating point number,
> and `int()` will convert a floating point number to an integer:
>
> ```
> print("string to float:", float("3.4"))
> print("float to int:", int(3.4))
> ```
>
> ```
> string to float: 3.4
> float to int: 3
> ```
>
> If the conversion doesn't make sense, however, an error message will occur.
>
> > <code-in-title>Python</code-in-title>
> > ```
> > print("string to float:", float("Hello world!"))
> > ```
> {: .code-in}
>
> > <code-out-title></code-out-title>
> > ```
> > Traceback (most recent call last):
> >   File "<stdin>", line 1, in <module>
> > ValueError: could not convert string to float: 'Hello world!'
> > ```
> {: .code-out}
>
>
> Given this information, what do you expect the following program to do?
>
> What does it actually do?
>
> Why do you think it does that?
>
> ```
> print("fractional string to int:", int("3.4"))
> ```
>
> > <solution-title></solution-title>
> > What do you expect this program to do? It would not be so unreasonable to expect the Python 3 `int` command to
> > convert the string "3.4" to 3.4 and an additional type conversion to 3. After all, Python 3 performs a lot of other
> > magic - isn't that part of its charm?
> >
> > ```
> > int("3.4")
> > ```
> > However, Python 3 throws an error.
> > ```
> > Traceback (most recent call last):
> >   File "<stdin>", line 1, in <module>
> > ValueError: invalid literal for int() with base 10: '3.4'
> > ```
> >
> > Why? To be consistent, possibly. If you ask Python to perform two consecutive
> > typecasts, you must convert it explicitly in code.
> >
> > ```
> > int(float("3.4"))
> > ```
> >
> > ```
> > 3
> > ```
> {: .solution}
{: .question}

```python
# Test out solutions here!
```

> <question-title>Arithmetic with Different Types</question-title>
>
> Which of the following will return the floating point number `2.0`?
> Note: there may be more than one right answer.
>
> ```
> first = 1.0
> second = "1"
> third = "1.1"
> ```
>
> 1. `first + float(second)`
> 2. `float(second) + float(third)`
> 3. `first + int(third)`
> 4. `first + int(float(third))`
> 5. `int(first) + int(float(third))`
> 6. `2.0 * second`
>
> > <solution-title></solution-title>
> >
> > Answer: 1 and 4 give exactly 2.0.
> > Answer 5 gives the value `2` which may be considered equivalent, but is not returning a float specifically.
> {: .solution}
{: .question}

```python
# Test out solutions here!
```
