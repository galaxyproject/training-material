---
layout: tutorial_hands_on

title: Python - Lists & Strings & Dictionaries
level: Introductory
requirements: []
follow_up_training: []
questions:
- "How can I store multiple values?"

objectives:
- "Explain why programs need collections of values."
- "Write programs that create flat lists, index them, slice them, and modify them through assignment and method calls."

time_estimation:  1H
key_points:
- "A list stores many values in a single structure."
- "Use an item's index to fetch it from a list."
- "Lists' values can be replaced by assigning to them."
- "Appending items to a list lengthens it."
- "Use `del` to remove items from a list entirely."
- "The empty list contains no values."
- "Lists may contain values of different types."
- "Character strings can be indexed like lists."
- "Character strings are immutable."
- "Indexing beyond the end of the collection is an error."

enable: false
subtopic: python
contributors:
  - carpentries
  - hexylena
  - dirowa

priority: 1
notebook:
  language: python
---



> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Lists

Doing calculations with a hundred variables called `pressure_001`, `pressure_002`, etc. would be at least as slow as doing them by hand. Using a *list* to store many values together solves that problems. Lists are surrounded by square brackets: `[`, `]`, with values separated by commas:

```python
pressures = [0.273, 0.275, 0.277, 0.275, 0.276]
print('pressures:', pressures)
print('length:', len(pressures))
```

## Indexing

You can use an item's index to fetch it from a list.

```python
print('zeroth item of pressures:', pressures[0])
print('fourth item of pressures:', pressures[4])
```

## Replacement

Lists' values can be changed or replaced by assigning a new value to the position in the list.

```
pressures[0] = 0.265
print('pressures is now:', pressures)
```

Note how the first item has changed from `0.273`

## Appending

Appending items to a list lengthens it. You can do `list_name.append()` to add items to the end of a list.

```python
primes = [2, 3, 5]
print('primes is initially:', primes)
primes.append(7)
print('primes has become:', primes)
```


`.append()` is a *method* of lists. It's like a function, but tied to a particular object. You use `object_name.method_name` to call methods, which deliberately resembles the way we refer to things in a library.

We will meet other methods of lists as we go along. Use `help(list)` for a preview. `extend` is similar to `append`, but it allows you to combine two lists.  For example:

```python
teen_primes = [11, 13, 17, 19]
middle_aged_primes = [37, 41, 43, 47]
print('primes is currently:', primes)
primes.extend(teen_primes)
print('primes has now become:', primes)
primes.append(middle_aged_primes)
print('primes has finally become:', primes)
```

Note that while `extend` maintains the "flat" structure of the list, appending a list to a list makes the result two-dimensional - the last element in `primes` is a list, not an integer.

This starts to become a more complicated data structure, and we'll use more of these later. A list containing both integers and a list can be called a "hetereogenous" list, since it has multiple different data types. This is relatively uncommon, most of the lists you'll encounter will have a single data type inside of them. Sometimes you'll see a list of lists, which can be used to store positions, like a chessboard.

## List Indices

In computer science and programming we number the positions within a list starting from `0`, rather than from `1`.

```python
# Position  0         1          2            3           4
weekdays = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday']
print(weekdays[0])
print(weekdays[4])
print(weekdays[3])
```

But if you try an access a position that is outside of the list, you'll get an error

```python
print(weekdays[9])
```

returns a `IndexError: list index out of range`.

> ### {% icon tip %} Tip: Reading Error Messages
> So how do you read this?
> ```
> 1 | ---------------------------------------------------------------------------
> 2 | IndexError                                Traceback (most recent call last)
> 3 | /tmp/ipykernel_648319/137030145.py in <module>
> 4 | ----> 1 print(weekdays[9])
> 5 |
> 6 | IndexError: list index out of range
> ```
>
> 1. This is just a line of `-`s as a separator
> 2. `IndexError`, here Jupyter/CoCalc/etc are trying to be helpful and highlight the error for us. This is the important bit of information!
> 3. This is the path to where the code is, Jupyter/CoCalc/etc create temporary files to execute your code.
> 4. Here an arrow points to the line number where something has broken. 1 shows that it's the first line within the cell, and it points to the print statement. Really it's pointing at the `weekdays[9]` within the print statement.
> 5. Blank
> 6. This is where we normally look for the **most important part of the Traceback**. The error message. An `IndexError`, namely that the list index (9) is out of the range of possible values (the length of the list.)
{: .question}

However, sometimes you want to access the very end of a list! You can either start at the beginning and count along to find the last item or second to last item, or you can use Negative Indices

```python
# Position  0         1          2            3           4
weekdays = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday']
# Position  -5        -4         -3           -2          -1

print(weekdays[-1])
print(weekdays[-2])
print(weekdays[-4])
```

## Removing Items.

You can use `del` to remove items from a list entirely. We use `del list_name[index]` to remove an element from a list (in the example, 9 is not a prime number) and thus shorten it. `del` is not a function or a method, but a statement in the language.

```python
primes = [2, 3, 5, 7, 9]
print('primes before removing last item:', primes)
del primes[4]
print('primes after removing last item:', primes)
```

## Empty Lists

The empty list contains no values. When you want to make a new list, use `[]` on its own to represent a list that doesn't contain any values. This is helpful as a starting point for collecting values, which we'll see soon.

## Heterogeneous Lists

Lists may contain values of different types. A single list may contain numbers, strings, and anything else.

```
goals = [1, 'Create lists.', 2, 'Extract items from lists.', 3, 'Modify lists.']
```

# Strings are like Lists

Text is often called a "string" in the programming world. Strings of text like `name = "Helena"` or `patient_id = "19237zud830"` are very similar conceptually to lists. Except instead of being a list of numbers, they're a lists of characters.

In a number of older programming languages, strings are indeed arrays of numbers internally. However python hides a lot of that complexity from us, so we can just work with text.

Still, many of the operations you use on lists, can also be used on strings as well! Strings can be indexed like lists so you can get single elements from lists.

```python
element = 'carbon'
print('zeroth character:', element[0])
print('third character:', element[3])
```

## Strings are immutable.

Strings, however, cannot be modified, you can't change a single letter in a string. Things that cannot be modified after creation are called *immutable* or sometimes *frozen*, compared to things which can be modified which are called *mutable*.
Python considers the string to be a single value with parts, not a collection of values.

```python
element[0] = 'C'
```

Lists and strings are both *collections*.

## Bounds

You cannot access values beyond the end of the list, this will result in an error. Python reports an `IndexError` if we attempt to access a value that doesn't exist. This is a kind of **runtime error**, as it cannot be detected as the code is parsed. Imagine if you had a script which let you read in a file, depending on how many lines were in the file, whether index 90 was valid or invalid, would depend on how big your file was.

```python
print('99th element of element is:', element[99])
```

# Tuples

Tuples are like lists, except after you create them, they're immutable. You can create a tuple with `(1, 2, 3, ...)` instead of `[1, 2, 3]`

```python
t = (1, 2, 3)
print(type(t))
# You can access indicies just like lists
print(t[0])
# But you can't modify them
# t.push(1) # This will fail. They're immutable!
```

# Sets

Sets are similar to lists, except:

- They're unique
- You cannot access individual indices within the set or slice them (next section)

You won't need sets very often, but once you get comfortable with python there will be times you might want to use them. And often you'll see code like:

```python
data = ['a', 'a', 'b', 'c', 'c', 'c']
unique_data = list(set(data))
```

in order to get the unique values in the list. Here we convert it to a set (unique), before converting it back to a list, so it behaves like we except. Sets are created with `{1, 2, 3, ...}` instead of `[1, 2, 3, ...]`

# Slicing & Dicing

All of the data types we've talked about today can be sliced, and this will be a key part of using lists.

```python
elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F']
# Instead of accessing a single element
print(elements[0])
# We'll access a range
print(elements[0:4])
```

Accessing only a portion of a list is commonly used, say if you have a list of FastQ files from paired end sequencing, perhaps you want two of them at a time. You could access those with `[0:2]`.

```python
# You don't need to start at 0
print(elements[6:8])
```

```python
# But your end should be bigger than your start.
# What do you think this will return?
# Make a guess before you run it
print(elements[6:5])
```

## The end of the list

But how do you access the end of the list? There are a couple options.

```python
elements[len(elements)-5:]
```

That is a really cumbersome way to access lists, so Python makes your life easier. If you need to access somewhere off of the end of the list, say, the last two items, you can just use negative values!

```python
elements[-5:]
```

Notice that we also haven't filled out the end value. If you don't supply an end value, Python will default to going to the end of the list. Likewise, if you don't provide a start value, Python will use `0` as the start by default, until whatever end value you subscribe.


> ### {% icon question %} Question: Valid and Invalid Slices
> Which of these do you think will be valid? Which are invalid? Predict what they will return:
> ```
> elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F']
> elements[0:3]
> elements[:3]
> elements[-3:3]
> elements[-8:-3]
> elements[:]
> elements[0:20]
> elements['H':'Li']
> elements[1.5:]
> ```
>
> > ### {% icon solution %} Solution
> > All of these are valid except the last two.
> > 1. If you dont' fill in a position, Python will use the default. 0 for the left hand side of the `:`, and `len(elements)` for the right hand side.
> > 2. You can request a slice longer than your list (e.g. up to 20), but Python may not give you that many items back.
> > 3. List slicing can only be done with integers, not floats.
> {: .solution}
{: .question}

```
# Check your answers here!
```

## Stride

However, list slicing can be more complicated. You can additionally use a 'stride' parameter, which is how Python should strep through the list. To take every other element from a list:

```python
values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
print(values[0:12:2]) # every other value
print(values[1:12:2]) # every other value from the second value
print(values[::2]) # the start and end are optional
print(values[::3]) # every third value in the list.
```

So list slicing together is either `list[low:high]` or `list[low:high:stride]`, where low and high are optional if you just want to go to the end of the list.

# Type Conversion

Just list with converting `"1.5"` to an float with the `float()` function, or `3.1` to a string with `str()`, we can do the same with lists using the `list()` function, and sets with `set()`:

```python
# Convert text to a list
print(list("sometext"))
```

Converting a list back into text is likewise possible, but you need to use the special function `join`. Join is a function of a `str`, which accepts a list

```python
word = ['c', 'a', 'f', 'e']
print("-".join(word))
```

It takes the string you called it on, and uses that as a separator. Then for the list that you provide, it joins that together with the separator.

# Exercise Time

> ### {% icon question %} Question: Fill in the Blanks
>
> Fill in the blanks so that the program below produces the output shown.
>
> ```
> values = ____
> values.____(1)
> values.____(3)
> values.____(5)
> print('first time:', values)
> values = values[____]
> print('second time:', values)
> ```
>
> ```
> first time: [1, 3, 5]
> second time: [3, 5]
> ```
>
> > ### {% icon solution %} Solution
> > ```
> > values = []
> > values.append(1)
> > values.append(3)
> > values.append(5)
> > print('first time:', values)
> > values = values[1:]
> > print('second time:', values)
> > ```
> {: .solution}
{: .question}

```
# Fill in the blanks here!

values = ____
values.____(1)
values.____(3)
values.____(5)
print('first time:', values) # Should print [1, 3, 5]
values = values[____]
print('second time:', values) # [should print 3, 5]
```

> ## How Large is a Slice?
>
> If `start` and `stop` are both non-negative integers,
> how long is the list `values[start:stop]`?
>
> > ## Solution
> > The list `values[start:stop]` has up to `stop - start` elements.  For example,
> > `values[1:4]` has the 3 elements `values[1]`, `values[2]`, and `values[3]`.
> > Why 'up to'? As we saw in [episode 2]({{ page.root }}/02-variables/),
> > if `stop` is greater than the total length of the list `values`,
> > we will still get a list back but it will be shorter than expected.
> {: .solution}
{: .question}

> ## From Strings to Lists and Back
>
> Given this:
>
> ```
> print('string to list:', list('tin'))
> print('list to string:', ''.join(['g', 'o', 'l', 'd']))
> ```
> ```
> string to list: ['t', 'i', 'n']
> list to string: gold
> ```
>
> 1.  What does `list('some string')` do?
> 2.  What does `'-'.join(['x', 'y', 'z'])` generate?
>
> > ## Solution
> > 1. [`list('some string')`](https://docs.python.org/3/library/stdtypes.html#list) converts a string into a list containing all of its characters.
> > 2. [`join`](https://docs.python.org/3/library/stdtypes.html#str.join) returns a string that is the _concatenation_
> >    of each string element in the list and adds the separator between each element in the list. This results in
> >    `x-y-z`. The separator between the elements is the string that provides this method.
> {: .solution}
{: .question}

> ## Working With the End
>
> What does the following program print?
>
> ```
> element = 'helium'
> print(element[-1])
> ```
>
> 1.  How does Python interpret a negative index?
> 2.  If a list or string has N elements,
>     what is the most negative index that can safely be used with it,
>     and what location does that index represent?
> 3.  If `values` is a list, what does `del values[-1]` do?
> 4.  How can you display all elements but the last one without changing `values`?
>     (Hint: you will need to combine slicing and negative indexing.)
>
> > ## Solution
> > The program prints `m`.
> > 1. Python interprets a negative index as starting from the end (as opposed to
> >    starting from the beginning).  The last element is `-1`.
> > 2. The last index that can safely be used with a list of N elements is element
> >    `-N`, which represents the first element.
> > 3. `del values[-1]` removes the last element from the list.
> > 4. `values[:-1]`
> {: .solution}
{: .question}

> ## Stepping Through a List
>
> What does the following program print?
>
> ```
> element = 'fluorine'
> print(element[::2])
> print(element[::-1])
> ```
>
> 1.  If we write a slice as `low:high:stride`, what does `stride` do?
> 2.  What expression would select all of the even-numbered items from a collection?
>
> > ## Solution
> > The program prints
> > ```
> > furn
> > eniroulf
> > ```
> > 1. `stride` is the step size of the slice.
> > 2. The slice `1::2` selects all even-numbered items from a collection: it starts
> >    with element `1` (which is the second element, since indexing starts at `0`),
> >    goes on until the end (since no `end` is given), and uses a step size of `2`
> >    (i.e., selects every second element).
> {: .solution}
{: .question}

> ## Slice Bounds
>
> What does the following program print?
>
> ```
> element = 'lithium'
> print(element[0:20])
> print(element[-1:3])
> ```
>
> > ## Solution
> > ```
> > lithium
> >
> > ```
> The first statement prints the whole string, since the slice goes beyond the total length of the string.
> The second statement returns an empty string, because the slice goes "out of bounds" of the string.
> {: .solution}
{: .question}

# Dictionaries

When you think of a Dictionary, you should think of a real life Dictionary, they map some key to a value. Like a term to it's definition

Key | Value
--- | ---
`Eichhörnchen` | Squirrel
`火锅` | Hot Pot

Or a Country to it's population

Key               | Value
---               | ---
South Sudan       | 492,970
Australia         | 411,667
Guinea            | 1,660,973
Morocco           | 573,895
Maldives          | 221,678
Wallis and Futuna | 1,126
Eswatini          | 94,874
Namibia           | 325,858
Turkmenistan      | 1,031,992

In Python we create a dictionary with `{}` and use `:` to separate keys and values. Turning the above list into a Python dictionary, it would look like:

```python
populations = {
  "South Sudan": 492970,
  "Australia": 411667,
  "Guinea": 1660973,
  "Morocco": 573895,
  "Maldives": 221678,
  "Wallis and Futuna": 1126,
  "Eswatini": 94874,
  "Namibia": 325858,
  "Turkmenistan": 1031992,
}
```

You can see a string (the country name) being used for the key, and then the number (an integer) as the value. (Would a float make sense? Why or why not?)

> ### {% icon tip %} Tip: Other Names
> They're also sometimes called associative arrays (because they're an array or list of values that associate a key to a value) or maps (because they map a key to a value), depending on what you're reading.
{: .tip}

## Methods

You can access both the keys, and the values

```python
print(populations.keys())
print(populations.values())
```

These will print out two list-like objects. They will become more useful in the future when we talk about looping over dictionaries and processing all of the values within.

## Accessing Values

Just like lists where you access by the position in the list

```python
print(populations["Namibia"])
```

And just like lists, if you try an access a key that isn't there or an index outside of the range of the list:

```
print(populations["Mars"])
```


> ### {% icon tip %} Tip: Dictionaries are faster than lists for looking up values
> Just like in real life, searching a dictionary for a specific term is quite fast. Often a lot faster than searching a list for a specific value.
>
> For those of you old enough to remember the paper version of a dictionary, you knew that As would be at the start and Zs at the end, and probably Ms around the middle. And if you were looking for a word like "Squirrel", you'd open the dictionary in the middle, maybe decide it was in the second half of the book, randomly choose a page in the second half, and you could keep deciding if it was "before" or "after" the current page, never even bothering to search the first half.
>
> Conceptually, compared with a list, you can't make this guess of if the item is in the first or second half. You need to search item by item, it would be like reading page by page until you get to Squirrel in the dictionary.
{: .tip}

## Modifying Dictionaries

Adding new values to a dictionary is easy, it's very similar to replacing a value in a list.

```
# For lists we did
x = ['x', 'y', 'z']
x[0] = 'a'
print(x)
```

For dictionaries, it's essentially the same, we access the 'place' in the dictionary just like we did with a list, and set it to a value

```
populations["Mars"] = 6 # robots
print(populations)
```

And similarly, removing items is the same as it was for lists:

```
del x[0] # Removes the first item
print(x)
```

And with dictionaries you delete by specifying which position/key you want to remove

```
del populations['Australia']
print(populations)
```

## Exercises

> ### {% icon question %} Question: DNA Complement
>
> DNA is usually in the form of dsDNA, a paired strand, where A maps to T and C maps to G and vice versa.
> But when we're working with DNA sequences in bioinformatics, we often only store one strand, because we can calculate the *complement* on the fly, when we need.
>
> Write a dictionary that lets you look up the letters A, C, T, and G and find their complements.
>
> > ### {% icon solution %} Solution
> > You need to have the complements of every base. If you just defined 'A' and 'C', how would you look up the complement when you want to translate a 'T' or a 'G'? It's not easily possible to look up a key by a value, only to search a key and find a value.
> > ```
> > translation = {
> > 'A': 'T',
> > 'T': 'A',
> > 'C': 'G',
> > 'G': 'C',
> > }
> > ```
> {: .solution}
{: .question}

> ### {% icon question %} Question: Modifying an array
>
> Fill in the blanks to make the execution correct:
>
> ```
> variants = {
>   'B.1.1.7': 26267,
>   'B.1.351': 439,
> }
> variants[_____] =  _____
> print(variants) # Should print {'B.1.1.7': 26267, 'B.1.351': 439, 'P.1': 384}
> __________
> print(variants) # Should print {'B.1.1.7': 26267, 'B.1.351': 439, 'P.1': 384, 'B.1.617.2': 43486}
> # Maybe we've exterminated B.1.1.7 and B.1.351, remove their numbers.
> del _______
> del _______
> print(variants[______]) # Should print 384
> print(variants[______]) # Should print 43486
> ```
>
> > ### {% icon solution %} Solution
> >
> > variants = {
> >   'B.1.1.7': 26267,
> >   'B.1.351': 439,
> > }
> > variants['P.1'] = 384
> > print(variants) # Should print {'B.1.1.7': 26267, 'B.1.351': 439, 'P.1': 384}
> > variants['B.1.617.2'] = 43486
> > print(variants) # Should print {'B.1.1.7': 26267, 'B.1.351': 439, 'P.1': 384, 'B.1.617.2': 43486}
> > # Maybe we've exterminated B.1.1.7 and B.1.351, remove their numbers.
> > del variants['B.1.1.7']
> > del variants['B.1.351']
> > print(variants['P.1']) # Should print 384
> > print(variants['B.1.617.2']) # Should print 43486
> > ```
> >
> {: .solution}
{: .question}


# Choosing the Right Data Type

Choosing the correct data type can sometimes require some thought, and even discussion with colleagues. And don't be afraid to search the internet for how other people have done it!


Data type       | Examples                  | When to use it | When **not** to use it
---             | ---                       | ---- | ----
Boolean (`bool`) | `True`, `False` | If there are only two possible states, true or false | If your data is not binary
Integer (`int`) | 1, 0, -1023, 42           | Countable, singular items. How many patients are there, how many events did you record, how many variants are there in the sequence | If doubling or halving the value would not make sense: do not use for e.g. patient IDs, or phone numbers. If these are integers you might accidentally do math on the value.
Float (`float`) | 123.49, 3.14159, -3.33334 | If you need more precision or partial values. Recording distance between places, height, mass, etc.
Strings (`str`) | 'patient_12312', 'Jane Doe', '火锅' | To store free text, identifiers, sequence IDs, etc. | If it's truly a numeric value you can do calculations with, like adding or subtracting or doing statistics.
List / Array (`list`) | `['A', 1, 3.4, ['Nested']]` | If you need to store a list of items, like sequences from a file. Especially if you're reading in a table of data from a file. | If you want to retrieve individual values, and there are clear identifiers it might be better as a dict.
Dictionary / Associative Array / map (`dict`) | `{"weight": 3.4, "age": 12, "name": "Fluffy"}` | When you have identifiers for your data, and want to look them up by that value. E.g. looking up sequences by an identifier, or data about students based on their name. Counting values. | If you just have a list of items without identifiers, it makes more sense to just use a list.
Tuple (`tuple`) | `(1, 2, 3)` | Probably not necessary, but if you want an immutable list | If you want to ever modify your list.
Set (`set`) | `{1, 2, 3}` | Use this to make a list of data unique | If you want to retain a specific order of your list (e.g. the order you read the data in) a set is inappropriate.

## Exercises

> ### {% icon question %} Question: Which Datatype
>
> 1. Name
> 2. Weight
> 3. Gender
> 4. Hair Colour
>
> > ### {% icon solution %} Solution
> >
> > 1. Here a string would be a good choice. (And probably just a single `name` string, rather than a `first` and `last` name, as not all humans have two names! And some have more than two.)
> > 2. A float is good type for storing weight. You could use an integer as well, but you would be losing the decimal places. In some applications that is ok, in others, you require the precision.
> > 3. This is a case where you should consider carefully the application, but `bool` is usally the wrong answer. Are you recording patient data? Is their expressed gender the correct variable or did you need sex? Chromosomal sex is also more complicated and cannot be stored with a true/false value, as people with [Kleinfelters](https://en.wikipedia.org/wiki/Klinefelter_syndrome) exist. A string can be an ok choice here, depending on the application.
> > 4. There is a limited vocabulary humans use to describe hair colour, so a string can be used, or a data type we haven't discussed! An `enum` is an `enumeration`, and when you have a limited set of values that are possible, you can use a `enum` to double check that whatever value is being used (or read from a file, or entered by a user) matches one of the "approved" values.
> >
> {: .solution}
{: .question}
