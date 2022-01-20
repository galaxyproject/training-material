---
layout: tutorial_hands_on

title: Python - Lists & Strings
level: Introductory
requirements: []
follow_up_training: []
questions:
- "How can I store multiple values?"

objectives:
- "Explain why programs need collections of values."
- "Write programs that create flat lists, index them, slice them, and modify them through assignment and method calls."

time_estimation:  30M
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

## Heterogenous Lists

Lists may contain values of different types. A single list may contain numbers, strings, and anything else.

```
goals = [1, 'Create lists.', 2, 'Extract items from lists.', 3, 'Modify lists.']
```

# Strings are like Lists

Strings of text like `name = "Helena"` or `patient_id = "19237zud830"` are very similar conceptually to lists. Except instead of being a list of numbers, they lists of characters.

In a number of older programing languages, strings are indeed arrays of numbers internally. However python hides a lot of that complexity from us, so we can just work with text.

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

# Sets

Sets are like lists, except:

1. They're unique
2. After you create them, they're immutable.

You generally do not need sets, but once you get comfortable with python there will be times you might want to use them. And commonly you'll see code like:

```python
data = ['a', 'a', 'b', 'c', 'c', 'c']
unique_data = set(data)
```

in order to get the unique values in the list.

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

That is a really cumbursome way to access lists, so Python makes your life easier. If you need to access somewhere off of the end of the list, say, the last two items, you can just use negative values!

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
