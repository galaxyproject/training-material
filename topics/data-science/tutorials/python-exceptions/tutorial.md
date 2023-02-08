---
layout: tutorial_hands_on

title: 'Python - Try & Except'
level: Introductory
requirements: []
follow_up_training: []
questions:
- How do I try to execute code, knowing it might fail?
- What are some situations where this is important?
- How can I write my own exceptions.

objectives:
- catch an exception
- raise your own exception

time_estimation:  20M
key_points:
- raise lets your raise your own `Exception`s
- This is mostly used by library authors (which you might become!)
- Use `try`/`except` to catch expected errors and work around them (if possible)
- finally lets you cleanup your temporary files, if you created some.

subtopic: python-modular
contributors:
  - hexylena
  - dirowa
  - bazante1
priority: 8
notebook:
  language: python
---

Try/except are a construct in Python used to catch a potential exception. Sometimes things go wrong in your code! Or in someone else's code in a module. Sometimes some errors might be expected like when you try and read a user supplied file, maybe it isn't available because they've specified the wrong path.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Raise

When you're writing code, sometimes there are errors you might need to handle:

- If you're calculating the mean of a set of numbers, and the list is empty. You don't want to divide by `len(numbers)` and trigger a Zero Division Error.
- If you're opening a file, there maybe a chance that the file is not found.

Using try/excepts allow you to:

- provide friendlier and more useful error messages
- handle expected error cases

For instance, returning to the `mean()` function example, what should `mean([])` return for a hypothetical mean function that calculates the average of a list of numbers. Should it return 0? It probably should return an error. Let's look at some example code:


```python
def mean(numbers):
    return sum(numbers) / len(numbers)

mean([1, 2, 3])
mean([])
```

This raises a `ZeroDivisionError`  but we can make this a more friendly error message by raising our own exception.

```python
def mean(numbers):
    if not numbers:
        raise ValueError("You are calculating the mean of an empty list, which is not possible.")
    return sum(numbers) / len(numbers)

mean([1, 2, 3])
mean([])
```

> <tip-title>Where do ValueError, ZeroDivisionError come from?</tip-title>
> There are loads of different types of exception codes! [The python documentation has a large list](https://docs.python.org/3/library/exceptions.html) of exceptions and some descriptions for when or why those exceptions might be raised.
{: .tip}

Now we get a much more useful error message from the function! Using `raise` is especially important for *library authors*, people writing python modules that we all use. If they provide useful error messages, it helps you as an end user understand what's happening.

## Try & Except

Let's look at how you would handle one of these exceptions, we'll continue with the `mean()` example above.

```python
numbers = []
try:
    result = mean(numbers)
    print(f"Result is {result}")
except ValueError:
    print("We cannot calculate the mean")
```

Here we use `try:` to setup a new block, and this code is **tried**, Python attempts to execute it. Below are one or more `except:` blocks which catch specific errors. Here we have specifically said we know a ValueError can happen, and decided to handle it.

Or for another example, accessing a user supplied file. Oftentimes users will call your program and supply a non-existent, or inacessible file. Here you can use multiple `except` blocks to catch all of those potential errors.

```python
user_supplied_file = 'does-not-exist.txt'
try:
    open(user_supplied_file, 'r')
except FileNotFoundError:
    print(f"The path you supplied ({user_supplied_file}) doesn't exist, please double check it!")
except PermissionError:
    print(f"You supplied a valid file, but it is unreadable. Try changing it's permissions with `chmod +r {user_supplied_file}`")
```

Failing to open a file raises a `FileNotFoundError` which indicates the file isn't available, and `PermissionError` indicates that a file is unreadable. However in practice, sometimes you'll see something like this:

```
# Bad!
try:
    doSomething()
except:
    print("error")
```

This is called a **bare exception**, and will catch any exception, compared with `except ValueError` which only catches value errors. People consider this generally a bad idea, termed *code smell*. (Because it smells (appears) bad!)

## Finally

The last portion of the `try:`/`except:` block is `finally:`, a third block which lets us do cleanup. It's often very nice to your users that if your program fails halfway through, that you cleanup after yourself.

```python
import os

try:
    with open('gene_query.fa', 'w') as handle:
        handle.write(">some fasta sequence we want to search against a database")

    # runQuery('gene_query.fa')

    # But here we have an error! Something goes wrong!
    1 / 0
except:
    # And now our results are invalid.
    print("Something went wrong! Oh no! Check your inputs.")
finally:
    # So we should cleanup this temporary file we created, so it doesn't cause
    # problems or distract the user from the results file.
    # This function will delete a file from your computer:
    os.remove('gene_query.fa')
```

## Fallback

Sometimes we can use `try`/`except` to have a fallback option. Consider the pseudocode below:

```
try:
    runBlast()
except BlastNotAvailable:
    try:
        runBLAT()
    except BLATnotAvailable:
        print("Neither Blast nor BLAT were available.")
```

Sometimes you have a fallback option, some other tool you can use in its place. When that's possible, you can use `try` and `except` to handle those cases and work around potential issues. But this isn't always the case, sometimes you just need to print your error message and stop executing code.
