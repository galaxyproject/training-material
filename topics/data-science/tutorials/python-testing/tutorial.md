---
layout: tutorial_hands_on

title: Python - Testing
level: Intermediate
requirements: []
follow_up_training: []
time_estimation:  45M
questions:
- "Does the code we develop work the way it should do?"
- "Can we (and others) verify these assertions for themselves?"
- "To what extent are we confident of the accuracy of results that appear in publications?"
objectives:
- "Explain the reasons why testing is important"
- "Describe the three main types of tests and what each are used for"
- "Implement and run unit tests to verify the correct behaviour of program functions"
key_points:
- "The three main types of automated tests are **unit tests**, **functional tests** and **regression tests**."
- "We can write unit tests to verify that functions generate expected output given a set of specific inputs."
- "It should be easy to add or change tests, understand and run them, and understand their results."
- "We can use a unit testing framework like `unittest` to structure and simplify the writing of tests."
- "We should test for expected errors in our code."
- "Testing program behaviour against both valid and invalid inputs is important and is known as **data validation**."

subtopic: python-modular
contributions:
  authorship:
  - carpentries
  - hexylena
  editing:
  - dirowa
  - bazante1
  funding:
  - carpentries
  - avans-atgm

priority: 10 
notebook:
  language: python
---

Here we will cover the basics of testing, an important part of software development. Testing lets you know that your code is correct in many situations that matter to you.

> <comment-title></comment-title>
>
> This tutorial is significantly based on [the Carpentries](https://carpentries.org) lesson ["Intermediate Research Software Development"](https://carpentries-incubator.github.io/python-intermediate-development/).
>
{: .comment}

Being able to demonstrate that a process generates the right results is important in any field of research, whether it's software generating those results or not. So when writing software we need to ask ourselves some key questions:

- Does the code we develop work the way it should do?
- Can we (and others) verify these assertions for themselves?
- Perhaps most importantly, to what extent are we confident of the accuracy of results that appear in publications?

If we are unable to demonstrate that our software fulfills these criteria, why would anyone use it? Having well-defined tests for our software are useful for this, but manually testing software can prove an expensive process.

Automation can help, and automation where possible is a good thing - it enables us to define a potentially complex process in a repeatable way that is far less prone to error than manual approaches. Once defined, automation can also save us a lot of effort, particularly in the long run. In this episode we'll look into techniques of automated testing to improve the predictability of a software change, make development more productive, and help us produce code that works as expected and produces desired results.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## What Is Software Testing?

For the sake of argument, if each line we write has a 99% chance of being right, then a 70-line program will be wrong more than half the time. We need to do better than that, which means we need to test our software to catch these mistakes.

We can and should extensively test our software manually, and manual testing is well-suited to testing aspects such as graphical user interfaces and reconciling visual outputs against inputs. However, even with a good test plan, manual testing is very time consuming and prone to error. Another style of testing is automated testing, where we write code that tests the functions of our software. Since computers are very good and efficient at automating repetitive tasks, we should take advantage of this wherever possible.

There are three main types of automated tests:

- **Unit tests** are tests for fairly small and specific units of functionality, e.g. determining that a particular function returns output as expected given specific inputs.
- **Functional or integration tests** work at a higher level, and test functional paths through your code, e.g. given some specific inputs, a set of interconnected functions across a number of modules (or the entire code) produce the expected result. These are particularly useful for exposing faults in how functional units interact.
- **Regression tests** make sure that your program's output hasn't changed, for example after making changes your code to add new functionality or fix a bug.
- **Property tests** are an advanced testing strategy to find corner cases.

For the purposes of this course, we'll focus on unit tests. But the principles and practices we'll talk about can be built on and applied to the other types of tests too.

## Example Codebase

```python
def rle_encode(input_string):
    count = 1
    prev = ""
    lst = []
    for character in input_string:
        if character != prev:
            if prev:
                entry = (prev, count)
                lst.append(entry)
            count = 1
            prev = character
        else:
            count += 1
    entry = (character, count)
    lst.append(entry)
    return "".join([f"{k}{v}" for k, v in lst])


def rle_decode(lst):
    q = ""
    for character, count in zip(lst[::2], lst[1::2]):
        q += character * int(count)
    return q
```

This is a simple run length encoding and decoding function. It does not handle a lot of corner cases or odd input, but it may do what we need: compressing long strings of repetitive texts. It could potentially provide good genomic data compression. (Spoilers, it generally doesn't. There are better compression algorithms.)

```python
print(rle_encode("aaabbc"))
```

Let's look at how we can test this cod.


## Writing Tests to Verify Correct Behaviour

### One Way to Do It?

One way to test our functions would be to write a series of checks or tests, each executing a function we want to test with known inputs against known valid results, and throw an error if we encounter a result that is incorrect

```python
test_input = "abba"
test_output = "a1b2a1"

assert rle_encode(test_input) == test_output
```

So we use the assert keyword - part of Python - to test that our calculated result is the same as our expected result. This function explicitly checks that the two values are the same, and throws an AssertionError if they are not.

Let's write some more test cases:

```python
test_input = "wwwwssssb"
test_output = "w5s4b1"
assert rle_encode(test_input) == test_output

test_round_trip = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
assert test_round_trip == rle_decode(rle_encode(test_round_trip))

test_input = "baabaa"
test_output = "b1a2b1a2"
assert rle_encode(test_input) == test_output
```

However, if we were to enter these in this order, we’ll find we get the following after the first test:

```
AssertionError                            Traceback (most recent call last)
<ipython-input-30-05489b2ef047> in <module>
      1 test_input = "wwwwssssb"
      2 test_output = "w5s4b1"
----> 3 assert rle_encode(test_input) == test_output

AssertionError: 
```

This tells us our function couldn't handle these inputs.

We could put these tests in a separate script to automate the running of these tests. But a Python script halts at the first failed assertion, so the second and third tests aren't run at all. It would be more helpful if we could get data from all of our tests every time they're run, since the more information we have, the faster we're likely to be able to track down bugs. It would also be helpful to have some kind of summary report: if our set of tests - known as a **test suite** - includes thirty or forty tests (as it well might for a complex function or library that's widely used), we'd like to know how many passed or failed.


Going back to our failed first test, what was the issue? As it turns out, the test itself was incorrect, and should have read:

```python
test_input = "wwwwssssb"
test_output = "w4s4b1"
assert rle_encode(test_input) == test_output
```

Which highlights an important point: as well as making sure our code is returning correct answers, we also need to ensure the tests themselves are also correct. Otherwise, we may go on to fix our code only to return an incorrect result that *appears* to be correct. So a good rule is to make tests simple enough to understand so we can reason about both the correctness of our tests as well as our code. Otherwise, our tests hold little value.

### Using a Testing Framework

Keeping these things in mind, here's a different approach that builds on the ideas we've seen so far but uses a **unit testing framework**. In such a framework we define our tests we want to run as functions, and the framework automatically runs each of these functions in turn, summarising the outputs. And unlike our previous approach, it will run every test regardless of any encountered test failures.

Most people don't enjoy writing tests, so if we want them to actually do it, it must be easy to:

- Add or change tests,
- Understand the tests that have already been written,
- Run those tests, and
- Understand those tests' results

Test results must also be reliable. If a testing tool says that code is working when it's not, or reports problems when there actually aren't any, people will lose faith in it and stop using it.

```python
from main import rle_encode, rle_decode
import unittest

class TestRunLengthEncoding(unittest.TestCase):
    def test_encode(self):
      i = "aaabbcccd"
      expected = "a3b2c3d1"
      actual = rle_encode(i)
      self.assertEqual(expected, actual)

    def test_decode(self):
        i = "a3b2c3d1"
        expected = "aaabbcccd"
        actual = rle_decode(i)
        self.assertEqual(expected, actual)

    def test_empty(self):
        test_input = "wwwwssssb"
        test_output = "w4s4b1"
        self.assertEqual(test_input, test_output)

    def test_roundtrip_short(self):
        test_input = "bbb"
        test_output = rle_decode(rle_encode(test_input))
        self.assertEqual(test_input, test_output)

    def test_roundtrip_long(self):
        test_input = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
        test_output = rle_decode(rle_encode(test_input))
        self.assertEqual(test_input, test_output)

if __name__ == '__main__':
    unittest.main()
```

So here, although we have specified our tests as separate functions in a class, they run the same assertions. Each of these test functions, in a general sense, are called **test cases** - these are a specification of:

- Inputs 
- Execution conditions - any imports we might need to do
- Testing procedure
- Expected outputs

And here, we're defining each of these things for a test case we can run independently that requires no manual intervention.

> <hands-on-title>Run this locally</hands-on-title>
> 1. Create a new folder somewhere
> 2. Save the RLE encode/decode functions as `main.py` in that folder
> 3. Save the unittests above as `test.py` in that folder
> 4. Run `python test.py`
{: .hands_on}

Going back to our list of requirements, how easy is it to run these tests? We can do this using a Python package called `unittest`. This is a built-in testing framework that allows you to write test cases using Python. You can use it to test things like Python functions, database operations, or even things like service APIs - essentially anything that has inputs and expected outputs. 

> <tip-title>What About Unit Testing in Other Languages?</tip-title>
>
> Other unit testing frameworks exist for Python, including Nose and pyunit, and the approach to unit testing can be translated to other languages as well, e.g. FRUIT for Fortran, JUnit for Java (the original unit testing framework), Catch for C++, etc.
{: .tip}


### What About Testing for Errors?

There are some cases where seeing an error is actually the correct behaviour, and Python allows us to test for exceptions. Add this test in `tests/test_models.py`:

```python
import unittest

class ExpectedFailureCase(unittest.TestCase):

    @unittest.expectedFailure
    def test_bad_input():
        from main import rle_decode

        # invalid input, we expect pairs of symbols and numbers.
        rle_decode("a1b")
```

Run all your tests as before.

> <tip-title>Why Should We Test Invalid Input Data?</tip-title>
>
> Testing the behaviour of inputs, both valid and invalid, is a really good idea and is known as *data validation*. Even if you are developing command line software that cannot be exploited by malicious data entry, testing behaviour against invalid inputs prevents generation of erroneous results that could lead to serious misinterpretation (as well as saving time and compute cycles which may be expensive for longer-running applications). It is generally best not to assume your user's inputs will always be rational.
>
{: .tip}

### Testing Frameworks

There are many, many testing frameworks for Python that solve similar problems. Look around and see what fits well with your project!

- [unittest](https://docs.python.org/3/library/unittest.html) is built into Python's stdlib.
- [pytest](https://docs.pytest.org/en/7.1.x/) is a heavily recommended framework.
- [nose](https://nose.readthedocs.io/en/latest/) is an alternative to pytest that does many of the same things.
- [hypothesis](https://hypothesis.readthedocs.io/en/latest/index.html) does property testing, the author's favourite.
- [doctest](https://docs.python.org/3/library/doctest.html) tests examples written directly in the documentation, which keeps test cases close to functions.
