---
layout: tutorial_hands_on

title: Python - Functions
level: Introductory
requirements:
 -
   type: internal
   topic_name: data-science
   tutorials:
     - python-math
follow_up_training: []
questions:
- How do I write functions in Python?
- What is a function?
- What do they look like?
- Fill in the missing part of a function
#- Discussion of the results
#- Write a new function that does a sepcific computation, building on previous portion (exercise)

objectives:
- Understand the structure of a "function" in order to be able to construct their own functions and predict which functions will not work.

time_estimation:  30M
key_points:
- Functions are foundational in Python
- Everything you do will require functions
- Functions keep your code DRY (don't repeat yourself), reducing your copying and pasting or rewriting the same block of code.
- Deciding what part of your code should, or should not be, a function is something that will come with practice.

subtopic: python-modular
contributors:
  - carpentries
  - hexylena
  - dirowa
  - bazante1

priority: 2
notebook:
  language: python
  pyolite: true

abbreviations:
  DRY: Don't Repeat Yourself
---

Functions are the basic unit of all work in Python! Absolutely everything you do uses functions. Conceptually, functions are super simple. Just like in maths, they take an input, do some transformation, and return an output. As an example, `f(x) = x + 2` is a function that calculates whatever value you request, plus two. But functions are foundational, so, you should understand them well before moving on.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


## What is a Function

Functions are a way to re-use some computation you want to do, multiple times. If you don't have a function, you need to re-write the calculation every time, so we use functions to collect those statements and make them easy to re-run. Additionally they let us "parameterise" some computation. Instead of computing the same value every time, we can template it out, we can re-run the computation with new inputs, and see new results.

> > <code-in-title>Math</code-in-title>
> > ```
> > # Define our function
> > f(x) = 3 * x
> > # Compute some value
> > f(3) # is 9
> > ```
> {: .code-in}
>
> > <code-out-title>Python</code-out-title>
> > ```
> > # Define our function
> > def f(x):
> >    return 3 * x
> > # Compute some value
> > f(3)
> > ```
> {: .code-out}
>
{: .code-2col}

We've talked about mathematical functions before, but now we'll talk about more programing-related functions


## Create Functions

Human beings can only keep a few items in working memory at a time. Breaking down larger/more complicated pieces of code in functions helps in understanding and using it. A function can be re-used. Write one time, use many times. (Known as staying {DRY} in the industry.)

Knowing what Americans mean when they talk about temperatures and weather can be difficult, but we can wrap the temperature conversion calculation ($$^{\circ}\text{C} = (^{\circ}\text{F} - 32) * \dfrac{5}{9}$$) up as a function that we can easily re-use.

```python
def fahr_to_celsius(temp):
    return ((temp - 32) * (5/9))
```

![The above function fahr to celsius is shown except annotated. def is labelled "def statement", fahr_to_celsius is noted as the function name. Inside parentheses is temp and an arrow shows it is called parameter names. The next line which is indented is annotated as the function body which has a return statement and the calculation from above.](../../images/python-basics/Figure7_Functions.png)

The function definition opens with the keyword `def` followed by the name of the function `fahr_to_celsius` and a parenthesized list of parameter names `temp`. The body of the function — the statements that are executed when it runs — is indented below the definition line. The body concludes with a `return` keyword followed by the return value.

When we call the function, the values we pass to it are assigned to those variables so that we can use them inside the function. Inside the function, we use a return statement to send a result back to whoever asked for it.

Let’s try running our function.

```python
fahr_to_celsius(32)
print(f"freezing point of water: {fahr_to_celsius(32)}C")
print(f"boiling point of water: {fahr_to_celsius(212)}C")
```

> <tip-title>Formatting Strings</tip-title>
> There are several ways to print out a few values in python. We'd recommend you to use `f-strings` as it's the cleanest and most modern way to do it.
>
> `f`-strings start with an `f` (very descriptive name eh?). Within the text between the single or double quotes (`'`/`"`) you can use curly braces to refer to variables or python code which will be placed there in the string.
>
> ```
> a = 10
> b = f"Here is the value of a: {a}"
> print(b)
> print(f"Here is the value of a: {5 + 5}")
> print(f"Here is the value of a: {function_that_returns_10()}")
> ```
>
> All of those would print out `Here is the value of a: 10`.
>
> f-strings can be a lot fancier for formatting decimal places, but we don't need that for now. Just know:
>
> 1. Start with an `f`
> 2. Use braces to use the value of a variable, a function, or some python expression.
{: .tip}

We’ve successfully called the function that we defined, and we have access to the value that we returned.

Now that we’ve seen how to turn Fahrenheit into Celsius, we can also write the function to turn Celsius into Kelvin:

```python
def celsius_to_kelvin(temp_c):
    return temp_c + 273.15

print(f'freezing point of water in Kelvin: {celsius_to_kelvin(0.)}')
```

> <tip-title>What is `0.`</tip-title>
> That's a float! A `.` in a number makes it a float, rather than an integer.
{: .tip}

What about converting Fahrenheit to Kelvin? We could write out both formulae, but we don’t need to. Instead, we can *compose* the two functions we have already created:

```python
def fahr_to_kelvin(temp_f):
    temp_c = fahr_to_celsius(temp_f)
    temp_k = celsius_to_kelvin(temp_c)
    return temp_k

print(f'boiling point of water in Kelvin: {fahr_to_kelvin(212.0)}')
```

This is our first taste of how larger programs are built: we define basic operations, then combine them in ever-larger chunks to get the effect we want. Real-life functions will usually be larger than the ones shown here — typically half a dozen to a few dozen lines — but they shouldn’t ever be much longer than that, or the next person who reads it won’t be able to understand what’s going on.

### Documentation

Documenting your code is extremely, *extremely*, **extremely** important to do. We all forget what we're doing, it's only normal, so documenting what you're doing is key to being able to restart work later.

```python
def fahr_to_kelvin(temp_f):
    """
    Converts a temperature from Fahrenheit to Kelvin

    temp_f: the temperature in Fahrenheit

    returns the temperature in Celsius
    """
    temp_c = fahr_to_celsius(temp_f)
    temp_k = celsius_to_kelvin(temp_c)
    return temp_k

print(f'boiling point of water in Kelvin: {fahr_to_kelvin(212.0)}')
```

For a function this small, with such a descriptive name (`fahr_to_kelvin`) it feels quite obvious what the function should do, what inputs it takes, what outputs it produces. However, you will thank yourself in the future if you do this now. You may think you will remember what the code does, but, be kind to your future self who is busy and stressed and may not want to spend time reading the code over again to figure out what every single function does.

> <question-title>Converting statements to functions</question-title>
> A lot of what you'll do in programing is to turn a procedure that you want to do, into statements and a function.
>
> Fill in the missing portions of this function, two average numbers. Then use it to find the average of 32326 and 631
>
> > <solution-title></solution-title>
> > ```
> > def average2(a, b):
> >     c = (a + b) / 2
> >     return c
> > ```
> > We call it "average2" here because it will only average two numbers. It will not work for three numbers or a list of them.
> {: .solution}
{: .question}

```python
# Test out solutions here!
def average2(a, b):
    c =
    return c

print(average2(32326, 631))
```

> <question-title>A more complicated example</question-title>
>
> The formula for a 90° triangle can be expressed as: $$c = \sqrt{a^2 + b^2}$$
>
> 1. Write a function which takes `a` and `b`, and calculates `c`
> 2. Name this function "pythagorus"
> 3. Remember to import math, if you haven't already
>
> > <solution-title></solution-title>
> > ```
> > def pythagorus(a, b):
> >     c = math.sqrt(math.pow(a, 2) + math.pow(b, 2))
> >     return c
> > ```
> {: .solution}
{: .question}

```python
# Test out solutions here!


print(pythagorus(1234, 4321)) # Should return 4493.750883170984
```


### Variable Scope

In composing our temperature conversion functions, we created variables inside of those functions, `temp`, `temp_c`, `temp_f`, and `temp_k`. We refer to these variables as local variables because they no longer exist once the function is done executing. If we try to access their values outside of the function, we will encounter an error:

```python
print(f'Again, temperature in Kelvin was: {temp_k}')
```

If you want to reuse the temperature in Kelvin after you have calculated it with fahr_to_kelvin, you can store the result of the function call in a variable:

```python
temp_kelvin = fahr_to_kelvin(212.0)
print(f'temperature in Kelvin was: {temp_kelvin}')
```

Watch out for scope issues:

- Variables created inside a function will stay inside the function
- Variables created outside of the function can be accessible inside the function, but you should not do this!
- Ensure that variables are properly scoped will prevent later errors when working with modules or testing big projects.


```python
a = 1
b = 2.0

# Location 1

def some_generic_computation(x, y):
    c = x + y
    d = x * y
    e = c / d
    # Location 2
    return e

# Location 3
```

> <question-title>Scope</question-title>
> Given the above code, which variables are accessible at Locations 1, 2, and 3?
>
> > <solution-title></solution-title>
> > 1. a, b
> > 2. a and b are there but you shouldn't use these. x, y, c, d, e are also accessible.
> > 3. a, b.
> {: .solution}
{: .question}

### Defining Default parameters

If we usually want a function to work one way, but occasionally need it to do something else, we can allow people to pass a parameter when they need to but provide a default to make the normal case easier. The example below shows how Python matches values to parameters:

```python
def display(a=1, b=2, c=3):
    print(f'a: {a}, b: {b}, c: {c}')

# no parameters:
display()
# one parameter:
display(55)
# two parameters:
display(55, 66)
```

As this example shows, parameters are matched up from left to right, and any that haven’t been given a value explicitly get their default value. We can override this behavior by naming the value as we pass it in:

```python
# only setting the value of c
display(c=77)
```


> <question-title>Exercise: Signing a message</question-title>
> Let's test out a default argument. Imagine you are printing out a message, and at the bottom it should have a signature.
>
> Inputs:
> - `message`: a variable that is always provided to the function, it has no default.
> - `signature`: a variable that can be *optionally* provided, it should have a default like your name.
>
> You can accomplish this with *three print statements*:
> 1. Print the message
> 2. Print nothing (i.e. `print()`)
> 3. Print a signature variable.
>
> > <solution-title></solution-title>
> > ```
> > def myFunction(message, signature="Your name"):
> >     print(message)
> >     print()
> >     print(signature)
> > ```
> {: .solution}
{: .question}

```python
# Test things out here!
def myFunction # Fix this function!


# Here are some test cases, for you to check if your function works!
myFunction('This is a message')
myFunction('This is a message', signature='Jane Doe')
```
