---
title: Regular Expressions 101
area: tools
box_type: tip
layout: faq
contributors: [shiltemann]
---

Regular expressions are a standardized way of describing patterns in textual data. They can be extremely useful for tasks such as finding and replacing data. They can be a bit tricky to master, but learning even just a few of the basics can help you get the most out of Galaxy.

#### Finding

Below are just a few examples of basic expressions:

| Regular expression | Matches |
|--------------------|-------------|
| `abc`              | an occurrence of `abc` within your data |
| `(abc|def)`        | `abc` *or* `def` |
| `[abc]`            | a single character which is either `a`, `b`, or `c` |
| `[^abc]`           | a character that is NOT `a`, `b`, nor `c` |
| `[a-z]`            | any lowercase letter |
| `[a-zA-Z]`         | any letter (upper or lower case)|
| `[0-9]`            | numbers 0-9 |
| `\d`               | any digit (same as `[0-9]`) |
| `\D`               | any non-digit character |
| `\w`               | any alphanumeric character |
| `\W`               | any non-alphanumeric character |
| `\s`               | any whitespace |
| `\S`               | any non-whitespace character |
| `.`                | *any* character |
| `\.`
| `{x,y}`            | between *x* and *y* repetitions |
| `^`                | the beginning of the line |
| `$`                | the end of the line |

**Note:** you see that characters such as `*`, `?`, `.`, `+` etc have a special meaning in a regular expression. If you want to match on those characters, you can escape them with a backslash. So `\?` matches the question mark character exactly.

**Examples**

| Regular expression   | matches                |
|----------------------|------------------------|
| `\d{4}`              | 4 digits (e.g. a year) |
| `chr\d{1,2}`         | `chr` followed by 1 or 2 digits |
| `.*abc$`             | anything with `abc` at the end of the line |
| `^$`                 | empty line |



#### Replacing

Sometimes you need to capture the exact value you matched on, in order to use it in your replacement, we do this using capture groups `(...)`, which we can refer to using `\1`, `\2` etc for the first and second captured values.

| Regular expression      | Input        | Captures               |
|-------------------------|--------------|------------------------|
| `chr(\d{1,2})`          | `chr14`      | `\1 = 14`              |
| `(\d{2}) July (\d{4})`  | 24 July 1984 | `\1 = 24`, `\2 = 1984` |

An expression like `s/find/replacement/g` indicates a replacement expression, this will search (`s`) for any occurrence of `find`, and replace it with `replacement`. It will do this globally (`g`) which means it doesn't stop after the first match.

Example: `s/chr(\d{1,2})/CHR\1/g` will replace `chr14` with `CHR14` etc.

Note: In Galaxy, you are often asked to provide the find and replacement expressions separately, so you don't have to use the `s/../../g` structure.

There is a lot more you can do with regular expressions, and there are a few different flavours in different tools/programming languages, but these are the most important basics that will already allow you to do many of the tasks you might need in your analysis.

**Tip:** [RegexOne](https://regexone.com/) is a nice interactive tutorial to learn the basics of regular expressions.

**Tip:** [Regex101.com](https://regex101.com/) is a great resource for interactively testing and constructing your regular expressions, it even provides an explanation of a regular expression if you provide one.

**Tip:** [Cyrilex](https://extendsclass.com/regex-tester.html) is a visual regular expression tester.
