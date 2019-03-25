---
layout: tutorial_hands_on

title: RNA Seq Counts to Viz in R
zenodo_link: ""
questions:
- How can I manipulate data using R in Galaxy?
- How can I create neat visualizations of the data?
objectives:
- Know advantages of analyzing data in R
- Know advantages of using RStudio
- Create an RStudio project, and know the benefits of working within a project

- Why use R?
- Why use RStudio and how does it differ from R?
- What are the basic features of the R language?
- What are the most common objects in R?
- How do I get started with tabular data (e.g. spreadsheets) in R?
- What are some best practices for reading data into R?
- How do I save tabular data generated in R?
- How can I manipulate dataframes without repeating myself?
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- bebatut

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

R is one of the most widely-used and powerful programming languages in bioinformatics. R especially shines where a variety of statistical tools are required (e.g. RNA-Seq, population genomics, etc.) and in the generation of publication-quality graphs and figures. Rather than get into an R vs. Python debate (both are useful), keep in mind that many of the concepts you will learn apply to Python and other programming languages.

**This tutorial is significantly based on [the Carpentries](https://carpentries.org) ["Intro to R and RStudio for Genomics"](https://carpentrieslab.github.io/genomics-r-intro/) lesson**

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introducing R and RStudio IDE

In this section we will take you through the very first things you need to get R working.

## A Brief History of R

[R](https://en.wikipedia.org/wiki/R_(programming_language)) has been around since 1995, and was created by Ross Ihaka and Robert Gentleman at the University of Auckland, New Zealand. R is based off the S programming language developed at Bell Labs and was developed to teach intro statistics. See this [slide deck](https://www.stat.auckland.ac.nz/~ihaka/downloads/Massey.pdf) by Ross Ihaka for more info on the subject.

## Advantages of using R

At more than 20 years old, R is fairly mature and [growing in popularity](https://www.tiobe.com/tiobe-index/r/). However, programming isn’t a popularity contest. Here are key advantages of analyzing data in R:

 - **R is [open source](https://en.wikipedia.org/wiki/Open-source_software)**. This means R is free - an advantage if you are at an institution where you have to pay for your own MATLAB or SAS license. Open source, is important to your colleagues in parts of the world where expensive software in inaccessible. It also means that R is actively developed by a community (see [r-project.org](https://www.r-project.org/)), and there are regular updates.
 - **R is widely used**. Ok, maybe programming is a popularity contest. Because, R is used in many areas (not just bioinformatics), you are more likely to find help online when you need it. Chances are, almost any error message you run into, someone else has already experienced.
- **R is powerful**. R runs on multiple platforms (Windows/MacOS/Linux). It can work with much larger datasets than popular spreadsheet programs like Microsoft Excel, and because of its scripting capabilities is far more reproducible. Also, there are thousands of available software packages for science, including genomics and other areas of life science.

## Opening up RStudio

In these lessons, we will be making use of a software called [RStudio](https://www.rstudio.com/products/RStudio/), an [Integrated Development Environment (IDE)](https://en.wikipedia.org/wiki/Integrated_development_environment). RStudio, like most IDEs, provides a graphical interface to R, making it more user-friendly, and providing dozens of useful features. We will introduce additional benefits of using RStudio as you cover the lessons. In this case, we are specifically using [RStudio Server](https://www.rstudio.com/products/RStudio/#Server), a version of RStudio that can be accessed in your web browser. RStudio Server has the same features of the Desktop version of RStudio you could download as standalone software.

***TODO***: *Description on how to open RStudio on Galaxy*

You should now be looking at a page with the RStudio interface:

![rstudio default session](../../images/rna-seq-counts-to-viz-in-r/rstudio_session_default.png "RStudio default session")

## Creating your first R script

Now that we are ready to start exploring R, we will want to keep a record of the commands we are using. To do this we can create an R script:

Click the **File** menu and select **New File** and then **R Script**. Before we go any further, save your script by clicking the save/disk icon that is in the bar above the first line in the script editor, or click the **File** menu and select **save**. In the "Save File" window that opens, name your file **"genomics_r_basics"**. The new script **genomics_r_basics.R** should appear under "files" in the output pane. By convention, R scripts end with the file extension **.R**.

## Overview and customization of the RStudio layout

Here are the major windows (or panes) of the RStudio environment:

![rstudio default session](../../images/rna-seq-counts-to-viz-in-r/rstudio_session_4pane_layout.png "RStudio default session")

- **Source**: This pane is where you will write/view R scripts. Some outputs (such as if you view a dataset using `View()`) will appear as a tab here.
- **Console/Terminal**: This is actually where you see the execution of commands. This is the same display you would see if you were using R at the command line without RStudio. You can work interactively (i.e. enter R commands here), but for the most part we will run a script (or lines in a script) in the source pane and watch their execution and output here. The "Terminal" tab give you access to the BASH terminal (the Linux operating system, unrelated to R).
- **Environment/History**: Here, RStudio will show you what datasets and objects (variables) you have created and which are defined in memory. You can also see some properties of objects/datasets such as their type and dimensions. The "History" tab contains a history of the R commands you've executed R.
- **Files/Plots/Packages/Help**: This multipurpose pane will show you the contents of directories on your computer. You can also use the "Files" tab to navigate and set the working directory. The "Plots" tab will show the output of any plots generated. In "Packages" you will see what packages are actively loaded, or you can attach installed packages. "Help" will display help files for R functions and packages.

> ### {% icon tip %} Tip: Uploads and downloads in the cloud
>
> In the "Files" tab you can select a file and download it from your cloud
> instance (click the "more" button) to your local computer.
> Uploads are also possible.
{: .tip}

All of the panes in RStudio have configuration options. For example, you can minimize/maximize a pane, or by moving your mouse in the space between panes you can resize as needed. The most important customization options for
pane layout are in the **View** menu. Other options such as font sizes, colors/themes, and more are in the **Tools** menu under **Global Options**.

> ### {% icon comment %} You are working with R
> Although we won't be working with R at the terminal, there are lots of reasons
> to. For example, once you have written an RScript, you can run it at any Linux
> or Windows terminal without the need to start up RStudio. We don't want
> you to get confused - RStudio runs R, but R is not RStudio. For more on
> running an R Script at the terminal see this [Software Carpentry lesson](https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/).
{: .comment}

## Using functions in R, without needing to master them

A function in R (or any computing language) is a short program that takes some input and returns some output. The next sections will help you understand what is happening in any R script.

> ### {% icon question %} What do these functions do?
>
> Try the following functions by writing them in your script. See if you can
> guess what they do, and make sure to add comments to your script about your
> assumed purpose.
> - `dir()`
> - `sessionInfo()`
> - `date()`
> - `Sys.time()`
>
> > ### {% icon solution %} Solution
> >
> > - `dir()` # Lists files in the working directory
> > - `sessionInfo()` # Gives the version of R and additional info including
> >    on attached packages
> > - `date()` # Gives the current date
> > - `Sys.time()` # Gives the current time
> >
> > *Notice*: Commands are case sensitive!
> >
> {: .solution}
{: .question}

You have hopefully noticed a pattern - an R function has three key properties:
- Functions have a name (e.g. `dir`, `getwd`); note that functions are case   sensitive!
- Following the name, functions have a pair of `()`
- Inside the parentheses, a function may take 0 or more arguments

An argument may be a specific input for your function and/or may modify the function's behavior. For example the function `round()` will round a number with a decimal:

```console
# This will round a number to the nearest integer
> round(3.14)
[1] 3
```

## Getting help with function arguments

What if you wanted to round to one significant digit? `round()` can do this, but you may first need to read the help to find out how. To see the help (In R sometimes also called a "vignette") enter a `?` in front of the function name:

```console
?round()
```

The **Help** tab will show you information (often, too much information). You will slowly learn how to read and make sense of help files. Checking the **Usage** or **Examples** headings is often a good place to look first. If you look under **Arguments**, we also see what arguments we can pass to this function to modify its behavior. You can also see a function's argument using the `args()` function:

```console
args(round)
function (x, digits = 0)
NULL
```

`round()` takes two arguments, `x`, which is the number to be rounded, and a `digits` argument. The `=` sign indicates that a default (in this case 0) is already set. Since `x` is not set, `round()` requires we provide it, in contrast to `digits` where R will use the default value 0 unless you explicitly provide a different value. We can explicitly set the digits parameter when we call the function:

```console
round(3.14159, digits = 2)
[1] 3.14
```

Or, R accepts what we call "positional arguments", if you pass a function arguments separated by commas, R assumes that they are in the order you saw when we used `args()`. In the case below that means that `x` is 3.14159 and digits is 2.

```console
round(3.14159, 2)
[1] 3.14
```

Finally, what if you are using `?` to get help for a function in a package not installed on your system, such as when you are running a script which has dependencies.

```console
?geom_point()
```

will return an error:

> ### {% icon warning %}
> ```
> Error in .helpForCall(topicExpr, parent.frame()) :
>   no methods for ‘geom_point’ and no documentation for it as a function
> ```
{: .warning}

Use two question marks (i.e. `??geom_point()`) and R will return results from a search of the documentation for packages you have installed on your computer in the **Help** tab. Finally, if you think there should be a function, for example a statistical test, but you aren't sure what it is called in R, or what functions may be available, use the `help.search()` function.

> ### {% icon question %} Searching for R functions
>
> Use `help.search()` to find R functions for the following statistical
> functions. Remember to put your search query in quotes inside the function's
> parentheses.
>
> - Chi-Squared test
> - Student-t test
> - mixed linear model
>
> > ### {% icon solution %} Solution
> >
> >   While your search results may return several tests, we list a few you might
> >   find:
> > - Chi-Squared test: `stats::Chisquare`
> > - Student-t test: `stats::TDist`
> > - mixed linear model: `stats::lm.glm`
> >
> {: .solution}
{: .question}

We will discuss more on where to look for the libraries and packages that contain functions you want to use. For now, be aware that two important ones are [CRAN](https://cran.r-project.org/) - the main repository for R, and [Bioconductor](http://bioconductor.org/) - a popular repository for bioinformatics-related R packages.

## RStudio contextual help

Here is one last bonus we will mention about RStudio. It's difficult to remember all of the arguments and definitions associated with a given function. When you start typing the name of a function and hit the **Tab** key, RStudio will display functions and associated help:

![rstudio default session](../../images/rna-seq-counts-to-viz-in-r/studio_contexthelp1.png "RStudio default session")

Once you type a function, hitting the **Tab** inside the parentheses will show you the function's arguments and provide additional help for each of these arguments.

![rstudio default session](../../images/rna-seq-counts-to-viz-in-r/studio_contexthelp2.png "RStudio default session")



# R Basics

It comes first a description of the step: some background and some theory.
Some image can be added there to support the theory explanation:

The idea is to keep the theory description before quite simple to focus more on the practical part.

***TODO***: *Consider adding a detail box to expand the theory*

> ### {% icon details %} More details about the theory
>
> But to describe more details, it is possible to use the detail boxes which are expandable
>
{: .details}

A big step can have several subsections or sub steps:


## Sub-step with **My Tool**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **My Tool** {% icon tool %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: File
>    - *"Parameter"*: `a value`
>
>    ***TODO***: *Check parameter descriptions*
>
>    ***TODO***: *Consider adding a comment or tip box*
>
>    > ### {% icon comment %} Comment
>    >
>    > A comment about the tool or something else. This box can also be in the main text
>    {: .comment}
>
{: .hands_on}

***TODO***: *Consider adding a question to test the learners understanding of the previous exercise*

> ### {% icon question %} Questions
>
> 1. Question1?
> 2. Question2?
>
> > ### {% icon solution %} Solution
> >
> > 1. Answer for question1
> > 2. Answer for question2
> >
> {: .solution}
>
{: .question}


## Re-arrange

To create the template, each step of the workflow had its own subsection.

***TODO***: *Re-arrange the generated subsections into sections or other subsections.
Consider merging some hands-on boxes to have a meaningful flow of the analyses*

# R Basics continued - factors and data frames



# Aggregating and Analyzing Data with dplyr



# Data Visualization with ggplot2


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
