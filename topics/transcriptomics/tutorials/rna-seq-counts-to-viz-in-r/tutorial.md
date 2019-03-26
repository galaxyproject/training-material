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

## "The fantastic world of R awaits you" OR "Nobody wants to learn how to use R"

Before we begin this lesson, we want you to be clear on the goal of the workshop and these lessons. We believe that every learner can **achieve competency with R**. You have reached competency when you find that you are able to **use R to handle common analysis challenges in a reasonable amount of time** (which includes time needed to look at learning materials, search for answers online, and ask colleagues for help). As you spend more time using R (there is no substitute for regular use and practice) you will find yourself gaining competency and even expertise. The more familiar you get, the more complex the analyses you will be able to carry out, with less frustration, and in less time - the fantastic world of R awaits you!

## What these lessons will not teach you

Nobody wants to learn how to use R. People want to learn how to use R to analyze their own research questions! Ok, maybe some folks learn R for R's sake, but these lessons assume that you want to start analyzing genomic data as soon as possible. Given this, there are many valuable pieces of information about R that we simply won't have time to cover. Hopefully, we will clear the hurdle of giving you just enough knowledge to be dangerous, which can be a high bar in R! We suggest you look into the additional learning materials in the tip box below.

**Here are some R skills we will *not* cover in these lessons**

- How to create and work with R matrices and R lists
- How to create and work with loops and conditional statements, and the "apply"
  family of functions (which are super useful, read more [here](https://www.r-bloggers.com/r-tutorial-on-the-apply-family-of-functions/))
- How to do basic string manipulations (e.g. finding patterns in text using grep, replacing text)
- How to plot using the default R graphic tools (we *will* cover plot creation, but will do so using the popular plotting package `ggplot2`)
- How to use advanced R statistical functions

> ### {% icon tip %} Tip: Where to learn more
>
> The following are good resources for learning more about R. Some of them
> can be quite technical, but if you are a regular R user you may ultimately
> need this technical knowledge.
> - [R for Beginners](https://cran.r-project.org/doc/contrib/Paradis-rdebuts_en.pdf):
    By Emmanuel Paradis and a great starting point
> - [The R Manuals](https://cran.r-project.org/manuals.html): Maintained by the
    R project
> - [R contributed documentation](https://cran.r-project.org/other-docs.html):
    Also linked to the R project; importantly there are materials available in
    several languages
> - [R for Data Science](http://r4ds.had.co.nz/): A wonderful collection by
    noted R educators and developers Garrett Grolemund and Hadley Wickham
> - [Practical Data Science for Stats](https://peerj.com/collections/50-practicaldatascistats/):
    Not exclusively about R usage, but a nice collection of pre-prints on data science
    and applications for R
> - [Programming in R Software Carpentry lesson](https://software-carpentry.org/lessons/):
    There are several Software Carpentry lessons in R to choose from
> - [Data Camp Introduction to R](https://www.datacamp.com/courses/free-introduction-to-r):
    This is a fun online learning platform for Data Science, including R.
{: .tip}

## Creating objects in R

> ### {% icon comment %} Reminder
> At this point you should be coding along in the "**genomics_r_basics.R**"
> script we created in the last episode. Writing your commands in the script
> (and commenting it) will make it easier to record what you did and why.
{: .comment}

What might be called a variable in many languages is called an **object** in R.

**To create an object you need:**

- a name (e.g. 'a')
- a value (e.g. '1')
- the assignment operator ('<-')

In your script, "**genomics_r_basics.R**", using the R assignment operator '<-', assign '1' to the object 'a' as shown. Remember to leave a comment in the line above (using the '#') to explain what you are doing:

```console
# this line creates the object 'a' and assigns it the value '1'
a <- 1
```

Next, run this line of code in your script. You can run a line of code by hitting the **Run** button that is just above the first line of your script in the header of the Source pane or you can use the appropriate shortcut:

- Windows execution shortcut: <KBD>Ctrl</KBD>+<KBD>Enter</KBD>
- Mac execution shortcut: <KBD>Cmd(⌘)</KBD>+<KBD>Enter</KBD>

To run multiple lines of code, you can highlight all the line you wish to run and then hit **Run** or use the shortcut key combo listed above.

In the RStudio 'Console' you should see:

```console
a <- 1
>
```

The 'Console' will display lines of code run from a script and any outputs or status/warning/error messages (usually in red).

In the 'Environment' window you will also get a table:

|Values||
|------|-|
|a|1|

The 'Environment' window allows you to keep track of the objects you have created in R.

> ### {% icon question %} Exercise: Create some objects in R
>
> Create the following objects; give each object an appropriate name
> (your best guess at what name to use is fine):
>
> 1. Create an object that has the value of number of pairs of human chromosomes
> 2. Create an object that has a value of your favorite gene name
> 3. Create an object that has this URL as its value: "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-39/fasta/bacteria_5_collection/escherichia_coli_b_str_rel606/"
> 4. Create an object that has the value of the number of chromosomes in a
> diploid human cell
>
> > ### {% icon solution %} Solution
> >
> > Here as some possible answers to the challenge:
> > ```{r, purl = FALSE}
> > human_chr_number <- 23
> > gene_name <- 'pten'
> > ensemble_url <- 'ftp://ftp.ensemblgenomes.org/pub/bacteria/release-39/fasta/bacteria_5_collection/escherichia_coli_b_str_rel606/'
> > human_diploid_chr_num <-  2 * human_chr_number
> > ```
> >
> {: .solution}
{: .question}


## Naming objects in R

Here are some important details about naming objects in R.

- **Avoid spaces and special characters**: Object names cannot contain spaces or the minus sign (`-`). You can use '_' to make names more readable. You should avoid using special characters in your object name (e.g. ! @ # . , etc.). Also, object names cannot begin with a number.
- **Use short, easy-to-understand names**: You should avoid naming your objects using single letters (e.g. 'n', 'p', etc.). This is mostly to encourage you to use names that would make sense to anyone reading your code (a colleague, or even yourself a year from now). Also, avoiding excessively long names will make your code more readable.
- **Avoid commonly used names**: There are several names that may already have a definition in the R language (e.g. 'mean', 'min', 'max'). One clue that a name already has meaning is that if you start typing a name in RStudio and it gets a colored highlight or RStudio gives you a suggested autocompletion you have chosen a name that has a reserved meaning.
- **Use the recommended assignment operator**: In R, we use '<- ' as the preferred assignment operator. '=' works too, but is most commonly used in passing arguments to functions (more on functions later). There is a shortcut for the R assignment operator:
  - Windows execution shortcut: <KBD>Alt</KBD>+<KBD>-</KBD>
  - Mac execution shortcut: <KBD>Option</KBD>+<KBD>-</KBD>
There are a few more suggestions about naming and style you may want to learn more about as you write more R code. There are several "style guides" that have advice, and one to start with is the [tidyverse R style guide](http://style.tidyverse.org/index.html).

> ### {% icon tip %} Tip: Pay attention to warnings in the script console
>
> If you enter a line of code in your script that contains an error, RStudio
> may give you an error message and underline this mistake. Sometimes these
> messages are easy to understand, but often the messages may need some figuring
> out. Paying attention to these warnings will help you avoid mistakes. In the example below, our object name has a space, which
> is not allowed in R. The error message does not say this directly,
> but R is "not sure"
> about how to assign the name to "human_ chr_number" when the object name we want is "human_chr_number".
>
> ![rstudio script warning](../../images/rna-seq-counts-to-viz-in-r/rstudio_script_warning.png "rstudio script warning")
{: .tip}

## Reassigning object names or deleting objects

Once an object has a value, you can change that value by overwriting it. R will not give you a warning or error if you overwriting an object, which may or may not be a good thing depending on how you look at it.

```console
# gene_name has the value 'pten' or whatever value you used in the challenge.
# We will now assign the new value 'tp53'
gene_name <- 'tp53'
```

You can also remove an object from R's memory entirely. The `rm()` function will delete the object.

```console
# delete the object 'gene_name'
rm(gene_name)
```

If you run a line of code that has only an object name, R will normally display the contents of that object. In this case, we are told the object no longer exists.

> ### {% icon warning %} Danger: You can lose data!
> ```
> Error: object 'gene_name' not found
> ```
{: .warning}

## Understanding object data types (modes)

In R, **every object has two properties**:

- **Length**: How many distinct values are held in that object
- **Mode**: What is the classification (type) of that object.

We will get to the "length" property later in the lesson. The **"mode" property** **corresponds to the type of data an object represents**. The most common modes you will encounter in R are:

|Mode (abbreviation)|Type of data|
|----|------------|
|Numeric (num)| Numbers such floating point/decimals (1.0, 0.5, 3.14), there are also more specific numeric types (dbl - Double, int - Integer). These differences are not relevant for most beginners and pertain to how these values are stored in memory |
|Character (chr)|A sequence of letters/numbers in single '' or double " " quotes|
|Logical| Boolean values - TRUE or FALSE|

There are a few other modes (i.e. "complex", "raw" etc.) but these are the three we will work with in this lesson.

Data types are familiar in many programming languages, but also in natural language where we refer to them as the parts of speech, e.g. nouns, verbs, adverbs, etc. Once you know if a word - perhaps an unfamiliar one - is a noun, you can probably guess you can count it and make it plural if there is more than one (e.g. 1 [Tuatara](https://en.wikipedia.org/wiki/Tuatara), or 2 Tuataras). If something is a adjective, you can usually change it into an adverb by adding "-ly" (e.g. [jejune](https://www.merriam-webster.com/dictionary/jejune) vs. jejunely). Depending on the context, you may need to decide if a word is in one category or another (e.g "cut" may be a noun when it's on your finger, or a verb when you are preparing vegetables). These concepts have important analogies when working with R objects.

> ### {% icon question %} QuestionsExercise: Create objects and check their modes
>
> Create the following objects in R, then use the `mode()` function to verify
> their modes. Try to guess what the mode will be before you look at the solution
>
> 1. `chromosome_name <- 'chr02'`
> 2. `od_600_value <- 0.47`
> 3. `chr_position <- '1001701'`
> 4. `spock <- TRUE`
> 5. `pilot <- Earhart`
>
> > ### {% icon solution %} Solution
> >
> > ```console
> > chromosome_name <- 'chr02'
> > od_600_value <- 0.47
> > chr_position <- '1001701'
> > spock <- TRUE
> > pilot <- Earhart
> > [1] Error in eval(expr, envir, enclos): object 'Earhart' not found

> > ```
> >
> > ```console
> > mode(chromosome_name)
> > [1] "character"
> >
> > mode(od_600_value)
> > [1] "numeric"
> >
> > mode(chr_position)
> > [1] "character"
> >
> > mode(spock)
> > [1] "logical"
> >
> > mode(pilot)
> > Error in mode(pilot): object 'pilot' not found
> >
> > ```
> >
> {: .solution}
{: .question}

Notice from the solution that even if a series of numbers are given as a value R will consider them to be in the "character" mode if they are enclosed as single or double quotes. Also, notice that you cannot take a string of alphanumeric characters (e.g. Earhart) and assign as a value for an object. In this case, R looks for an object named `Earhart` but since there is no object, no assignment can be made. If `Earhart` did exist, then the mode of `pilot` would be whatever the mode of `Earhart` was originally. If we want to create an object called `pilot` that was the **name** "Earhart", we need to enclose `Earhart` in quotation marks.

```console
pilot <- "Earhart"
mode(pilot)

[1] "character"
```

## Mathematical and functional operations on objects

Once an object exists (which by definition also means it has a mode), R can appropriately manipulate that object. For example, objects of the numeric modes can be added, multiplied, divided, etc. R provides several mathematical (arithmetic) operators including:

|Operator|Description|
|--------|-----------|
|+|addition|
|-|subtraction|
|*|multiplication|
|/|division|
|^ or **|exponentiation|
|a%%b|modulus (returns the remainder after division)|

These can be used with literal numbers:

```console
> (1 + (5 ** 0.5))/2
[1] 1.618034
```

and importantly, can be used on any object that evaluates to (i.e. interpreted by R) a numeric object:

```console
> human_chr_number <- 23

# multiply the object 'human_chr_number' by 2
> human_chr_number * 2
[1] 46
```

> ### {% icon question %} Exercise: Compute the golden ratio
>
> One approximation of the golden ratio (φ) can be found by taking the sum of 1
> and the square root of 5, and dividing by 2 as in the example above. Compute
> the golden ratio to 3 digits of precision using the `sqrt()` and `round()`
> functions. Hint: remember the `round()` function can take 2 arguments.
>
> > ### {% icon solution %} Solution
> >
> > ```console
> > round((1 + sqrt(5))/2, digits = 3)
> > [1] 1.618
> > ```
>> Notice that you can place one function inside of another.
> >
> {: .solution}
{: .question}

## Vectors

Vectors are probably the most used commonly used object type in R.
**A vector is a collection of values that are all of the same type (numbers, characters, etc.)**.
One of the most common ways to create a vector is to use the `c()` function - the "concatenate" or "combine" function. Inside the function you may enter one or more values; for multiple values, separate each value with a comma:

```console
# Create the SNP gene name vector
snp_genes <- c("OXTR", "ACTN3", "AR", "OPRM1")
```

Vectors always have a **mode** and a **length**. You can check these with the `mode()` and `length()` functions respectively. Another useful function that gives both of these pieces of information is the `str()` (structure) function.

```console
# Check the mode, length, and structure of 'snp_genes'
> mode(snp_genes)
[1] "character"

> length(snp_genes)
[1] 4

> str(snp_genes)
chr [1:4] "OXTR" "ACTN3" "AR" "OPRM1"
```

Vectors are quite important in R. Another data type that we will work with later in this lesson, data frames, are collections of vectors. What we learn here about vectors will pay off even more when we start working with data frames.

## Creating and subsetting vectors

Let's create a few more vectors to play around with:

```console
# Some interesting human SNPs
# while accuracy is important, typos in the data won't hurt you here
> snps <- c('rs53576', 'rs1815739', 'rs6152', 'rs1799971')
> snp_chromosomes <- c('3', '11', 'X', '6')
> snp_positions <- c(8762685, 66560624, 67545785, 154039662)
```

Once we have vectors, one thing we may want to do is specifically retrieve one or more values from our vector. To do so, we use **bracket notation**. We type the name of the vector followed by square brackets. In those square brackets we place the index (e.g. a number) in that bracket as follows:

```console
# get the 3rd value in the snp_genes vector
> snp_genes[3]
[1] "AR"
```

In R, every item your vector is indexed, starting from the first item (1) through to the final number of items in your vector. You can also retrieve a range of numbers:

```console
# get the 1st through 3rd value in the snp_genes vector
> snp_genes[1:3]
[1] "OXTR"  "ACTN3" "AR"  
```

If you want to retrieve several (but not necessarily sequential) items from a vector, you pass a **vector of indices**; a vector that has the numbered positions you wish to retrieve.

```console
# get the 1st, 3rd, and 4th value in the snp_genes vector
> snp_genes[c(1, 3, 4)]
[1] "OXTR"  "AR"    "OPRM1"
```

There are additional (and perhaps less commonly used) ways of subsetting a vector (see [these
examples](https://thomasleeper.com/Rcourse/Tutorials/vectorindexing.html)). Also, several of these subsetting expressions can be combined:

```console
# get the 1st through the 3rd value, and 4th value in the snp_genes vector
# yes, this is a little silly in a vector of only 4 values.
> snp_genes[c(1:3,4)]
[1] "OXTR"  "ACTN3" "AR"    "OPRM1"
```

## Adding to, removing, or replacing values in existing vectors

Once you have an existing vector, you may want to add a new item to it. To do so, you can use the `c()` function again to add your new value:

```console
# add the gene 'CYP1A1' and 'APOA5' to our list of snp genes
# this overwrites our existing vector
> snp_genes <- c(snp_genes, "CYP1A1", "APOA5")
```

We can verify that "snp_genes" contains the new gene entry

```console
> snp_genes
[1] "OXTR"   "ACTN3"  "AR"     "OPRM1"  "CYP1A1" "APOA5"
```

Using a negative index will return a version a vector with that index's value removed:

```console
> snp_genes[-6]
[1] "OXTR"   "ACTN3"  "AR"     "OPRM1"  "CYP1A1"
```

We can remove that value from our vector by overwriting it with this expression:

```console
> snp_genes <- snp_genes[-6]
> snp_genes
[1] "OXTR"   "ACTN3"  "AR"     "OPRM1"  "CYP1A1"
```

We can also explicitly rename or add a value to our index using double bracket notation:

```console
> snp_genes[[7]]<- "APOA5"
> snp_genes
[1] "OXTR"   "ACTN3"  "AR"     "OPRM1"  "CYP1A1" NA       "APOA5"
```

Notice in the operation above that R inserts an `NA` value to extend our vector so that the gene "APOA5" is an index 7. This may be a good or not-so-good thing depending on how you use this.

> ### {% icon question %} Exercise: Examining and subsetting vectors
>
> Answer the following questions to test your knowledge of vectors
>
> Which of the following are true of vectors in R?
> A) All vectors have a mode **or** a length  
> B) All vectors have a mode **and** a length  
> C) Vectors may have different lengths  
> D) Items within a vector may be of different modes  
> E) You can use the `c()` to one or more items to an existing vector  
> F) You can use the `c()` to add a vector to an exiting vector  
>
> > ### {% icon solution %} Solution
> >
> > A) False - Vectors have both of these properties  
> > B) True  
> > C) True  
> > D) False - Vectors have only one mode (e.g. numeric, character); all items in  
> > a vector must be of this mode.
> > E) True  
> > F) True  
> >
> {: .solution}
{: .question}

## Logical Subsetting

There is one last set of cool subsetting capabilities we want to introduce. It is possible within R to retrieve items in a vector based on a logical evaluation or numerical comparison. For example, let's say we wanted get all of the SNPs in our vector of SNP positions that were greater than 100,000,000. We could index using the '>' (greater than) logical operator:

```console
> snp_positions[snp_positions > 100000000]
[1] 154039662
```

In the square brackets you place the name of the vector followed by the comparison operator and (in this case) a numeric value. Some of the most common logical operators you will use in R are:

  | Operator | Description              |
  |----------|--------------------------|
  | <        | less than                |
  | <=       | less than or equal to    |
  | >        | greater than             |
  | >=       | greater than or equal to |
  | ==       | exactly equal to         |
  | !=       | not equal to             |
  | !x       | not x                    |
  | a \| b   | a or b                   |
  | a & b    | a and b                  |


> ### {% icon comment %} The magic of programming
>
> The reason why the expression `snp_positions[snp_positions > 100000000]` works
> can be better understood if you examine what the expression "snp_positions > 100000000"
>evaluates to:
>
> ```console
> > snp_positions > 100000000
> [1] FALSE FALSE FALSE  TRUE
> ```
>
> The output above is a logical vector, the 4th element of which is TRUE. When
> you pass a logical vector as an index, R will return the true values:
>
> ```console
> > snp_positions[c(FALSE, FALSE, FALSE, TRUE)]
> [1] 154039662
> ```
>
> If you have never coded before, this type of situation starts to expose the
> "magic" of programming. We mentioned before that in the bracket notation you
> take your named vector followed by brackets which contain an index:
> **named_vector[index]**. The "magic" is that the index needs to *evaluate to*
> a number. So, even if it does not appear to be an integer (e.g. 1, 2, 3), as
> long as R can evaluate it, we will get a result. That our expression
> `snp_positions[snp_positions > 100000000]` evaluates to a number can be seen
> in the following situation. If you wanted to know which **index** (1, 2, 3, or
> 4) in our vector of SNP positions was the one that was greater than
> 100,000,000?
>
> We can use the `which()` function to return the indices of any item that
> evaluates as TRUE in our comparison:
>
> ```console
> > which(snp_positions > 100000000)
> [1] 4
> ```
>
> **Why this is important**
>
> Often in programming we will not know what inputs
> and values will be used when our code is executed. Rather than put in a
> pre-determined value (e.g 100000000) we can use an object that can take on
> whatever value we need. So for example:
>
> ```console
> > snp_marker_cutoff <- 100000000
> > snp_positions[snp_positions > snp_marker_cutoff]
> [1] 154039662
> ```
>
> Ultimately, it's putting together flexible, reusable code like this that gets
> at the "magic" of programming!
{: .comment}

## A few final vector tricks

Finally, there are a few other common retrieve or replace operations you may want to know about. First, you can check to see if any of the values of your vector are missing (i.e. are `NA`). Missing data will get a more detailed treatment later, but the `is.NA()` function will return a logical vector, with TRUE for any NA value:

```console
# current value of 'snp_genes':
# chr [1:7] "OXTR" "ACTN3" "AR" "OPRM1" "CYP1A1" NA "APOA5"
> is.na(snp_genes)
[1] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
```

Sometimes, you may wish to find out if a specific value (or several values) is present a vector. You can do this using the comparison operator `%in%`, which will return TRUE for any value in your collection that is in the vector you are searching:

```console
# current value of 'snp_genes':
# chr [1:7] "OXTR" "ACTN3" "AR" "OPRM1" "CYP1A1" NA "APOA5"
# test to see if "ACTN3" or "APO5A" is in the snp_genes vector
# if you are looking for more than one value, you must pass this as a vector
> c("ACTN3","APOA5") %in% snp_genes
[1] TRUE TRUE
```

> ### {% icon question %} Review Exercise 1
>
> What data types/modes are the following vectors?
>    a. `snps`  
>    b. `snp_chromosomes`  
>    c. `snp_positions`
>
> > ### {% icon solution %} Solution
> >
> > ```console
> > > typeof(snps)
> > [1] "character"
> > > typeof(snp_chromosomes)
> > [1] "character"
> > > typeof(snp_positions)
> > [1] "double"
> > ```
> >
> {: .solution}
{: .question}


> ### {% icon question %} Review Exercise 2
>
> Add the following values to the specified vectors:
>    a. To the `snps` vector add: 'rs662799'  
>    b. To the `snp_chromosomes` vector add: 11  
>    c. To the `snp_positions` vector add: 	116792991
>
> > ### {% icon solution %} Solution
> >
> > ```console
> > > snps <- c(snps, 'rs662799')
> > > snps
> > [1] "rs53576"   "rs1815739" "rs6152"    "rs1799971" "rs662799"
> > > snp_chromosomes <- c(snp_chromosomes, "11") # did you use quotes?
> > > snp_chromosomes
> > [1] "3"  "11" "X"  "6"  "11"
> > > snp_positions <- c(snp_positions, 116792991)
> > > snp_positions
> > [1]   8762685  66560624  67545785 154039662 116792991
> > ```
> >
> {: .solution}
{: .question}



> ### {% icon question %} Review Exercise 3
>
> Make the following change to the `snp_genes` vector:
>
> Hint: Your vector should look like this in 'Environment':
> `chr [1:7] "OXTR" "ACTN3" "AR" "OPRM1" "CYP1A1" NA "APOA5"`. If not
> recreate the vector by running this expression:
> `snp_genes <- c("OXTR", "ACTN3", "AR", "OPRM1", "CYP1A1", NA, "APOA5")`
>
>    a. Create a new version of `snp_genes` that does not contain CYP1A1 and then  
>    b. Add 2 NA values to the end of `snp_genes`
>
> > ### {% icon solution %} Solution
> >
> > ```console
> > > snp_genes <- snp_genes[-5]
> > > snp_genes <- c(snp_genes, NA, NA)
> > > snp_genes
> > [1] "OXTR"  "ACTN3" "AR"    "OPRM1" NA      "APOA5" NA      NA    
> > ```
> >
> {: .solution}
{: .question}


> ### {% icon question %} Review Exercise 4
>
> Using indexing, create a new vector named `combined` that contains:
>    - The the 1st value in `snp_genes`
>    - The 1st value in `snps`
>    - The 1st value in `snp_chromosomes`
>    - The 1st value in `snp_positions`
>
> > ### {% icon solution %} Solution
> >
> > ```console
> > > combined <- c(snp_genes[1], snps[1], snp_chromosomes[1], snp_positions[1])
> > > combined
> > [1] "OXTR"    "rs53576" "3"       "8762685"
> > ```
> >
> {: .solution}
{: .question}


> ### {% icon question %} Review Exercise 5
>
> What type of data is `combined`?
>
> > ### {% icon solution %} Solution
> >
> > ```console
> > > typeof(combined)
> > [1] "character"
> > ```
> >
> {: .solution}
{: .question}


## Bonus material: Lists

Lists are quite useful in R, but we won't be using them in the genomics lessons. That said, you may come across lists in the way that some bioinformatics programs may store and/or return data to you. One of the key attributes of a list is that, unlike a vector, a list may contain data of more than one mode. Learn more about creating and using lists using this [nice tutorial](http://r4ds.had.co.nz/lists.html). In this one example, we will create a named list and show you how to retrieve items from the list.

```console
# Create a named list using the 'list' function and our SNP examples
# Note, for easy reading we have placed each item in the list on a separate line
# Nothing special about this, you can do this for any multiline commands
# To run this command, make sure the entire command (all 4 lines) are highlighted
# before running
# Note also, as we are doing all this inside the list() function use of the
# '=' sign is good style
> snp_data <- list(genes = snp_genes,
                 refference_snp = snps,
                 chromosome = snp_chromosomes,
                 position = snp_positions)
# Examine the structure of the list
> str(snp_data)

List of 4
 $ genes         : chr [1:8] "OXTR" "ACTN3" "AR" "OPRM1" ...
 $ refference_snp: chr [1:5] "rs53576" "rs1815739" "rs6152" "rs1799971" ...
 $ chromosome    : chr [1:5] "3" "11" "X" "6" ...
 $ position      : num [1:5] 8.76e+06 6.66e+07 6.75e+07 1.54e+08 1.17e+08

```

To get all the values for the `position` object in the list, we use the `$` notation:

```console
# return all the values of position object
> snp_data$position
[1]   8762685  66560624  67545785 154039662 116792991
```

To get the first value in the `position` object, use the `[]` notation to index:

```console
# return first value of the position object
> snp_data$position[1]
[1] 8762685
```


# R Basics continued - factors and data frames



# Aggregating and Analyzing Data with dplyr



# Data Visualization with ggplot2


# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
