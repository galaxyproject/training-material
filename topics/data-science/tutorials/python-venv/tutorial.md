---
layout: tutorial_hands_on

title: "Virtual Environments For Software Development"

level: Intermediate
requirements: []
follow_up_training: []
time_estimation:  30M

questions:
- "What are virtual environments in software development and why you should use them?"
- "How can we manage Python virtual environments and external (third-party) libraries?"
objectives:
- "Set up a Python virtual environment for our software project using `venv` and `pip`."
- "Run our software from the command line."

key_points:
- "Virtual environments keep Python versions and dependencies required by different projects separate."
- "A virtual environment is itself a directory structure."
- "Use `venv` to create and manage Python virtual environments."
- "Use `pip` to install and manage Python external (third-party) libraries."
- "`pip` allows you to declare all dependencies for a project in a separate
file (by convention called `requirements.txt`) which can be shared with collaborators/users and used to replicate a virtual environment."
- "Use `pip3 freeze > requirements.txt` to take snapshot of your project's dependencies."
- "Use `pip3 install -r requirements.txt` to replicate someone else's virtual environment on your machine from
the `requirements.txt` file."

subtopic: python-modular
contributions:
  authorship:
  - carpentries
  - hexylena
  funding:
  - carpentries
  - avans-atgm

priority: 10
notebook:
  language: bash 

---

"Virtual Environments" allow you to easily manage your installed Python packages and prevent conflicts between different project's dependencies. In general most modern projects should use `conda` for dependency management, but `venv` can be convenient for Python-only projects.

> <comment-title></comment-title>
>
> This tutorial is significantly based on [the Carpentries](https://carpentries.org) lesson ["Intermediate Research Software Development"](https://carpentries-incubator.github.io/python-intermediate-development/).
>
{: .comment}


If you have a python project you are using, you will often see something like
following two lines somewhere at the top.

```python
from matplotlib import pyplot as plt
import numpy as np
```

This means that our code requires two *external libraries* (also called third-party packages or dependencies) -
`numpy` and `matplotlib`.
Python applications often use external libraries that don’t come as part of the standard Python distribution. This means
that you will have to use a *package manager* tool to install them on your system.
Applications will also sometimes need a
specific version of an external library (e.g. because they require that a particular
bug has been fixed in a newer version of the library), or a specific version of Python interpreter.
This means that each Python application you work with may require a different setup and a set of dependencies so it
is important to be able to keep these configurations separate to avoid confusion between projects.
The solution for this problem is to create a self-contained *virtual
environment* per project, which contains a particular version of Python installation plus a number of
additional external libraries.

Virtual environments are not just a feature of Python - all modern programming languages use them to isolate code
of a specific project and make it easier to develop, run, test and share code with others. In this tutorial, we learn how
to set up a virtual environment to develop our code and manage our external dependencies.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Virtual Environments

So what exactly are virtual environments, and why use them?

A Python virtual environment is an **isolated working copy** of a specific version of
Python interpreter together with specific versions of a number of external libraries installed into that
virtual environment. A virtual environment is simply a *directory with a particular
structure* which includes links to and enables multiple side-by-side installations of
different Python interpreters or different versions of the same external library to coexist on your machine and only one to be selected for each of our projects. This allows you to work on a particular
project without worrying about affecting other projects on your machine.

As more external libraries are added to your Python project over time, you can add them to
its specific virtual environment and avoid a great deal of confusion by having separate (smaller) virtual environments
for each project rather than one huge global environment with potential package version clashes. Another big motivator
for using virtual environments is that they make sharing your code with others much easier (as we will see shortly).
Here are some typical scenarios where the usage of virtual environments is highly recommended (almost unavoidable):

- You have an older project that only works under Python 2. You do not have the time to migrate the project to Python 3
  or it may not even be possible as some of the third party dependencies are not available under Python 3. You have to
  start another project under Python 3. The best way to do this on a single machine is to set up two separate Python virtual
  environments.
- One of your Python 3 projects is locked to use a particular older version of a third party dependency. You cannot use the
  latest version of the
  dependency as it breaks things in your project. In a separate branch of your project, you want to try and fix problems
  introduced by the new version of the dependency without affecting the working version of your project. You need to set up
  a separate virtual environment for your branch to 'isolate' your code while testing the new feature.

You do not have to worry too much about specific versions of external libraries that your project depends on most of the time.
Virtual environments enable you to always use the latest available version without specifying it explicitly.
They also enable you to use a specific older version of a package for your project, should you need to.

> <tip-title>A Specific Python or Package Version is Only Ever Installed Once</tip-title>
> Note that you will not have a separate Python or package installations for each of your projects - they will only
ever be installed once on your system but will be referenced from different virtual environments.
{: .tip}

### Managing Python Virtual Environments

There are several commonly used command line tools for managing Python virtual environments:
- `venv`, available by default from the standard `Python` distribution from `Python 3.3+`
- `virtualenv`, needs to be installed separately but supports both `Python 2.7+` and `Python 3.3+`versions
- `pipenv`, created to fix certain shortcomings of `virtualenv`
- `conda`, package and environment management system (also included as part of the Anaconda Python distribution often used by the scientific community)
- `poetry`, a modern Python packaging tool which handles virtual environments automatically

While there are pros and cons for using each of the above, all will do the job of managing Python
virtual environments for you and it may be a matter of personal preference which one you go for.
In this course, we will use `venv` to create and manage our
virtual environment (which is the preferred way for Python 3.3+).

Until you encounter the needs of a project which goes beyond what is available
in the Python ecosystem, e.g. when you depend on external packages like htslib
or bioinformatics tools that are simply not distributed as part of PyPI, then
`venv` is a good choice to get started with.


### Managing Python Packages

Part of managing your (virtual) working environment involves installing, updating and removing external packages
on your system. The Python package manager tool `pip` is most commonly used for this - it interacts
 and obtains the packages from the central repository called [Python Package Index (PyPI)](https://pypi.org/).
`pip` can now be used with all Python distributions (including Anaconda).

> <tip-title>A Note on Anaconda and `conda`</tip-title>
> Anaconda is an open source Python
> distribution commonly used for scientific programming - it conveniently installs Python, package and environment management `conda`, and a 
> number of commonly used scientific computing packages so you do not have to obtain them separately. 
> `conda` is an independent command line tool (available separately from the Anaconda distribution too) with dual functionality: (1) it is a package manager that helps you find Python packages from
> remote package repositories and install them on your system, and (2) it is also a virtual environment manager. So, you can use `conda` for both tasks instead of using `venv` and `pip`.
{: .tip}

### Many Tools for the Job

Installing and managing Python distributions, external libraries and virtual environments is, well,
complex. There is an abundance of tools for each task, each with its advantages and disadvantages, and there are different
ways to achieve the same effect (and even different ways to install the same tool!).
Note that each Python distribution comes with its own version of
`pip` - and if you have several Python versions installed you have to be extra careful to use the correct `pip` to
manage external packages for that Python version.

`venv` and `pip` are considered the *de facto* standards for virtual environment and package management for Python 3.
However, the advantages of using Anaconda and `conda` are that you get (most of the) packages needed for
scientific code development included with the distribution. If you are only collaborating with others who are also using
Anaconda, you may find that `conda` satisfies all your needs. It is good, however, to be aware of all these tools,
and use them accordingly. As you become more familiar with them you will realise that equivalent tools work in a similar
way even though the command syntax may be different (and that there are equivalent tools for other programming languages
too to which your knowledge can be ported).

![Python environment hell XKCD comic  showing boxes like pip, easy_install, homebrew 2.7, anaconda, homebrew 3.6, /usr/local/Cellar, ~/python/, and a chaotic mess of arrows moving between them all. At the bottom is the text: My python environment has become so degraded that my laptop has been declared a superfund site. (A superfund site is generally an environmental disaster area.)](../../images/xkcd/python_environment.png "Python Environment Hell from XKCD 1987 (CC-BY-NC 2.5)")

Let us have a look at how we can create and manage virtual environments from the command line using `venv` and manage packages using `pip`.

### Creating a `venv` Environment
Creating a virtual environment with `venv` is done by executing the following command:

```bash
python3 -m venv /path/to/new/virtual/environment
```

where `/path/to/new/virtual/environment` is a path to a directory where you want to place it - conventionally within
your software project so they are co-located.
This will create the target directory for the virtual environment (and any parent directories that don’t exist already).

For our project, let's create a virtual environment called `venv` off the project root:

```bash
python3 -m venv venv
```

If you list the contents of the newly created `venv` directory, on a Mac or Linux system 
(slightly different on Windows as explained below) you should see something like:

```bash
ls -l venv
```

So, running the `python3 -m venv venv` command created the target directory called `venv`
containing:

- `pyvenv.cfg` configuration file with a home key pointing to the Python installation from which the command was run,
- `bin` subdirectory (called `Scripts` on Windows) containing a symlink of the Python interpreter binary used to create the
environment and the standard Python library,
- `lib/pythonX.Y/site-packages` subdirectory (called `Lib\site-packages` on Windows) to contain its own independent set of installed Python packages isolated from other projects,
- various other configuration and supporting files and subdirectories.

> <tip-title>Naming Virtual Environments</tip-title>
> What is a good name to use for a virtual environment? Using "venv" or ".venv" as the
name for an environment and storing it within the project's directory seems to be the recommended way -
this way when you come across such a subdirectory within a software project,
by convention you know it contains its virtual environment details.
A slight downside is that all different virtual environments
on your machine then use the same name and the current one is determined by the context of the path
you are currently located in. A (non-conventional) alternative is to
use your project name for the name of the virtual environment, with the downside that there is nothing to indicate
that such a directory contains a virtual environment. In our case, we have settled to use the name "venv" since it is 
not a hidden directory and we want it to be displayed by the command line when listing directory contents (hence, 
no need for the "." in its name that would, by convention, make it hidden). In the future,
you will decide what naming convention works best for you. Here are some references for each of the naming conventions:
- [The Hitchhiker's Guide to Python](https://docs.python-guide.org/dev/virtualenvs/) notes that "venv" is the general convention used globally
- [The Python Documentation](https://docs.python.org/3/library/venv.html) indicates that ".venv" is common
- ["venv" vs ".venv" discussion](https://discuss.python.org/t/trying-to-come-up-with-a-default-directory-name-for-virtual-environments/3750)
{: .tip}

Once you’ve created a virtual environment, you will need to activate it:

```bash
source venv/bin/activate
```

Activating the virtual environment will change your command line’s prompt to show what virtual environment
you are currently using (indicated by its name in round brackets at the start of the prompt),
and modify the environment so that running Python will get you the particular
version of Python configured in your virtual environment.

You can verify you are using your virtual environment's version of Python by checking the path using `which`:
```bash
which python3
```

When you’re done working on your project, you can exit the environment with:

```bash
deactivate
```

If you've just done the `deactivate`, ensure you reactivate the environment ready for the next part:

```bash
source venv/bin/activate
```

> <tip-title>Python Within A Virtual Environment</tip-title>
> 
> Within a virtual environment, commands `python` and `pip` will refer to the version of Python you created the environment with. If you create a virtual environment with `python3 -m venv venv`, `python` will refer to `python3` and `pip` will refer to `pip3`. 
>
> On some machines with Python 2 installed, `python` command may refer to the copy of Python 2 installed outside of the virtual environment instead, which can cause confusion. You can always check which version of Python you are using in your virtual environment with the command `which python` to be absolutely sure. We continue using `python3` and `pip3` in this material to avoid confusion for those users, but commands `python` and `pip` may work for you as expected.
{: .tip}

### Installing External Libraries in an Environment with `pip`

We noticed earlier that our code depends on two *external libraries* - `numpy` and `matplotlib`. In order
for the code to run on your machine, you need to
install these two dependencies into your virtual environment.

To install the latest version of a package with `pip` you use pip's `install` command and specify the package’s name, e.g.:

```bash
pip3 install numpy
pip3 install matplotlib
```

or like this to install multiple packages at once for short:

```bash
pip3 install numpy matplotlib
```

> <tip-title>How About `python3 -m pip install`?</tip-title>
> Why are we not using `pip` as an argument to `python3` command, in the same way we did with `venv` 
> (i.e. `python3 -m venv`)? `python3 -m pip install` should be used according to the 
> [official Pip documentation](https://pip.pypa.io/en/stable/user_guide/#running-pip); other official documentation
> still seems to have a mixture of usages. Core Python developer Brett Cannon offers a 
> [more detailed explanation](https://snarky.ca/why-you-should-use-python-m-pip/) of edge cases when the two options may produce 
> different results and recommends `python3 -m pip install`. We kept the old-style command (`pip3 install`) as it seems more 
> prevalent among developers at the moment - but it may be a convention that will soon change and certainly something you should consider. 
{: .tip}

If you run the `pip3 install` command on a package that is already installed, `pip` will notice this and do nothing.

To install a specific version of a Python package give the package name followed by `==` and the version number, e.g.
`pip3 install numpy==1.21.1`.

To specify a minimum version of a Python package, you can
do `pip3 install numpy>=1.20`.

To upgrade a package to the latest version, e.g. `pip3 install --upgrade numpy`.

To display information about a particular installed package do:

```bash
pip3 show numpy
```

To list all packages installed with `pip` (in your current virtual environment):

```bash
pip3 list
```

To uninstall a package installed in the virtual environment do: `pip3 uninstall package-name`.
You can also supply a list of packages to uninstall at the same time.

### Exporting/Importing an Environment with `pip`

You are collaborating on a project with a team so, naturally, you will want to share your environment with your
collaborators so they can easily 'clone' your software project with all of its dependencies and everyone
can replicate equivalent virtual environments on their machines. `pip` has a handy way of exporting,
saving and sharing virtual environments.

To export your active environment - use `pip freeze` command to
produce a list of packages installed in the virtual environment.
A common convention is to put this list in a `requirements.txt` file:

```bash
pip3 freeze > requirements.txt
cat requirements.txt
```

The first of the above commands will create a `requirements.txt` file in your current directory.
The `requirements.txt` file can then be committed to a version control system and
get shipped as part of your software and shared with collaborators and/or users. They can then replicate your environment and
install all the necessary packages from the project root as follows:

```bash
pip3 install -r requirements.txt
```

As your project grows - you may need to update your environment for a variety of reasons. For example, one of your project's dependencies has
just released a new version (dependency version number update), you need an additional package for data analysis
(adding a new dependency) or you have found a better package and no longer need the older package (adding a new and
removing an old dependency). What you need to do in this case (apart from installing the new and removing the
packages that are no longer needed from your virtual environment) is update the contents of the `requirements.txt` file
accordingly by re-issuing `pip freeze` command and propagate the updated `requirements.txt` file to your collaborators
via your code sharing platform (e.g. GitHub).

> <tip-title>Official Documentation</tip-title>
> For a full list of options and commands, consult the [official `venv` documentation](https://docs.python.org/3/library/venv.html)
> and the [Installing Python Modules with `pip` guide](https://docs.python.org/3/installing/index.html#installing-index). Also check out the guide ["Installing packages using `pip` and virtual environments"](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#installing-packages-using-pip-and-virtual-environments).
{: .tip}

## Running Python Scripts From Command Line
Congratulations! Your environment is now activated and set up to run your script
from the command line.
