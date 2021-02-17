---
layout: tutorial_hands_on
logo: "GTN"
questions:
- Who might want to use the ToolFactory for quick tools from scripts
- How can I get the ToolFactory working since it should not be on a public server?
- How does a bioinformatician get from an Integrated Environment script to a proper Galaxy tool?
title: "Introduction to Quick tools from scripts using a Galaxy tool generator"
type: tutorial_hands_on
contributors:
  - fubar2
---

# ToolFactory intended audience, warnings and installation options.

## Intended Audience

> The ToolFactory is for bioinformaticians.
> It provides a quick route to turn simple scripts developed using Interactive Environment into real Galaxy tools.
> It is a code generator so is limited in many ways compared to the recommended Galaxy developer tool development software.
> It has no conditionals. Yet. Pull requests accepted.

---

## Not useful unless you already have a script that does something useful on the command line

- The ToolFactory is an unusual Galaxy tool.
- It makes new Galaxy tools. It is not designed for ordinary scientist users.
- It is something a bioinformatician who is already comfortable with scripting languages on a linux command line might want to run.

Although it is an ordinary looking Galaxy tool, it automates much of the work needed to prepare a new Galaxy
tool using information provided on the ToolFactory form. *However, the ToolFactory does not do any of the work needed to prepare a working script.*
Galaxy is far too clumsy as an IDE to be used for that purpose.

The user can use it to wrap any simple script, but it must already be proven to run correctly on the command line with some small test input samples and indeed,
this is exactly what the ToolFactory does best.
Normal scientist Galaxy users are unlikely to need it unless they are also capable of confidently scripting for themselves.
---
## Warning !

>- *Please do not install it on a public facing server*
>- Although it will only run for administrative users, it allows unlimited scripting and that is never a good idea on a public facing machine. Please install it locally as described below.
>- For this reason, the training materials can't make use of existing public Galaxy infrastructure like most of the GTN material. Fortunately, there are a number of local installation alternatives available, depending on how you prefer to work.

---
## The quickest and simplest way

- Make a new (potentially throw away) directory for the Planemo installation - e.g. `mkdir tftute` and `cd tftute`
- [click here to download a bash script to build a local (potentially throw-away) copy of the ToolFactory] (https://zenodo...)
- Save the script in the new directory and run it.
- It will download a planemo fork and a separate local copy of Galaxy. You can change the location of the galaxy directory
and remove the code to download one if you wish to save time.
- It will take some time - so watch the Hello World demonstration while you wait.
---
## Using Docker

- The Docker container provided with this topic is very different from most GTN Docker containers. It does not include this tutorial.
It runs planemo tool_factory for you and exposes it on port 9090 so you can do all the same things as you
can with a local venv described above - but a little slower and isolated in a container.

- There is a more complex but integrated solution using the ToolFactory docker container. It provides an integrated toolshed and allows
tools to be installed and used in the Galaxy used to run the ToolFactory. The big advantage is that it can be persisted as shown
in the documentation for docker-galaxy-stable upon which it is based.

---

# Learning how to use the ToolFactory


---
Please follow our
[tutorial to learn how to fill the slides]({{ site.baseurl }}/topics/contributing/tutorials/create-new-tutorial-slides/slides.html)
