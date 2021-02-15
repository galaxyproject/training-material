---
layout: tutorial_hands_on
logo: "GTN"

title: Tool building with the ToolFactory in Galaxy
objectives:
  - Learn to build simple tools quickly using a Galaxy tool - the ToolFactory
questions:
  - What does the ToolFactory (TF) do?
  - When might it be useful?
  - What resources are needed to generate a new Galaxy tool?
  - Limits and types of tools the TF can make.
time_estimation: 20M
key_points:
  - Integrated environments are very handy for bioinformaticians
  - Easy way to develop new scripts
  - The ToolFactory can quickly turn these into Galaxy tools
  - Generated tools have a test, are lint free and indistinguishable from hand written tools.
  - No conditionals. Easy to manage in scripts but some Conda dependencies won't like empty default values.
contributors:
  - fubar2
---
> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
>
{: .agenda}

---

# Integrated Environments and tools in Galaxy: bridging the gap for simple scripts

- Bioinformaticians like the Galaxy IE feature. In a persistable integrated environment notebook, they can quickly write and test scripts that do useful things for some analysis.
They are very popular and are a great way to develop code that can run in Galaxy, free from the constraint of running pre-installed Galaxy tools.
- Graphical notebook environments are far more flexible and generalised than any possible form driven Galaxy tool. Like tools, they can be persisted and shared. Technically skilled users can run them easily enough, but many users who
are only familiar with the Galaxy tool form interface may find them confusing. They are not available for user workflows.
- Getting from an IE script to a Galaxy tool requires getting the script working on a command line, then writing a tool wrapper. The recommended, fully featured tool wrapper development tools such as Planemo and VCS with the new Galaxy
Language server, are covered in other training sessions. They are complex and have substantial learning curves. They are not limited once the necessary skills and knowledge are acquired, but it
can take time. That time may not be easy for a busy bioinformatician to find.
> - There is a gap for newcomers to Galaxy who have not yet had time to come up to speed with the developer tools. They can quickly get things working using integrated environments and can develop scripts that can be used in tools, but
making a new tool for the first time manually requires some substantial effort. A tool generator can automate most of the work of wrapping simple scripts.
> - This tutorial introduces a Galaxy tool that is a tool generator, enabling a bioinformatician to quickly convert any simple working script into a *new* real tool. There are limits compared to the developer tools but for many
simple cases, the generator works well. A Galaxy tool is used to generate the new tool, so a bioinformatician already familiar with the Galaxy forms GUI requires relatively little training. The
generation of the tool is entirely form driven and happens inside Galaxy. The Galaxy GUI is horrible when forms get very large so although the ToolFactory can probably cope with a very large number
of form fields, the bioinformatician might not. Simple tools with a few parameters are easy and quick.
- Scripts can easily be changed to accommodate some of the ToolFactory's limits. Some conda dependencies need conditional parameters because they fail if given empty values on the command line.
This is not a limit for a bioinformatican IE developed script, because it can cope with empty parameters being always passed. The ToolFactory can will generate a positional or argparse style command line
that can usually be made to work with a script. Many conda packages can probably be wrapped using the ToolFactory, but
scripts, particularly those known to work on the command line, are particularly suitable.

---


# The ToolFactory. A Galaxy tool to generate new simple tools from scripts

- Form driven simple tools from scripts in Galaxy, using a Galaxy tool. Any simple script. Important that it works on the command line first. The ToolFactory is too clumsy for debugging scripts. Much easier to use command line
development tools and get it working correctly.

- Script (e.g.) initially developed in an IE needs to be changed. It must accept command line parameters somehow. For bash, positional parameters work well. For scripting languages, argparse style parameters are preferable to avoid
mixing parameters up. The ToolFactory can provide both. Debugging will require a set of small sample input files and these will be needed to build the new tool.

- Sample datasets are then uploaded to the history and the ToolFactory tool is run, revealing the long and complicated form that needs to be completed.

- Specify dependencies, I/O and parameters. Paste optional script.

- Generate a new tool with a test, ready for any toolshed or server.

- **Simple tools only** - limited compared to Galaxy developer tools

- No conditionals - ok for scripts where can interpret default '' as None if not set by user

- Problematic for many packages that may not cope - a bash caller might work

---

### The ToolFactory is a *Galaxy tool*. It generates new Galaxy tools

- Wrap many Conda (or system) executable or interpreter+script

- Standard Galaxy tool XML with a test based on the samples

- Auto-generated from the details supplied on a Galaxy form

- Packaged ready for upload to a toolshed

- Could probably wrap any IE script


---


### When the TF might be useful

- Busy bioinformatician supporting Galaxy users

- Not (yet) familiar with the recommended Galaxy developer tools

- Users demanding functionality not (yet!) available in the toolshed

- IE sessions yield working scripts ready for production

- Scripts needed (quickly!) for repeatable user workflows.



---

### TF options: Planemo

- planemo tool\_factory … is very handy but **not persistent**

- PR not accepted yet so need to use a fork.

- Clone fubar2/planemo; make a venv; python setup.py install

- Tool archives unpacked topath can be loaded with **--extra\_tools path**



---

### TF options: Docker

- TF docker - has toolshed, workflow and example tools built in.

- Can be **persistent** if local volumes are mounted in image.

- Clone fubar2/toolfactory-galaxy-docker.

- Edit & run startclean.sh

- FROM quay.io/bgruening/galaxy:latest


---


### Galaxy ToolFactory: Process overview
![](../../images/wrapper_big_picture_1.png)

---

### Galaxy ToolFactory: Tool form for hello world

![](../../images/ToolFactory_hello1form.png)

---
```xml
<tool name="hello_toolshed" id="hello_toolshed" version="0.01">
  <!--Source in git at: https://github.com/fubar2/toolfactory-->
  <!--Created by planemo@galaxyproject.org at 22/01/2021 13:48:27 using the Galaxy Tool Factory.-->
  <description>Says hello</description>
  <stdio>
    <exit_code range="1:" level="fatal"/>
  </stdio>
  <version_command><![CDATA[echo "0.01"]]></version_command>
  <command><![CDATA[bash
$runme
"$sayhelloto"
>
$Hello_output]]></command>
  <configfiles>
    <configfile name="runme"><![CDATA[
echo "Hello $1"
]]></configfile>
  </configfiles>
  <inputs>
    <param label="Say hello to" help="" value="Sailor!" type="text" name="sayhelloto" argument="sayhelloto"/>
  </inputs>
  <outputs>
    <data name="Hello_output" format="txt" label="Hello_output" hidden="false"/>
  </outputs>
  <tests>
    <test>
      <output name="Hello_output" value="Hello_output_sample" compare="diff" lines_diff="0"/>
      <param name="sayhelloto" value=""/>
    </test>
  </tests>
  <help><![CDATA[

**What it Does**

ToolFactory demonstration - hello world in Galaxy



------


Script::

    echo "Hello $1"

]]></help>
  <citations>
    <citation type="doi">10.1093/bioinformatics/bts573</citation>
  </citations>
</tool>
```

### Practical limits and usage patterns

- Total all parameters (i/o files; user settings;..) = 10 or fewer. Perhaps 20.

- No conditionals (yet)

- Bash/sed/awk... system utilities with positional parameters and/or scripts

- Python/R/…(anything in conda) scripts with argparse or positional parameters

- Any conda package(s) with sane command line structure should work.

- Tutorials: Hello world!,...,BWA/samtools Planemo “advanced” example.



---

### Next steps

- Watch “hello World” and other simple tool example videos

- View the planemo TF example workflow demonstration.

- Run it yourself!

- Prepare data and a script for your own purposes and try the TF

- If it works for you, please show your colleagues

- Otherwise, pull requests and issues welcome @ fubar2/toolfactory
