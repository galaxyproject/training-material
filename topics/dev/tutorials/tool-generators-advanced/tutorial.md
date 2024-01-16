---

layout: tutorial_hands_on
title: "ToolFactory: Generating Tools From More Complex Scripts"
key_points:
  - The ToolFactory Appliance includes an automated form driven tool generator.
  - It produces new tools as Toolshed ready archives
  - It was designed for scientists and developers who routinely write scripts for their analyses.
  - It can quickly turn a working command line script into a working Galaxy tool.
  - It generates tools from information entered on a Galaxy form in the familiar UI.
  - The new tool is installed in the appliance after generation so can be used and explored immediately.
  - New tools can be adjusted as needed by re-running the job that generated the tool and updating the form settings to suit.
  - Adding a test to the toolshed archive takes a minute or so and is only needed when the new tool is ready to export for sharing.
  - Experimenting with the ToolFactory samples will help acquire insight and skills useful for manual tool building.

objectives:
 - Further develop your ToolFactory Skills

questions:
 - What else can I do with the ToolFactory?

subtopic: tooldev
time_estimation: 1H

requirements:
  - type: "internal"
    topic_name: introduction
    tutorials:
      - galaxy-intro-short
      - galaxy-intro-101-everyone
  - type: "internal"
    topic_name: dev
    tutorials:
      - tool-integration
      - interactive-environments
      - tool-generators

contributors:
  - fubar2

---

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

---

# Background and a user's guide to this training material

Galaxy users who write and share scripts useful for scientific analyses are likely to be reading this material, perhaps after seeing the "Hello Galaxy"
demonstration. It was written to help you find out about the capabilities and limits of the ToolFactory by experimenting with it yourself.
It is hoped that this advanced tutorial will introduce some features that potentially make the ToolFactory useful in your work.

This training material is unlike most other GTN tutorials. There is no specific tool building curriculum on offer because it is hard to know how
that might be orgranised. The introductory `Hello Galaxy` demonstrations and
hands-on are the only formal attempt at explaining any tool's structure and construction.

The reason is that tools are mostly scientific domain specific, while the ToolFactory is applicable to any domain where Galaxy might be useful.
Developers and the tools they might be interested in creating are so diverse
that it was decided to provide the infrastructure and environment to encourage "learning by doing".
The result of editing and regenerating a tool is immediately available in the ToolFactory for evaluation to encourage self-guided exploration.
Developers are free to explore their own specific requirements.

Simple demonstration tools are provided. The are both the functional documentation, and useful models to adjust and explore.
They each address a specific feature (such as the `plotter` illustrating R parameter passing and collectionscollections),
and they illustrate some of the styles of tool building supported by the ToolFactory. The `sed` example does not use a script. All the others do.
The ToolFactory is likely to be more useful for scripts than directly driving Conda packages. Many Conda package interfaces have requirements that the ToolFactory cannot
currently address. Those will probably always require a skilled tool document developer. Conda dependencies requiring minimal conditional parameter
logic for the construction of a command line or for adjusting outputs, such as the `sed` example, may be more tractable.

This training material will introduce some ToolFactory features and will show how some of the limitations can be worked around.
The introductory tutorial walked through some simple tools. This advanced tutorial shows some examples but the focus is on self-guided
exploration of the sample tools to see how they might be used in your work.


## The ToolFactory is for scientists and developers who write command line scripts in their work.

The ToolFactory is distributed as a docker image or as a self-configuring clone of the Galaxy server code. This is intended to make it easy to create new tools, analyses and workflows can be developed to help speed uptake in areas of data-intensive science where
few Galaxy tools are currently available.

The ToolFactory automates much of the work needed to prepare a new Galaxy tool using information provided by the script writer
on a Galaxy tool form. It can generate XML to wrap any simple script that runs correctly on the linux command line with some small test input samples. This is potentially
handy for developers new to Galaxy, and for Galaxy users who are capable of correctly scripting on the command line for themselves, because those working scripts can
be wrapped and tested to make toolshed ready, shareable Galaxy tools, or if too trivial to be worth sharing, used on the host desktop Galaxy for real analyses.

Generated tools are immediately installed and ready to run. This provides instant feedback for the developer in an integrated tool
development environment for simple scripts. Jobs that generate tools can always be rerun using {% icon galaxy-refresh %}. The form reappears as it was
when the tool was generated. Data and user configurable parameters can be added or removed. Text that appears when the tool is run, such as user parameter or data input
labels and help, can be updated by editing the form. When the updated tool is generated, installed and run as a tool, the tool form will include all the changes, ready to run.


The development server is a self-configuring, disposable platform for a data intensive scientist to learn how Galaxy tools work, in preparation for the more flexible and efficient
Galaxy project supported manual tool building infrastructure. At the same time, they can wrap and refine new tools from existing scripts
on their own workstations, preparing them for sharing and deployment to production servers. It also may be useful as a private sandbox for learning about tools and
Galaxy administration or experimenting with code development.


> <tip-title>Under the hood:</tip-title>
>  - It uses [galaxyml](https://github.com/hexylena/galaxyxml) to generate the tool XML from ToolFactory form settings.
{: .tip}

---

## Limits and scope

- The ToolFactory can generate, install and run new tools from scripts on your desktop.
    - It is a fully functional development Galaxy instance, so it can import any existing tool from a toolshed and can handle as much data as your desktop disks will fit.
    - It takes 30 minutes to install locally or a few minutes as a docker image, and can be completely removed when no longer useful.
    - The user is entirely responsible for securing and backing up all their work in the docker version. The desktop version is a normal persistent Galaxy development server.
    - It is not recommended for production and is *not supported by the Galaxy team*.
    - A well endowed modern workstation or laptop with plenty of cores, disk and RAM is needed. Older commodity hardware may struggle.
- The ToolFactory works best wrapping simple R/Bash/Python and other interpreted scripts, with a few user supplied parameters and a few I/O history files.
- Conditional parameters, XML macros, output filters and many other advanced features are not currently supported. Manual coding will be needed.
- Scripts are easier than some Conda packages
    - This is because the tool builder can modify the code to respond to default empty parameters as if they had not been passed.
    - This may not sound like much, but as a result, advanced tool building elements such as conditionals and related tricks requiring manual coding, can sometimes be avoided.
    - In contrast, some Conda dependencies or combinations will require output filters
or other complex tool XML constructs outside the `<command>` section, that are not easy to generate automatically and can only be coded in the document.
    - While some simple requirements may be manageable, complex ones will not be suitable for the ToolFactory.
- Compared to the more usual linux shell and text editor, the ToolFactory appliance is a rather clumsy way to debug scripts.
    - Starting a new ToolFactory tool with a know good command line and data is strongly recommended.
      - You will know exactly what to expect from the tool test for a first sanity check.
    - Corrolary: Unless there is a working script that needs to be wrapped into a toolshed-ready Galaxy tool, the ToolFactory is of little use
      - other than for learning about Galaxy tools.
- The ToolFactory Appliance is for developers and code-writing scientists not yet familiar with the more flexible and complex manual tools, and who need to wrap scripts that are simple
enough for the ToolFactory.
    - Compared to the more flexible manual Galaxy tool development software, there is far less to learn to get up to speed with a form driven,
automated code generator in a tailored, readily deployed appliance.
    - The cost of this convenience is that ToolFactory is limited to a subset of simple script and package wrappers.
    - Much can be learned in the process of trying things out in the ToolFactory that will help acquire skills for manually building tools.

---

# Installation (Copied from the introductory tutorial)

> <warning-title>Security advisory!</warning-title>
>- *Please do not install the ToolFactory on any public server*
>- Configured as a default server, it has none of the usual additional security layers required for a safe public Galaxy server.
>- Although it will only run for administrative users, it allows unlimited scripting and exposes unacceptable risks.
>- Please install it locally.
>- For this reason, the training materials can't make use of existing public Galaxy infrastructure like most of the GTN material.
{: .warning}


> <hands-on-title>Installing and managing a Galaxy ToolFactory development server</hands-on-title>
> # Logging in as an administrator
> 
> Once you have a working installation running as described below, the server should be ready in 20-30 seconds or so, at [http://localhost:8080](http://localhost:8080).
> *The ToolFactory will only execute for administrative users as a minimal security barrier.* 
> When the webserver starts, immediately login using the administrative user email *toolfactory@galaxy.org* with the password *ChangeMe!* which of course you should change.
>   
> # Docker
> Note: *Nothing is persistent in the image*. Useful work must be manually exported and saved.
> For following the GTN tutorial and for test driving the ToolFactory for the first time, the docker version is recommended.
> Non-persistent means it does not remember anything after you stop the container. The next time it starts it will be a fresh installation.
> You can save your work by exporting the histories and tool tarballs you want to keep. You need the history to rerun the job that generated a tarball, so the history is the most important thing to preserve if you make a useful tool.
> The image must be pulled first, then run with port 8080 open for the Galaxy server.
> 
> ```
> docker pull quay.io/fubar2/galaxy_toolfactory:latest
> docker run -d -p 8080:8080 quay.io/fubar2/galaxy_toolfactory:latest
> ```
> Check the docker logs until gunicorn is ready to serve or wait 
> about 20-30 seconds, then browse to [http://localhost:8080](http://localhost:8080)
> If a Galaxy server appears, proceed with the login instructions above and you should see a history containing all the example tools. There is also a workflow that can reproduce that history if you set all the inputs to the right input datasets in that history.
>  
> # Local workstation development Galaxy server installation
> 
> A persistent desktop development ToolFactory server can be built by cloning the [bootstrap github repository](https://github.com/fubar2/galaxy_tf_overlay), and using the included *localtf.sh* script, to build a complete, new development server with the ToolFactory installed and ready to run.
>
> From a convenient directory, download the overlay configuration code repository, then
> run the *localtf.sh* setup script from that cloned repository root directory:
>
> ```
> git clone https://github.com/fubar2/galaxy_tf_overlay.git
> cd galaxy_tf_overlay
> sh ./localtf.sh
> ```
> Running *localtf.sh* will create a new directory, *galaxytf*, in the parent directory of the *galaxy_tf_overlay* clone.
> The script will configure a fresh Galaxy 23.0 server, with the ToolFactory installed, in that new directory.
> This takes 20 minutes or more to complete since the client must be built once.
>
> The resulting development server directory will occupy ~9GB of disk space, so be sure your machine has plenty of room.
> It will be based in a single directory, *galaxytf* in the same directory as the galaxy_tf_overlay repository was cloned into.
> That is where the script should be run as shown above.
>
> Rerunning the *localtf.sh* script will *destroy the entire galaxytf directory* - all ~9GB, and create a clean new installation.
> It should only need to be run once in the life of the development server.
>
> Remove the *galaxytf* directory to remove the entire development server when it is no longer needed. Save all your tools and histories,
> because the jobs in a history can be used to update a tool easily, and a history can be imported into a fresh development instance
> when needed.
>
>
> Once installation is complete:
>  * start the server from the *galaxytf* directory with *sh run.sh*. The logs will be displayed.
>  * ^c (control+c) will stop it from the console.
>  * In 23.0 that is equivalent to *.venv/bin/galaxyctl start* and *.venv/bin/galaxyctl stop*.
>
>
{: .hands_on}


The trivial `Hello World!` tool example is readily extended to suit many situations where a tool is needed quickly for a workflow. Try adding another parameter.
For example, the planemo `lint` tool example (described below) can be derived by adding a history toolshed archive as input, plus a few more lines of bash script.
In practice, it's a flexible basis for generating many simple tools.


# ToolFactory tips and tricks illustrated by some of the examples.

## Before you begin a new tool

- Make sure it runs correctly on a command line with your sample inputs and default parameter settings. That is how the test will be generated.
- You cannot specify any inputs on the form without providing samples, and those samples must run correctly, with the supplied defaults.
- Easiest to make sure those samples are in the history before you begin.
- Be well prepared. The ToolFactory cannot do that preparation for  you.

## STDIN and STDOUT

- Demonstration tools often capture output from a bash script using the special STDOUT designation for output files
- This can save sending the output path as a parameter to the script or executable
- STDIN is also available as a special designation for history inputs if the script takes input from STDIN when it runs.
- The option `Tool reads selected input file from STDIN and writes STDOUT with no parameters` in the parameter passing model
selection list will generate a simple filter tool, executing a script or Conda dependency with input on STDIN and output on STDOUT
for those rare situations where that's all you need. No i/o or other parameters for the user to set. Used in the tacrev demonstration tool.

## Repeats in the ToolFactory form theoretically permit *any* number of parameters.

- Inputs, outputs and additional user parameters are all presented within repeat elements on the form.
- For any given tool, repeats allow an unlimited number of each of those tool elements - from a technical perspective.
- The practical limit is reached when your patience runs out as you add them in the current Galaxy forms UI.
- The form gets very long, very quickly.
- A handful is easily manageable.

## Repeated parameters in generated tools allow the end user to supply multiple values

- Repeats *on the generated tool form* are supported for input and user edited parameters
- The generated tool form has the usual "Add..." button associated with the parameter, so the tool user can add any number of them.
- The script must be able to parse and deal correctly with multiple instances of the same parameter name.
- A Python sample script using argparse with `action="append"` parameters to deal with potentially any number of multiples of each parameter command
line will echo all the repeated parameters is shown in the example shown in the example below.
- Repeats do not make sense for `positional parameters` because their number is unpredictable. The ToolFactory will ignore repeats in positional mode. A warning is issued in the log.
- Repeats on `output parameters` are not supported - use a `collection` described below when an unpredictable number of output files are required.
- The sample shown below has one input data parameter and one text parameter. The `repeat` option is selected for both on the ToolFactory form.
- In the generated tool XML wrapper, these are embedded within `<repeat>` tags in the `<inputs>` and also in the `<tests>` sections automatically.
- Note that repeats are limited to *single parameters* - groups of parameters cannot currently be repeated. They will require manual tool writing.
- The end user sees the form with repeatable fields as shown at the end of the detail below.
- The demonstration returns whatever the user chose to repeat.


> <details-title>Repeats demonstration generated tool XML</details-title>
>
>
> ```xml
> <tool name="repeats_demo" id="repeats_demo" version="2.00">
>   <!--Source in git at: https://github.com/fubar2/galaxy_tf_overlay-->
>   <!--Created by admin@galaxy.org at 29/05/2021 07:46:08 using the Galaxy Tool Factory.-->
>   <description>Repeated parameter demonstration</description>
>   <requirements/>
>   <stdio>
>     <exit_code range="1:" level="fatal"/>
>   </stdio>
>   <version_command><![CDATA[echo "2.00"]]></version_command>
>   <command><![CDATA[python
> $runme
> #for $rep in $R_mi:
> --mi "$rep.mi"
> #end for
> #for $rep in $R_mp:
> --mp "$rep.mp"
> #end for
> >
> $repeats_out]]></command>
>   <configfiles>
>     <configfile name="runme"><![CDATA[#raw
>
> import argparse
> parser = argparse.ArgumentParser()
> a = parser.add_argument
> a("--mi", action="append")
> a("--mp", action="append")
> args = parser.parse_args()
> if args.mi:
>    print(" file and ".join(args.mi))
> if args.mp:
>    print(" string and ".join(args.mp))
> if not (args.mi or args.mp):
>    print('Nothing was selected')
>
> #end raw]]></configfile>
>   </configfiles>
>   <inputs>
>     <repeat name="R_mi" title="Add as many Multiple input files from your history - as many as you like as needed">
>       <param name="mi" type="data" optional="false" label="Multiple input files from your history - as many as you like" help="" format="html,txt,xml" multiple="false"/>
>     </repeat>
>     <repeat name="R_mp" title="Add as many Multiple user supplied text strings - as many different ones as you like as needed">
>       <param name="mp" type="text" value="Multiple user supplied text strings - as many different ones as you like" label="Multiple user supplied text strings - as many different ones as you like" help=""/>
>     </repeat>
>   </inputs>
>   <outputs>
>     <data name="repeats_out" format="txt" label="repeats_out" hidden="false"/>
>   </outputs>
>   <tests>
>     <test>
>       <output name="repeats_out" value="repeats_out_sample" compare="diff" lines_diff="6"/>
>       <repeat name="R_mi">
>         <param name="mi" value="mi_sample"/>
>       </repeat>
>       <repeat name="R_mp">
>         <param name="mp" value="Multiple user supplied text strings - as many different ones as you like"/>
>       </repeat>
>     </test>
>   </tests>
>   <help><![CDATA[
>
> **What it Does**
>
> Simple python Argparse sample to echo repeated user selections - how to use repeated inputs and user parameters.
>
> Unpredictable or messy "repeated" outputs can use a collection if they are not useful downstream but otherwise require manual wrapping - see the GTN advanced tutorial.
>
>
>
> ------
>
>
> Script::
>
>     import argparse
>     parser = argparse.ArgumentParser()
>     a = parser.add_argument
>     a("--mi", action="append")
>     a("--mp", action="append")
>     args = parser.parse_args()
>     if args.mi:
>        print(" file and ".join(args.mi))
>     if args.mp:
>        print(" string and ".join(args.mp))
>     if not (args.mi or args.mp):
>        print('Nothing was selected')
>
> ]]></help>
>   <citations>
>     <citation type="doi">10.1093/bioinformatics/bts573</citation>
>   </citations>
> </tool>
>
>
> ```
>
>  - The user sees the following form after adding 3 repeats for each of the two available items
>
>  ![User can add as many repeats as they want for repeated input file and parameter values](../../images/toolfactory_repeats_sample_form.png)
{: .details}

## ToolFactory `collection` outputs are handy for hiding dozens of miscellaneous tool outputs in a single history item

- If the script writes a large number of different output file types (images, reports..) that are not in themselves useful for downstream analyses.
- A Collection can be used to hide them in an organised list as a single new history item.
- Collections are special kinds of Galaxy history items that can present a directory full of files arranged inside one single new history item.
- When viewed, that history item presents them all as a list, each a typical, viewable Galaxy history item.
- They don't clutter the user's history because they are all hidden in the collection.
- The plotter example uses an Rscript.
- It generates as many pairs of random plots as the parameter supplied requires.
- The script sends them into the the collection that appears in the history after the job runs.
- The user's history shows only one new item after it runs - it must be viewed to see all the individual contents.

> <warning-title>The default generated test for output collections always passes because it doesn't test anything.</warning-title>
>
>    - Supplying a test over-ride is recommended for collections.
>    - Example code is shown on the sample tool's form and in the original example code below - removed from the current sample.
>    - For a real test, one or more expected <element.../> tags must be provided so the test really does test something.
>    - Add your own file names to the sample to make a real test for your own collections.
>    - Otherwise, the generated test will be empty, because the code generator knows nothing about what's going to appear in the collection at the time the `<test>` code is generated.
>    - Introspecting arbitrary scripts to reliably populate the test with actual file names?. Not likely any time soon.
>
{: .warning}

> <details-title>`plotter` collection output demonstration tool form, generated XML and outputs</details-title>
>
>  - The ToolFactory form for the plotter example tool
>  is configured as shown below, from "rerunning" the plotter job from the sample history.
>
> ![ToolFactory form configuration of the output collection in the plotter example](../../images/toolfactory_plotter_demo_form.png)
>
>  The Rscript is contained in a configfile so`#` is escaped - this is automatic.
>
> ```xml
> <tool name="plotter" id="plotter" version="0.01">
>   <!--Source in git at: https://github.com/fubar2/tgalaxy_tf_overlay-->
>   <!--Created by admin@galaxy.org at 24/01/2021 05:02:33 using the Galaxy Tool Factory.-->
>   <description>ToolFactory collection demonstration - random plots</description>
>   <requirements>
>     <requirement version="" type="package">r-base</requirement>
>   </requirements>
>   <stdio>
>     <exit_code range="1:" level="fatal"/>
>   </stdio>
>   <version_command><![CDATA[echo "0.01"]]></version_command>
>   <command><![CDATA[Rscript
> $runme
> "$nplot"]]></command>
>   <configfiles>
>     <configfile name="runme"><![CDATA[
> \# demo
> args = commandArgs(trailingOnly=TRUE)
> if (length(args)==0) {
>    n_plots = 3
> } else {
>    n_plots = as.integer(args[1]) }
> dir.create('plots')
> for (i in 1:n_plots) {
>     foo = runif(100)
>     bar = rnorm(100)
>     bar = foo + 0.05*bar
>     pdf(paste('plots/yet',i,"anotherplot.pdf",sep='_'))
>     plot(foo,bar,main=paste("Foo by Bar plot \#",i),col="maroon", pch=3,cex=0.6)
>     dev.off()
>     foo = data.frame(a=runif(100),b=runif(100),c=runif(100),d=runif(100),e=runif(100),f=runif(100))
>     bar = as.matrix(foo)
>     pdf(paste('plots/yet',i,"anotherheatmap.pdf",sep='_'))
>     heatmap(bar,main='Random Heatmap')
>     dev.off()
> }
>
> ]]></configfile>
>   </configfiles>
>   <inputs>
>     <param label="Number of random plots pairs to draw" help="" value="3" type="text" name="nplot" argument="nplot"/>
>   </inputs>
>   <outputs>
>     <collection name="plots" type="list" label="Plots">
>       <discover_datasets pattern="__name_and_ext__" directory="plots" visible="false"/>
>     </collection>
>   </outputs>
>
>
>   <tests>
>     <test>
>       <param name="nplot" value="3" />
>       <output_collection name="plots" type="list">
>      <element file="yet_1_anotherplot_sample" name="yet_1_anotherplot" ftype="pdf" compare="sim_size" delta_frac="0.05"/>
>     </output_collection>
>  </test>
>   </tests>
>
>
>
>   <help><![CDATA[
>
> **What it Does**
>
> Draws as many random plot pairs as you need
>
>
>
> ------
>
>
> Script::
>
>     # demo
>     args = commandArgs(trailingOnly=TRUE)
>     if (length(args)==0) {
>        n_plots = 3
>     } else {
>        n_plots = as.integer(args[1]) }
>     dir.create('plots')
>     for (i in 1:n_plots) {
>         foo = runif(100)
>         bar = rnorm(100)
>         bar = foo + 0.05*bar
>         pdf(paste('plots/yet',i,"anotherplot.pdf",sep='_'))
>         plot(foo,bar,main=paste("Foo by Bar plot #",i),col="maroon", pch=3,cex=0.6)
>         dev.off()
>         foo = data.frame(a=runif(100),b=runif(100),c=runif(100),d=runif(100),e=runif(100),f=runif(100))
>         bar = as.matrix(foo)
>         pdf(paste('plots/yet',i,"anotherheatmap.pdf",sep='_'))
>         heatmap(bar,main='Random Heatmap')
>         dev.off()
>     }
>
> ]]></help>
>   <citations>
>     <citation type="doi">10.1093/bioinformatics/bts573</citation>
>   </citations>
> </tool>
> ```
>  After requesting 25 pairs of plots from the sample tool, a collection appears in the history and is shown below.
>  One of them is displayed by clicking the "eye" icon.
>  Collections are ideal for messy analysis reporting outputs such as images, pdfs and other material that is not useful as an input to a downstream tool.
>  It is material that the user will want kept together, so a single history item is ideal to avoid unnecessary clutter.
>  As shown above, the script only has to write the files to a directory.
>  Note that the test is over-ridden in the ToolFactory form to generate this tool.
>  Without an over-ride to suit your script, an empty test will be generated that will not do any testing - see warning above.
>
>
> ![Plotter collection opens to show all the plots as separate items to be viewed or downloaded](../../images/toolfactory_plotter_sample_output.png)
{: .details}

## Selects as user supplied parameters

- Additional parameter types are selected from a drop down list. It includes text, numeric and select parameters.
- Selects offer a list of prespecified options to the user.
- There is a repeat on the ToolFactory form to collect options as pairs of names/values.
- It is clumsy and suitable only for a limited number of options.
- Galaxyxml generates appropriate select parameters on the generated tool as shown in the select demonstration tool.

> <details-title>`select_test` select field demonstration tool generated XML</details-title>
>
> The ToolFactory form section for user configurable command line settings is
>  configured as shown here for the select demonstration
>
> ![Form to configure an additional parameter as a select](../../images/toolfactory_select_demo_form.png)
>
> The generated XML is shown below.
>
> ```xml
> <tool name="select_test" id="select_test" version="0.01">
>   <!--Source in git at: https://github.com/fubar2/galaxy_tf_overlay-->
>   <!--Created by admin@galaxy.org at 24/01/2021 05:03:21 using the Galaxy Tool Factory.-->
>   <description>ToolFactory select demonstration</description>
>   <stdio>
>     <exit_code range="1:" level="fatal"/>
>   </stdio>
>   <version_command><![CDATA[echo "0.01"]]></version_command>
>   <command><![CDATA[bash
> $runme
> "$choose"
> >
> $select_out]]></command>
>   <configfiles>
>     <configfile name="runme"><![CDATA[
> echo "You chose \$1"
> ]]></configfile>
>   </configfiles>
>   <inputs>
>     <param label="Choose" help="" type="select" name="choose" argument="choose">
>       <option value="won">one</option>
>       <option value="too">two</option>
>       <option value="free">three</option>
>     </param>
>   </inputs>
>   <outputs>
>     <data name="select_out" format="txt" label="select_out" hidden="false"/>
>   </outputs>
>   <tests>
>     <test>
>       <output name="select_out" value="select_out_sample" compare="diff" lines_diff="0"/>
>       <param name="choose" value="won"/>
>     </test>
>   </tests>
>   <help><![CDATA[
>
> **What it Does**
>
> Echoes your selection
>
>
>
> ------
>
>
> Script::
>
>     echo "You chose $1"
>
> ]]></help>
>   <citations>
>     <citation type="doi">10.1093/bioinformatics/bts573</citation>
>   </citations>
> </tool>
> ```
>  The generated tool form from the select demonstration shows the three options and returns the one selected.
>
> ![Generated form seen by users of the select demonstration tool](../../images/toolfactory_select_test_tool.png)
{: .details}



## The ToolFactory can wrap some Conda packages correctly using a simple wrapper script to work around limitations.

- Many complex tools require manual coding to make them useable.
  - Often complex tools hide data being moved between multiple dependencies such as samtools in the bwa tool.
  - These are often not suitable for the simple automated tool generator in the ToolFactory.
     - Sometimes a suitable script may be able to perform the complications needed so a tool can be created.
     - A lot depends on the ingenuity and preferences of the developer.
     - For specialist use and complex tools, the project developer supported tool building infrastructure is recommended.

> <details-title>`planemo lint` demonstration tool generated XML</details-title>
> The labels for outputs have been edited to include the tool name in this sample - this is not possible at present in the ToolFactory.
>
> ```xml
> <tool name="planemo_lint" id="planemo_lint" version="0.01">
>   <!--Source in git at: https://github.com/fubar2/galaxy_tf_overlay->
>   <!--Created by admin@galaxy.org at 18/05/2021 01:16:17 using the Galaxy Tool Factory.-->
>   <description>Runs Planemo lint on any ToolFactory xml history file</description>
>   <requirements>
>     <requirement version="0.74.3" type="package">planemo</requirement>
>     <requirement version="3.7" type="package">python</requirement>
>     <requirement type="package" version="4.6.3">lxml</requirement>
>   </requirements>
>   <stdio>
>     <exit_code range="1:" level="fatal"/>
>   </stdio>
>   <version_command><![CDATA[echo "0.01"]]></version_command>
>   <command><![CDATA[python
> $runme
> $ToolFactory_XML_to_be_linted
> >
> $lint_output]]></command>
>   <configfiles>
>     <configfile name="runme"><![CDATA[#raw
>
> import lxml.etree as ET
> import os
> import subprocess
> import sys
>
> def main():
>     assert len(sys.argv) >= 2, 'Must have input xml on command line'
>     xmlin = sys.argv[1]
>     tree = ET.parse(xmlin)
>     root = tree.getroot()
>     toolname = root.get('id')
>     toolxml = os.path.join('/export/galaxy/tools/TFtools', toolname, '%s.xml' % toolname)
>     cl = ['planemo', 'lint', toolxml]
>     print('Running', cl)
>     p = subprocess.run(cl, shell=False)
>     if p.returncode > 0:
>          print('Planemo lint call returned error %d')
>     else:
>          print('Lint report ends')
> main()
>
>
> #end raw]]></configfile>
>   </configfiles>
>   <inputs>
>     <param name="ToolFactory_XML_to_be_linted" type="data" optional="false" label="ToolFactory XML to be linted" help="" format="xml" multiple="false"/>
>   </inputs>
>   <outputs>
>     <data name="lint_output" format="txt" label="${ToolFactory_XML_to_be_linted.name}_lint_output" hidden="false"/>
>   </outputs>
>   <tests>
>     <test>
>       <output name="lint_output" value="lint_output_sample" compare="diff" lines_diff="5"/>
>       <param name="ToolFactory_XML_to_be_linted" value="ToolFactory_XML_to_be_linted_sample"/>
>     </test>
>   </tests>
>   <help><![CDATA[
>
> **What it Does**
>
> ToolFactory demonstration script using bash to run planemo lint from a history XML representing a tool.
>
>
>
> ------
>
>
> Script::
>
>     import lxml.etree as ET
>     import os
>     import subprocess
>     import sys
>     def main():
>         assert len(sys.argv) >= 2, 'Must have input xml on command line'
>         xmlin = sys.argv[1]
>         tree = ET.parse(xmlin)
>         root = tree.getroot()
>         toolname = root.get('id')
>         toolxml = os.path.join('/export/galaxy/tools/TFtools', toolname, '%s.xml' % toolname)
>         cl = ['planemo', 'lint', toolxml]
>         print('Running', cl)
>         p = subprocess.run(cl, shell=False)
>         if p.returncode > 0:
>              print('Planemo lint call returned error %d')
>         else:
>              print('Lint report ends')
>     main()
> #end raw
> ]]></help>
>   <citations>
>     <citation type="doi">10.1093/bioinformatics/bts573</citation>
>   </citations>
> </tool>
>
>
>
> ```
{: .details}


## Other languages such as Lisp and Prolog scripts in Galaxy.

- Any interpreted language in Conda can be used, as long as it will take a command line script path or read from STDIN.
- Trivial samples using Lisp and Prolog are included in the demonstration history.
- They may be old, but they show how it is easy to use any Conda scripting language package in Galaxy for long running applications.
- The Prolog tool has an inbuilt script. Substitute the sample script for real code and add inputs, user-configurable parameters and outputs to produce a tool wrapping a Prolog script if you ever need one in Galaxy.
- The Lisp tool will try to execute the selected input text file, hoping it contains a Lisp program like the `hello world` example supplied. This is a terrible idea for a public server and is shown only as a skeletal example of what is possible, not what is sensible.

> <details-title>`Prolog` and `Lisp` demonstration tools</details-title>
> ```xml
> <tool name="prolog_demo" id="prolog_demo" version="0.01">
>   <!--Source in git at: https://github.com/fubar2/galaxy_tf_overlay-->
>   <!--Created by admin@galaxy.org at 04/04/2021 17:21:32 using the Galaxy Tool Factory.-->
>   <description>Runs a prolog script</description>
>   <requirements>
>     <requirement version="" type="package">swi-prolog</requirement>
>   </requirements>
>   <stdio>
>     <exit_code range="1:" level="fatal"/>
>   </stdio>
>   <version_command><![CDATA[echo "0.01"]]></version_command>
>   <command><![CDATA[swipl
> -q
> -g
> main
> -s
>
> $runme
> >
> $prolog_out]]></command>
>   <configfiles>
>     <configfile name="runme"><![CDATA[
> parent(pam,bob).
> parent(tom,bob).
> parent(tom,liz).
> parent(bob,ann).
> parent(bob,pat).
> parent(pat,jim).
>
> main :-
>     parent(X,jim),
>     format('~a is the parent of jim~n', [X]),
>     halt.
> ]]></configfile>
>   </configfiles>
>   <inputs/>
>   <outputs>
>     <data name="prolog_out" format="txt" label="prolog_out" hidden="false"/>
>   </outputs>
>   <tests>
>     <test>
>       <output name="prolog_out" value="prolog_out_sample" compare="diff" lines_diff="0"/>
>     </test>
>   </tests>
>   <help><![CDATA[
>
> **What it Does**
>
>
>
> Prolog demonstration in the ToolFactory
>
>
>
> ------
>
>
> Script::
>
>     parent(pam,bob).
>     parent(tom,bob).
>     parent(tom,liz).
>     parent(bob,ann).
>     parent(bob,pat).
>     parent(pat,jim).
>     main :-
>         parent(X,jim),
>         format('~a is the parent of jim~n', [X]),
>         halt.
>
> ]]></help>
>   <citations>
>     <citation type="doi">10.1093/bioinformatics/bts573</citation>
>   </citations>
> </tool>
> ```
>
>  Lisp too!
>  In this case, the tool takes a lisp script as input
>
>  A lisp expression was typed into the upload tool text box and saved as hello_lisp.txt ```(write-line "Hello, ToolFactory does Lisp!")```.
>  It or any other lisp script can be used as input to the SBCL Conda dependency.
>
> ```xml
> <tool name="lisp_demo" id="lisp_demo" version="0.01">
>   <!--Source in git at: https://github.com/fubar2/galaxy_tf_overlay-->
>   <!--Created by admin@galaxy.org at 04/04/2021 16:38:47 using the Galaxy Tool Factory.-->
>   <description>Runs SBCL hello world demonstration</description>
>   <requirements>
>     <requirement version="" type="package">sbcl</requirement>
>   </requirements>
>   <stdio>
>     <exit_code range="1:" level="fatal"/>
>   </stdio>
>   <version_command><![CDATA[echo "0.01"]]></version_command>
>   <command><![CDATA[sbcl
> --script
> $script
> >
> $lisp_out]]></command>
>   <inputs>
>     <param optional="false" label="SBCL lisp script to execute" help="" format="txt" multiple="false" type="data" name="script" argument="--script"/>
>   </inputs>
>   <outputs>
>     <data name="lisp_out" format="txt" label="lisp_out" hidden="false"/>
>   </outputs>
>   <tests>
>     <test>
>       <output name="lisp_out" value="lisp_out_sample" compare="diff" lines_diff="0"/>
>       <param name="script" value="script_sample"/>
>     </test>
>   </tests>
>   <help><![CDATA[
>
> **What it Does**
>
> Lisp in the ToolFactory
>
>  ]]></help>
>   <citations>
>     <citation type="doi">10.1093/bioinformatics/bts573</citation>
>   </citations>
> </tool>
> ```
{: .details}


## Command over-ride and test over-ride

- If necessary, the automatically generated `<command>` section and/or the `<test>` section can be over-written by a developer supplied text.
- This is a hybrid solution.
  - Where the generated command needs to be embellished, logic can be placed directly in the tool document this way
  - Only useful for developers comfortable with writing Galaxy XML namespace template code.
- There are two simple BWA wrappers based on a Planemo documentation advanced example
- One uses a command over-ride pasted into the appropriate text box on the ToolFactory form
- This was based on the one shown in the example it is copied from.
- It allows templating - `${reference_fasta_filename}` is replaced by the current value of the `reference_fasta_filename` parameter on the tool form.
- The pasted over-ride completely replaces the galaxyxml generated ones.

> <details-title>`bwa_test_command_override` sample - the command override</details-title>
>  ToolFactory command over-ride section adapted from the Planemo BWA example.
>
> ```
>  ## Build reference
> #set $reference_fasta_filename = "localref.fa"
> ln -s "${ref_file}" "${reference_fasta_filename}" ;
> bwa index -a is "${reference_fasta_filename}" ;
> bwa mem -t "2" -v 1 "${reference_fasta_filename}" "${fastq_input1}"  | samtools view -Sb - > temporary_bam_file.bam ;
> samtools sort -o "${bwa_test_commover_bam_output}" temporary_bam_file.bam
> ```
{: .details}

- There is another sample bwa_test tool that achieves the same results using a bash script.
- It does not need a command over-ride but is more typing because three positional parameters are named for readability.
- Bash is probably more familiar to many ToolFactory users than mako/cheetah templating.
- The effects of templating the command line can usually be achieved in a script, at the expense of needing to script the handling of parameters. It's just logic.
- The `<test>` section allows writing more useful tests than the automatically generated one, or for a real test when using collections as discussed above.
- At present, other embedded logic such as filters in output tags cannot be over-ridden.
  - These and other extensions are welcomed as pull requests if there is a need for them.
  - But if a tool needs them, it probably needs an experienced tool developer.

> <details-title>`bwa_test_toolfactory_positional_bash` sample alternative.</details-title>
> ToolFactory form bash script to replace above command over-ride section:
>
> ```
> REFFILE=$1
> FASTQ=$2
> BAMOUT=$3
> ln -s "$REFFILE" "refalias"
> bwa index -a is "refalias"
> bwa mem -t "2" -v 1 "refalias" "$FASTQ"  | samtools view -Sb - > temporary_bam_file.bam
> samtools sort -o "$BAMOUT" temporary_bam_file.bam
> ```
{: .details}


- Most users may never need to use the command over-ride option to access templating for wrapping their own scripts where they control the command line interpretation.
- Where the generated test is not sufficient, a hand written one can be substituted. Not needed for the simple BWA example.
- Test or command over-rides are likely to be edge cases more suited to the alternate, more powerful tools.

---

# Limits and workarounds

- The ToolFactory is an automated, form based code generator.
- A generator can replace manual editing by a skilled developer only in relatively constrained, simple cases.
- These are common enough in the daily work of most data intensive scientific fields to make a tool generator worth keeping handy.
- For simple scripts and appropriate Conda packages, it's potentially very useful.
- It is not hard to imagine using a Python wrapper to finesse some aspects of more complex tools just as bash was used in the `planemo lint` example.
- Logic in the `<command>` section can probably always be replaced by equivalent code in a script at the cost of extra work compared to templating.
- Other aspects of tool logic such as output filters based on other parameter values can only be implemented in the wrapper document, not in a tool script.
   - The ToolFactory relies on galaxyxml, so those kinds of extensions to galaxyxml will permit extensions to the ToolFactory. Pull requests are welcomed.
- The ToolFactory appliance is a convenient and efficient way to create and maintain simple Galaxy tools from simple working scripts.
- Tools can have command-override and test-override pasted in as in one of the BWA samples.
  - This can work around some of the current limitations on automation.
  - If the package requires that kind of complexity, it might be better to prepare the wrapper manually.

---

# Notes on common issues

## My Rscript generates a strange R error on STDOUT about an invalid operation on a closure called 'args' ?

- Did your code declare the `args vector` with something like `args = commandArgs(trailingOnly=TRUE)` before it tried to access args[1] ?
- See the plotter tool for a sample

## I want to use a collection for outputs but it always passes the test even when the script fails. What gives?

- Collections are tricky for generating tests.
  - The contents appear only after the tool has been run and even then may vary with settings.
- A manual test override is currently the only way to test collections properly.
- Automation is hard. If you can help, pull requests are welcomed.
- Until it's automated, please take a look at the plotter sample.
- It is recommended that you modify the test over-ride that appears in that sample form. Substitute one
or more of the file names you expect to see after the collection is filled by your new tool for the `<element.../>` used in the plotter sample's tool test.

## Only one ToolFactory job runs at a time. Why doesn't the server allow more than one at once?

- When a new dependency is being installed in the Conda repository, there is no locking to prevent a second process from overwriting or otherwise
interfering with it's own independent repository update. The result is not pretty.
- Allowing two tests to run at once has proven to be unstable so the job queue is currently limited to one.

## The tool I just generated is not in the tool menu?

- Did you refresh your browser's Galaxy page? Click on the "Analyse data" or "home" or the left side of the masthead to do that, then open the "ToolFactory Generated Tools" submenu.
It should be there if the XML for the new tool appeared without complaint. Please export and send me a copy of the history so I can see what's going on?


# Appendices - material of potential interest for those interested in details

## ToolFactory functional test workflows

- There are two workflows supplied in the Appliance.
- One will make and install all the test tools and runs pretty fast.
- The other will do the same but also test each one.
  - It is the Appliance functional test and a good stress test for your installation if you need one.
  - It will take a long time to run because there is a deliberate bottleneck as only one planemo process runs at a time.
  - This is necessary to avoid damage to the Planemo Conda installation.
    - When two or more processes try to install new dependencies, Bad Things Happen.
    - Conda does a pretty good job of recovering if it can but it is not pretty. Best avoided.
  - Planemo takes a while because it is doing a lot. Planemo test builds a fresh Galaxy first, and it is run twice - first to generate the test outputs then to do the real test.


## Why is the form so complex?

> <details-title>Summary: details needed and how they are used to generate a new tool</details-title>
>
> #### What information is needed to generate a tool ?
>
> - The code generator requires enough detail to be able to create the appropriate command line
> template to call the script or executable and pass the required file paths and other settings correctly.
> - Small input samples and default settings are used to construct a test for the newly generated tool. These should be known to work with the script, having been used to debug
> the script on the command line.
> - Upload the samples to the current history before
> starting a new tool in the ToolFactory. No tool will be generated without sample inputs. This test becomes part of the XML and of the toolshed archive.
> - The outputs from running the script during the first planemo run become sample outputs to be compared with test outputs in the archive.
>
> - In addition to an ID and name, a tool may have any combination of:
>
>     - Multiple dependencies. Conda is currently supported.
>     - Interpreters such as python, r-base and perl are typically used for ToolFactory tools.
>     - System utilities such as bash and sed can be used.
>     - Argparse (named) or positional (ordered) style parameter passing at tool execution time depending on the script requirements. Positional works well for bash scripts with only a handful of parameters. Argparse is preferred for clarity.
>     - Unlimited individual input data files to be selected from the user's history.
>     - Unlimited individual output files to be written to the user's history, paths determined at tool execution.
>     - Unlimited additional command line parameters that the user can control on the new tool form.
>     - an (optional) script to execute. Running a script to call an executable using parameters passed from the user can be useful to overcome some limitations of the ToolFactory for more complex tools.
>
> - Many of these generate parameter input boxes and history data selects on the new tool form.
> - Metadata about command line formatting together with text strings for the form seen by the user are needed.
>
> - Many of these are components of the generated command line template.
> - This can be seen in the new tool XML. Galaxy file paths for the script are only determined at generated tool execution. The generated template ensures that these are correct.
>
> #### The Galaxy UI imposes additional limits
>
> - The ToolFactory has limited flexibility and works best for simple tools.
> - Even then, the form becomes complicated as more parameters are added.
> - Tools can have unlimited numbers of some items, including input files, output files, citations and user parameters.
> - Each one has half a dozen metadata or text details. Galaxy form repeats are used for those.
> - As more repeats are added, the Galaxy UI becomes increasingly unwieldy.
> - In theory, the Toolfactory can potentially generate very complicated tools with large numbers if inputs, outputs and user modifiable parameters.
> - Great patience would be required.
> - That is why manual methods are likely more productive for complicated requirements.
>
>
{: .details}


## Motivation

- Uptake of Galaxy in new quantitative scientific fields requiring complex computing for analyses is arguably rate-limited early on by the availability of domain specific tools.
   - many scientists routinely write their own analysis code in quantitative disciplines - probably more commonly than in biological domains where Galaxy started.
   - lowering the barriers to those scientists generating their own new tools may speed up adoption of Galaxy in their domains.
- Galaxy tool wrapping has a well established and growing range of project supported infrastructure.
- Much complex Galaxy tool logic is embedded in the tool document namespace.
  - Some document logic cannot be easily automated. This includes output filters and other constructs that are hard to generate automatically.
  - Being in the namespace saves substantial effort compared to the alternative of writing the command  section logic in a script, as provided by the ToolFactory.
  - That extra effort is needed to pass those parameters on the command line and then to parse them in the script so the logic can be implemented.
- Manual templating is far more efficient of developer time and effort, particularly as conditional parameter and tool complexities grow.
  - It is widely preferred by dedicated developers.
  - It is the only way to satisfy some complex tool requirements.
    - An example is an output filter making an output datatype conditional on another parameter - that requires manual code in the output parameter.
  - The syntax requires some getting used to for newcomers.
  - Errors are not always pleasant or convenient to debug in the template.
- The ToolFactory may make it easier for some developers new to Galaxy to begin creating the tools needed for scientists from their new domain.
  - some may prefer a GUI tool form.
  - some may prefer the persistent IDE aspect.
  - some may feel more at home in their favourite scripting language.
  - any scripting language should be capable of implementing equivalent logic.
- Negatives include limited complexity and lower efficiency for complex tools.
  - Specialised manual development tools can make the process far more efficient.
  - ToolFactory form becomes increasingly unwieldy as complexity grows with many parameters and files.
    - sections help some but it requires determination.
- May be a convenient way to learn by building simple tools before diving into the project supported infrastructure.
- It extends the range of options available for creating new tools for new scientific domains.
- Developers can choose the method that suits them best for each new task.
  - For simple tools, the ToolFactory provides a convenient, pop-up, persistent but easily disposable integrated development environment.

---

# Your turn! Please help improve this community resource!
- tell Ross ([ross.lazarus@gmail.com](mailto:ross.lazarus@gmail.com)) what went well and what could be improved for the benefit of future students
- This tutorial has had almost no testing yet.
- Please use the fork at [fubar2/training-material](https://github.com/fubar2/training-material) for issues and pull requests for the nonce.
- The ToolFactory has had little testing in the present form so expect many bugs
- If the ToolFactory works well for you, please tell your friends.
- If you find bugs, please tell me by raising an issue at [fubar2/galaxy_tf_overlay](https://github.com/fubar2/galaxy_tf_overlay), or better, a pull request with the fix :)


# Acknowledgements

This tutorial is based on the work of thousands of individual contributers to the Galaxy project over the last 18 years or so.
Thanks all! It has been a lot of fun.

Special thanks to:

- {% include _includes/contributor-badge.html id="hexylena" %} for
    - review and contribution to the tutorial and associated code.
    - the vision of instant installation of generated tools for developer feedback.
    - elegantly generated lint-free XML provided by [galaxyml code](https://github.com/hexylena/galaxyxml)
- {% include _includes/contributor-badge.html id="mvdbeek" %} for thoughtful comments on the role of the ToolFactory that helped motivate the tutorial.
