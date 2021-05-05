---

layout: tutorial_hands_on
title: "ToolFactory: Generating Tools From More Complex Scripts"
key_points:
  - The ToolFactory Appliance is an automated form driven tool generator
  - It runs in a flavour of docker-galaxy-stable and produces new tools as Toolshed ready archives
  - It was designed for scientists and developers who routinely write scripts for their analyses.
  - It can quickly turn a working command line script into a toolshed-ready archive.
  - It generates tools from information entered on a Galaxy form in the familiar UI.
  - The new tool is installed in the appliance after generation so can be used and tested immediately.
  - Adding a test to the toolshed archive takes a few minutes and is only needed when the new tool is ready for sharing.

objectives:
 - Further develop your ToolFactory Skills

questions:
 - What else can I do with the ToolFactory?

time_estimation: 1H

requirements:
  - type: "internal"
    topic_name: introduction
    tutorials:
      - galaxy-intro-short
      - galaxy-intro-101-everyone
    topic_name: dev
    tutorials:
      - tool-integration
      - interactive-environments
      - tool-generators

contributors:
  - fubar2
  - hexylena

---

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

## The ToolFactory Appliance: A pop-up MYOT (make your own tools) Galaxy for scientists who write command line scripts in their work.

The ToolFactory automates much of the work needed to prepare a new Galaxy tool using information provided by the script writer,
on the ToolFactory form. The ToolFactory can wrap any simple script that runs correctly on the linux command line with some small test input samples. This is potentially
handy for developers new to Galaxy, and for Galaxy users who are capable of correctly scripting on the command line for themselves.

Untested tools are immediately available providing instant feedback for the developer. Wording in labels and help text can be edited by re-running the job. The
newly generated version will have all the text updated so it can be checked again.

The Appliance Galaxy server is fully functional. Tools can be generated or added from the Toolshed to create a
tailored Galaxy for any kind of analysis.

The appliance is an ideal way for any data intensive scientist to quickly develop and refine new tools on their workstations,
ready for deployment in production and sharing.


> ### {% icon tip %} Under the hood:
>
>  - It uses [galaxyml](https://github.com/hexylena/galaxyxml) to generate the tool XML from ToolFactory form settings.
>  - It uses [Planemo](https://github.com/galaxyproject/planemo) to generate the test outputs and then again to test newly generated code
>  - The appliance in Docker is built on [docker-galaxy-stable](https://github.com/bgruening/docker-galaxy-stable/compose) infrastructure
{: .tip}



## Limits and scope

- The Appliance can generate, install and run new tools from scripts on your desktop.
    - However, it is a fully functional Galaxy instance, so it can also import any existing
tool from a Toolshed and [can easily be connected to a cluster for real work](https://github.com/bgruening/docker-galaxy-stable/compose).
- The ToolFactory works best wrapping simple R/Bash/Python and other interpreted scripts, with a few user supplied parameters and a few I/O history files.
- Scripts are easier than some Conda packages
    - This is because the tool builder can modify the code to respond to default empty parameters as if they had not been passed.
    - This may not sound like much, but as a result, advanced tool building elements such as conditionals and related tricks requiring manual coding, can sometimes be avoided.
    - In contrast, some Conda dependencies or combinations will require XML conditionals
or other complex tool XML constructs that are not easy to generate automatically.
    - While some simple requirements may be manageable, complex ones will not be suitable for the ToolFactory.
- Compared to the more usual linux shell and text editor, the ToolFactory appliance is a rather clumsy way to debug scripts.
    - **Starting a new ToolFactory tool with a know good command line and data** is strongly recommended. You will know exactly what to expect from the tool test for a first sanity check.
    - Corrolary: Unless there is a working script that needs to be wrapped into a toolshed-ready Galaxy tool, the ToolFactory is of little use.
- The ToolFactory Appliance is for developers and code-writing scientists not yet familiar with the more flexible and complex manual tools, and who need to wrap scripts that are simple
enough for the ToolFactory.
    - Compared to the more flexible manual Galaxy tool development software, there is far less to learn to get up to speed with a form driven,
automated code generator in a tailored, readily deployed appliance.
    - The cost of this convenience is that ToolFactory is limited to a limited subset of simple script and package wrappers.


# 2. Getting your hands on a ToolFactory Appliance for some hands-on training.

- If you found the introductory material relevant to your own needs, you may wish to start the DIY/hands-on part of the tutorial that follows.
- Install your own ToolFactory Appliance as described below.
- Start exploring the provided samples to figure out if and how it might help your work.
- Tutorial material that follows **can only be completed with a working ToolFactory**.

>#### Active Tutorial content follows
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Installation

> ### {% icon warning %} Security advisory!
>- *Please do not install the ToolFactory on any public server*
>- Although it will only run for administrative users, it allows unlimited scripting and that exposes unwise security weakness for any public facing machine.
>- In fact, Galaxy is very good at isolating tools to stop them doing mischief. But that's no reason to chance your arm. They keep inventing better mice.
>- Please install it locally as described below.
>- For this reason, the training materials can't make use of existing public Galaxy infrastructure like most of the GTN material.
{: .warning}

# Running the ToolFactory

> ### {% icon hands_on %} Hands-on: Launching the Container
>>
>> 1. [Install Docker](https://docs.docker.com/engine/install/) following the appropriate instructions for your platform. Then `pip3 install docker-compose`.
>>
>> 2. Go to [the ToolFactory appliance github repository](https://github.com/fubar2/toolfactory-galaxy-server)
>>
>> 3. Clone it or download the zip and unzip it somewhere handy - such as `~/toolfactory-galaxy-server-main`
>>
>> 4. Change to the compose directory - `cd ~/toolfactory-galaxy-server-main/compose`
>>
>>Something like this should get it started:
>>
>>```
>>git clone https://github.com/fubar2/toolfactory-galaxy-server
>>cd toolfactory-galaxy-server/compose
>>docker-compose pull
>>docker-compose up -d
>>```
>>
>>
>>
>>    > ### {% icon code-in %} Input: Bash
>>    > ```bash
>>    > wget https://github.com/fubar2/toolfactory-galaxy-server/archive/refs/heads/main.zip
>>    > unzip main.zip
>>    > cd toolfactory-galaxy-server-main/compose
>>    > docker-compose pull
>>    > docker-compose up -d
>>    > ```
>>
>>Your appliance should be running with a local Galaxy on  [port 8080 of your workstation](http://localhost:8080) after a fair bit of activity.
>>
>> -  Out of the box login is 'admin@galaxy.org' and the password is 'password'
>>    - This is obviously insecure but convenient and easily changed at first login.
>>    - Or more permanently in the docker-compose.yml if you prefer.
>>
>>- The container `/export` directory is mounted locally at `compose/export` .
>>
>>## Demonstration tools are the functional documentation
>>
>>- At first login to a new Appliance, you will find the demonstration history ready to explore if you waited a minute for all the Conda activity to die down
>>- It takes a minute or two to import because the dependencies for the ToolFactory must first be installed.
>>- If it's not there, you can import it manually from Zenodo as described in the Welcome page text.

>>
>>- To explore an example, open the toolshed archive by clicking on the name, and select the `rerun` button from the expanded view
>>    - The form that generated that tool will appear for you to examine
>>    - Edit the form - add parameters and change the script to suit - and rerun to create an *updated* tool. The history has previous versions.
>>    - Change the tool ID to change the tool name.
>>
>>## To safely shut the appliance down
>> - If the console was not detached using the --detach/-d flag
>>      - `<ctrl><c>` in the console will gracefully shut the server down - takes time but your work will be preserved.
>> - If the -d flag was used,
>>    - `docker-compose down` from the same directory it was started `.../compose`, should shut it down nicely
>>    {: .code-in}
>>    > ### {% icon tip %} Tip: Patience!
>>    > When you run the ToolFactory for the first time inside the container and whenever you run a new tool with new dependencies, it will require some time to build the conda environment.
>>    > Check for Conda or other processes if things seem stuck.
>>    {: .tip}
>>
>>
{: .hands_on}

> ### {% icon hands_on %} Hands-on: Brief Guide to ToolFactory Appliance Operation
>>
>>## Generating new tools - what happens when you press `execute` on a valid ToolFactory form?
>>
>> - The form is processed and a new tool generated.
>> - The new tool is installed to the Appliance.
>>    - The tool generation process takes a few seconds.
>>    - The `Home` or `Analysis` tab should be selected so the screen is refreshed after building.
>>         - Otherwise the new tool menu will not be loaded so the newly generated tool will not be there
>> - Choose the names thoughtfully and be warned: there are no checks on tool names
>> - Any existing installed tool with the same name will be overwritten permanently.
>> - The history will retain all the generating jobs if you accidentally overwrite a tool.
>> - Rerun the job and adjust the form. Rinse and repeat until ready.
>>
>> - Note that the generated tool has not been run to generate test outputs, so the archive is not complete although the installed tool may work fine.
>>
>> - To generate a "proper" tested toolshed archive, use the ToolFactory generated `planemo_test` tool in the ToolFactory section of the tools menu.
>> - It will run Planemo to generate outputs, then run a real test and return a proper toolshed archive and test reports.
>>     - An archive containing the tool with proper test will be returned with the planemo report in the history.
>>     - The archive can be downloaded and shared in the usual ways. It is a normal Galaxy tool that wraps the supplied script
>>
{: .hands_on}

----

## Import ToolFactory functional documentation - the demonstration tools.

- Congratulations on getting this far and acquiring a local instance of the ToolFactory
- There is a sample history built in that shows some sample tools.
- Note that this is installed the first time a new Appliance is run. It takes a minute - wait until all the installation processes have settled down before first logging in.
- If there is an empty history when you first log in, check the histories after a minute - if still not there, follow the Welcome page instructions to install it manually.
- You can examine how these were generated by using the Galaxy job redo button on either of the two history items sharing the tool name.
- This will show you the fully completed ToolFactory form used to generate the sample
- You can edit the form and regenerate a new tool with your changes incorporated.
- When the job completes, a new edited version of that tool will be installed and ready to run.
- Refresh the Galaxy panels to find it in the new "ToolFactory Generated Tools" section.
- Change the tool name to generate a different tool - the Appliance will overwrite any existing tool of the same name without warning.


> ### {% icon announcement %} Note!
> - This is the **first step** recommended
> - Sample ToolFactory tools are the best way to learn how the ToolFactory works and how you might adapt the variations shown in your own work.
> - They provide functional documentation to help you become comfortable using the ToolFactory.
{: .announcement}

- Viewing the samples history, you will see a large number of pairs of history items sharing the same tool name, and some data files used for all the building and tool testing.
- There is an untested toolshed archive for each tool, and a collection. Opening it (click on the name) will reveal the generated tool XML and the future test input file samples.
- Both these items for each tool have a {% icon galaxy-refresh %} rerun button.
    - Click that button and the ToolFactory form that generated the sample tool will appear.
    - You can see how the tool was built using the ToolFactory's limited capacities.
    - Most of them are trivial of course.
    - They are models that show how different kinds of scripts and tasks can be accomplished

<sup id='section3'>*</sup>
# 3. Hands-on: Learning to use the ToolFactory

> ### {% icon tip %} Using an Appliance involves dependency installation that may cause long pauses...
>>- There will be delays as any new dependencies are installed for the first time
>>        - the first ToolFactory run after first starting a new Appliance will involve Conda installing the ToolFactory dependencies before running the job.
>>        - the first time any new tool with a new dependency is run, Conda must install it locally taking a variable amount of time depending on complexity.
>>        - the first time the planemo_test tool is run, there will be a 10+ minute delay as Conda grinds away.
>>        - Check for Conda and other running processes before assuming it has frozen.
{: .tip}


> ### {% icon hands_on %} Exploring the sample tools by regenerating their original ToolFactory forms
>
> * With the ToolFactory working and the sample history active as described above
> * Select any of the generated toolshed archive history items.
> * This should open the item details up, so you can see the circular "redo" button
> * Click that button - the ToolFactory form that generated that tool will appear.
> * Examine the form settings used to generate the tool.
> * Try changing names or prompts. Add new parameters or inputs/outputs; press `execute`; check the new version of the tool
> * For example, change the default for the Hello example to `Galaxy Training Network` and generate an updated version.
{: .hands_on}


The best way to explore the kinds of tasks that can be achieved with simple scripts is to take a look at each sample tool. Note how the various
options have been configured and what kinds of scripts this could be used for in your work. The example script can be swapped out for another one known to work and additional
new parameters added to suit, to extend the toy examples and create tools of use to your users. Change the tool name on the newly edited form, press `execute` and
rerun the job to generate a new toolshed archive and test report collection.

The trivial `Hello World!` tool example is readily extended to suit many situations where a tool is needed quickly for a workflow. Try adding another parameter.
For example, the planemo `lint` tool example (described below) can be derived by adding a history toolshed archive as input, plus a few more lines of bash script.
In practice, it's a flexible basis for generating many simple tools.

> ### {% icon details %} Summary: details needed and how they are used to generate a new tool
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
>     - The Appliance server always exposes them to tools, but this may/should not happen on production servers.
>         - It is recommended that they be added as dependencies after testing locally before testing and export of production ready toolshed archives.
>         - They can be added to the tool form like other Conda dependencies
>         - This will ensure that they are available when the script runs. Bash, sed and so on are all available in Conda.
>         - Versions may not matter as much as other packages. Latest will be used by default but a version can be specified to ensure reproducibility.
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


#### Workflow used to create the tools in the demonstration history

The workflow at https://zenodo.org/record/4686436/files/TFdemo_wf_april13_planemo.ga?download=1 was used to create the sample history examples. It requires some data files as inputs. They
can be copied from the imported history to a new history and the workflow can be imported and run from there. After connecting appropriate data sets to the different inputs, it will
re-create all the samples for anyone wanting to see them run.

## ToolFactory tips and tricks illustrated by some of the examples.

#### Before you begin a new tool

- Make sure it runs correctly on a command line with your sample inputs and default parameter settings. That is how the test will be generated.
- You cannot specify any inputs on the form without providing samples, and those samples must run correctly, with the supplied defaults.
- Easiest to make sure those samples are in the history before you begin.
- Be well prepared. The ToolFactory cannot do that preparation for  you.

#### STDIN and STDOUT

- Demonstration tools often capture output from a bash script using the special STDOUT designation for output files
- This can save sending the output path as a parameter to the script or executable
- STDIN is also available as a special designation for history inputs if the script takes input from STDIN when it runs.
- The option `Tool reads selected input file from STDIN and writes STDOUT with no parameters` in the parameter passing model
selection list will generate a simple filter tool, executing a script or Conda dependency with input on STDIN and output on STDOUT
for those rare situations where that's all you need. No i/o or other parameters for the user to set. Used in the tacrev demonstration tool.

#### Repeats in the ToolFactory form theoretically permit *any* number of parameters.

- Inputs, outputs and additional user parameters are all presented within repeat elements on the form.
- For any given tool, repeats allow an unlimited number of each of those tool elements - from a technical perspective.
- The practical limit is reached when your patience runs out as you add them in the current Galaxy forms UI.
- The form gets very long, very quickly.
- A handful is easily manageable.

#### Repeated parameters in generated tools allow the end user to supply multiple values

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
- The demonstration trivially returns whatever the user chose to repeat.


> ### {% icon details %} Repeats demonstration tool XML
> >
> >
> >```xml
> >  <tool name="tool1" id="tool1" version="0.01">
> >    <!--Source in git at: https://github.com/fubar2/toolfactory-->
> >    <!--Created by planemo@galaxyproject.org at 07/04/2021 14:55:09 using the Galaxy Tool Factory.-->
> >    <description>test repeats</description>
> >    <expand macro="requirements"/>
> >    <stdio>
> >      <exit_code range="1:" level="fatal"/>
> >    </stdio>
> >    <expand macro="stdio"/>
> >    <version_command><![CDATA[echo "0.01"]]></version_command>
> >    <command><![CDATA[python
> >  $runme
> >   #for $rep in $R_repeat:
> >  --repeat $rep.repeat
> >  #end for
> >  >
> >  $rep_out]]></command>
> >    <configfiles>
> >      <configfile name="runme"><![CDATA[
> >  import argparse
> >  parser = argparse.ArgumentParser()
> >  a = parser.add_argument
> >  a('--repeat',default=[],action="append")
> >  args = parser.parse_args()
> >  s = ' and '.join(args.repeat)
> >  print(s)
> >  ]]></configfile>
> >    </configfiles>
> >    <inputs>
> >      <repeat name="R_repeat" title="Add as many Things to pass as needed">
> >        <param name="repeat" type="text" value="add lots of repeats" label="Things to pass" help=""/>
> >      </repeat>
> >    </inputs>
> >    <outputs>
> >      <data name="rep_out" format="txt" label="rep_out" hidden="false"/>
> >    </outputs>
> >    <tests>
> >      <test>
> >        <output name="rep_out" value="rep_out_sample" compare="diff" lines_diff="0"/>
> >        <repeat name="R_repeat">
> >          <param name="repeat" value="add lots of repeats"/>
> >        </repeat>
> >      </test>
> >    </tests>
> >    <help><![CDATA[
> >
> >  **What it Does**
> >
> >
> >
> >  ------
> >
> >
> >  Script::
> >
> >      import argparse
> >      parser = argparse.ArgumentParser()
> >      a = parser.add_argument
> >      a('--repeat',default=[],action="append")
> >      args = parser.parse_args()
> >      s = ' and '.join(args.repeat)
> >      print(s)
> >
> >  ]]></help>
> >    <citations>
> >      <citation type="doi">10.1093/bioinformatics/bts573</citation>
> >    </citations>
> >  </tool>
> >```
> > - The user sees the following form after adding 3 repeats for each of the two available items
> >
> > ![Form to configure an additional parameter as a select](../../images/toolfactory_repeats_sample_form.png)
{: .details}

#### ToolFactory `collection` outputs are handy for hiding dozens of miscellaneous tool outputs in a single history item

- If the script writes a large number of different output file types (images, reports..) that are not in themselves useful for downstream analyses.
- A Collection can be used to hide them in an organised list as a single new history item.
- Collections are special kinds of Galaxy history items that can present a directory full of files arranged inside one single new history item.
- When viewed, that history item presents them all as a list, each a typical, viewable Galaxy history item.
- They don't clutter the user's history because they are all hidden in the collection.
- The plotter example uses an Rscript.
- It generates as many pairs of random plots as the parameter supplied requires.
- The script sends them into the the collection that appears in the history after the job runs.
- The user's history shows only one new item after it runs - it must be viewed to see all the individual contents.

> ### {% icon warning %} The default generated test for output collections always passes because it doesn't test anything.
>
>    - Supplying a test over-ride is recommended for collections.
>    - Example code is shown on the sample tool's form and in the original example code below - removed from the current sample.
>    - For a real test, one or more expected <element.../> tags must be provided so the test really does test something.
>    - Add your own file names to the sample to make a real test for your own collections.
>    - Otherwise, the generated test will be empty, because that is what the code generator knows about what's in the collection when the `<test>` code is generated. Nothing.
>    - Introspecting arbitrary scripts to reliably populate the test with actual file names?. Not likely any time soon.
>    - Introspecting the Planemo test result to write the test and then retest might be possible but is not planned.
>
{: .warning}

> ### {% icon details %} `plotter` collection output demonstration tool form, generated XML and outputs
> >
> > - The ToolFactory form for the plotter example tool
> > is configured as shown below, from "rerunning" the plotter job from the sample history.
> >
> >![ToolFactory form configuration of the output collection in the plotter example](../../images/toolfactory_plotter_demo_form.png)
> >
> > The Rscript is contained in a configfile so`#` is escaped - this is automatic.
> >
> >```xml
> ><tool name="plotter" id="plotter" version="0.01">
> >  <!--Source in git at: https://github.com/fubar2/toolfactory-->
> >  <!--Created by admin@galaxy.org at 24/01/2021 05:02:33 using the Galaxy Tool Factory.-->
> >  <description>ToolFactory collection demonstration - random plots</description>
> >  <requirements>
> >    <requirement version="" type="package">r-base</requirement>
> >  </requirements>
> >  <stdio>
> >    <exit_code range="1:" level="fatal"/>
> >  </stdio>
> >  <version_command><![CDATA[echo "0.01"]]></version_command>
> >  <command><![CDATA[Rscript
> >$runme
> >"$nplot"]]></command>
> >  <configfiles>
> >    <configfile name="runme"><![CDATA[
> >\# demo
> >args = commandArgs(trailingOnly=TRUE)
> >if (length(args)==0) {
> >   n_plots = 3
> >} else {
> >   n_plots = as.integer(args[1]) }
> >dir.create('plots')
> >for (i in 1:n_plots) {
> >    foo = runif(100)
> >    bar = rnorm(100)
> >    bar = foo + 0.05*bar
> >    pdf(paste('plots/yet',i,"anotherplot.pdf",sep='_'))
> >    plot(foo,bar,main=paste("Foo by Bar plot \#",i),col="maroon", pch=3,cex=0.6)
> >    dev.off()
> >    foo = data.frame(a=runif(100),b=runif(100),c=runif(100),d=runif(100),e=runif(100),f=runif(100))
> >    bar = as.matrix(foo)
> >    pdf(paste('plots/yet',i,"anotherheatmap.pdf",sep='_'))
> >    heatmap(bar,main='Random Heatmap')
> >    dev.off()
> >}
> >
> >]]></configfile>
> >  </configfiles>
> >  <inputs>
> >    <param label="Number of random plots pairs to draw" help="" value="3" type="text" name="nplot" argument="nplot"/>
> >  </inputs>
> >  <outputs>
> >    <collection name="plots" type="list" label="Plots">
> >      <discover_datasets pattern="__name_and_ext__" directory="plots" visible="false"/>
> >    </collection>
> >  </outputs>
> >
> >
> >  <tests>
> >    <test>
> >      <param name="nplot" value="3" />
> >      <output_collection name="plots" type="list">
> >     <element file="yet_1_anotherplot_sample" name="yet_1_anotherplot" ftype="pdf" compare="sim_size" delta_frac="0.05"/>
> >    </output_collection>
> > </test>
> >  </tests>
>>
>>
> >
> >  <help><![CDATA[
> >
> >**What it Does**
> >
> >Draws as many random plot pairs as you need
> >
> >
> >
> >------
> >
> >
> >Script::
> >
> >    # demo
> >    args = commandArgs(trailingOnly=TRUE)
> >    if (length(args)==0) {
> >       n_plots = 3
> >    } else {
> >       n_plots = as.integer(args[1]) }
> >    dir.create('plots')
> >    for (i in 1:n_plots) {
> >        foo = runif(100)
> >        bar = rnorm(100)
> >        bar = foo + 0.05*bar
> >        pdf(paste('plots/yet',i,"anotherplot.pdf",sep='_'))
> >        plot(foo,bar,main=paste("Foo by Bar plot #",i),col="maroon", pch=3,cex=0.6)
> >        dev.off()
> >        foo = data.frame(a=runif(100),b=runif(100),c=runif(100),d=runif(100),e=runif(100),f=runif(100))
> >        bar = as.matrix(foo)
> >        pdf(paste('plots/yet',i,"anotherheatmap.pdf",sep='_'))
> >        heatmap(bar,main='Random Heatmap')
> >        dev.off()
> >    }
> >
> >]]></help>
> >  <citations>
> >    <citation type="doi">10.1093/bioinformatics/bts573</citation>
> >  </citations>
> ></tool>
>>```
>> After requesting 25 pairs of plots from the sample tool, a collection appears in the history and is shown below.
>> One of them is displayed by clicking the "eye" icon.
>> Collections are ideal for messy analysis reporting outputs such as images, pdfs and other material that is not useful as an input to a downstream tool.
>> It is material that the user will want kept together, so a single history item is ideal to avoid unnecessary clutter.
>> As shown above, the script only has to write the files to a directory.
>> Note that the test is over-ridden in the ToolFactory form to generate this tool.
>> Without an over-ride to suit your script, an empty test will be generated that will not do any testing - see warning above.
>>
>>
>>![Plotter collection opens to show all the plots as separate items to be viewed or downloaded](../../images/toolfactory_plotter_sample_output.png)
{: .details}

#### Selects as user supplied parameters

- Additional parameter types are selected from a drop down list. It includes text, numeric and select parameters.
- Selects offer a list of prespecified options to the user.
- There is a repeat on the ToolFactory form to collect options as pairs of names/values.
- It is clumsy and suitable only for a limited number of options.
- Galaxyxml generates appropriate select parameters on the generated tool as shown in the select demonstration tool.

> ### {% icon details %} `select_test` select field demonstration tool generated XML
>>
>>The ToolFactory form section for user configurable command line settings is
>> configured as shown here for the select demonstration
>>
>>![Form to configure an additional parameter as a select](../../images/toolfactory_select_demo_form.png)
> >
> >The generated XML is shown below.
> >
> >```xml
> ><tool name="select_test" id="select_test" version="0.01">
> >  <!--Source in git at: https://github.com/fubar2/toolfactory-->
> >  <!--Created by admin@galaxy.org at 24/01/2021 05:03:21 using the Galaxy Tool Factory.-->
> >  <description>ToolFactory select demonstration</description>
> >  <stdio>
> >    <exit_code range="1:" level="fatal"/>
> >  </stdio>
> >  <version_command><![CDATA[echo "0.01"]]></version_command>
> >  <command><![CDATA[bash
> >$runme
> >"$choose"
> >>
> >$select_out]]></command>
> >  <configfiles>
> >    <configfile name="runme"><![CDATA[
> >echo "You chose \$1"
> >]]></configfile>
> >  </configfiles>
> >  <inputs>
> >    <param label="Choose" help="" type="select" name="choose" argument="choose">
> >      <option value="won">one</option>
> >      <option value="too">two</option>
> >      <option value="free">three</option>
> >    </param>
> >  </inputs>
> >  <outputs>
> >    <data name="select_out" format="txt" label="select_out" hidden="false"/>
> >  </outputs>
> >  <tests>
> >    <test>
> >      <output name="select_out" value="select_out_sample" compare="diff" lines_diff="0"/>
> >      <param name="choose" value="won"/>
> >    </test>
> >  </tests>
> >  <help><![CDATA[
> >
> >**What it Does**
> >
> >Echoes your selection
> >
> >
> >
> >------
> >
> >
> >Script::
> >
> >    echo "You chose $1"
> >
> >]]></help>
> >  <citations>
> >    <citation type="doi">10.1093/bioinformatics/bts573</citation>
> >  </citations>
> ></tool>
> >```
> > The generated tool form from the select demonstration shows the three options and returns the one selected.
> >
> >![Generated form seen by users of the select demonstration tool](../../images/toolfactory_select_test_tool.png)
{: .details}



#### The ToolFactory can wrap some Conda packages correctly using a simple wrapper script to work around limitations.

- There are two demonstration tools that use Planemo as a Conda dependency
- One runs `planemo test...` and the other `planemo lint....` on toolshed archives in a history
- The linter XML is available below. It's a variant of the Hello example in using a bash script.
- Instead of echo "Hello $1", it takes an input toolshed archive and writes the planemo lint output to STDOUT with this script pasted into the box:

```
cp \$1 foo.tar
tar -xvf foo.tar
TOOLNAME=`find . -name "*.xml"`
planemo lint $TOOLNAME >> $2
```

- The ToolFactory makes exposing these Planemo functions as Galaxy tools fairly easy. Planemo test will not work.
- Similarly tractable Conda dependencies are also potential candidates for being quickly wrapped as tools

> ### {% icon details %} `planemo lint` demonstration tool generated XML
> >```xml
> ><tool name="planemo_lint" id="planemo_lint" version="0.01">
> >  <!--Source in git at: https://github.com/fubar2/toolfactory-->
> >  <!--Created by planemo@galaxyproject.org at 08/01/2021 17:34:35 using the Galaxy Tool Factory.-->
> >  <description>Lints a ToolFactory or other xml using planemo</description>
> >  <requirements>
> >    <requirement version="0.74.1" type="package">planemo</requirement>
> >  </requirements>
> >  <stdio>
> >    <exit_code range="1:" level="fatal"/>
> >  </stdio>
> >  <version_command><![CDATA[echo "0.01"]]></version_command>
> >  <command><![CDATA[bash
> >$runme
> >$input1
> >$lint_output]]></command>
> >  <configfiles>
> >    <configfile name="runme"><![CDATA[
> >cp \$1 foo.tar
> >tar -xvf foo.tar
> >TOOLNAME=`find . -name "*.xml"`
> >planemo lint \$TOOLNAME >> \$2
> >]]></configfile>
> >  </configfiles>
> >  <inputs>
> >    <param optional="false" label="Toolshed archive to be linted" help="" format="tgz" multiple="false" type="data" name="input1" argument="input1"/>
> >  </inputs>
> >  <outputs>
> >    <data name="lint_output" format="txt" label="lint_output" hidden="false"/>
> >  </outputs>
> >  <tests>
> >    <test>
> >      <output name="lint_output" value="lint_output_sample" compare="diff" lines_diff="5"/>
> >      <param name="input1" value="input1_sample"/>
> >    </test>
> >  </tests>
> >  <help><![CDATA[
> >
> >*What it Does**
> >
> >planemo lint
> >
> >-----
> >
> >Script::
> >
> >  tar -xvf $1
> >  TOOLNAME=`find . -name "*.xml"`
> >  echo "$$$$$TOOLNAME = $TOOLNAME" > $2
> >  planemo lint $TOOLNAME >> $2
> >
> >]></help>
> ><citations>
> >  <citation type="doi">10.1093/bioinformatics/bts573</citation>
> ></citations>
> >/tool>
>>
>>```
{: .details}


#### Other languages such as Lisp and Prolog scripts in Galaxy.

- Any interpreted language in Conda can be used, as long as it will take a command line script path or read from STDIN.
- Trivial samples using Lisp and Prolog are included in the demonstration history.
- They may be old, but they show how it is easy to use any Conda scripting language package in Galaxy for long running applications.
- The Prolog tool has an inbuilt script. Substitute the sample script for real code and add inputs, user-configurable parameters and outputs to produce a tool wrapping a Prolog script if you ever need one in Galaxy.
- The Lisp tool will try to execute the selected input text file, hoping it contains a Lisp program like the `hello world` example supplied. This is a terrible idea for a public server and is shown only as a skeletal example of what is possible, not what is sensible.

> ### {% icon details %} `Prolog` and `Lisp` demonstration tools
> >```xml
> ><tool name="prolog_demo" id="prolog_demo" version="0.01">
> >  <!--Source in git at: https://github.com/fubar2/toolfactory-->
> >  <!--Created by admin@galaxy.org at 04/04/2021 17:21:32 using the Galaxy Tool Factory.-->
> >  <description>Runs a prolog script</description>
> >  <requirements>
> >    <requirement version="" type="package">swi-prolog</requirement>
> >  </requirements>
> >  <stdio>
> >    <exit_code range="1:" level="fatal"/>
> >  </stdio>
> >  <version_command><![CDATA[echo "0.01"]]></version_command>
> >  <command><![CDATA[swipl
> >-q
> >-g
> >main
> >-s
> >
> >$runme
> >>
> >$prolog_out]]></command>
> >  <configfiles>
> >    <configfile name="runme"><![CDATA[
> >parent(pam,bob).
> >parent(tom,bob).
> >parent(tom,liz).
> >parent(bob,ann).
> >parent(bob,pat).
> >parent(pat,jim).
> >
> >main :-
> >    parent(X,jim),
> >    format('~a is the parent of jim~n', [X]),
> >    halt.
> >]]></configfile>
> >  </configfiles>
> >  <inputs/>
> >  <outputs>
> >    <data name="prolog_out" format="txt" label="prolog_out" hidden="false"/>
> >  </outputs>
> >  <tests>
> >    <test>
> >      <output name="prolog_out" value="prolog_out_sample" compare="diff" lines_diff="0"/>
> >    </test>
> >  </tests>
> >  <help><![CDATA[
> >
> >**What it Does**
> >
> >
> >
> >Prolog demonstration in the ToolFactory
> >
> >
> >
> >------
> >
> >
> >Script::
> >
> >    parent(pam,bob).
> >    parent(tom,bob).
> >    parent(tom,liz).
> >    parent(bob,ann).
> >    parent(bob,pat).
> >    parent(pat,jim).
> >    main :-
> >        parent(X,jim),
> >        format('~a is the parent of jim~n', [X]),
> >        halt.
> >
> >]]></help>
> >  <citations>
> >    <citation type="doi">10.1093/bioinformatics/bts573</citation>
> >  </citations>
> ></tool>
> >```
> >
> > Lisp too!
> > In this case, the tool takes a lisp script as input
> >
> > A lisp expression was typed into the upload tool text box and saved as hello_lisp.txt ```(write-line "Hello, ToolFactory does Lisp!")```.
> > It or any other lisp script can be used as input to the SBCL Conda dependency.
> >
> >```xml
> ><tool name="lisp_demo" id="lisp_demo" version="0.01">
> >  <!--Source in git at: https://github.com/fubar2/toolfactory-->
> >  <!--Created by admin@galaxy.org at 04/04/2021 16:38:47 using the Galaxy Tool Factory.-->
> >  <description>Runs SBCL hello world demonstration</description>
> >  <requirements>
> >    <requirement version="" type="package">sbcl</requirement>
> >  </requirements>
> >  <stdio>
> >    <exit_code range="1:" level="fatal"/>
> >  </stdio>
> >  <version_command><![CDATA[echo "0.01"]]></version_command>
> >  <command><![CDATA[sbcl
> >--script
> >$script
> >>
> >$lisp_out]]></command>
> >  <inputs>
> >    <param optional="false" label="SBCL lisp script to execute" help="" format="txt" multiple="false" type="data" name="script" argument="--script"/>
> >  </inputs>
> >  <outputs>
> >    <data name="lisp_out" format="txt" label="lisp_out" hidden="false"/>
> >  </outputs>
> >  <tests>
> >    <test>
> >      <output name="lisp_out" value="lisp_out_sample" compare="diff" lines_diff="0"/>
> >      <param name="script" value="script_sample"/>
> >    </test>
> >  </tests>
> >  <help><![CDATA[
> >
> >**What it Does**
> >
> >Lisp in the ToolFactory
> >
> > ]]></help>
> >  <citations>
> >    <citation type="doi">10.1093/bioinformatics/bts573</citation>
> >  </citations>
> ></tool>
> >```
{: .details}


#### Command over-ride and test over-ride

- There are two simple BWA wrappers based on a Planemo documentation advanced example
- One uses a command over-ride pasted into the appropriate text box on the ToolFactory form
- This was based on the one shown in the example it is copied from.
- It allows templating - `$foo` is interpreted by mako.
- The pasted over-ride completely replaces the galaxyxml generated ones.

> ### {% icon details %} `bwa_test_command_override` sample - the command override
>> ToolFactory command over-ride section adapted from the Planemo BWA example.
>>
>>```
>> ## Build reference
>>#set $reference_fasta_filename = "localref.fa"
>>ln -s "${ref_file}" "${reference_fasta_filename}" ;
>>bwa index -a is "${reference_fasta_filename}" ;
>>bwa mem -t "\${GALAXY_SLOTS:-4}" -v 1 "${reference_fasta_filename}" "${fastq_input1}"  | samtools view -Sb - > temporary_bam_file.bam ;
>>samtools sort -o "${bwa_test_commover_bam_output}" temporary_bam_file.bam
>>```
{: .details}

- There is another sample bwa_test tool that achieves the same results using a bash script.
- It does not need a command over-ride but is more typing because three positional parameters are named for readability.
- Bash is probably more familiar to many ToolFactory users than mako templating.
- The effects of templating the command line can usually be achieved using bash or Python at the expense of needing to script the handling of parameters.

> ### {% icon details %} `bwa_test_toolfactory_positional_bash` sample alternative.
>> ToolFactory form bash script to replace above command over-ride section:
>>
>>```
>>REFFILE=$1
>>FASTQ=$2
>>BAMOUT=$3
>>rm -f "refalias"
>>ln -s "$REFFILE" "refalias"
>>bwa index -a is "refalias"
>>bwa mem -t "2"  -v 1 "refalias" "$FASTQ"  > tempsam
>>samtools view -Sb tempsam > temporary_bam_file.bam
>>samtools sort -o "$BAMOUT" temporary_bam_file.bam
>>```
{: .details}


- Most users may never need to use the command over-ride option to access templating for wrapping their own scripts where they control the command line interpretation.
- Where the generated test is not sufficient, a hand written one can be substituted. Not needed for the simple BWA example.
- Test or command over-rides are likely to be edge cases more suited to the alternate, more powerful tools.


## Limits and workarounds

- The ToolFactory is an automated, form based code generator.
- A generator can replace manual editing by a skilled developer only in relatively constrained, simple cases.
- These are common enough in the daily work of most data intensive scientific fields to make a tool generator potentially worth keeping handy.
- For simple scripts and appropriate Conda packages, it's potentially very useful.
- It is not hard to imagine using a Python wrapper to finesse more complex tools just as bash was used in the `planemo lint` example.
- The ToolFactory appliance is a convenient and efficient way to create and maintain Galaxy tools from working scripts.
- Tools can have command-override and test-override pasted in as in one of the BWA samples. This can solve some of the limitations. However, if the package requires that kind of complexity, it might be better to prepare the wrapper manually.


## Notes on some commonly reported issues

#### First Appliance job I submitted container remains grey or running for a long time - is it broken?

- Check with `top` or your system monitor - if Conda is running, things are working but it's slow the first time a dependency is installed.
- The first run generally takes a while to install all the needed dependencies.
- Subsequent runs should start immediately
- Installing new Conda dependencies also takes time so tools that have new Conda packages will take longer to generate as they must be installed before the tool can be tested.
- In general, a complete ToolFactory job usually takes less than a minute - planemo has to build and tear down a new Galaxy for generating test results and then
again for testing properly. Longer if the tool has Conda dependencies....

#### My Rscript generates a strange R error on STDOUT about an invalid operation on a closure called 'args' ?

- Does your code create the args vector at the start of the script with something like `args = commandArgs(trailingOnly=TRUE)` before Rscript tries to access args[1] ?

#### I want to use a collection for outputs but it always passes the test even when the script fails. What gives?

- Unfortunately, collections are tricky to generate automated tests for. The contents are never known until the tool has been run.
- A manual test override is currently the only way to test collections properly.
- Automation is hard. If you can help, pull requests are welcomed.
- Until it's automated, please take a look at the plotter sample.
- It is recommended that you modify the test over-ride that appears in that sample form. Substitute one
or more of the file names you expect to see after the collection is filled by your new tool for the `<element.../>` used in the plotter sample's tool test.



# Your turn! Please help improve this community resource!
- tell Ross (ross.lazarus@gmail.com) what went well and what could be improved for the benefit of future students
- This tutorial has had almost no testing yet.
- It is in pre-review.
- A PR will soon appear once we get the showstoppers sorted.
- Please use the fork at https://github.com/fubar2/training-material for issues and pull requests for the nonce.
- The ToolFactory has had little testing in the present form so expect many bugs
- If the ToolFactory works well for you, please tell your friends.
- If you find bugs, please tell me by raising an issue at https://github.com/fubar2/toolfactory, or better, a pull request with the fix :)
