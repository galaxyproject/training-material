---

layout: tutorial_hands_on
title: "ToolFactory: Generating Tools From Simple Scripts"
key_points:
  - The ToolFactory is an automated Galaxy tool generator for scientists and developers who routinely write command line scripts.
  - It can turn a working command line script into a proper Galaxy tool with a test in a few minutes.
  - It automatically generates simple, complete Galaxy tools from information provided by filling in a normal Galaxy form in the familiar UI.
  - It is installed as a self-configuring development code clone.
  - The Galaxy server is suited any modern Linux workstation or high-end laptop.
  - Useful for learning about system administration or framework code, and for developing tools - all on your own private pop-up server.

objectives:
 - Learn why you might want to use the ToolFactory development server
 - Watch a video demonstration and explore the generated code - Hello Galaxy Training Network!
 - Run it locally in Docker.
 - Install and explore the simple examples provided.
 - Modify and re-generate them to see how the changes affect the tool
 - Generate new simple Galaxy tools using your own scripts

questions:
 - What options exist for new-to-Galaxy developers to convert functioning command line scripts into Galaxy tools?
 - Can any command line script I've written be wrapped as a Galaxy tool?
 - How can I get the ToolFactory working locally, since you tell me it should never be exposed on a public server?
 - What is the difference between a hand-crafted tool and a ToolFactory generated one?

time_estimation: 1H
subtopic: tooldev

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

follow_up_training:
  - type: "internal"
    topic_name: dev
    tutorials:
      - tool-generators-advanced

contributors:
  - fubar2
  - hexylena
---

The Toolfactory and these tutorials are for developers and researchers learning about Galaxy, who routinely develop their own analysis scripts using
bash, Python, Perl, Rscript or other common scientific scripting languages. The tutorials show a convenient way to bridge the gap between a
working command line script and a new tool that "wraps" that script so users can use it like any other tool in Galaxy.

The ToolFactory is a Galaxy tool. Tools are constructed through the normal Galaxy interface when it is run. It is distributed as a self-installing
configuration on a freshly cloned copy of the Galaxy source code. Generated tools are *immediately installed* and ready to run so you can see
what the end user will see. Jobs can be re-run to edit and update generated tools, so Galaxy becomes an integrated development environment for Galaxy tools.

A [`Hello Galaxy!` demonstration](https://youtu.be//DK1eKz5TRs4) using the ToolFactory is available if you'd like to see a walk-through of some of
the hands-on material in this tutorial. You can see whether it looks useful for your work and decide whether to read the material below.

This first tutorial is a slow introduction. For some developers, it may be too slow and the second tutorial may be a better place to start. This one steps
in some detail through the process of using the ToolFactory
to generate `Hello World!` style simple demonstration Galaxy tools.

The reader will soon learn if it might be adapted to their work. If so, an [Advanced ToolFactory tutorial]({% link topics/dev/tutorials/tool-generators-advanced/tutorial.md %}) is
available if the material here is relevant to your needs and you would like to learn more details about the different kinds of tools and features the ToolFactory offers.

Experienced galaxy tool developers already have specialised tools and training so may not gain much from this material.
Users new to Galaxy from other scientific disciplines not yet familiar with the manual tool development process,
may find the ToolFactory appliance useful for familiarising themselves with tool development in Galaxy.


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Tools, tool wrappers and the ToolFactory in Galaxy.

Tools are the fundamental building blocks for analyses in Galaxy. Thousands are available in the [project tool library](https://toolshed.g2.bx.psu.edu) from
many different kinds of science. Galaxy
tools can be created from almost any Linux command line software packages or interpreter scripts. What happens
at tool execution is decoupled from the framework itself, so the framework is indifferent to scientific discipline or coding language.
This is one reason for rapid uptake in new data intensive scientific areas. In many cases, growing sets of discipline specific tools are attracting
whole new communities of scientists into the Galaxy community.

Tool execution is tightly constrained and secured. The only things that can be changed
before execution are tool form exposed user-controlled settings, and data input selections. Everything else is fixed.
This is ideal for non-programmer Galaxy users who rely on prepared tools for their work. Tools are
utterly predictable, to the point of being about as reproducible as any useful complex computing component is likely to be.
This is one of the strengths of the Galaxy framework for users requiring reproducible and shareable scientific analyses.

## Galaxy Interactive Tools and Interactive Environments (GxIT/GIE)

In contrast to tools, GIE allow unconstrained scripting in a Galaxy environment. They offer complete freedom from the constraints of existing tools for
appropriately skilled researchers and developers, because they allow code to run inside Galaxy that is not available in any existing tool.
These can be shared and reused and GIE can now be run in workflows.

If a shareable tool that performs exactly the same functions as a GIE script is needed, the code can be extracted and turned into a parameterised command line script.
Any functional script can be turned into a typical Galaxy tool.

## Pathways from scripts to tools

Any command line script that runs correctly with sample test data input files and parameter settings can be turned into a Galaxy tool.
These may be derived from notebook scripts that have been consolidated and rewritten to take positional or argparse style parameters and tested
in a shell using small data input files. Alternatively, skilled users can develop scripts and test them using small input data files on the command line
without using Galaxy.

The Galaxy developers support extensive software infrastructure for manually creating new tools including Planemo and the new Galaxy language server.
These are complex and powerful with substantial
learning curves but can be used to turn almost any software package into a tool.

For those new to Galaxy, in many simple cases, it may be possible to generate a new tool "wrapping" a script in a few minutes, using an XML code generator for tool wrappers
in a specialised Galaxy tool for developers.

---

# The ToolFactory development server

The ToolFactory implements an automated, form driven XML document generator, and an installer for newly generated tools so you can try them straight away
in Galaxy. The ToolFactory can be "popped up" as a docker container conveniently but all work must be exported and saved before shutting down because no changes are persisted in the docker image. It can be installed locally by cloning and running the setup script from the git repository - ready to run in about 20 minutes as a fully functional, local, throw-away Galaxy development server.

> <tip-title>The ToolFactory installs and configures a fresh development server clone.</tip-title>
> - Private local desktop Galaxy server or docker image
>     - ideal for tinkering and experimentation
>     - learning how the Galaxy server works and
>     - developing new tools for new kinds of scientists.
> - Any Galaxy tool from the toolshed can be installed and used.
> - Simple scripts can have tool wrapper XML generated, and installed for testing.
> - They appear as a shareable Toolshed ready archive in the history.
> - Newly generated tools appear in the tool menu after a refresh, and can be viewed and used as the user will see them.
> - Tool generation jobs can be rerun using the {% icon galaxy-refresh %} button on the history item after editing the form to make changes to the tool the user will see in Galaxy.
> - The development server is a self-installing clone of the Galaxy code, configured with the ToolFactory and sample tools to explore.
>    - It can be backed up and persisted for as long as required, or it can be treated as a throw-away instance and deleted when no longer needed.
>    - There is almost zero technical friction. Only time is required to initialise the server.
>    - On a modern workstation or well-endowed laptop, it will perform well out of the box.
>    - It is suitable only for development in a private deployment.
>    - Please, never expose as a public server.
{: .tip }

The server was developed for programmers who need scripts they write turned in to new Galaxy tools for their own use and if sufficiently useful,
for others to share. Any user comfortable with scientific or general scripting languages on a Linux command
line may find it useful if they ever need a Galaxy tool that wraps a working script. Linux command line utilities and scripting language interpreters supported by Conda can be used.
Some Conda packages can also be used without a script, but the focus is on scripts.

Generated tools pass Planemo lint and are functionally indistinguishable from equivalent
manually written tools. The tested toolshed archives contain a test based on the test data provided at tool generation.

Working demonstration scripts are provided that use bash, Python, Rscript, Lisp, Prolog, Perl, sed, BWA and samtools, as described below.
Many demonstrate ToolFactory features and all can be updated and changed easily, supporting learning by experimenting.
More useful tools can be developed using more complex scripts and as many inputs, outputs
and user supplied parameters as needed by that script. Note that many tool complexities are not easily automated, so the
XML generator provides only limited features. Tools using those limited features may still be useful in many situations but a specialised
tool developer will be needed for many requirements.

If you are a scientist/programmer or software developer new to Galaxy and new to the dark arts of Galaxy tool building, this tutorial may be of help.
It introduces an automated way to convert any useful script into a toolshed ready tool, quickly *inside* Galaxy.


> <tip-title>Alternative ways to generate and see tools:</tip-title>
> - The [Galaxy Language Server](https://github.com/galaxyproject/galaxy-language-server)
>   - Undergoing rapid development.
>   - Specialised semi-automated tool building environment with VCS bindings.
> - Planemo can manually [generate and serve tool XML](https://planemo.readthedocs.io/en/latest/writing_standalone.html) with an optional test.
>    - Recommended for developers who will focus on building Galaxy tools on the command line.
>    - Outstanding documentation.
>    - Widely used. Requires relatively little time to figure out - Galaxy tool syntax takes longer.
>    - No GUI, although can serve tools on a web port.
>    - Command line only. Can create archives with additional steps.
>    - Need to pass all i/o and parameter details at once on the command line.
>    - Takes longer to learn to use and less accessible to many users than a form driven GUI might be.
>    - Manual XML editing required for selects and collections.
>    - See the recommended next steps at the end of this tutorial for Planemo related training.
> - Many Unix utilities (sed, awk...) are already available as IUC tools.
>    - They are `generic` in the sense that a script is supplied at **run time** by the tool user.
>    - This is possible with the ToolFactory, but for reproducible workflows, a specific script can be permanently built-in to make a reproducible tool.
>        - The user supplies only those I/O and parameter values specified on the ToolFactory form.
>        - Nothing else can be changed - just like with most Galaxy tools.
> - Choose whichever one fits best for the task at hand.
{: .tip }

---

# `Hello World!` with the ToolFactory.

A `Hello World!` Galaxy tool is a good place to start, just like any other new programming environment.
It requires planning and preparation. The ToolFactory can automate the generation of a wrapper, but the developer must supply a working script and
configure the inputs, outputs, user supplied parameters and metadata for the tool to be useful.

## Planning the new tool

A very simple bash script can be used to say "hello" but we make it a little more like a real
Galaxy tool by adding a text box so the user may designate whatever they want to add after that
such as "Hello, Galaxy Training Network".

Save the following sample as `hello.sh`:

 > <code-in-title>Starting bash script: Hello World</code-in-title>
 > ```bash
 > #!/bin/bash
 > echo "Hello $1!"
 > ```
 {: .code-in}

Test it on the command line by running:

`bash hello.sh ToolFactory`

In this case, `ToolFactory` is the first command line parameter. `Hello ToolFactory!` should appear as the output.

Once the script works and produces the expected outputs, the next step is to plan how the generated tool form should look to the user when run as a Galaxy tool.
In this case, a single text string is needed from the user.

Tool definition involves configuring the major sections of the ToolFactory form for the new tool.

Five categories needed to generate a script:

- Conda dependency requirements
- History data inputs
- History outputs
- User controlled parameters
- Developer supplied code/scripts to embed

For the `hello` tool case:

Conda dependency requirements:

- There are no dependencies usually because bash is available and version is not important.
- For completeness, it could be included as a Conda package. It's your tool.

History data inputs:

- This tool requires no history input files.

History outputs - data and collections:
- It produces one text output file.

User controlled parameters:
- The tool form should show a single input text field for the user to supply.

Developer supplied code/scripts to embed:
- These are optional - many Conda packages will not need them.
- Most simple use-cases will involve developer supplied, known working code.
- Executing the tool is expected to write the decorated string to a new history item.
- This must be known to work with suitable inputs on the command line. 
  - Broken code == broken tool.
- Such as `echo "Hello $1!"` as a Bash script
- There are other advanced features, such as command over-rides where code can be embedded

At this point, the plan for this new tool is:

- The Galaxy user should see a helpfully labelled text input field on the tool form, and the usual tool `execute` button.
- When the tool executes, that text should be passed to the script running under bash, as the first positional parameter.
- The script output should appear as a new output file in the history.
- It should contain the expected decorated input text.
- Galaxy tools need a test.
   - A simple test would be to supply a default value for the text string, run the tool and check that the output is correct.

## Putting the plan into action using the ToolFactory

The form collects all the information needed for a new Galaxy tool. It is long and complex as a result, particularly with many repeated form elements for more complex tools.
Much of what is collected is used to construct a command line for the script when the generated tool runs.
Other information such as the name and dependencies are needed to construct the relevant
sections of the generated XML file in the toolshed archive. The ToolFactory form configured to generate the `Hello` example can be viewed below.

> <details-title>Detail to explore: Annotated ToolFactory form for the `Hello World` example</details-title>
> ![First part of the form](../../images/ToolFactory_hello1form.png "The first part of the form collects the new tool name and dependencies to be installed. In this case, no Conda dependency is used. bash can be specified as a conda dependency, but it is not very version dependent and usually available. Reproducibility is not an issue for this trivial example. When it is, specify the dependencies and their versions here and the generated tool will always use them. If run in a shell, the bash script <code>echo "Hello $1"</code> in the text box will emit a string that includes the first command line parameter - such as "Hello Galaxy Training Network" This will be collected from STDOUT (configured below) into a new history output file (named and configured below). Positional parameters are chosen so the first parameter on the command line will be emitted when the script runs.")
>
> ![Second part of the form](../../images/ToolFactory_hello2form.png "The second section shows the new generated history output. It uses the special name <code>STDOUT</code> so the tool will take whatever the bash script writes and create a new text file called <code>Hello_output</code>. When the test is generated, the pass criterion is that the default value <code>Galaxy Training Network</code> should appear as the message in <code>hello_output</code>. no difference. Other criteria including <code>sim_size</code> are available for the test applied to each output file. There is no limit (other than your patience) to the number of new generated history outputs. Note that this example has no history input files. Again, any number of these can be specified on the form using the repeat.")
>
> ![Third part of the form](../../images/ToolFactory_hello3form.png "The third section shows the form settings for the user supplied parameter to be passed in to the bash script on the command line. It will be the first positional parameter because the ordinal position is 1. Argparse parameters are shown in other samples. The help and label text for each input file and user defined parameter will appear on the generated tool form for the user so make them informative. This is where you can change the default from "World" to "Galaxy Training Network" on the sample provided and regenerate it to make a new tool later in the tutorial.")
>
{: .details }


The generated tool XML appears in the history after the ToolFactory is executed and the tool itself is installed in the `Local Tools` submenu.
Text on the form is specified in the XML and it all comes from the ToolFactory form.

> <details-title>Detail to explore: Generated XML and tool form</details-title>
>
> [Galaxy XML documentation is here](https://docs.galaxyproject.org/en/latest/dev/schema.html)
>
> Note how text from the form appears in the generated tool XML
>
> ```xml
> <tool name="hello_toolshed" id="hello_toolshed" version="0.01">
>   <!--Source in git at: https://github.com/fubar2/galaxy_tf_overlay-->
>   <!--Created by planemo@galaxyproject.org at 22/01/2021 13:48:27 using the Galaxy Tool Factory.-->
>   <description>Says hello</description>
>   <stdio>
>     <exit_code range="1:" level="fatal"/>
>   </stdio>
>  <version_command><![CDATA[echo "0.01"]]></version_command>
>  <command><![CDATA[bash
>  $runme
>  "$sayhelloto" > $Hello_output]]>
>  </command>
>   <configfiles>
>     <configfile name="runme"><![CDATA[
>  echo "Hello $1"
>  ]]></configfile>
>   </configfiles>
>   <inputs>
>     <param label="Say hello to" help="" value="Galaxy Training Network!!" type="text" name="sayhelloto" argument="sayhelloto"/>
>   </inputs>
>   <outputs>
>     <data name="Hello_output" format="txt" label="Hello_output" hidden="false"/>
>   </outputs>
>   <tests>
>     <test>
>       <output name="Hello_output" value="Hello_output_sample" compare="diff" lines_diff="0"/>
>       <param name="sayhelloto" value="Galaxy Training Network!!"/>
>     </test>
>   </tests>
>   <help><![CDATA[
>
> **What it Does**
>
> ToolFactory demonstration - hello world in Galaxy
>
>
>
> ------
>
>
> Script::
>
>     echo "Hello $1"
>
> ]]></help>
>   <citations>
>     <citation type="doi">10.1093/bioinformatics/bts573</citation>
>   </citations>
> </tool>
> ```
>
> Which, when seen loaded into Galaxy looks like an ordinary tool:
> ![Generated form seen by the new tool user](../../images/toolfactory_hello_demo_form.png "The form displayed when the generated Hello tool is executed is below. The user sees a text box to enter any string. When executed, it will be echoed to a new history file called <code>Hello_output</code>")
> When a user runs the tool and enters some text in the text box, the decorated output will appear in a new history text dataset.
> This may not seem very exciting but it provides a useful pattern that can easily be adapted and extended.
> The script could do something far more interesting and could take unlimited input datasets, user configurable parameters and can produce as many outputs as needed in the history.
>
{: .details}


> <tip-title>If this is confusing</tip-title>
>
> If you are not yet familiar with the basics of Galaxy tools covered in the
> [tool integration training material]({% link topics/dev/tutorials/tool-integration/slides.html %}), the example form and XML
> will be confusing. You may gain more by reviewing the introductory part of that material, and then coming back here?
> It's a lot to learn and it is complicated. While a form driven code generator can hide much of the complexity of generating the code,
> the user must supply valid inputs for the code to be useful.
>
{: .tip}

## Extending this trivial example

This demonstrates a script based model, that can be extended to do more useful things with more complex scripts. More useful tools will ingest one or more user supplied input files, emit more complex outputs such as collections, and allow more user controlled parameters. 
These ToolFactory features are illustrated in the examples, and discussed in the advanced Tutorial.

Bash was used for this demonstration. Literally any scripting language, or any other useful package available in Conda, can be made available for tool execution by adding the appropriate dependency to the ToolFactory form.

> <comment-title>ToolFactory limitations and scope</comment-title>
>  
> The ToolFactory includes a very simple, automated, form driven XML code generator. Code automation is limited to the most common and easily implemented features. For example, conditionals are not available, so many complicated tools cannot be automatically generated with this tool. Send code.
> - It works best wrapping simple R/Bash/Python and other interpreted scripts, with few user supplied parameters and a few input and output files.
> - Scripts are easier than some Conda packages
>   - Scripts can often be modified to suit any ToolFactory limitations.
>   - Advanced tool components such as conditional logic and related tricks, requiring manual coding, can sometimes be worked around.
>         - where those features are needed, a skilled developer will be required.
>   - Many Conda dependencies will require XML conditionals or other tool XML constructs that are not easy to generate automatically.
>   - While some simple requirements may be manageable, complex ones will not be suitable for the ToolFactory.
> - Compared to the more usual shell and a text editor, The ToolFactory in Galaxy is a slow and clumsy way to debug your scripts.
> - **Starting a new ToolFactory tool with a know good command line and data** is strongly recommended.
>      - You will know exactly what to expect from the tool test for a first sanity check.
> - Corrolary: Unless there is a working script that needs to be wrapped into a toolshed-ready Galaxy tool, the ToolFactory is of little use.
{: .comment}

---

# Installation Options


> <warning-title>Security advisory!</warning-title>
>- *Please do not install the ToolFactory on any public server*
>- Configured as a default server, it has none of the usual additional security layers required for a safe public Galaxy server.
>- Although it will only run for administrative users, it allows unlimited scripting and exposes unwanted risk.
>- Install it locally and do not expose to the public internet.
>- For this reason, the training materials can't make use of existing public Galaxy infrastructure like most of the GTN material.
{: .warning}


> <hands-on-title>Installing and managing a Galaxy ToolFactory development server</hands-on-title>
> 
> # Logging in as an administrator to a new ToolFactory server
> 
> Once you have a working installation running, as described below, the server should be ready after 20-30 seconds, at [http://localhost:8080](http://localhost:8080).
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
> This takes twenty minutes or more to complete since the client must be built once. Visualisation plugins are not built to save some time.
>
> The resulting development server directory will occupy ~9GB of disk space, so be sure your machine has plenty of room.
> It will be based in a single directory, *galaxytf* in the same level as the cloned galaxy_tf_overlay repository directory, where the script should be run as shown above.
>
> Rerunning the *localtf.sh* script will *destroy the entire galaxytf directory* - all ~9GB, and create a clean new installation.
> It should only need to be run once in the life of the development server.
>
> Remove the *galaxytf* directory to remove the entire development server when it is no longer needed. 
> Before doing that, save any useful histories, because the jobs in a ToolFactory history can be used to update a tool, because that history can be imported into a fresh development instance when needed and the job rerun with the necessary adjustments to the form.
>
>
> Once installation is complete:
>  * start the server from the *galaxytf* directory with *sh run.sh*. The logs will be displayed.
>  * ^c (control+c) will stop it from the console.
>  * In 23.0, *.venv/bin/galaxyctl start* and *.venv/bin/galaxyctl stop* should work.
>
>
{: .hands_on}


## Exploring the ToolFactory

- The best way to understand what can be done, is to look at the sample tools in the default administrator initial history.
- As you explore the forms for each sample tool, you can see how the various options have been configured and what kinds of scripts or Conda packages this could be used for in your work.
- The example script can be swapped out for another one known to work and additional new parameters added to suit, to extend the toy examples and create tools of use to your users.
- Change the tool name when you do this on the newly edited form, then press `execute`
  - The new wrapper XML will appear
  - The new tool will be installed in the `Local Tools` submenu.
- If the tool name is not changed before re-generating a tool, the original installed tool will be updated with the new configuration. The old job can still be rerun from the history if necessary. Galaxy can be a clumsy but not entirely useless integrated development environment.

---

# Hello World!

> <hands-on-title>Building the Hello World example</hands-on-title>
>
> 1. Run {% tool ToolFactory %} with the following parameters:
>    - "Dependencies, optional script and script interpreter"
>      - *"Interpreter for the script"*: `bash`
>      - *"Script for executable above to interpret"*: `echo "Hello $1"`
>    - "Data file input, output and settings forming the executable or script command line"
>      - *"Command line parameter passing method to use"*: `positional`
>      - "Input and output files"
>        - {% icon param-repeat %} *"Insert Outputs"*
>          - *"Name for this output to appear in new history"*: `Hello_output`
>          - *"Select the datatype for this output"*: `txt`
>          - *"Position"*: `STDOUT`
>      - "Arguments"
>        - {% icon param-repeat %} *"Insert Command Line Parameters"*
>          - *"Choose the name for this parameter - MUST not be blank!"*: `say_hello_to`
>          - *"Enter this parameter's label for the form"*: `Say hello to`
>          - *"Positional ordinal \| argparse argument name"*: `1`
>
> 2. Execute
>
> 3. Explore the outputs - do they match what you expected?
>
> 4. Refresh the page - click the home icon (or the "Analysis" tab) - to see the new tool in the `Local Tools` section of the tools menu.
>
> 5. Run the tool that has been added - Select the new tool and examine the form. Check that all the changes are as they should be.
>
{: .hands_on}

## The Development Cycle

1. Test on the command line and confirm it produces correct output with defaults and test data.
1. In the development server, start a new history
1. Upload all input samples used on the command line if any, for use in the tool test.
1. Open the ToolFactory tool form.
    1. Define the tool metadata, dependencies, interpreter and optionally, paste the script.
    1. Add the required history inputs using the small samples as examples.
    1. Specify all the output files to be created in the user's history.
    1. Add any user adjustable command line parameters such as text fields.
    1. Look at the samples to see how the ToolFactory form can be used.
1. Execute the tool when the form is completed.
1. When the job is complete, refresh the page (Home icon or Analysis tab). The new tool will be found in the `Local Tools` section, ready to run.
1. Run the new tool and check that it does what you expect, or re-generate after adjusting the form settings as needed.
1. If it needs any changes, open the XML history item created when the tool was generated and use the {% icon galaxy-refresh %} rerun button to
recreate the ToolFactory form as it was when you last ran it. Adjust as needed and use the tool form`Execute` button to run the ToolFactory again with updated settings.
    1. Rinse, repeat.
1. Warning: generating a tool with an existing tool id such as `mytool` will overwrite the installed version of any previously generated tool with id "mytool".
    1. Persisted jobs in user histories allow previous versions to be recreated to restore older versions.

Galaxy can be used as an Integrated Development Environment for tools - clunky but oddly satisfying.
Note this is distinct from debugging the script - that is not at all satisfying in Galaxy unless you like waiting for jobs to finish.

A shell is much better for that.

![Galaxy as an IDE for tools with the ToolFactory](../../images/ToolFactory_big_picture.png "Galaxy can be used as a tool development environment for users who can write their own scripts as shown in this process overview slide.")

## Hello World: Continued

> <hands-on-title>Modifying the Hello World example</hands-on-title>
>
> 1. Rerun the output of your previous job, and make the following changes
>
>    - "Dependencies, optional script and script interpreter"
>      - *"Script for executable above to interpret"*: `echo "Hello $1"; echo "Goodbye $2";`
>    - "Data file input, output and settings forming the executable or script command line"
>      - "Arguments"
>        - Add a second {% icon param-repeat %} *"Insert Command Line Parameters"*
>          - *"Choose the name for this parameter - MUST not be blank!"*: `say_bye_to`
>          - *"Enter this parameter's label for the form"*: `Say bye to`
>          - *"Positional ordinal \| argparse argument name"*: `2`
>
{: .hands_on}

## Hello Collections!

> <hands-on-title>Building a File Splitter</hands-on-title>
>
> 1. Run {% tool ToolFactory %} with the following parameters:
>    - *"New tool ID and title for outputs"*: `file_splitter`
>    - "Dependencies, optional script and script interpreter"
>      - *"Interpreter for the script"*: `bash`
>      - *"Script for executable above to interpret"*:
>
>        ```bash
>        mkdir -p outputs/;
>        split --lines=$2 --additional-suffix=.txt $1 outputs/
>        ```
>
>    - "Data file input, output and settings forming the executable or script command line"
>      - *"Command line parameter passing method to use"*: `positional`
>      - "Input and output files"
>        - {% icon param-repeat %} *"Insert Inputs"*
>          - *"Select an input file from your history"*: Choose any XML file from your history, we'll use this as an example
>          - *"Select the datatype for this output"*: `txt`
>          - *"This will become the user prompt for the form so please make it informative"*: `File to split`
>          - *"Positional: ordinal integer. Argparse: argument name. STDIN if the executable/script expects it"*: `1`
>        - {% icon param-repeat %} *"Insert Output Collections"*
>          - *"Select the kind of collection for this output"*: `List`
>          - *"Label for this collection"*: `File Parts`
>      - "Arguments"
>        - {% icon param-repeat %} *"Insert Command Line Parameters"*
>          - *"Choose the name for this parameter - MUST not be blank!"*: `lines`
>          - *"Select the type for this parameter"*: `Integer`
>            - *"Enter this parameter's default integer value"*: `4`
>          - *"Enter this parameter's label for the form"*: `Number of lines in each split file`
>          - *"Positional ordinal \| argparse argument name"*: `2`
>
> 2. Execute
>
> 3. Refresh the Galaxy Page and locate your `collections_test` tool
>
> 4. Run the tool, selecting any text file in your history, e.g. the XML output from the ToolFactory that created this tool.
>
> 5. You should see a collection filled with files named `aa` to `ao` (or so), each with 4 lines from your file.
{: .hands_on}

This is presented as a motivating example. It is up to you to imagine what else you might be able to accomplish with collection outputs!
You could extend this example to split a `.fastq` file which could speed up your processing, or you could use collections to store extra images or plots produced by your tool.
They are also a convenient way to keep outputs together that are unlikely to be used downstream in a workflow, such as reports and images
that a user might want to be able to easily view if they want without cluttering up the history.

> <warning-title>Collection Testing</warning-title>
> The default generated test for output collections always passes because it doesn't test anything.
> Supplying a test over-ride is recommended for collections.
> For a real test, one or more expected <element.../> tags must be provided so the test really does test something.
{: .warning}

## Done!

> <hands-on-title>To safely shut the server down</hands-on-title>
>
> 1. Type <kbd>Ctrl-C</kbd> in the terminal where you ran `sh run.sh`.
>
{: .hands_on}

### Limits, workarounds and value proposition

- The ToolFactory Appliance is a slightly clumsy but useable way to create, test and maintain Galaxy tools in a web GUI.
- The ToolFactory tool is an automated code generator installed in the appliance tool menu
- No generator can replace manual editing by a skilled developer other than in constrained simple cases.
- These are common enough in the daily work of most data intensive scientific fields to make a tool generator potentially worth keeping handy.
- For simple scripts and appropriate Conda packages, a professional Galaxy tool developer can probably do it quickly by hand, but those skills take time to acquire
and are not widely available, particularly in scientific fields coming to Galaxy.
- Tools can have command-override and test-override pasted in as in one of the BWA samples.
   - This can solve some of the limitations but if it is needed, it might be better to prepare the wrapper manually if a skilled developer is available.
   - Any logic in the `<command>` section can probably always be replaced by equivalent code in a script at the cost of time and effort compared to templating.
   - Other aspects of tool logic such as output filters based on other parameter values can only be implemented in the wrapper document, not in a tool script.
- The ToolFactory can help new scientists and developers to quickly get some simple tools working for their colleagues while awaiting help with the complex ones.

# Next Steps

Expand your knowledge further with the [Advanced ToolFactory tutorial]({% link topics/dev/tutorials/tool-generators-advanced/tutorial.md %})

# Acknowledgements

This tutorial is based on the work of thousands of individual contributers to the Galaxy project over the last 15 years or so.
Thanks all! It has been a lot of fun.

Special thanks to:

- {% include _includes/contributor-badge.html id="hexylena" %} for
    - review and contribution to the tutorial and associated code.
    - the vision of instant installation of generated tools for developer feedback.
    - elegantly generated lint-free XML provided by [galaxyml code](https://github.com/hexylena/galaxyxml)
- {% include _includes/contributor-badge.html id="mvdbeek" %} for thoughtful comments on the role of the ToolFactory that helped motivate the tutorial.
