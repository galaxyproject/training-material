---

layout: tutorial_hands_on
logo: "GTN"
objectives:
- Learn why you might want to try the ToolFactory
- See a video demonstration - Hello Galaxy Training Network!
- Learn how to run it locally
- Explore the simple examples provided
- Modify and re-generate them to see how the changes affect the tool
- Try the ToolFactory on your own script!

questions:
 - How does a bioinformatician get from a functioning Integrated Environment to a `real` workflow compatible, shareable Galaxy tool?
 - Who might want to use the ToolFactory for quick tools from scripts?
 - How can I get the ToolFactory working since it should not be on a public server?

title: "Using a Galaxy tool generator for tools from simple scripts"
type: tutorial_hands_on
key_points:
  - "The ToolFactory can turn working command line scripts into proper Galaxy tools"
  - "A Galaxy tool to generate complete Galaxy tools by filling in a Galaxy form and clicking `Execute`"
  - "Designed for bioinformaticians developing scripts in Galaxy using IEs."
  - "Dedicated Galaxy tool developers use more powerful tools without limits"
  - "The ToolFactory is a code generator so covers only simple situations. These are common in many small scale, bespoke Galaxy analyses."
  - "Please do not upload trivial tools to the main toolshed!"

requirements:
 -
    type: "internal"
    topic_name: introduction
    tutorials:
      - galaxy-intro-short
      - galaxy-intro-101-everyone

follow_up_training:
 -
    type: "internal"
    topic_name: dev
    tutorials:
      - tool-integration
      - interactive-environments
contributors:
  - fubar2

---

> ### {% icon tip %} Before starting, check that this training will be useful for *your* work in Galaxy?
>
> * Read the *Brief Introduction* section below.
> * Most non-programmer scientists using Galaxy do not need this training.
> * The ToolFactory is designed for bioinformaticians and researchers who develop scripts using interactive environments in Galaxy
> * Galaxy tool developers who already have the tools they need do not need this training
> * It is particularly useful for developers coming to Galaxy from other scientific fields
> * Watch the `Hello Galaxy Training Network!` tool generation demonstration video to see if it looks interesting for your work
{: .tip }


# A Galaxy tool that generates simple Galaxy tools

#### Background

- Galaxy Integrated Environments are very popular and useful for skilled researchers and developers because they allow interactive
development in scripting languages such as Python or R in Galaxy.
- Scripts developed in IE's can be turned into command line scripts and tested using suitable small input data sets.
- Once the script is working on the command line, **the ToolFactory provides a quick route to a real Galaxy tool**.
- Any scripting language executable supported by Conda can be used such as R or Python, or system utilities like Bash.
- It is an automated, form driven code generator that runs in Galaxy.
- It is ideal for turning **simple** scripts from Interactive Environments into tools for workflows.
- Developed for bioinformaticians needing to produce "real" Galaxy tools for their users from scripts they have developed in Galaxy IEs.
- It is much easier to learn to use, but consequently has limits in terms of the complexity it can deal with compared to the recommended Galaxy developer
tool development software.
- Training is available for those in the "Development in Galaxy" section of the GTN - for example, "Tool development and integration into Galaxy" linked in the follow-up
training recommended at the end of this document

---

#### Overview

- Turns scripts that run correctly on the command line into real Galaxy tools
- It is not designed for scientist users, unless they also write working scripts.
- It makes new Galaxy tools from scripts. Developed in Galaxy by bioinformaticians for example.
- A bioinformatician who is comfortable with scripting languages on a linux command line might want to run it.
- Produces new Galaxy tools that wrap the supplied script.
- They pass Planemo lint, and are no different from manually written tools.
- They contain a test based on the test data provided at tool generation.

---

# The ToolFactory

#### A form driven Galaxy tool generator for bioinformaticians

- The ToolFactory is a Galaxy tool and can be found in the ToolShed.
- It runs in Galaxy like any other tool
- It automates much of the work needed to prepare a new Galaxy tool using information provided by the script writer, on the ToolFactory form.

> ### {% icon comment %} Note!
> - *The ToolFactory does not do any of the hard work needed to prepare a script to run correctly on the command line.*
> - *Galaxy is far too clumsy as an IDE to be used for that purpose.*
{: .comment}

- It can wrap any simple script that runs correctly on the command line with some small test input samples.
- This is exactly what the ToolFactory does best.
- Normal scientist Galaxy users are unlikely to need it unless they are also capable of confidently scripting for themselves.

---

#### *Simple* scripts

- Ideal for simple R/Bash/Python/.... scripts with a few user supplied parameters and a few i/o history files. The script can easily be modified to respond to default
empty parameters as if they had not been passed so conditionals and related tricks are not needed.
- For many Conda dependencies, wrappers need conditionals and other tool XML constructs that are not easy to generate automatically so while some may
be manageable, complex ones will often not be possible.
- Galaxy developer tools for building tools include Planemo and the new Galaxy Language Server in VCS. These are far more flexible than the ToolFactory
- They are recommended for full time tool developers in Galaxy willing to learn to use them and needing the flexibility and power.
- *The ToolFactory is for developers and bioinformaticians not yet familiar with those far more flexible tools. Sometimes the scripts they want to wrap are simple enough for the ToolFactory.*


---

> ### {% icon warning %} Security advisory!
>- *Please do not install the ToolFactory on a public facing server*
>- Although it will only run for administrative users, it allows unlimited scripting and that is never a good idea on a public facing machine. Please install it locally as described below.
>- For this reason, the training materials can't make use of existing public Galaxy infrastructure like most of the GTN material. Fortunately, there are a number of local installation alternatives available, depending on how you prefer to work.
{: .warning}

---

# *Hello Galaxy Training Network* sample

ToolFactory demonstration video

> ### {% icon tip %} Annotated toolFactory form sections for the Hello demonstration
>>>![](../images/ToolFactory_hello1form.png)
>
> - **The first section of the completed form collects the new tool name and dependencies.**
> - In this case, no Conda dependency is used although bash could be specified it is usually available on the command line.
> - The script pasted into the text box simply echos the first command line parameter.
> - Positional parameters are chosen so the first parameter on the command line will be emitted when the script runs.
>
>
>>> ![](../images/ToolFactory_hello2form.png)
> - **The second section shows the new generated history output.**
> - It uses the special name `STDOUT` - the tool will take whatever the bash script writes and create a new text file called `hello_output`
> - When the test is generated, the pass criterion is that the default value `Galaxy Training Network` should appear as the message in `hello_output`
with no difference. Other criteria including `sim_size` are available for the test applied to each output file.
>
>
>>> ![](../images/ToolFactory_hello3form.png)
> - **The third section shows the user supplied parameter to be passed in to the bash script on the command line**
> - It will be the first positional parameter because the ordinal position is 1
>
>
>
>>> ![](../images/ToolFactory_hello4form.png)
> - **The fourth section controls ToolFactory actions and optional outputs**
> - If you supply appropriate API keys, the ToolFactory can upload the newly generated tool to a toolshed. Optionally it can be installed
back into the Galaxy server specified.
> - *This is potentially annoying and dangerous if you have API keys you can misuse - so please be mindful*
>
{: .tip }

>### {% icon tip %} Generated XML Hello Galaxy Training Network sample
>```xml
><tool name="hello_toolshed" id="hello_toolshed" version="0.01">
>  <!--Source in git at: https://github.com/fubar2/toolfactory-->
>  <!--Created by planemo@galaxyproject.org at 22/01/2021 13:48:27 using the Galaxy Tool Factory.-->
>  <description>Says hello</description>
>  <stdio>
>    <exit_code range="1:" level="fatal"/>
>  </stdio>
> <version_command><![CDATA[echo "0.01"]]></version_command>
> <command><![CDATA[bash
> $runme
> ">$sayhelloto"
> $Hello_output]]></command>
>  <configfiles>
>    <configfile name="runme"><![CDATA[
> echo "Hello $1"
> ]]></configfile>
>  </configfiles>
>  <inputs>
>    <param label="Say hello to" help="" value="Sailor!" type="text" name="sayhelloto" argument="sayhelloto"/>
>  </inputs>
>  <outputs>
>    <data name="Hello_output" format="txt" label="Hello_output" hidden="false"/>
>  </outputs>
>  <tests>
>    <test>
>      <output name="Hello_output" value="Hello_output_sample" compare="diff" lines_diff="0"/>
>      <param name="sayhelloto" value="Galaxy Training Network!!"/>
>    </test>
>  </tests>
>  <help><![CDATA[
>
>**What it Does**
>
>ToolFactory demonstration - hello world in Galaxy
>
>
>
>------
>
>
>Script::
>
>    echo "Hello $1"
>
>]]></help>
>  <citations>
>    <citation type="doi">10.1093/bioinformatics/bts573</citation>
>  </citations>
></tool>
>```
{: .tip }

---

# Run the ToolFactory locally

- Depending on your preferences, the best way to install your own ToolFactory should be chosen from one of the options below.
- The sections after this assume that you have successfully installed it, so please choose a method


#### Install into an existing local development Galaxy

See [the tutorial on installing tools from the toolshed](https://galaxyproject.org/admin/tools/add-tool-from-toolshed-tutorial)

> ### {% icon hands_on %} Hands-on: Steps to get the ToolFactory installed in a private development Galaxy
>
> 1. Log in to Galaxy as an administrative user
> 2. Select the`Admin` tab from the top bar in Galaxy;
> 3. Under the `Tool Management` option, select `Install and Uninstall - Search and install new tools and other Galaxy utilities from the Tool Shed. See the tutorial.`
> 2. Make sure the Main toolshed is selected so the entire category list is displayed. Choose the `Tool Generators` link.
> 3. Select `tool_factory_2 updated version of the tool factory`
> 3. Select `Install` for the first listed version - revision 119 at present.
> 4. Wait a few minutes - it takes some time for Conda to install all the dependencies
{: .hands_on}


- If you have a local disposable development Galaxy on your desktop or laptop, this is the easiest way to run a ToolFactory
- Please do not install on a public facing or production Galaxy server
- https://toolshed.g2.bx.psu.edu/view/fubar/tool_factory_2/f267ed4a3026
- Once installed, it appears on the tool menu for all users, **but only local administrative users can successfully execute it**
- It will fail with an explanation for non-administrative users.

---

#### Install a potentially disposable ToolFactory inside Planemo in a virtual environment

- Make a new (potentially throw away) directory for the Planemo installation - e.g. `mkdir tftute` and `cd tftute`
- [click here to download a bash script to build a local (potentially throw-away) copy of the ToolFactory] (https://zenodo...)
- Save the script in the new directory and run it.
- It will create a virtual environment, download a fork of planemo and a separate local copy of galaxy-dev.
Edit the location of the galaxy directory `$GALDIR` and remove the code to download galaxy-dev if you wish to save time and space.
- It will take some time - so watch the Hello World demonstration while you wait.

> ### {% icon tip %} Sample script to install a local disposable ToolFactory in a planemo virtual environment
> > ### {% icon code-in %} Input: topics/tool-builders/docker/maketf.sh
> > ```bash
> > GALDIR="galaxy-central"
> > PDIR="planemo"
> > git clone --recursive https://github.com/fubar2/planemo.git $PDIR
> > rm -rf $PDIR/docs
> > mkdir -p $GALDIR
> > curl -L -s https://github.com/galaxyproject/galaxy/archive/dev.tar.gz | tar xzf - --strip-components=1 -C $GALDIR
> > cp $PDIR/planemo_ext/welcome.html $GALDIR/static/welcome.html
> > cp $PDIR/planemo_ext/welcome.html $GALDIR/static/welcome.html.sample
> > mkdir -p $PDIR/mytools
> > cd $PDIR
> > python3 -m venv .venv
> > . .venv/bin/activate
> > python3 setup.py build
> > python3 setup.py install
> > planemo conda_init --conda_prefix ./con
> > planemo tool_factory --galaxy_root $GALDIR --port 9090 --host 0.0.0.0 --conda_prefix $PDIR/con
> > ```
> {: .code-in}
{: .tip}

---

#### Build a simple Docker container - Training Docker

- The Docker script provided with this topic builds a different Galaxy from most GTN Docker containers. It does not include this tutorial.
It runs planemo tool_factory for you and exposes it on port 9090 so you can do all the same things as you
can with a local venv described above - but a little slower and isolated in a container.

> ### {% icon tip %} Sample Dockerfile to build a simple version of the ToolFactory in Planemo
> > ### {% icon code-in %} Input: topics/tool-builders/docker/Dockerfile
> > ```docker
> ># Galaxy - Using Galaxy tools to generate new Galaxy tools
> >#
> ># To build the docker image, go to root of the training repo and
> >#    docker build -t tool-generators -f topics/tool-generators/docker/Dockerfile .
> ># Take a break. Takes a while!
> ># To run image, make an export directory where you want to run it regularly and then
> >#    docker run -p "9090:9090" -v export:/export/  -t tool-generators
> ># ToolFactory planemo will be available on localhost:9090
> ># Training material is not yet installed - not sure how to do that in Planemo ?
> >
> >FROM ubuntu:latest
> >
> >MAINTAINER Ross Lazarus
> >ENV TARGDIR "/galaxy-central"
> >ENV PDIR "/planemo"
> >RUN apt update -y -qq && apt install -y -qq python3-dev gcc python3-pip build-essential python3-venv python3-wheel nano curl wget git python3-setuptools gnupg curl mercurial \
> >&& python3 -m pip install --upgrade pip \
> >&& curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - \
> >&& apt upgrade -y \
> >&& mkdir -p $TARGDIR \
> >&& curl -L -s https://github.com/galaxyproject/galaxy/archive/dev.tar.gz | tar xzf - --strip-components=1 -C $TARGDIR \
> >&& git clone --recursive https://github.com/fubar2/planemo.git $PDIR \
> >&& cd $PDIR \
> >&& mkdir mytools \
> >&& rm -rf $PDIR/doc \
> >&& python3 setup.py build \
> >&& python3 setup.py install \
> >&& planemo conda_init --conda_prefix $PDIR/con \
> >&& hg clone https://fubar@toolshed.g2.bx.psu.edu/repos/fubar/tacrev  /planemo/tacrev \
> >&& planemo test --galaxy_root $TARGDIR /planemo/tacrev \
> >&& cp $TARGDIR/config/datatypes_conf.xml.sample $TARGDIR/config/datatypes_conf.xml \
> >&& sed -i 's/<\/registration>/<datatype extension="tgz" type="galaxy.datatypes.binary:Binary" subclass="true" mimetype="multipart\/x-gzip" display_in_upload="true"\/> <\/registration>/' $TARGDIR/config/datatypes_conf.xml \
> >&& sed -i 's/<datatype extension="html"/<datatype extension="html" display_in_upload="true"/' $TARGDIR/config/datatypes_conf.xml \
> >&& apt-get clean && apt-get purge \
> >&&  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
> >ADD topics/tool-generators/docker/welcome.html $TARGDIR/static/welcome.html.sample
> >ADD topics/tool-generators/docker/welcome.html $TARGDIR/static/welcome.html
> >
> >
> >ENV GALAXY_CONFIG_BRAND "ToolFactory in Planemo"
> >EXPOSE 9090
> >ENTRYPOINT ["/usr/local/bin/planemo" ,"tool_factory", "--galaxy_root" ,"/galaxy-central", "--port", "9090", "--host", "0.0.0.0", "--conda_prefix", "/planemo/con"]
```
> {: .code-in}
{: .tip}

---

#### Install the ToolFactory docker container with integrated toolshed

- There is a more complex but integrated solution using the [ToolFactory docker container](https://github.com/fubar2/toolfactory-galaxy-docker).
- It provides an integrated toolshed and allows tools to be installed and used in the Galaxy used to run the ToolFactory.
- Like installation in a local Galaxy server, the docker container can be persisted as shown in the documentation for docker-galaxy-stable upon which it is based.

---

# Import the demonstration history

> ### {% icon comment %} Note!
> - This is the **first step** recommended after any of the installation options above until you are comfortable using the ToolFactory
> - It will give access to some sample ToolFactory tools that can be used to learn how the ToolFactory works
> - This is already done in the [ToolFactory docker container](https://github.com/fubar2/toolfactory-galaxy-docker)
{: .comment}


On the welcome page of the ToolFactory in Planemo, there is a [zenodo link](https://zenodo.org/record/4542837/files/planemo_demohistory_jan23.tar.gz?download=1).
Copy it and paste it into the URL box on the screen for importing a remote history.

> ### {% icon hands_on %} Hands-on: Steps to get the URL box to paste into for history upload of the Zenodo link
>
> 1. Select the`User` tab from the top bar in Galaxy;
> 2. Select `Histories`
> 3. Select `Import`
{: .hands_on}


It will take a few minutes to import. Wait patiently and when it's complete select the link to view histories and switch to the new one.
You will see a large number of pairs of history items and 4 data files used for testing.
Each pair comprises a toolshed ready archive containing a generated tool and a test, and a collection including a Planemo test report, the tool XML and a job log.
The archive history object has a circular "redo" button. Click that button and the ToolFactory form that generated the sample tool will appear. You can see how the tool was
built using the ToolFactory's limited capacities. Most of them are trivial of course. They are meant to be models rather than real examples.

---

# Learning to use the ToolFactory


> ### {% icon hands_on %} Exploring the sample tool forms
>
> * Install the ToolFactory safely to suit your needs as described below
> * Import the sample history
> * Select any of the generated toolshed archive history items so you can see the circular "redo" button
> * Click "redo"
> * Examine the form settings used to generate the tool
> * Try changing things - add new parameters or inputs/outputs; press `execute`; check the new version of the tool
{: .hands_on}



Probably the best way is to take a look at each sample tool by rerunning the job. Think about how the options have been configured and try extending it or using the form
as the basis for a new tool - but remember to change the tool name before you press execute to rerun the job and generate a new toolshed archive and report.

The Hello tool is a model for any simple bash script requiring only a few parameters and is easily extended to many situations
where a tool is needed quickly for a workflow. Try adding another parameter

---

# Testing newly generated tools

#### Using the ToolFactory installed directly into a development instance

> ### {% icon hands_on %} Hands-on: Loading new generated tools in a normal development Galaxy
>
> - You need a local toolshed that the Galaxy server is configured to talk to in `config/tool_sheds_conf.xml`
> - From the Galaxy root, `sh run_tool_shed.sh` should start one on `localhost:9009`
> - Make a new repository and upload the new toolshed archive
> - In the `Admin` menu, search for the new tool in your local toolshed and then choose `install`
> - Refresh the tool menu when installation is complete
{: .hands_on}

- You need a local toolshed that the Galaxy server is configured to talk to in `config/tool_sheds_conf.xml`
- From the Galaxy root, `sh run_tool_shed.sh` should start one on `localhost:9009`
- Make a new repository and upload the new toolshed archive
- In the `Admin` menu, search for the new tool in your local toolshed and then choose `install`
- Refresh the tool menu when installation is complete


#### Using Planemo in a virtual environment or GTN docker container

> ### {% icon hands_on %} Hands-on: Loading new generated tools in a venv or docker ToolFactory in Planemo
>
> - Download the toolshed archive from the Galaxy history where you generated the tool.
> - Unpack the archive into the planemo/mytools directory. It should appear as a single directory containing a test-data subdirectory and the new tool xml.
> - Stop planemo - `^c` will do it somewhat messily
> - Restart planemo with an additional parameter `--extra_tools planemo/mytools/`
> - The new tool should be ready to test on the tool menu
{: .hands_on}


---

# Uploading generated archives to toolsheds


> ### {% icon warning %} Trivial tools do not belong in the public toolsheds!
>- *Please do not upload trivial tools to the main ToolShed!*
>- The ToolFactory provides an option to upload a newly built tool to a toolshed. This was designed for the persistent docker option but works in the
planemo tool_factory.
>- Please do not abuse it by adding trivial tools to confuse users looking for useful tools.
{: .warning}



- The [ToolFactory docker container](https://github.com/fubar2/toolfactory-galaxy-docker) includes a local toolshed
- This allows new tools to be automatically installed back into the Galaxy running the ToolFactory.
- Please run your own local toolshed for trivial tools that are so simple or specialised that they are not likely to ever be useful for other scientists
- Uploading them to the main ToolShed is unlikely to help anyone

---

# Conclusion
