---

layout: tutorial_hands_on
title: "ToolFactory: Generating Tools From More Complex Scripts"
key_points:
  - The ToolFactory is a fully automated Galaxy tool generator for scientists and developers who routinely write command line scripts.
  - It can turn a working command line script into a proper Galaxy tool with a test in a few minutes.
  - It automatically generates simple, complete Galaxy tools from information provided by filling in a normal Galaxy form in the familiar UI.

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

## The ToolFactory supports users who routinely write command line scripts in their work.

The ToolFactory can be found in the main ToolShed under the `tool-generators` category. It automates much of the work needed to prepare a
new Galaxy tool using information provided by the script writer,
on the ToolFactory form. The ToolFactory can wrap any simple script that runs correctly on the linux command line with some small test input samples. This is potentially
handy for developers new to Galaxy, and for Galaxy users who are capable of correctly scripting on the command line for themselves.


> ### {% icon tip %} Under the hood:
>
>  - It uses [galaxyml](https://github.com/hexylena/galaxyxml) to generate the tool XML from ToolFactory form settings.
>  - It uses [Planemo](https://github.com/galaxyproject/planemo) to generate the test outputs
>  - Then again to test newly generated code
{: .tip}



## Limits and scope

- It works best wrapping simple R/Bash/Python and other interpreted scripts with a few user supplied parameters and a few i/o history files.
- Scripts are easier than some Conda packages because they can easily be modified to respond to default empty parameters as if they had not been passed. As a result, advanced tool building elements
such as conditionals and related tricks requiring manual coding, can often be avoided.
- On the other hand, many Conda dependencies will require XML conditionals
or other tool XML constructs that are not easy to generate automatically. While some simple requirements may be manageable, complex ones will not be suitable for the ToolFactory.
- Compared to the more usual shell and a text editor, The ToolFactory in Galaxy is a slow and clumsy way to debugging scripts. More than a minute per cycle because`planemo test` is run twice, building and tearing down a Galaxy each time.
- **Starting a new ToolFactory tool with a know good command line and data** is strongly recommended. You will know exactly what to expect from the tool test for a first sanity check.
- Corrolary: Unless there is a working script that needs to be wrapped into a toolshed-ready Galaxy tool, the ToolFactory is of little use.


**The ToolFactory is for developers and informaticians not yet familiar with those far more flexible tools.**
**Scripts they need to wrap are frequently simple enough for the ToolFactory.**

Compared to other Galaxy tool development software, there is far less to learn in order to get up to speed when using a form driven, automated code generator. The
cost of this convenience is that ToolFactory is limited to automated generation of a large but limited subset of simple script and package wrappers.


# 2. Getting your hands on a ToolFactory for some hands-on training.

#### Run the ToolFactory locally and adapt the sample tools

- If you found the introductory material presented so far relevant to your own needs, you may wish to start the DIY/hands-on part of the tutorial that follows
- Set up your own working ToolFactory, install the samples in a history and then start exploring it and figuring out how it might help your work.
- Depending on your preferences, install your own ToolFactory from one of the options described below.
- The sections after this can **only be completed with a working ToolFactory**.
- Work done in a Planemo ToolFactory will not be `persistent`. For any serious use, this is a problem although saving histories or converted workflows can be used
to manually persist the ToolFactory configuration for each new tool. Some options involve `persistent` Galaxy servers and these are much more useful for building and
more importantly, updating tools. It can be disappointing to learn that the history recording all your ToolFactory work is no longer available the next time you start working.
- Non persistent options are **recommended only for testing or teaching - not production**.
- In all cases, the first time they are run and the first time a tool is built, most versions take 10 minutes or so  - there's a lot that needs to be installed.
Check for Conda and other running processes before assuming it has frozen.

>#### Active Tutorial content follows
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Installation

> ### {% icon warning %} Security advisory!
>- *Please do not install the ToolFactory on a public server*
>- Although it will only run for administrative users, it allows unlimited scripting and that is a high security risk opportunity for any public facing machine.
>- In fact, Galaxy is very good at isolating tools to stop them doing mischief. But that's no reason to chance your arm. They keep inventing better mice.
>- Please install it locally as described below.
>- For this reason, the training materials can't make use of existing public Galaxy infrastructure like most of the GTN material.
{: .warning}

# Running the ToolFactory

> ### {% icon hands_on %} Hands-on: Launching the Container
>
> 1. [Install Docker and Docker Compose](https://docs.docker.com/engine/install/) following the appropriate instructions for your platform
>
> 2. Go to [the ToolFactory appliance github repository](https://github.com/fubar2/toolfactory-galaxy-server)
>
> 3. Clone it or download the zip and unzip it somewhere handy - such as `~/toolfactory-galaxy-server-main`
>
> 4. Change to the compose directory - `cd ~/toolfactory-galaxy-server-main/compose`
>
> 5. `docker-compose up`
>
>
>    > ### {% icon code-in %} Input: Bash
>    > ```bash
>    > wget https://github.com/fubar2/toolfactory-galaxy-server/archive/refs/heads/main.zip
>    > unzip main.zip
>    > cd toolfactory-galaxy-server-main/compose
>    > docker-compose up
>    > ```
>    {: .code-in}
>
>    > ### {% icon tip %} Tip: Patience!
>    > When you run the ToolFactory for the first time inside the container and whenever you run a new tool with new dependencies, it will require some time to build the conda environment.
>    > Check for Conda or other processes if things seem stuck.
>    {: .tip}
>
>
> 6. Wait a minute or until process activity dies down
>
> 7. Browse to [port 8080 on your local machine - http://localhost:8080](http://localhost:8080)
>
{: .hands_on}

## Import ToolFactory functional documentation - the demonstration tools.

- Congratulations on getting this far and acquiring a local instance of the ToolFactory
- There is a history you should import that shows some sample tools
- You can examine how these were generated by using the Galaxy job redo button.
- This will show you the fully completed ToolFactory form used to generate the sample
- You can edit the form and regenerate a new tool with your changes incorporated.


> ### {% icon announcement %} Note!
> - This is the **first step** recommended after any of the installation options above until you are comfortable using the ToolFactory
> - It provides access to the sample ToolFactory tools.
> - They are the best way to learn how the ToolFactory works and how you might adapt the variations shown in your own work.
> - It provides functional documentation and is not needed once you are comfortable using the ToolFactory.
> - It is pre-installed in the [ToolFactory docker container](https://github.com/fubar2/toolfactory-galaxy-docker)
{: .announcement}


- Use this url `https://zenodo.org/record/4729971/files/TF_demo_history_April30.tar.gz?download=1`
-[zenodo link](https://zenodo.org/record/4729971/files/TF_demo_history_April30.tar.gz?download=1).
- Copy it and paste it into the URL box on the screen for importing a remote history.
- The link is also on the welcome page of the virtualenv Planemo installation described above.

> ### {% icon hands_on %} Hands-on: Steps to use that URL to import the history
>
> 1. Select the`User` tab from the top bar in Galaxy;
> 2. Select `Histories`
> 3. Select `Import`
> 4. Paste the URL into the URL field and press `import`
{: .hands_on}


- It will take a few minutes to import.
- Get up and have a stretch for a minute.
- When it's complete, select the link to view histories and choose `switch` from the drop down arrow on the new history to make the imported one your current history.
- Viewing the new history, you will see a large number of history items and 5 data files used for testing.
- Each item is a collection. Opening it will reveal a toolshed ready archive containing a generated tool and the generated tool XML.
- These items have a {% icon galaxy-refresh %} rerun button. Click that button and the ToolFactory form that generated the sample tool will appear. You can see how the tool was
built using the ToolFactory's limited capacities. Most of them are trivial of course. They are meant to be models rather than useful examples.

<sup id='section3'>*</sup>
# 3. Hands-on: Learning to use the ToolFactory

> ### {% icon comment %}First time use involves a long pause in some installations
> - The first job takes longer in some installation scenarios because the ToolFactory dependencies are installed before the tool can run.
{: .comment}


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

Consider the trivial `Hello World!` tool example. It is readily extended to suit many situations where a tool is needed quickly for a workflow. Try adding another parameter.
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
>     - Multiple dependencies. Conda is currently supported. System utilities can be used assuming the target server exposes them to tools, or they can be provided as Conda dependencies to ensure they will always be available
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
> ><tool name="reps" id="reps" version="0.01">
> >  <!--Source in git at: https://github.com/fubar2/toolfactory-->
> >  <!--Created by planemo@galaxyproject.org at 11/04/2021 10:18:15 using the Galaxy Tool Factory.-->
> >  <description></description>
> >  <requirements/>
> >  <stdio>
> >    <exit_code range="1:" level="fatal"/>
> >  </stdio>
> >  <version_command><![CDATA[echo "0.01"]]></version_command>
> >  <command><![CDATA[python
> >$runme
> >#for $rep in $R_mi:
> >--mi "$rep.mi"
> >#end for
> >#for $rep in $R_mp:
> >--mp "$rep.mp"
> >#end for
> >>
> >$repout]]></command>
> >  <configfiles>
> >    <configfile name="runme"><![CDATA[
> >import argparse
> >parser = argparse.ArgumentParser()
> >a = parser.add_argument
> >a("--mi", action="append")
> >a("--mp", action="append")
> >args = parser.parse_args()
> >print(" and ".join(args.mi))
> >print(" and".join(args.mp))
> >]]></configfile>
> >  </configfiles>
> >  <inputs>
> >    <repeat name="R_mi" title="Add as many Input lots of inputs as needed">
> >      <param name="mi" type="data" optional="false" label="Input lots of inputs" help="" format="ab1,affybatch,agilentbrukeryep.d.tar,agilentmasshunter.d.tar,analyze75,arff,asn1,
> > asn1-binary,augustus,axt,bam,bcf,bed,bed12,bed6,bedgraph,bedstrict,bigbed,bigwig,biom1,biom2,blastxml,blib,brukerbaf.d.tar,brukertdf.d.tar,bus,cel,chira.sqlite,chrint,cisml,ckpt,cmap,
> > cml,consensusxml,cool,cps,cpt,cram,csfasta,csv,ct,cuffdiff.sqlite,cxb,d3_hierarchy,daa,dada2_dada,dada2_errorrates,dada2_mergepairs,dada2_sequencetable,dada2_uniques,dbn,dbnsfp.tabular,
> > dcd,deeptools_compute_matrix_archive,deeptools_coverage_matrix,dlib,drf,dta,dta2d,dzi,edr,edta,eland,elandmulti,elib,embl,encodepeak,eset,excel.xls,expression.json,fai,fast5.tar,
> > fast5.tar.bz2,fast5.tar.gz,fasta,fasta.gz,fastg,fastq,fastq.bz2,fastq.gz,fastqcssanger,fastqcssanger.bz2,fastqcssanger.gz,fastqillumina,fastqillumina.bz2,fastqillumina.gz,fastqsanger,
> > fastqsanger.bz2,fastqsanger.gz,fastqsolexa,fastqsolexa.bz2,fastqsolexa.gz,featurexml,ffdata,ffindex,flv,fped,fphe,fps,fqtoc,gafa.sqlite,gal,gatk_dbsnp,gatk_interval,gatk_recal,gatk_report,
> > gatk_tranche,gemini.sqlite,genbank,genbank.gz,geojson,gfa1,gfa2,gff,gff3,gff3.bz2,gff3.gz,gii,gii.gz,gpr,grd,grd.tgz,gro,gtf,h5,h5ad,hdt,hhr,hlf,hmm2,hmm3,icm,idat,ideaspre,idpdb,idxml,
> > imgt.json,imzml,inchi,intermine_tabular,interval,ipynb,isa-json,isa-tab,itp,jellyfish,json,jsonld,kallisto.idx,kroenik,lav,ldindep,len,loom,lped,maf,malist,mascotxml,maskinfo-asn1,
> > maskinfo-asn1-binary,mcool,mdp,memepsp,memexml,meryldb,mgf,mkv,mol,mol2,mothur.accnos,mothur.align,mothur.align.check,mothur.align.report,mothur.axes,mothur.cons.taxonomy,
> > mothur.count_table,mothur.design,mothur.dist,mothur.filter,mothur.filtered.masked.quan,mothur.filtered.quan,mothur.freq,mothur.groups,mothur.list,mothur.lower.dist,mothur.map,
> > mothur.masked.quan,mothur.names,mothur.oligos,mothur.otu,mothur.otu.corr,mothur.otulabels,mothur.pair.dist,mothur.quan,mothur.rabund,mothur.rdp.taxonomy,
> > mothur.ref.taxonomy,mothur.relabund,mothur.sabund,mothur.seq.taxonomy,mothur.sff.flow,mothur.shared,mothur.square.dist,mothur.summary,mothur.tax.summary,mothur.tre,mp3,mp4,
> > mpg,mrm,ms2,msh,msp,mtx,mz.sqlite,mz5,mzdata,mzid,mzml,mzq,mztab,mztab2,mzxml,n3,ncbitaxonomy.sqlite,ndx,neostore.zip,netcdf,newick,nex,nhx,nii1,nii1.gz,nii2,nii2.gz,nmrml,
> > nt,obfs,obo,odgi,ome.tiff,osw,owl,oxlicg,oxligl,oxling,oxliss,oxlist,oxlits,paf,paf.gz,par,paramxml,parquet,pbed,pdb,pdbqt,pdf,peff,peplist,peptideshaker_archive,pepxml,pepxml.tsv,phylip,
> > phyloxml,picard_interval_list,pileup,plyascii,plybinary,png,postgresql,pphe,pqp,pqr,probam,probed,protobuf2,protobuf3,protxml,protxml.tsv,psms,pssm-asn1,ptkscmp,qcml,qname_sorted.bam,
> > qual454,qualillumina,qualsolexa,qualsolid,rdata,rdata.camera.negative,rdata.camera.positive,rdata.camera.quick,rdata.msnbase.raw,rdata.sce,rdata.xcms.fillpeaks,rdata.xcms.findchrompeaks,
> > rdata.xcms.group,rdata.xcms.raw,rdata.xcms.retcor,rdf,rdock_as,rma6,rna_eps,sam,sbml,scf,scidx,sdf,searchgui_archive,sf3,sff,shp,sif,smat,smi,snaphmm,snpeffdb,snpmatrix,snpsiftdbnsfp,
> > snptest,spalndba,spalndbnp,spec.xml,splib,splib_noindex,sqlite,sqmass,sra,sra_manifest.tabular,stl,stockholm,tabular,tabular.gz,tandem,tar,taxonomy,tck,textgrid,tgz,thermo.raw,tiff,top,tpr,
> > trackhub,trafoxml,traml,trk,trr,tsv,ttl,twobit,txt,uniprotxml,unsorted.bam,vcf,vcf_bgzip,vel,velvet,vg,vtkascii,vtkbinary,watersmasslynx.raw.tar,wav,wiff,wiff.tar,wig,xg,xgmml,xlsx,xmfa,xml,
> > xquest.xml,xtc,xvg,zip" multiple="false"/>
> >    </repeat>
> >    <repeat name="R_mp" title="Add as many things to add as needed">
> >      <param name="mp" type="text" value="Add lots of things" label="things to add" help=""/>
> >    </repeat>
> >  </inputs>
> >  <outputs>
> >    <data name="repout" format="txt" label="repout" hidden="false"/>
> >  </outputs>
> >  <tests>
> >    <test>
> >      <output name="repout" value="repout_sample" compare="sim_size" delta_frac="0.2"/>
> >      <repeat name="R_mi">
> >        <param name="mi" value="mi_sample"/>
> >      </repeat>
> >      <repeat name="R_mp">
> >        <param name="mp" value="Add lots of things"/>
> >      </repeat>
> >    </test>
> >  </tests>
> >  <help><![CDATA[
> >
> >**What it Does**
> >
> >
> >
> >------
> >
> >
> >Script::
> >
> >    import argparse
> >    parser = argparse.ArgumentParser()
> >    a = parser.add_argument
> >    a("--mi", action="append")
> >    a("--mp", action="append")
> >    args = parser.parse_args()
> >    print(" and ".join(args.mi))
> >    print(" and".join(args.mp))
> >
> >]]></help>
> >  <citations>
> >    <citation type="doi">10.1093/bioinformatics/bts573</citation>
> >  </citations>
> ></tool>
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

- The ToolFactory is an automated code generator.
- No generator can replace manual editing by a skilled developer other than in constrained, simple cases.
- These are common enough in the daily work of most data intensive scientific fields to make a tool generator potentially worth keeping handy.
- For simple scripts and appropriate Conda packages, it's potentially very useful.
- It is not hard to imagine using a Python wrapper to finesse more complex tools just as bash was used in the `planemo lint` example.
- The ToolFactory, if in a persistent form, is a slightly clumsy but useable way to create and maintain Galaxy tools.
- Tools can have command-override and test-override pasted in as in one of the BWA samples. This can solve some of the limitations. However, if the package requires that kind of complexity, it might be better to prepare the wrapper manually.


## Notes on some commonly reported issues

#### First job I submitted in Planemo ToolFactory or the Docker container remains grey or running for a long time - is it broken?

- Check with top or your system monitor - if Conda is running, things are working but it's slow the first time.
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
