---
title: How to Wrap Your Pipeline
description: Pre-existing pipelines can be wrapped in a variety of ways
area: tool-dev
box_type: tip
layout: faq
contributors: [hexylena]
---

{% include _includes/cyoa-choices.html title="Is this for wide distribution or to help a colleague internally?" disambiguation="distribute" option1="Internal only" option2="We want the world to use this in Galaxy" default="__hidden__" %}

## Pipeline Source

<div class="Internal-only">

{% include _includes/cyoa-choices.html disambiguation="package" title="Which workflow engine is used for this pipeline?" option1="Nextflow" option2="Snakemake" option3="Bash" default="__hidden__" %}

<div class="Nextflow" markdown="1">

You're probably a busy researcher! You may not want to spend a lot of time wrapping this. There is an *easy way out* for you. It's **not best practice**, but maybe you just need to get it done quickly

```xml
<tool name="Nextflow Pipeline: Sarek" id="internal.nextflow.sarek" version="1.0">
  <description>runs a specific nextflow pipeline</description>
    <requirements>
      <requirement type="package" version="24.04.4">nextflow</requirement>
    </requirements>
  <stdio>
    <exit_code range="1:" level="fatal"/>
  </stdio>
  <command><![CDATA[
  nextflow run nf-core/sarek \
   --input '$input' \
   --outdir outdir
  ]]></command>
  <inputs>
    <param label="Input File" type="data" format="csv" name="input"/>
  </inputs>
  <outputs>
    <discover_datasets pattern="(?P&lt;designation&gt;.+)\.(?P&lt;ext&gt;.*)" directory="outdir" />
  </outputs>
  <help><![CDATA[]]></help>
</tool>
```

</div>

<div class="Snakemake" markdown="1">

You're probably a busy researcher! You may not want to spend a lot of time wrapping this. There is an *easy way out* for you. It's **not best practice**, but maybe you just need to get it done quickly

```xml
<tool name="Snakemake Pipeline" id="internal.snakemake.pipeline" version="1.0" profile="21.05">
  <description>runs a specific pipeline</description>
    <requirements>
      <requirement type="package" version="8.20.5">snakemake</requirement>
    </requirements>
  <stdio>
    <exit_code range="1:" level="fatal"/>
  </stdio>
  <command><![CDATA[
  cp '$__tool_directory__/Snakefile' . ;
  snakemake
  ]]></command>
  <inputs>
    <param label="Input File" type="data" format="csv" name="input"/>
  </inputs>
  <outputs>
    <discover_datasets pattern="(?P&lt;designation&gt;.+)\.(?P&lt;ext&gt;.*)" directory="outdir" />
  </outputs>
  <help><![CDATA[]]></help>
</tool>
```

</div>

<div class="Bash" markdown="1">

You're probably a busy researcher! You may not want to spend a lot of time wrapping this. There is an *easy way out* for you. It's **not best practice**, but maybe you just need to get it done quickly

```xml
<tool name="Bash Pipeline" id="internal.bash.pipeline" version="1.0">
  <description>runs a specific pipeline</description>
    <requirements>
      <requirement type="package" version="1.0.0">my-conda-dependency</requirement>
    </requirements>
  <stdio>
    <exit_code range="1:" level="fatal"/>
  </stdio>
  <command><![CDATA[
  # Literally just paste your bash script here
  # And add some '$input' and 'outdir' variables
  ]]></command>
  <inputs>
    <param label="Input File" type="data" format="csv" name="input"/>
  </inputs>
  <outputs>
    <discover_datasets pattern="(?P&lt;designation&gt;.+)\.(?P&lt;ext&gt;.*)" directory="outdir" />
  </outputs>
  <help><![CDATA[]]></help>
</tool>
```

<div class="Bash Snakemake Nextflow" markdown="1">

Obviously this will not take optimal advantage of your cluster, you'll have to set the memory for this "tool" to the maximum value that the entire pipeline will need, but it'll get the job done.

Of course if you have time or want the world to use this pipeline in Galaxy, you can go back and follow the <button class="btn btn-primary" onclick="cyoaChoice('We-want-the-world-to-use-this-in-Galaxy', 'gtn-cyoadistribute');">Best Practice workflow</button>

</div>

</div>
</div>


<div class="We-want-the-world-to-use-this-in-Galaxy" markdown="1">

Ah you want the world to use this! Ok, that probably means you'll want to make it a best practice workflow, which means decomposing each and every individual step into the corresponding Galaxy tool.

This may mean updating existing tools to use new flags, or creating new tools entirely. 

## Tools

- tools
    - are ALL of those tools in Galaxy?
    - are SOME in Galaxy?
    - are NONE in Galaxy?


We'll start off with what you have

</div>
