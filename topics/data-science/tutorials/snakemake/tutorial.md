---
layout: tutorial_hands_on

title: Make & Snakemake
level: Introductory
zenodo_link: ""
requirements: []
follow_up_training: []

questions:
- What is Make and Snakemake
- What is a Makefile and a Snakefile
- How do these relate
- Why is Snakemake better for scientific research and how can I use it
objectives:
- Write a snakefile that does some simple computations
time_estimation: 1H
key_points:
- Make and Snakemake are ways to write pipelines in a declarative format
- Instead of writing individual actions you wish to take, you describe how to produce each file you need, and the system executes those steps only if they're needed
- These systems can significantly speed up your pipelines by providing automatic detection of inputs, and not re-creating them if they don't need
- These systems both also provide very easy access to paralellisation without complexity.
subtopic: sciwms
contributors:
  - hexylena
  - bazante1
abbreviations:
  SciWMS: Scientific Workflow Management System
---

Here you will learn to write both Make and Snakemake workflows. We teach two workflow engines because Snakemake uses a lot of the concepts of Make, and these concepts are somewhat complex and a very different way of thinking than you might be used to with workflow design.

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

We're going to go through this tutorial as "the evolution of a pipeline" from a simple bash script that we might run as individual commands on the command line, all the way up to a Snakemake workflow. This should give you some perspective of why the systems exist, and what benefits each system brings.

# Bash

We've set up a simple bash pipeline. It downloads some read files from a website, decompresses them, builds an index for a `genome.fa` we already have, and then does the alignment against that genome.

```bash
# Downloading our datasets
wget https://zenodo.org/api/files/TODO/GCA_000017985.1_ASM1798v1_genomic.fna.gz
wget https://zenodo.org/api/files/TODO/SRR2584866_1.fq.gz
wget https://zenodo.org/api/files/TODO/SRR2584866_2.fq.gz
wget https://zenodo.org/api/files/TODO/SRR2589044_1.fq.gz
wget https://zenodo.org/api/files/TODO/SRR2589044_2.fq.gz

# Generate FastQC Report
fastqc *.fq

# Run Trimmomatic to trim the bad reads out.
trimmomatic PE SRR2589044_1.fq SRR2589044_2.fq \
               SRR2589044_1.trim.fq SRR2589044_1un.trim.fq \
               SRR2589044_2.trim.fq SRR2589044_2un.trim.fq \
               SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

trimmomatic PE SRR2584866_1.fq SRR2589044_2.fq \
               SRR2584866_1.trim.fq SRR2589044_1un.trim.fq \
               SRR2584866_2.trim.fq SRR2589044_2un.trim.fq \
               SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

# Generate Updated FastQC Report
fastqc *.trim.fq

# Generate the genome index
bwa index GCA_000017985.1_ASM1798v1_genomic.fna

# Align reads to the genome
bwa mem GCA_000017985.1_ASM1798v1_genomic.fna SRR2589044_1.trim.fq SRR2589044_2.trim.fq | samtools sort -O bam -o SRR2589044.bam
bwa mem GCA_000017985.1_ASM1798v1_genomic.fna SRR2584866_1.trim.fq SRR2584866_2.trim.fq | samtools sort -O bam -o SRR2584866.bam
```

This is a fine start to our analysis pipeline, and might simply be summarising
the commands we executed interactively at the command line.

# Make

Make is a build automation tool that has been around since 1976. It lets you describe build processes, writing down the individual processes which will occur, and then providing you a single entrypoint to run your workflow, a lot like how {SciWMS} work! But importantly it is declarative, rather than imperative which is a big change if you're familiar with programming languages like bash or python.

Instead of defining step by step all of the steps that should be executed like the bash example above,
you instead write more general rules on how to produce individual files. We'll look at them one-by-one.

Whereas bash scripts you write something like:
- `command input` and it implicitly creates the output file automatically, often by adding a suffix
- `command input > output` redirecting output of the command to a file
- `command input -O output` where we explicitly state where the output is

In Make you will write:

<pre class="highlight"><code>
<span class="nb">output</span>: <span class="kt">input</span>
	<span class="s2">command [and here the correct variant from above]</span>
</code></pre>

The important difference being that Make is aware of the inputs and outputs of each step, so it knows what input files need to exist before it can execute this step, and it knows precisely which output files you will generate. It needs this information so that when you have multiple rules, it can decide which dependencies need to be executed.

## Line-By-Line Comparison

<table class="table">
<thead><tr><th>Bash</th><th>Make</th></tr></thead>
<tbody>
	<tr>
		<td><pre class="highlight"><code>
<span class="s2">wget https://.../GCA_000017985.1_ASM1798v1_genomic.fna.gz
wget https://.../SRR2584866_1.fq.gz</span>
</code></pre></td>
		<td><pre class="highlight"><code>
<span class="nb">%.gz</span>:
	<span class="s2">wget https://.../$@ -O $@</span>
</code></pre></td>
	</tr>
	<tr>
		<td>Make the directory and download each read file one by one</td>
		<td>Here is a generic rule which can be used any time you need to download a file for this project. <code>$@</code> is used as the name of the output file and in the templated source url we will download from. If you ran <code>make SRR2584866_1.fq.gz</code> it would template out the <code>wget</code> command and run it.</td>
	</tr>
	<tr>
		<td><pre class="highlight"><code>
<span class="s2">fastqc *.fq</span>
</code></pre></td>
		<td><pre class="highlight"><code>
<span class="nb">%.fastqc.html</span>: <span class="kt">%.fq</span>
	<span class="s2">fastqc $<</span>
</code></pre></td>
	</tr>
	<tr>
		<td>Generate ALL of the FastQC reports</td>
		<td>Here is a rule to generate a single FastQC report from a single FastQ file</td>
	</tr>

	<tr>
		<td><pre class="highlight"><code>
<span class="s2">trimmomatic PE</span> <span class="kt">SRR2589044_1.fq SRR2589044_2.fq</span> \
               <span class="nb">SRR2589044_1.trim.fq SRR2589044_1un.trim.fq \
               SRR2589044_2.trim.fq SRR2589044_2un.trim.fq</span> \
               <span class="s2">SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15</span>

<span class="s2">trimmomatic PE</span> <span class="kt">SRR2584866_1.fq SRR2584866_2.fq</span> \
               <span class="nb">SRR2584866_1.trim.fq SRR2584866_1un.trim.fq \
               SRR2584866_2.trim.fq SRR2584866_2un.trim.fq</span> \
               <span class="s2">SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15</span>
</code></pre></td>
		<td><pre class="highlight"><code>

<span class="nb">%_1.trim.fq %_2.trim.fq %_1un.trim.fq %_2un.trim.fq</span>: <span class="kt">%_1.fq %_2.fq</span>
	<span class="s2">trimmomatic PE $^ \
		$(shell basename $(word 1,$^) .fq).trim.fq \
		$(shell basename $(word 1,$^) .fq)un.trim.fq \
		$(shell basename $(word 2,$^) .fq).trim.fq \
		$(shell basename $(word 2,$^) .fq)un.trim.fq \
		SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15</span>
</code></pre></td>
	</tr>
	<tr>
		<td>Run these two individual trimmomatic commands with these hardcoded filenames. If you need to run a new file make sure to carefully change all of the uses of it!</td>
		<td>A more generic rule to generate the trimmed fastq files from whichever pair of sequence identifiers were provided. By using <code>%</code> we also get some safety in this function that we didn't have in the bash version! We know for sure those are identical values, and no one accidentally made a small typo in any one of the <b>6</b> times the same identifier was repeated. <u>But it comes at the cost</u> of using some complex Make statements like <code>$(shell)</code> which lets us execute arbitrary shell commands inside our pipeline. <br/> You'll notice that what we're really doing is extracting the inputs and outputs from the command and making them much more uniform</td> </tr>

</tbody>
</table>

> ### {% icon comment %} Don't remember this
> There is a lot of Make specific things going on in the above examples. You do not need to know it! We just wanted to provide working, correct examples. This tutorial is focused on Snakemake, so read these for context and understanding of why we use Snakemake, not reading them to learn Makefiles.
{: .comment}

## Why Make?

So now comes the question, why Make? Why would you want to write the rules in this declarative way, rather than the imperative way that is so much easier and requires so much less work? **Speed!** Why would you want to write these generic rules that require learning a lot of Make's syntax? **Reusability**.

Speed
:  Every time you want to run your declarative pipeline, you either need to run it from start to finish every time (re-downloading sequence data, re-building genome indicies, re-mapping reads) or to write a lot of bash code to check if those files exist and only downloading/indexing/alinging as-needed.
Reusability
:  That pipeline can only download those specific files, unless you write additional code to template out the name of the reads and genome you want to align. By writing more generic rules, as soon as someone gives you new datasets, you can immediately start processing those with your re-usable pipeline.

However with `make`, you've written generic rules which can be used to download any fastq files. And best of all, make can check if the files already exist, and if they do, it skips the step and goes on to the next step, rather than re-creating the file. Make does this by checking the

Aspect | Bash | Make
--- |--- | ---
Language style | Individual steps performed line-by-line | Generic rules are written which say how to do the operation, but not on which data
Partial re-run | You must run the entire script every time (or write extra code) | Only the missing files are created
How it is invoked | `bash script.sh` | `make`
Running with different data? | Edit the script to replace the identifiers, or support templated identifiers. | `make read1111.sam read2222.sam read3333.same read4444.sam`
Paralellisation? | None by default, you must edit the script to add it. | `make -j 8` runs each build step on one thread, with 8 threads available, until all tasks are finished.
Filename | Anything ending in `.sh` | `Makefile` is the default name, and you should name your makefile this, unless you want people to have to type `make -f other-file.mk`
Dependencies | Up to you to manage | Up to you to manage
Multiple output files per step | n/a | [Very tricky to get completely right](https://www.gnu.org/software/automake/manual/html_node/Multiple-Outputs.html), Makefiles really expect one rule produces one output file (and can only check e.g. update times of a single file.)
Cluster/HPC Friendliness | Everything is manual | Everything is manual

## Backwards

A full makefile example

```
all: SRR2589044.bam

# Here we've hardcoded the genome name because it's less likely to change for a
# single pipeline than the individual data files are.
%.bam: %_1.trim.fq %_1.trim.fq GCA_000017985.1_ASM1798v1_genomic.fna.bwt
	bwa mem GCA_000017985.1_ASM1798v1_genomic.fna $(word 1,$^) $(word 2,$^) | \
		samtools sort -O bam -o $@

# This indexing step however will work for any possible
%.fna.bwt: %.fna
	bwa index $<

# This handles ALL fastqc reporting, so we don't have to do it in two sections,
# but, unless we ask for these reports they won't be generated.
%.fastqc.html: %.fq
	fastqc $<

# This rule violates some Makefile internal expectations by having multiple
# outputs which is not handled well by all implementations
%_1.trim.fq %_2.trim.fq %_1un.trim.fq %_2un.trim.fq: %_1.fq %_2.fq
	trimmomatic PE $^ \
		$(shell basename $(word 1,$^) .fq).trim.fq \
		$(shell basename $(word 1,$^) .fq)un.trim.fq \
		$(shell basename $(word 2,$^) .fq).trim.fq \
		$(shell basename $(word 2,$^) .fq)un.trim.fq \
		SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

# And here's finally our download step
%:
	wget https://.../$(shell basename $@) \
		-O $@
```

You'll notice a couple things about this above example:

1. There is a new rule called `all`. In a makefile, by default, the very first rule is executed when you run `make`. If you want to execute other rules you can, but it defaults to the first one. By convention, many people name it `all`.
2. Inside there we've also written a file we'd like created, `SRR2589044.bam`, which doesn't exist yet. Make sees this as a dependency to finishing the (empty) all rule, and then goes on to figure out how to create it.
3. It is written backwards, we've started with what we want to output, and for each line, we figured out what we needed for that, and wrote a rule on how to create it. This is relatively common in makefiles.

You'll also notice some weird additional things we've had to do like `$(word 2,$^)` to get the second input file to that rule, this is really kind of ugly and hard to understand and is a great motivation for learning Snakemake which helps address these issues.

Make will read the above makefile like so:
<TODO>

Reading the above you should be able to imagine a tree of tasks that Make is creating internally:

<TODO>

This is called a "Directed Acyclic Graph", it is a graph of nodes (tasks that need to be executed), with connections between nodes that have a direction (this task depends on outputs of that task), and there are no cycles (no outputs depend on inputs.) These are very common in {SciWMS}s because they make computation faster. Instead of executing step-by-step, we can build this graph of what work needs to be done, and then starting with the leaves of the graph (the end nodes without other dependencies) we can start executing and removing them.

This also how we can really easily parallelise workflows: because we know the dependencies of each step, we know which can be executed right now, and we can execute them in parallel because we know for sure they do not depend on each other.

# Snakemake

> The Snakemake workflow management system is a tool to create reproducible and scalable data analyses. Workflows are described via a human readable, Python based language. They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment. [source](https://snakemake.readthedocs.io/en/stable/index.html)
{: .quote}

Snakemake addresses a lot of the issues with make for use in scientific contexts: clearer pipelines and dependencies. We did not talk about it in the previous section, but where did `bowtie2` and `bowtie2-build` come from? How did those get installed? What versions are they? None of that information is included in the Makefile

Snakemake rules are a bit more complex, in Snakemake you will write rules that follow this form:

<pre class="highlight"><code>
rule {name}
	<span class="kt">input:
		"something",
		"something-else"</span>
	<span class="nb">output:
		"output-1",
		"output-2"</span>
	conda:
		"envs/mapping.yaml"
	<span class="s2">shell:
		"cat {input} > {output}"</span>
</code></pre>

## Line-By-Line Comparison

<table class="table">
<thead><tr><th>Make</th><th>Snakemake</th></tr></thead>
<tbody>
	<tr>
		<td><pre class="highlight"><code>
<span class="nb">%.fq.gz</span>:
	<span class="s2">wget https://ncbi.example.org/$@</span>
</code></pre></td>
		<td><pre class="highlight"><code>
rule download:
	<span class="nb">output:
		"{sample}.fq.gz"</span>
	<span class="s2">shell:
		"wget https://ncbi.example.org/{sample}.fq.gz -O {output}"</span>
</code></pre></td>
	</tr>
	<tr>
		<td>Generic download rule, the <code>$@</code> and <code>%</code> used are a bit opaque, you need to know what they mean to understand how the rule works.</td>
		<td>This is much more explicit, the outputs are listed and <code>{sample}</code> is used as the variable to be templated out, a lot like you might recognise from Python's <code>format</code> function or <code>f""</code> strings.</td>
	</tr>


	<tr>
		<td><pre class="highlight"><code>
<span class="nb">%.fastqc.html</span>: <span class="kt">%.fq</span>
	<span class="s2">fastqc $<</span>
</code></pre></td>
		<td><pre class="highlight"><code>
rule fastqc:
	<span class="kt">input:
		"{sample}.fq"</span>
	<span class="nb">output:
		"{sample}.fastqc.html"</span>
	conda:
		"envs/fastqc.yaml"
	<span class="s2">shell:
		"fastqc {sample}"</span>
</code></pre></td>
	</tr>
	<tr>
		<td>Here is a rule to generate a single FastQC report from a single FastQ file</td>
		<td>Essentially the same, but now we've also added a Conda environment in which our job will run. This makes dependency management a lot simpler. The <code>envs/fastqc.yaml</code> file contains a list of the necessary dependencies at the correct version and if the environment for it does not exist, Snakemake will create it.</td>
	</tr>

	<tr>
		<td><pre class="highlight"><code>
<span class="nb">%_1.trim.fq %_2.trim.fq %_1un.trim.fq %_2un.trim.fq</span>: <span class="kt">%_1.fq %_2.fq</span>
	<span class="s2">trimmomatic PE $^ \
		$(shell basename $(word 1,$^) .fq).trim.fq \
		$(shell basename $(word 1,$^) .fq)un.trim.fq \
		$(shell basename $(word 2,$^) .fq).trim.fq \
		$(shell basename $(word 2,$^) .fq)un.trim.fq \
		SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15</span>
</code></pre></td>
		<td><pre class="highlight"><code>
rule trimmomatic:
	<span class="kt">input:
		r1="{sample}_1.fq",
		r2="{sample}_2.fq"</span>
	<span class="nb">output:
		o1="{sample}_1.trim.fq",
		o2="{sample}_2.trim.fq",
		o1un="{sample}_1un.trim.fq",
		o2un="{sample}_2un.trim.fq"</span>
	conda:
		"envs/trimmomatic.yaml"
	<span class="s2">shell:
		"trimmomatic PE "
		"{input.r1} {input.r2} "
		"{output.o1} {output.o1un} "
		"{output.o2} {output.o2un} "
		"SLIDINGWINDOW:4:20 MINLEN:25 "
		"ILLUMINACLIP:NexteraPE-PE.fa:2:40:15"</span>
</code></pre></td>
	</tr>

	<tr>
		<td>Here we take our very complicated and hard to understand Make rule (shell? basename? word?) with ugly and potentially quite broken multiple output syntax</td>
		<td>And turn it into a much more readable and clear Snakemake step!</td>
	</tr>

</tbody>
</table>

Now that you have seen a few rules, let's write the rest:

> ### {% icon question %} Question
>
> How would you write the following task in Snakemake?
>
> The command is
> <pre class="highlight"><code>
> <span class="s2">bwa index</span> <span class="kt">GCA_000017985.1_ASM1798v1_genomic.fna</span></code></pre>
> and it creates `GCA_000017985.1_ASM1798v1_genomic.fna.bwt`
>
> > ### {% icon solution %} Solution
> >
> > <pre class="highlight"><code>
> > rule indexgenome:
> > 	<span class="kt">input:
> > 		"GCA_000017985.1_ASM1798v1_genomic.fna"</span>
> > 	<span class="nb">output:
> > 		"GCA_000017985.1_ASM1798v1_genomic.fna.bwt"</span>
> > 	conda:
> > 		"envs/bwa.yaml"
> > 	<span class="s2">shell:
> > 		"bwa index {input}"</span>
> > </code></pre></td>
> {: .solution}
{: .question}

> ### {% icon question %} Question
>
> The command is `bwa mem GCA_000017985.1_ASM1798v1_genomic.fna SRR2584866_1.trim.fq SRR2584866_2.trim.fq | samtools sort -O bam -o SRR2584866.bam`
>
> 1. What are the inputs?
> 2. What are the outputs?
> 3.  How would you write the following task in Snakemake?
>
>
> > ### {% icon solution %} Solution
> >
> > 1. `GCA_000017985.1_ASM1798v1_genomic.fna.bwt`, the index file (but beware it does not get passed in as-is, the indexing tool expects just the `GCA_000017985.1_ASM1798v1_genomic.fna` portion.)
> >
> >    Also we have our two sequence files `SRR2584866_1.trim.fq SRR2584866_2.trim.fq`
> >
> > 2. Our output is `SRR2584866.bam`
> >
> > 3. There are a couple of options we have here, we can supply both the `.fna` and the `.fna.bwt` file as inputs (not strictly true, we don't need the fasta file) and then just use the `fna` file in the command line, or we can pass in just the `.fna.bwt` file and try and calculate the `.fna` version that is expected as the index name. We will show the second option as it is more complicated.
> > <pre class="highlight"><code>
> > rule align:
> > 	<span class="kt">input:
> > 		r1="{sample}_1.trim.fq",
> > 		r2="{sample}_1.trim.fq",
> > 		index="{genome}.fna.bwt"</span>
> > 	<span class="nb">output:
> > 		"{genome}/{sample}.bam"</span>
> > 	conda:
> > 		"envs/bwa.yaml"
> > 	<span class="s2">shell:
> > 		"bwa mem {wildcard.genome}.fna {input.r1} {input.r2} | "
> > 		"samtools sort -O bam -o {output}"</span>
> > </code></pre></td>
> >
> > Here we used a number of features to accomplish what we need, and we'll now go through them. First is [Wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards) which can be used to take a portion of the output name or a portion of the input name and to re-use that in the command line. Here we declared that the first part of the index name up to `.fna.bwt` was going to be the `genome` wildcard.
> >
> > Importantly, we also used this in our output. What would have happened if we didn't? It would be unresolvable! We would run `snakemake ... output.bam` and it would say "I don't know what value genome should be set to", so we need to have that value somewhere in our output filename in order to be able to figure that out.
> >
> > That isn't the only way to solve that problem, we could also hardcode this or write it in a [config file](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) that is used by snakemake.
> {: .solution}
{: .question}

## Best Practices

But this was our very first attempt at a workflow, so what might a best practice workflow look like?

1. Conda for reproducibility
2. Saves log files to ensure we can see what went wrong if anything.

> ### {% icon code-in %} `envs/bwa.yaml`
> ```yaml
> channels:
>   - bioconda
>   - conda-forge
> dependencies:
>   - samtools =1.9
> ```
{: .code-in}



## Why Snakemake

So now comes the question, why Snakemake? **Better for science**. While it is quite similar to good old `make`, Snakemake adds several features that are important for science like dependency management with Conda/Docker/Singularity, and better execution on HPCs and Clusters.

Aspect | Make | Snakemake
--- |--- | ---
Language style | Generic rules are written which say how to do the operation, but not on which data | Same
Partial re-run | Only the missing files are created | Same
How it is invoked | `make` | `snakemake`
Running with different data? | `make read1111.sam read2222.sam` | `snakemake read1111.sam read2222.sam`
Paralellisation? | `make -j 8` | `snakemake --cores 8`
Filename | `Makefile` | `Snakefile`
Dependencies | Up to you to manage | Built-in dependency management with Conda
Multiple output files per step | A bit tricky | Incredibly easy
Cluster/HPC Friendliness | Everything is manual | Very good support
