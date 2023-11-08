---
layout: tutorial_hands_on

title: 'Generating a single cell matrix using Alevin (bash + R)'
subtopic: single-cell-CS-code
priority: 1
zenodo_link: 

questions:
  - I have some single cell FASTQ files I want to analyse. Where do I start?
  - How to generate a single cell matrix using command line?

objectives:
  - Generate a cellxgene matrix for droplet-based single cell sequencing data
  - Interpret quality control (QC) plots to make informed decisions on cell thresholds
  - Find relevant information in GTF files for the particulars of their study, and include this in data matrix metadata

time_estimation: 1H

key_points:
  - Create a scanpy-accessible AnnData object from FASTQ files, including relevant gene metadata

requirements:
-
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin

follow_up_training:
  -
    type: "internal"
    topic_name: single-cell
    tutorials:
        - scrna-case_alevin-combine-datasets

tags:
- single-cell
- 10x
- paper-replication
- jupyter-notebook
- interactive-tools


contributions:
  authorship:
    - wee-snufkin
    - nomadscientist

  funding:
    - 

notebook:
  language: bash
  snippet: topics/single-cell/tutorials/alevin-commandline/preamble.md
---


# Setting up the environment 

Alevin is a tool integrated with the salmon software, so first we need to get salmon. You can install salmon using bioconda, but in this tutorial we will show an alternative method - downloading the pre-compiled binaries from the [releases page](https://github.com/COMBINE-lab/salmon/releases).

```bash
wget -nv https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
```

Once you've downloaded a specific binary (here we're using version 1.10.0), just extract it like so:

```bash
tar -xvzf salmon-1.10.0_linux_x86_64.tar.gz
```


We're going to use Alevin {% cite article-Alevin %} for demonstration purposes, but we do not endorse one method over another.

## Get Data

We continue working on the same example data - a very small subset of the reads in a mouse dataset of fetal growth restriction {% cite Bacon2018 %} (see the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6945/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6945/)). For the purposes of this tutorial, the datasets have been subsampled to only 50k reads (around 1% of the original files). Those are two fastq files - one with transcripts and the onther one with cell barcodes. You can download the files by running the code below:

```bash
wget -nv https://zenodo.org/
wget -nv https://zenodo.org/
```

<!---
begin
might add question if it's formatted nicely 
how would you know which file is which
-->

 > <question-title></question-title>
>
> Test rendering
>
> > <solution-title></solution-title>
> >
> > is it ok?
> >
> {: .solution}
>
{: .question}

<!---
end
might add question if it's formatted nicely 
how would you know which file is which
-->

Additionally, to map your reads, you will need a transcriptome to align against (a FASTA) as well as the gene information for each transcript (a gtf) file. These files are included in the data import step below. 

```bash
wget -nv https://zenodo.org/
wget -nv https://zenodo.org/
```

<!---
begin
leave as normal text or make it into a tip or extract into data ingest tutorial?
-->

You can also download these for your species of interest [from Ensembl](https://www.ensembl.org/info/data/ftp/index.html). Once you find the cDNA FASTA file you are interested in, right click on the link and choose "Copy link address" and paste it along the command `wget -nv`, then extract it using `tar`. Here is the example how to do it:

```bash
# Getting FASTA file
wget -nv https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/
tar  
```
Do exactly the same to get the GTF file:

```bash
# Getting GTF file
wget -nv https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus
tar  
```

<!---
end
leave as normal text or make it into a tip or extract into data ingest tutorial?
-->

Why do we need FASTA and GTF files? 
To generate gene-level quantifications based on transcriptome quantification, Alevin and similar tools require a conversion between transcript and gene identifiers. We can derive a transcript-gene conversion from the gene annotations available in genome resources such as Ensembl. The transcripts in such a list need to match the ones we will use later to build a binary transcriptome index. If you were using spike-ins, you'd need to add these to the transcriptome and the transcript-gene mapping.

We will use the murine reference annotation as retrieved from Ensembl in GTF format. This annotation contains gene, exon, transcript and all sorts of other information on the sequences. We will use these to generate the transcript-gene mapping by passing that information to a tool that extracts just the transcript identifiers we need.


## Generate a transcript to gene map and filtered FASTA

You can have a look at the Terminal tab again. Has the package `atlas-gene-annotation-manipulation` been installed yet? If yes, you can execute the code cell below and while it's running, I'll explain all the parameters we set here. 

```bash
gtf2featureAnnotation.R -g gtf.gff -c fasta.fasta -d "transcript_id" -t "transcript" -f "transcript_id" -o map_code -l "transcript_id,gene_id" -r -e filtered_fasta_code 
```

In essence, [gtf2featureAnnotation.R script](https://github.com/ebi-gene-expression-group/atlas-gene-annotation-manipulation) takes a GTF annotation file and creates a table of annotation by feature, optionally filtering a cDNA file supplied at the same time. Therefore the first parameter `-g` stands for "gtf-file" and requires a path to a valid GTF file. Then `-c` takes a cDNA file for extracting meta info and/or filtering - that's our FASTA! Where --parse-cdnas (that's our `-c`) is specified, we need to specify, using `-d`, which field should be used to compare to identfiers from the FASTA. We set that to "transcript_id" - feel free to inspect the GTF file to explore other attributes. We pass the same value in `-f`, meaning first-field, ie. the name of the field to place first in output table. To specify which other fields to retain in the output table, we provide comma-separated list of those fields, and since we're only interested in transcript to gene map, we put those two names ("transcript_id,gene_id") into `-l`. `-t` stands for the feature type to use, and in our case we're using "transcript". Guess what `-o` is! Indeed, that's the output annotation table - here we specify the file path of our transcript to gene map. We will also have another output denoted by `-e` and that's the path to a filtered FASTA. Finally, we also put `-r` which is there only to suppress header on output. Summarising, output will be a an annotation table, and a FASTA-format cDNAs file with unannotated transcripts removed.

Why filtered FASTA?
Sometimes it's important that there are no transcripts in a FASTA-format transcriptome that cannot be matched to a transcript/gene mapping. Salmon, for example,  used to produce errors when this mismatch was present. We can synchronise the cDNA file by removing mismatches as we have done above.


## Generate a transcriptome index

We will use Salmon in mapping-based mode, so first we have to build a salmon index for our transcriptome. We will run the salmon indexer as so:

```bash
salmon-latest_linux_x86_64/bin/salmon index -t filtered_fasta_code -i salmon_index_code -k 31
```

Where `-t` stands for our filtered FASTA file, and `-i` is the output the mapping-based index. To build it, the funciton is using an auxiliary k-mer hash over k-mers of length 31. While the mapping algorithms will make used of arbitrarily long matches between the query and reference, the k size selected here will act as the minimum acceptable length for a valid match. Thus, a smaller value of k may slightly improve sensitivity. We find that a k of 31 seems to work well for reads of 75bp or longer, but you might consider a smaller k if you plan to deal with shorter reads. Also, a shorter value of k may improve sensitivity even more when using selective alignment (enabled via the –validateMappings flag). So, if you are seeing a smaller mapping rate than you might expect, consider building the index with a slightly smaller k.


> <details-title>What is the index?</details-title>
>
> To be able to search a transcriptome quickly, salmon needs to convert the text (FASTA) format sequences into something it can search quickly, called an 'index'. The index is in a binary rather than human-readable format, but allows fast lookup by Alevin. Because the types of biological and technical sequences we need to include in the index can vary between experiments, and because we often want to use the most up-to-date reference sequences from Ensembl or NCBI, we can end up re-making the indices quite often. 
>
{: .details}

<!---
reference salmon
-->

<!---
check if we need decoy --decoys decoys.txt 
-->

## Use Alevin

Time to use Alevin now! Alevin works under the same indexing scheme (as salmon) for the reference, and consumes the set of FASTA/Q files(s) containing the Cellular Barcode(CB) + Unique Molecule identifier (UMI) in one read file and the read sequence in the other. Given just the transcriptome and the raw read files, alevin generates a cell-by-gene count matrix (in a fraction of the time compared to other tools).

> <question-title></question-title>
>
> How does Alevin work in detail?
>
> > <solution-title></solution-title>
> >
> > Alevin works in two phases. In the first phase it quickly parses the read file containing the CB and UMI information to generate the frequency distribution of all the observed CBs, and creates a lightweight data-structure for fast-look up and correction of the CB. In the second round, alevin utilizes the read-sequences contained in the files to map the reads to the transcriptome, identify potential PCR/sequencing errors in the UMIs, and performs hybrid de-duplication while accounting for UMI collisions. Finally, a post-abundance estimation CB whitelisting procedure is done and a cell-by-gene count matrix is generated.
> >
> {: .solution}
>
{: .question}

All the required input parameters are described in [the documentation](https://salmon.readthedocs.io/en/latest/alevin.html), but for the ease of use, they are presented below as well:

- `-l`: library type (same as salmon), we recommend using ISR for both Drop-seq and 10x-v2 chemistry.

- `-1`: CB+UMI file(s), alevin requires the path to the FASTQ file containing CB+UMI raw sequences to be given under this command line flag. Alevin also supports parsing of data from multiple files as long as the order is the same as in -2 flag. That's our barcodes_701.fastq file.

- `-2`: Read-sequence file(s), alevin requires the path to the FASTQ file containing raw read-sequences to be given under this command line flag. Alevin also supports parsing of data from multiple files as long as the order is the same as in -1 flag. That's our transcript_701.fastq file.

- `--dropseq` / `--chromium` / `--chromiumV3`: the protocol, this flag tells the type of single-cell protocol of the input sequencing-library. This is a study using the Drop-seq chemistry, so we specify that in the flag.

- `-i`: index, file containing the salmon index of the reference transcriptome, as generated by salmon index command.

- `-p`: number of threads, the number of threads which can be used by alevin to perform the quantification, by default alevin utilizes all the available threads in the system, although we recommend using ~10 threads which in our testing gave the best memory-time trade-off.

- `-o`: output, path to folder where the output gene-count matrix (along with other meta-data) would be dumped. We simply call it alevin_output_code

- `--tgMap`: transcript to gene map file, a tsv (tab-separated) file — with no header, containing two columns mapping of each transcript present in the reference to the corresponding gene (the first column is a transcript and the second is the corresponding gene). In our case, that's map_code generated by using gtf2featureAnnotation.R function. 

- `--freqThreshold` - minimum frequency for a barcode to be considered. We've chosen 3 as this will only remove cell barcodes with a frequency of less than 3, a low bar to pass but useful way of avoiding processing a bunch of almost certainly empty barcodes.

- `--keepCBFraction` - fraction of cellular barcodes to keep. We're using 1 to quantify all!

- `--dumpFeatures` - if activated, alevin dumps all the features used by the CB classification and their counts at each cell level. It’s generally used in pair with other command line flags.

We have also added some additional parameters (`--freqThreshold`, `--keepCBFraction`) and their values are derived from the [Alevin Galaxy tutorial]({% link topics/single-cell/tutorials/scrna-case_alevin/tutorial.md %}) after QC. 

Once all the above requirement are satisfied, Alevin can be run using the following command:

```bash
salmon-latest_linux_x86_64/bin/salmon alevin -l ISR -1 barcodes_701.fastq -2 transcript_701.fastq --dropseq  -i salmon_index_code -p 10 -o alevin_output_code --tgMap map_code --freqThreshold 3 --keepCBFraction 1 --dumpFeatures
```


This tool will take a while to run. Alevin produces many file outputs, not all of which we'll use. You can refer to the [Alevin documentation](https://salmon.readthedocs.io/en/latest/alevin.html) if you're curious what they all are, you can look through all the different files to find different parameters such as the mapping rate, but we'll just pass the whole output folder directory for downstream analysis. 


> <warning-title>Process stopping</warning-title>
>  
> The command above will display the log of the process and will say "Analyzed X cells (Y% of all)". For some reason, running Alevin may sometimes cause problems in Jupyter Notebook and this process will stop and not go to completion. This is the reason why we use hugely subsampled dataset here - bigger ones couldn't be fully analysed (they worked fine locally though). The dataset used in this tutorial shouldn't make any issues when you're using Jupyter notebook through galaxy.eu, however might not work properly on galaxy.org. If you're accessing Jupyter notebook via galaxy.eu and alevin process stopped, just restart the kernel and that should help.
> 
{: .warning}


<!---
check if we can get alevinQC to work - paste the info from the other tutorial?
-->


# Basic QC

The question we're looking to answer here, is: "do we mostly have a single cell per droplet"? That's what experimenters are normally aiming for, but it's not entirely straightforward to get exactly one cell per droplet. Sometimes almost no cells make it into droplets, other times we have too many cells in each droplet. At a minimum, we should easily be able to distinguish droplets with cells from those without.

Now, the image generated here (400k) isn't the most informative - but you are dealing with a fraction of the reads! If you run the total sample (so identical steps above, but with significantly more time!) you'd get the image below.

![raw droplet barcode plots-total](../../images/scrna-casestudy/wab-raw_barcodes-total.png "Total sample - 32,579,453 reads - raw")

This is our own formulation of the barcode plot based on a [discussion](https://github.com/COMBINE-lab/salmon/issues/362#issuecomment-490160480) we had with community members. The left hand plots with the smooth lines are the main plots, showing the UMI counts for individual cell barcodes ranked from high to low. We expect a sharp drop-off between cell-containing droplets and ones that are empty or contain only cell debris. Now, this data is not an ideal dataset, so for perspective, in an ideal world with a very clean 10x run, data will look a bit more like the following taken from the lung atlas (see the [study in Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6653/results/tsne) and the [project submission](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6653/)).

![raw droplet barcode plots - lung atlas](../../images/scrna-casestudy/wab-lung-atlas-barcodes-raw.png "Pretty data - raw")

In that plot, you can see the clearer 'knee' bend, showing the cut-off between empty droplets and cell-containing droplets.

The right hand plots are density plots from the first one, and the thresholds are generated either using [dropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) or by the method described in the discussion mentioned above. We could use any of these thresholds to select cells, assuming that anything with fewer counts is not a valid cell. By default, Alevin does something similar, and we can learn something about that by plotting just the barcodes Alevin retains.

In experiments with relatively simple characteristics, this 'knee detection' method works relatively well. But some populations (such as our sample!) present difficulties. For instance, sub-populations of small cells may not be distinguished from empty droplets based purely on counts by barcode. Some libraries produce multiple 'knees' for multiple sub-populations. The [emptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y) method has become a popular way of dealing with this. emptyDrops still retains barcodes with very high counts, but also adds in barcodes that can be statistically distinguished from the ambient profiles, even if total counts are similar. In order to ultimately run emptyDrops (or indeed, whatever tool you like that accomplishes biologically relevant thresholding), we first need to re-run Alevin, but prevent it from applying its own less than ideal thresholds.

To use emptyDrops effectively, we need to go back and re-run Alevin, stopping it from applying it's own thresholds. Click the re-run icon {% icon galaxy-refresh %} on any Alevin output in your history, because almost every parameter is the same as before, except you need to change the following:

Alevin outputs MTX format, which we can pass to the dropletUtils package and run emptyDrops. Unfortunately the matrix is in the wrong orientation for tools expecting files like those produced by 10X software (which dropletUtils does). We need to 'transform' the matrix such that cells are in columns and genes are in rows.

Alevin outputs MTX format, which we can pass to the dropletUtils package and run emptyDrops. Unfortunately the matrix is in the wrong orientation for tools expecting files like those produced by 10X software (which dropletUtils does). We need to 'transform' the matrix such that cells are in columns and genes are in rows.


> <hands-on-title>Generate gene information</hands-on-title>
>
> 1. {% tool [GTF2GeneList](toolshed.g2.bx.psu.edu/repos/ebi-gxa/gtf2gene_list/_ensembl_gtf2gene_list/1.52.0+galaxy0) %} with the following parameters:
>    - *"Feature type for which to derive annotation"*: `gene`
>    - *"Field to place first in output table"*: `gene_id`
>    - *"Suppress header line in output?"*: `Yes`
>    - *"Comma-separated list of field names to extract from the GTF (default: use all fields)"*: `gene_id,gene_name,mito`
>    - *"Append version to transcript identifiers?"*: `Yes`
>    - *"Flag mitochondrial features?"*: `Yes` - note, this will auto-fill a bunch of acronyms for searching in the GTF for mitochondrial associated genes. This is good!
>    - *"Filter a FASTA-format cDNA file to match annotations?"*: `No` - we don't need to, we're done with the FASTA!
> 2. Check that the output file type is `tabular`. If not, change the file type by clicking the 'Edit attributes'{% icon galaxy-pencil %} on the dataset in the history (as if you were renaming the file.) Then click `Datatypes` and type in `tabular`. Click `Change datatype`.)
> 2. Rename {% icon galaxy-pencil %} the annotation table to `Gene Information`
{: .hands_on}

Inspect {% icon galaxy-eye %} the **Gene Information** object in the history. Now you have made a new key for gene_id, with gene name and a column of mitochondrial information (false = not mitochondrial, true = mitochondrial). We need to add this information into the salmonKallistoMtxTo10x output 'Gene table'. But we need to keep 'Gene table' in the same order, since it is referenced in the 'Matrix table' by row.

> <hands-on-title>Combine MTX Gene Table with Gene Information</hands-on-title>
>
> 1. {% tool [Join two Datasets](join1) %} with the following parameters:
>    - *"Join"*: `Gene Table`
>    - *"Using column"*: `Column: 1`
>    - *"with"*: `Gene Information`
>    - *"and column"*: `Column: 1`
>    - *"Keep lines of first input that do not join with second input"*: `Yes`
>    - *"Keep lines of first input that are incomplete"*: `Yes`
>    - *"Fill empty columns"*: `No`
>    - *"Keep the header lines"*: `No`
>
>
>    If you inspect {% icon galaxy-eye %} the object, you'll see we have joined these tables and now have quite a few gene_id repeats. Let's take those out, while keeping the order of the original 'Gene Table'.
>
>
> 2.  {% tool [Cut columns from a table](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c1,c4,c5`
>    - *"Delimited by"*: `Tab`
>    - *"From"*: output of **Join two Datasets** {% icon tool %}
>
> 3. Rename output `Annotated Gene Table`
{: .hands_on}

Inspect {% icon galaxy-eye %} your `Annotated Gene Table`. That's more like it! You now have `gene_id`, `gene_name`, and `mito`. Now let's get back to your journey to emptyDrops and sophisticated thresholding of empty droplets!

# emptyDrops

emptyDrops {% cite article-emptyDrops %} works with a specific form of R object called a SingleCellExperiment. We need to convert our transformed MTX files into that form, using the DropletUtils Read10x tool:

> <hands-on-title>Converting to SingleCellExperiment format</hands-on-title>
>
> 1. {% tool [DropletUtils Read10x](toolshed.g2.bx.psu.edu/repos/ebi-gxa/dropletutils_read_10x/dropletutils_read_10x/1.0.4+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Expression matrix in sparse matrix format (.mtx)"*: `Matrix table`
>    - {% icon param-file %} *"Gene Table"*: `Annotated Gene Table`
>    - {% icon param-file %} *"Barcode/cell table"*: `Barcode table`
>    - *"Should metadata file be added?"*: `No`
>
> 2. Rename {% icon galaxy-pencil %} output: `SCE Object`
{: .hands_on}

Fantastic! Now that our matrix is combined into an object, specifically the SingleCellExperiment format, we can now run emptyDrops! Let's get rid of those background droplets containing no cells!

