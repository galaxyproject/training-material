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

# Get Data

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


# Generate a transcript to gene map and filtered FASTA

You can have a look at the Terminal tab again. Has the package `atlas-gene-annotation-manipulation` been installed yet? If yes, you can execute the code cell below and while it's running, I'll explain all the parameters we set here. 

```bash
gtf2featureAnnotation.R -g gtf.gff -c fasta.fasta -d "transcript_id" -t "transcript" -f "transcript_id" -o map_code -l "transcript_id,gene_id" -r -e filtered_fasta_code 
```

In essence, [gtf2featureAnnotation.R script](https://github.com/ebi-gene-expression-group/atlas-gene-annotation-manipulation) takes a GTF annotation file and creates a table of annotation by feature, optionally filtering a cDNA file supplied at the same time. Therefore the first parameter `-g` stands for "gtf-file" and requires a path to a valid GTF file. Then `-c` takes a cDNA file for extracting meta info and/or filtering - that's our FASTA! Where --parse-cdnas (that's our `-c`) is specified, we need to specify, using `-d`, which field should be used to compare to identfiers from the FASTA. We set that to "transcript_id" - feel free to inspect the GTF file to explore other attributes. We pass the same value in `-f`, meaning first-field, ie. the name of the field to place first in output table. To specify which other fields to retain in the output table, we provide comma-separated list of those fields, and since we're only interested in transcript to gene map, we put those two names ("transcript_id,gene_id") into `-l`. `-t` stands for the feature type to use, and in our case we're using "transcript". Guess what `-o` is! Indeed, that's the output annotation table - here we specify the file path of our transcript to gene map. We will also have another output denoted by `-e` and that's the path to a filtered FASTA. Finally, we also put `-r` which is there only to suppress header on output. Summarising, output will be a an annotation table, and a FASTA-format cDNAs file with unannotated transcripts removed.

Why filtered FASTA?
Sometimes it's important that there are no transcripts in a FASTA-format transcriptome that cannot be matched to a transcript/gene mapping. Salmon, for example,  used to produce errors when this mismatch was present. We can synchronise the cDNA file by removing mismatches as we have done above.


# Generate a transcriptome index

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

# Use Alevin

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

We have also added some additional parameters (`--freqThreshold`, `--keepCBFraction`) and their values are derived from the [Alevin Galaxy tutorial]({% link topics/single-cell/tutorials/scrna-case_alevin/tutorial.md %}) after QC to stop Alevin from applying its own thresholds.

Once all the above requirement are satisfied, Alevin can be run using the following command:

```bash
salmon-latest_linux_x86_64/bin/salmon alevin -l ISR -1 barcodes_701.fastq -2 transcript_701.fastq --dropseq  -i salmon_index_code -p 10 -o alevin_output_code --tgMap map_code --freqThreshold 3 --keepCBFraction 1 --dumpFeatures
```


This tool will take a while to run. Alevin produces many file outputs, not all of which we'll use. You can refer to the [Alevin documentation](https://salmon.readthedocs.io/en/latest/alevin.html) if you're curious what they all are, you can look through all the different files to find information such as the mapping rate, but we'll just pass the whole output folder directory for downstream analysis. 


> <warning-title>Process stopping</warning-title>
>  
> The command above will display the log of the process and will say "Analyzed X cells (Y% of all)". For some reason, running Alevin may sometimes cause problems in Jupyter Notebook and this process will stop and not go to completion. This is the reason why we use hugely subsampled dataset here - bigger ones couldn't be fully analysed (they worked fine locally though). The dataset used in this tutorial shouldn't make any issues when you're using Jupyter notebook through galaxy.eu, however might not work properly on galaxy.org. If you're accessing Jupyter notebook via galaxy.eu and alevin process stopped, just restart the kernel and that should help.
> 
{: .warning}


<!---
check if we can get alevinQC to work - paste the info from the other tutorial?
-->

# Alevin output to SummarizedExperiment 

Let's change gear a little bit. We've done the work in bash, and now we're switching to R to complete the processing. To do so, you have to change Kernel to R (either click on `Kernel` -> `Change Kernel...` in the upper left corner of your JupyterLab or click on the displayed current kernel in the upper right corner and change it).
![Figure showing the JupyterLab interface with an arrow pointing to the left corner, showing the option `Kernel` -> `Change Kernel...` and another arrow pointing to the right corner, showing the icon of the current kernel. The pop-up window asks which kernel should be chosen instead.](../../images//switch_kernel.jpg "Two ways of switching kernel.")

Now load the library that we have previously installed in terminal:

```r
library(tximeta)
```

The [tximeta package](https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html) REF (Love et al. 2020) is used for import of transcript-level quantification data into R/Bioconductor and requires that the entire output of alevin is present and unmodified. 

First, let's specify the path to the quants_mat.gz file: 

```r
path <- 'alevin_output/alevin/quants_mat.gz'
```
We will specify the following arguments when running *tximeta*:
- 'coldata' a data.frame with at least two columns:
  - files - character, paths of quantification files
  - names - character, sample names
- 'type' - what quantifier was used (can be 'salomon', 'alevin', etc.)

With that we can create a dataframe and pass it to tximeta to create SummarizedExperiment object.

```r
coldata <- data.frame(files = path, names="sample701")
alevin_se <- tximeta(coldata, type = "alevin")
```

Inspect the created object:
```r
alevin_se
```

As you can see, *rowData names* and *colData names* are still empty. Before we add some metadata,  we will first identify barcodes that correspond to non-empty droplets. 

# Identify barcodes that correspond to non-empty droplets 

Some sub-populations of small cells may not be distinguished from empty droplets based purely on counts by barcode. Some libraries produce multiple ‘knees’ (see the [Alevin Galaxy tutorial]() for multiple sub-populations. The [emptyDrops]() method has become a popular way of dealing with this. emptyDrops still retains barcodes with very high counts, but also adds in barcodes that can be statistically distinguished from the ambient profiles, even if total counts are similar. 

```r
library(DropletUtils)               # load the library and required packages
```

emptyDrops takes multiple arguments that you can read about in the [documentation](https://rdrr.io/github/MarioniLab/DropletUtils/man/emptyDrops.html). However, in this case, we will only specify the following arguments:

- `m` -	A numeric matrix-like object - usually a dgTMatrix or dgCMatrix - containing droplet data prior to any filtering or cell calling. Columns represent barcoded droplets, rows represent genes.
- `lower` - A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
- `niters` - An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations.
- `retain` - A numeric scalar specifying the threshold for the total UMI count above which all barcodes are assumed to contain cells.

Let's then extract the matrix from our `alevin_se` object. It's stored in *assays* -> *counts*. 

```r
matrix_alevin <- assays(alevin_se)$counts
```

And now run emptyDrops:
```r
# Identify likely cell-containing droplets
out <- emptyDrops(matrix_alevin, lower = 100, niters = 1000, retain = 20)
out
```
<!---
comment on those values
-->

False discovery rate - ???
```r
is.cell <- out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)
```

We got rid of the background droplets containing no cells, so now we will filter the matrix that we passed on to emptyDrops, so that it corresponds to the remaining cells. 

```r
emptied_matrix <- matrix_alevin[,which(is.cell),drop=FALSE]          # filter the matrix
dim(emptied_matrix)                                                  # check the dimension of the filtered matrix
```

From here, we can move on to adding cell metadata.

# Adding cell metadata

The genes IDs are stored in *colnames*. Let's exctract them into a separate object:
```r
barcode <- colnames(alevin_se)
```

Now, we can simply add those barcodes into *colData names* which stores cell metadata. To do this, we will create a column called `barcode` in *colData* and pass the stored values into there.

```r
colData(alevin_se)$barcode <- barcode
```

As we saw above, the dimension of the filtered matrix is A x B. It means that there are X cells and Y genes. We will now extract those cells from the filtered matrix. 

```r
retained_cells <- colnames(emptied_matrix)
retained_cells
```

Now, we can simply add those barcodes into *rowData names* which stores gene metadata. To do this, we will create a column called `gene_ID` in *rowData* and pass the stored values into there.

 
# Adding gene metadata 

As you saw above, the genes IDs are stored in *rownames*. Let's exctract them into a separate object:

```r
gene_ID <- rownames(alevin_se)
```

Now, we can simply add those genes IDs into *rowData names* which stores gene metadata. To do this, we will create a column called `gene_ID` in *rowData* and pass the stored values into there.

```r
rowData(alevin)$gene_ID <- gene_ID
```

## Adding genes symbols based on their IDs

Since gene symbols are much more informative than only gene IDs, we will add them to our metadata. We will base this annotation on Ensembl - the genome database – with the use of the library BioMart. We will use the archive Genome assembly GRCm38 to get the gene names. Please note that the updated version (GRCm39) is available, but some of the gene IDs are not in that EnsEMBL database. The code below is written in a way that it will work for the updated dataset too, but will produce ‘NA’ where the corresponding gene name couldn’t be found.

```r
# get relevant gene names
library("biomaRt")                                      # load the BioMart library
ensembl.ids <- gene_ID                               
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")    # connect to a specified BioMart database and dataset hosted by Ensembl
ensembl_m = useMart("ensembl", dataset="mmusculus_gene_ensembl", host='https://nov2020.archive.ensembl.org') 	

# The line above connects to a specified BioMart database and dataset within this database.
# In our case we choose the mus musculus database and to get the desired Genome assembly GRCm38,
# we specify the host with this archive. If you want to use the most recent version of the dataset, just run:
# ensembl_m = useMart("ensembl", dataset="mmusculus_gene_ensembl")
```
```r
genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
               filters = 'ensembl_gene_id',
               values = ensembl.ids,
               mart = ensembl_m)

# The line above retrieves the specified attributes from the connected BioMart database;
# 'ensembl_gene_id' are genes IDs,
# 'external_gene_name' are the genes symbols that we want to get for our values stored in ‘ensembl.ids’.
```
```r
# see the resulting data
head(genes)                          
```
```r
# replace IDs for gene names
gene_names <- ensembl.ids	 
count = 1 	 
for (geneID in gene_names)
{
 index <- which(genes==geneID)    # finds an index of geneID in the genes object created by getBM()
 if (length(index)==0)            # condition in case if there is no corresponding gene name in the chosen dataset
  {
    gene_names[count] <- 'NA'
  }
  else
  {
    gene_names[count] <- genes$external_gene_name[index] 	# replaces gene ID by the corresponding gene name based on the found geneID’s index
  }
 count = count + 1                # increased count so that every element in gene_names is replaced
}
```
```r
# add the gene names into rowData in a new column gene_name
rowData(alevin_se)$gene_name <- gene_names
```
```r
# see the changes
rowData(alevin_se)                  
```

If you are working on your own data and it’s not mouse data, you can check available datasets for other species and just use relevant dataset in `useMart()` function.
```r
listDatasets(mart)                # available datasets
```

> <warning-title>Ensembl connection problems</warning-title>
> Sometimes you may encounter some connection issues with Ensembl. To improve performance Ensembl provides several mirrors of their site distributed around the globe. When you use the default settings for useEnsembl() your queries will be directed to your closest mirror geographically. In theory this should give you the best performance, however this is not always the case in practice. For example, if the nearest mirror is experiencing many queries from other users it may perform poorly for you. In such cases, the other mirrors should be chosen automatically.
>
{: .warning}


<!---
add mito annotation
-->



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

