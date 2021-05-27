---
layout: tutorial_hands_on

title: "SRA Aligned Read Formats to Speed Up SARS-CoV-2 data Analysis"
questions:
- How can I search SRA SARS-CoV-2 metadata from within Galaxy?
- How can I import SRA aligned read files and extract the data in my format of choice?
- How can I import vcf files into Galaxy that have been generated for these Runs?
objectives:
- Learn about SRA aligned read format and vcf files for Runs containing SARS-CoV-2 content
- Understand how to search the metadata for these Runs to find your dataset of interest and then import that data in your preferred format
time_estimation: "30m"
key_points: []
contributors:
  - jontrow
  - adelaiderhodes
---

Galaxy tutorial to be presented by Jon Trow at GCC 2021

# Description of the Workshop

Traditionally, after a list of run accessions has been filtered on the NCBI website, the accessions are used to download and extract fastq using the SRA toolkit to enter into the next steps of the workflow. A newer compressed data type, generated from raw submitted data containing SARS-CoV-2 sequence, is also accessible to Galaxy users from SRA in the Cloud.


SRA Aligned Read Format (SARF) provides further output options beyond basic fastq format, for example:

1. contigs created from the raw reads in the run
2. reads aligned back to the contigs
3. reads with placeholder quality scores
4. VCF files can also be downloaded for these records relative to the SARS-CoV-2 RefSeq record

- These formats can speed up workflows such as assembly and variant calling.
- This data format is still referenced by the Run accession and accessed using the sra toolkit.
- This workshop describes the SARF data objects along with associated searchable metadata, and demonstrates a few ways to enter them into traditional workflows.


> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Introduction

The aim of this tutorial is to introduce you to some of SRA's new SARS-CoV-2 cloud resources and data formats, then show you how to filter for Runs of interest to you and access that data in your format of choice in Galaxy to use in your analysis pipeline.

## SRA Aligned Read Format

All data submitted to SRA is scanned with our SARS-CoV-2 Detection Tool which uses a Kmer-based approach to identify Runs with Coronaviridae content. The initial scope of the project is limited to those runs deposited in SRA with at least 100 hits for SARS-CoV-2 via the SARS-CoV-2 Detection Tool, a read length of at least 75, and generated using the Illumina platform.

1. For these Runs, Saute was used to assemble contigs via guided assembly, with the SARS-CoV-2 refseq genomic sequence (NC_045512.2) used as the guide.

2. If contigs were successfully assembled, reads were mapped back to the contigs and coverage calculated. These contigs with the reads mapped back and with quality scores removed (to keep the object size small) are the aligned read format files.

3. The SRA toolkit can be used to dump just the contigs in fasta format, the reads aligned to the contigs in sam format or the raw reads in fastq format with placeholder quality scores.

4. The contigs were also assessed via megablast against the nucleotide blast database and the results made available for search.

5. In addition, to support investigation of viral evolution during the pandemic and after the introduction of vaccines, variants are identified relative to the SARS-CoV-2 RefSeq record for each processed run using BCFTools.

The SRA aligned reads, the VCF files, the results of these analyses (such as BLAST and VIGOR3 annotation), and the associated BioSample and sequencing library metadata are available for free access from cloud providers.


> ### {% icon comment %} Comment
>
> These data can be dumped in `sam` format using the `sam-dump` tool in the SRA Toolkit, but this function doesn't work within Galaxy yet.
> We hope to include that functionality in a future update.
>
{: .comment}


# Finding Runs of Interest for SARS-CoV-2 Submissions to SRA

Metadata for SARS-CoV-2 submissions to the SRA includes submitted sample and library information, BLAST results, descriptive contig statistics, and variation and annotation information.  These metadata are updated daily and made available to query in the cloud using [Google's BigQuery](https://www.ncbi.nlm.nih.gov/sra/docs/sra-bigquery/) or [Amazon's Athena](https://www.ncbi.nlm.nih.gov/sra/docs/sra-athena/) services. However, the raw underlying information is also provided as a group of json files that can be **downloaded for free from the Open Data Platform** without logging in to the cloud. These json files can be imported to Galaxy and queried there to find **Runs of interest**.

> ### {% icon comment %} Comment
>
> Some of these tables include complex data fields (array of values) that don't have a clean analogue in a classic SQL database or table and these can't be easily queried in Galaxy currently. If you require access to [cloud tables](https://www.ncbi.nlm.nih.gov/sra/docs/aligned-metadata-tables/) or fields not available in Galaxy we recommend accessing those natively in [BigQuery](https://www.ncbi.nlm.nih.gov/sra/docs/sra-bigquery/) or [Athena](https://www.ncbi.nlm.nih.gov/sra/docs/sra-athena/).
>
{: .comment}


We will import the JSON files into Galaxy to query them directory, however the files are split up for efficient querying in the cloud and updated daily, so we first need to get the most up-to-date list of files so we can import those to Galaxy. We'll just be using a couple of tables in this training, but the other tables can be imported in the same way, using the index files and awk commands below.

> ### {% icon tip %} Metadata tables for SARF
>
> The file index for these tables can be found here:
>
> * [https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/contigs.filelist](https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/contigs.filelist)
> * [https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/annotated_variations.filelist](https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/annotated_variations.filelist)
> * [https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/blastn.filelist](https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/blastn.filelist)
> * [https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/metadata.filelist](https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/metadata.filelist)
> * [https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/peptides.filelist](https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/peptides.filelist)
> * [https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/tax_analysis.filelist](https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/tax_analysis.filelist)
> * [https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/variations.filelist](https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/variations.filelist)
>
{: .tip}


> ### {% icon comment %} Comment
>
> These metadata files are updated daily around 5:30pm EST. If you try to access the data around this time but encounter an error, trying again a short while later should resolve the issue.  [Time Zone Converter](https://www.thetimezoneconverter.com)
>
{: .comment}


> ### {% icon hands_on %} Hands-on 1: Loading SRA Aligned Read Format (SARF) Object Metadata into Galaxy
>
>This step needs to be repeated at the beginning of an analysis to refresh the metadata to the latest daily version.
>
> 1. Go to your Galaxy instance of choice such as one of the [usegalaxy.org](https://usegalaxy.org/), [usegalaxy.eu](https://usegalaxy.eu), [usegalaxy.org.au](https://usegalaxy.org.au) or any other. (This tutorial uses usegalaxy.org).
> 2. If your history is not already empty, then start a new history (see [here](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/history/tutorial.html) for more on Galaxy histories)
> 3. Click the `Upload Data` button, then click `Paste/Fetch data`, and copy/paste the cloud address for the `contigs.filelist` table into the provide box:
>
>    ```markdown
>     https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/contigs.filelist
>
>    ```
>
>    ![Paste Fetch Command ](../../images/sarf/1_paste_fetch.PNG "Paste Fetch Command")
>
> 4. Then click the `start` button followed by the `close` button. This download job should now appear to the right in your history bar.
> 5. From the history bar you can click on the `edit` icon to update the name, in this case we'll call it `contigs.json.list` and click the `save` button.
> This file contains the filenames for all json files that comprise the current 'contigs' table. It is necessary to import the latest version of this file because the metadata is updated daily.
>
> **Next we will convert this list of filenames to the HTTP URLs for easy import into Galaxy.**
>
> 6. Select the `text reformatting with awk` tool.
> 7. For `file to process` select `contigs.json.list`.
> 8. In the `AWK program` box enter:
>
>    ```markdown
>    {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/contigs/"$0}
>    ```
>
>    ![Text Reformatting](../../images/sarf/contigs_new.PNG "Text Reformating with AWK")
>
> 9. Click the `execute` button, then use the `edit` button on the history item to rename the output to `contig_urls`.
>
>    > ### {% icon tip %} AWK commands for the other metadata tables
>    > We are not going to use all of these tables in this demo.
>    > ```markdown
>    > {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/contigs/"$0}
>    > {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/annotated_variations/"$0}
>    > {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/blastn/"$0}
>    > {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/metadata/"$0}
>    > {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/peptides/"$0}
>    > {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/tax_analysis/"$0}
>    > {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/variations/"$0}
>    > ```
>    {: .tip}
>
> **Now we'll import the json files into Galaxy and construct a table for filtering.**
>
> 10. Click the `Upload Data` button, then select the `Rule-based` tab at the top.
>     - `Upload data as > Datasets`
>     - `Load tabular data from > History Dataset`
>     - `History dataset > contig_urls`
> 11. Click the `build` button, then in the `Build Rules for Uploading Collections`
>     - `+Rules > add/modify column definitions`
>     - `+Add definition > URL > Column A > apply`
>
> 12. Click the `upload` button.
>
> **Once those download jobs have all turned green in the history list, we'll concatenate these into a single file.**
>
>
> 13. Select the `Concatenate multiple datasets tail-to-head` tool
> 14. Click on the `multiple datasets` icon, then
>     - Select the `jsonlines.gz` files we just downloaded and click `execute`.
>     - Rename this history item to `contigs.json`
>
> **The last step in loading the data from the json files is to convert the entries to tabular format using the `JQ` tool**
>
> 15. Using the `Search Tools` box on the upper left side of the Galaxy board, enter `json` and click on `JQ`, it should be near the top of suggested items.
> 16. Select the `contigs.json` file that was just created in the `json input` box
> 17. To import the entire json file without filters, enter ` [.[]]` in the `jq filter` box
> 18. Select `yes` under `convert output to tabular`, then click the `execute` button.
> 19. Rename this file `contigs.tsv` using the `edit` icon on the tables history item
>
{: .hands_on}


# Query SARF Metadata

Now that a table has been generated, we will query the table to find the runs of interest.  It is a good idea to save the table for future queries on the same dataset.  Rerunning the import steps above without filtration will provide a different set of metadata each day.

> ### {% icon hands_on %} Hands-on 2: Query the SRA Metadata Table using SQLite
> Next we'll query this metadata using the `Query tabular` tool to get a list of all Runs containing contigs of greater than 20,000 nucleotides and average coverage of at least 100X.
> 1. In the Search Tools Box in the upper left corner of the Galaxy board, type `sql` to find the tool `Query Tabular` for the next steps.
> 2. Click `insert database table`, for `Tabular Dataset for Table` select `contigs.tsv`
> 3. Click `table options`,
>    * for `specify name for table` enter `SARS_contigs`,
>    * for `Specify Column Names` enter:
>    ```markdown
>    name,run,coverage,tax_id,hits,length,md5
>
>    ```
>
> > ### {% icon tip %} Column Headers for the Other Metadata Tables
> > We are not going to bring in the other metadata tables in this tutorial. Here is a list of column headers for contigs and the other tables.
> > You can find full definitions for these columns here:
> >
> > [https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-examples/](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-examples/)
> >
> > [https://www.ncbi.nlm.nih.gov/sra/docs/aligned-metadata-tables/](https://www.ncbi.nlm.nih.gov/sra/docs/aligned-metadata-tables/)
> >
> > **contigs**
> > ```
> > name,run,coverage,tax_id,hits,length,md5
> > ```
> > **annotated_variations**
> > ```
> > run,chrom,pos,id,ref,alt,qual,filter,info,format,sample_a,ac,an,bqb,dp,dp4,dp4_1,dp4_2,dp4_3,dp4_4,idv,imf,mq,mq0f,mqb,mqsb,rpb,sgb,vdb,g_gt,g_pl,g_pl_1,g_pl_2,g_dp,g_ad,g_ad_1,g_ad_2,protein_position,ref_codon,alt_codon,ref_aa,alt_aa,protein_name,protein_length,variation
> > ```
> > **blastn**
> > ```
> > acc,qacc,staxid,sacc,slen,length,bitscore,score,pident,sskingdom,evalue,ssciname
> > ```
> > **metadata**
> > ```
> > acc,assay_type,center_name,consent,experiment,sample_name,instrument,librarylayout,libraryselection,librarysource,platform,sample_acc,biosample,organism,sra_study,releasedate,bioproject,mbytes,loaddate,avgspotlen,mbases,insertsize,library_name,biosamplemodel_sam,collection_date_sam,geo_loc_name_country_calc,geo_loc_name_country_continent_calc,ena_first_public_run,ena_last_update_run,sample_name_sam,datastore_filetype,datastore_provider,datastore_region,attributes,jattr
> > ```
> > **peptides**
> > ```
> > name,contig,mat_peptide,run,location,gene,product,ref_db,ref_id,sequence
> > ```
> > **tax_analysis**
> > ```
> > run,contig,tax_id,rank,name,total_count,self_count,ilevel,ileft,iright
> > ```
> > **variations**
> > ```
> > run,chrom,pos,id,ref,alt,qual,filter,info
> > ```
> >
> >
> {: .tip}
>
> 4. In the `SQL Query to generate tabular output` box enter:
>
>    ```SQL
>     select distinct run from SARS_contigs where length > 20000 and coverage > 100 and run like '%SRR%' order by run asc limit 10
>    ```
>
> 5. Select `no` for `include query result column headers`
> 6. Click the `execute` button and rename the output file to `Run_list`
>
>    ![SQLite Input](../../images/sarf/3_sql_contig.PNG "Query Table with SQLite")
>
{: .hands_on}

> ### {% icon tip %} Download Fastq with Quality Scores
> If you would like to dump the raw, underlying data in fastq format with the original quality scores, you can stop here and use the `Faster Download and Extract Reads in FASTQ` tool.
> 1. select input type > `list of SRA accessions, one per line`
> 2. sra accession list > `Run_list`
> 3. Click `execute`
{: .tip}

> ### {% icon tip %} Importing a list of SRR from Athena or BigQuery
> If you opted to conduct your metadata search in the cloud using AWS Athena or GCP BigQuery instead of importing the json file to Galaxy, you can save a list of your Run accessions from that search result and import that file as the 'Run_list' to proceed with the rest of this tutorial.
{: .tip}


# Importing SARFs of Interest

Now that we have assembled a list of Runs that have contigs we are interested in, we'll construct the path to the **SARFS** in the cloud and import those to Galaxy so we can work with them.

> ### {% icon hands_on %} Hands-on 3: Importing SARFs of Interest
> 1. Select the `text reformatting with awk` tool.
> 2. For `file to process` select `Run_list`.
> 3. In the `AWK program` box enter:
>
>    ```bash
>    {print "https://sra-pub-sars-cov2.s3.amazonaws.com/RAO/"$0"/"$0".realign"}
>    ```
>
> 4. Click the `execute` button, rename the output to `sarf_path`.
>
> **Bringing in the SARF contig files for Runs of Interest**
>
> 5. Click the `upload data` button,
>    - select the `Rule based` tab,
>    - `upload data as` > `datasets`,
>    - `load tabular data from` > `history dataset`,
>    - select `sarf_path` from the popup box.
> 6. In the `Build Rules for Uploading Collections`:
>    - `+Rules > add/modify column definitions`
>    - `>+Add definition > URL > Column A > apply`
> 7. Click the `upload` button.
>
> ![SQLite Input](../../images/sarf/sarf_download.PNG "Query Table with SQLite")
>
> ![SQLite Input 2](../../images/sarf/sarf_download_2.PNG "Query Table with SQLite")
>
> > ### {% icon comment %} Comment
> >
> > Please note that there can be some lag in availability of SARF/VCF files in the cloud (particularly for newly submitted data). So it's possible to get a download error for a file that isn't yet present in the cloud. In these cases waiting ~24 hours will generally resolve the issue and allow you to access the file.
> >
> {: .comment}
>
> **Finally we will use the SRA toolkit to dump the contigs in fasta format.**
>
> 8. Select the `Download and Extract Reads in FASTA/Q` tool.
>
>    `Select input type` > SRA Archive in current history
>    sra archive > multiple datasets > select each .realign object
>
> 9. Click `advanced options`, for `Table name within cSRA object` enter `REFERENCE`
> 10. Click `execute`.
>
> ![Load Contigs in Fasta Format](../../images/sarf/fastq-dump.PNG "Load Contigs in Fastq Format")
>
> The resulting dataset includes the contigs generated from these Runs.
>
> > ### {% icon tip %} Fastq format option
> > If you prefer to dump the raw reads in fastq format with placeholder quality scores, leave the `Table name within cSRA object` field blank.
> >
> > This will generate fastq with quality placeholders `?`, as SARFs do not contain quality scores.
> {: .tip}
{: .hands_on}

# Finding VCFs for SARS-Cov-2 Runs

Metadata about these Runs, including submitted sample and library information, BLAST results, descriptive contig statistics, and variation and annotation information are available to query in the cloud using  [Google's BigQuery](https://www.ncbi.nlm.nih.gov/sra/docs/sra-bigquery/) or [Amazon's Athena](https://www.ncbi.nlm.nih.gov/sra/docs/sra-athena/) services. However, the raw underlying information is also provided as a group of json files that can be **downloaded for free from the Open Data Platform** without logging in to the cloud. The VCF metadata can be imported to Galaxy and queried there to find VCFs of interest.

> ### {% icon comment %} Comment
>
> Some of these tables include complex data fields (array of values) that don't have a clean analogue in a classic SQL database or table and these can't be easily queried in Galaxy currently. If you require access to [cloud tables](https://www.ncbi.nlm.nih.gov/sra/docs/aligned-metadata-tables/) or fields not available in Galaxy recommend accessing those natively in [BigQuery](https://www.ncbi.nlm.nih.gov/sra/docs/sra-bigquery/) or [Athena](https://www.ncbi.nlm.nih.gov/sra/docs/sra-athena/).
>
{: .comment}


> ### {% icon hands_on %} Hands-on 3: Importing VCFs of Interest
> First, we'll import the annotated variation directory index into Galaxy.
> 1. Click the `Upload Data` button, then click `Paste/Fetch data`, and copy/paste the cloud address for the annotated variation table index file into the provide box:
>
>    ```
>    https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/annotated_variations.filelist
>    ```
> 2. Click the `start` button followed by the `close` button. This download job should now appear to the right in your history bar.
> 3. From the history bar you can click on the `edit` icon to update the name, in this case we'll call it `annotated_variation.json.list` and click the `save` button.
>    This files contains the filenames for all json files that comprise the current `annotated variation` table. It is necessary to always import the latest version  because these dumps are updated daily (and the filenames contain the current date).
>
> **Next we will convert this list of filenames to the HTTP URLs for easy import into Galaxy.**
>
> 4. Select the `text reformatting with awk` tool.
> 5. For `file to process` select `annotated_variation.json.list`.
> 6. In the `AWK program` box enter:
>
>    ```
>    {print "https://storage.googleapis.com/nih-sequence-read-archive/SARS_COV_2/annotated_variations/"$0}
>    ```
>
> 7. Click the `execute` button, then use the `edit` button the history item to rename the output to `variation_urls`.
>
> Now we'll import the individual json files that comprise the annotated variation table into Galaxy.
>
> 8. Click the `Upload Data` button, then select the `Rule-based` tab at the top.
>    - `Upload data as > Datasets`
>    - `Load tabular data from > History Dataset`
>    - `History dataset > variation_urls`
> 9. Click the `build` button, then in the `Build Rules for Uploading Collections`
>    - `+Rules > add/modify column definitions`
>    - `+Add definition > URL > Column A > apply`
>
> 10. Click the `upload` button.
>
> **Once those download jobs have all turned green in the history list, we'll concatenate these into a single file.**
>
>
> 11. Select the `Concatenate multiple datasets tail-to-head` tool
> 12. Click on the `multiple datasets` icon, then
>    - Select the `jsonlines.gz` files we just downloaded and click `execute`.
>    - Rename this history item to `variation.json`
>
> Next we'll run the `JQ` tool on this file to convert it to a tabular format.
> 13. Type `json` in the `search tools` box and `JQ` should be one of the top suggested items, click on the tool.
> 14. For `json input` you'll select the `variation.json` file we just created.
> 15. We want to import the entire json file, so in the `jq filter` box enter: [.[]]
> 16. Select `yes` under `convert output to tabular`, then click the `execute` button.
> 17. Rename this file `variation.tsv` using the `edit` icon on its history item.
>
**Query for VCFs of Interest**
>Next we'll query this metadata using the `Query tabular` tool to get a list of all Runs containing a called `E484K` variant.
>
> 1. Click `insert database table`, select `variation.tsv`
> 2. Click `table options`, under `specify name for table` enter `sars_variation`, under `Specify Column Names` enter:
>
>    ```markdown
>    run,chrom,pos,id,ref,alt,qual,filter,info,format,sample_a,ac,an,bqb,dp,dp4,dp4_1,dp4_2,dp4_3,dp4_4,idv,imf,mq,mq0f,mqb,mqsb,rpb,sgb,vdb,g_gt,g_pl,g_pl_1,g_pl_2,g_dp,g_ad,g_ad_1,g_ad_2,protein_position,ref_codon,alt_codon,ref_aa,alt_aa,protein_name,protein_length,variation
>    ```
>
> 3. In the `SQL Query to generate tabular output` box enter:
>
>    ```SQL
>    SELECT distinct run FROM sars_variation where variation = `E484K` and run like `%SRR%` order by run asc limit 10
>    ```
> 4. Select `no` for `include query result column headers`
> 5. Click the `execute` button and rename the output file to `vcf_list`
>
> Now that we have assembled a list of Runs that have VCFs we are interested in, we'll construct the path to the VCF files in the cloud and import those to Galaxy so we can work with them.
> 5. Select the `text reformatting with awk` tool.
> 6. For `file to process` select `vcf_list`.
> 7. In the `AWK program` box enter:
>
>    ```
>    {print "https://sra-pub-sars-cov2.s3.amazonaws.com/VCF/"$0"/"$0".vcf"}
>    ```
> 8. Click the `execute` button, rename the output to `vcf_path`.
>
> **Next we will import these VCF files**
>
> 9. click the `upload data` button, select the `Rule based` tab,
>    - `upload data as > Datasets`,
>    - `load tabular data from > History Dataset`,
>    - `History Dataset > vcf_path`
> 10. In the `Build Rules for Uploading Collections`:
>    - `+Rules > add/modify column definitions`
>    - `+Add definition > URL > Column A > apply`
> 11. Click the `upload` button.
>
> **Use VCFs in Another Galaxy Tool**
> Once you have imported the VCF files, you can use them in your standard pipeline- here we will annotate them with SnpEff.
>
> 12. Select the `SnpEff eff: annotate variants for SARS-CoV-2`
>
>     > ### {% icon comment %} Comment
>     >
>     > Please note that there are 2 Snepff tools, please choose the one for SARS-CoV-2
>     >
>     {: .comment}
>
> 13. Under `Sequence changes` click the `multiple datasets` option, then select the vcf files that were downloaded.
> 14. Click `execute`.
>
{: .hands_on}


# Other NCBI Resources
* Getting started with SRA in the Cloud:  https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud/
* NCBI SARS-CoV-2 Resources: https://www.ncbi.nlm.nih.gov/sars-cov-2/
* STAT tool: https://www.ncbi.nlm.nih.gov/sra/docs/sra-taxonomy-analysis-tool/
* SRA Aligned Read Format: https://www.ncbi.nlm.nih.gov/sra/docs/sra-aligned-read-format/
* VCF generation: https://www.ncbi.nlm.nih.gov/sra/docs/sars-cov-2-variant-calling/

If you have questions or feedback you can email the SRA helpdesk: sra@ncbi.nlm.nih.gov

