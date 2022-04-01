---
title: How can i load data?
area: introduction
layout: faq
box_type: tip
contributors: [jennaj, bernandez]
---



## Video Tutorials

* Quick help: [**Manupulating NGS data with Galaxy: Getting Data In**](/tutorials/ngs/#getting-ngs-data-in)
* Full tutorial: [**Uploading data into Galaxy**](/tutorials/upload/)
* Dataset Collections, including creation during Upload: [**Processing many samples at once with collections**](https://galaxyproject.org/tutorials/collections/)
* FTP/FTPS tutorial: [**FTP Upload**](/ftp-upload/)

## GTN Tutorials

Most tutorials from the [Galaxy Training Network (GTN)](https://training.galaxyproject.org/) include one or more data upload steps.

* Start with this tutorial, section "Getting Data In", for a quick overview: [NGS Logistics](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-ngs-data-managment/tutorial.html)
* Advanced **Upload** methods are covered in the tutorial topic group [Data Manipulation](https://training.galaxyproject.org/training-material/topics/galaxy-data-manipulation/)
* All Tutorials: https://training.galaxyproject.org/

## Supplimental Step-by-Step

Data is loaded using the tools in the _**Get Data**_ tool group. Some access specific data provider sites that will load data back into your Galaxy history. To directly load your own local data or data from another source, use the tool _**Get Data → Upload File**_ (also accessible from the top of the left tool panel, as seen in the graphics below). Want to practice import/export functions with small sample data? Import the [Upload sample data history here](https://usegalaxy.org/u/usinggalaxy/h/upload).

* Each file loaded creates one dataset in the history.
* The maximum size limit is 50G (uncompressed).
* Most individual file compression formats are supported, but multi-file archives are not (`.tar`, `.zip`).
* If a `.tar` archive contains multiple datasets, only the first dataset inside the archive will upload.
* If your compressed data does not load correctly, try loading an uncompressed version.

### Get Data: Upload

1. Load by "browsing" for a local file. Best for small [datasets](/learn/managing-datasets). Some servers will support load data that is 2 GB or larger.If you are having problems with this method, try FTP.
2. Load using an **HTTP URL** or **FTP URL**.
3. Load a few lines of plain text.
4. Load using FTP. Either line command or with a desktop client.

### Get Data: EBI-SRA

* Search for your data directly in the tool and use the Galaxy links
* There are a few links, so [which data do I load?](/support/ebi-sra-data-load/)
* Be sure to check your sequence data for correct quality score formats and the metadata "datatype" assignment.
* [How to format fastq data for tools that require .fastqsanger format?](/support/fastqsanger/)
* [Understanding compressed fastq data (fastq.gz)](/support/compressed-fastq/)

### Get Data: Upload tool used with FTP

* Load the data using line command FTP or a client. [More help...](/ftp-upload/)
* Note that the FTP server name is specific to the Galaxy you are working on. This is by default the URL of the server.
    * For the public Galaxy Main instance at [http://usegalaxy.org](http://usegalaxy.org) the FTP server name to use is **usegalaxy.org**.
    * For a default local (with FTP enabled, see next) the FTP server name to use is **localhost:8080**. If the server URL was modified, use that custom URL.
* If on another server, the FTP server name will appear in the **Upload** tool pop-up window (see graphics below). When using a local Galaxy server, be certain to [configure your instance for FTP](/admin/config/upload-via-ftp/) first.
* Use your email and password for the same instance as your credentials to log in and save the data to your account.
* Once the data is loaded (confirm through FTP client), use the **Upload** tool to load the data into a History.
* **`FTPS` was enabled for all transfers to [http://usegalaxy.org](http://usegalaxy.org) on July 19, 2017**. If you are having trouble connecting the first time after this date, verifying the server certificate is required when using an [FTP client](/ftp-upload/#upload-from-client).
* Working at Galaxy EU https://usegalaxy.eu? Read the server-specific FTP help here: https://galaxyproject.eu/ftp/
* More help for FTP is at Galaxy Help https://help.galaxyproject.org. Search with the keyword "ftp". Example post: https://help.galaxyproject.org/t/ftp-help-guides-tutorials-and-troubleshooting/3449


### Upload tool location

![Upload tool location](/images/screenshots/Upload.png "Upload tool location")

### Upload tool option to move FTP datasets into a History

If you DO NOT see any files, **load data using FTP first**, then come back to the _Upload_ tool.

Make sure that the FTP transfer is completed before moving files into a History.

If you already moved a partial or truncated dataset into a History, you will need to FTP upload the data again.


### Upload tips

* **Data [quota](/main/#user-data-and-job-quotas) is at limit**, so _no new data can be loaded_. Disk usage and quotas are reported at _**User → Preferences**_ when logged in.
* **Password protected data** will require a special URL format. Ask the data source. Double check that it is _publicly accessible_.
* Use _**[FTP or FTPS](/ftp-upload/)**_, not _**SFTP**_. Check with local admin if not sure.
* **No HTML content.** The loading error generated may state this. Remove HTML fields from your dataset before loading into Galaxy or omit HTML fields from the query if importing from a data source (such as Biomart).
* Compression types **.gz/.gzip, .bz/.bzip, .bz2/.bzip2, and _single-file_ .tar and .zip are (usually) supported -- but if your tar/zip data does not load -- download the data locally, unpack the archive, and upload the data directly.**
* Only the **first file in any compressed archive** will load as a **[dataset](/learn/managing-datasets)**.
* Data must be **< 50G** (uncompressed) to be successfully uploaded and added as a dataset to a history, from any source.
* **Is the problem the dataset format or the assigned datatype?** Can this be corrected by editing the datatype or converting formats? See [Learn/Managing Datasets](/learn/managing-datasets/) for help or watch the screencast above for a how-to example.
* **Problems in the first step working with your loaded data?** It may not have _uploaded_ completely. If you used an FTP client, the transfer message will indicate if a load was successful or not and can often restart interrupted loads. This makes FTP a great choice for slower connections, even when loading small files.
