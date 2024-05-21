---
layout: tutorial_hands_on

title: "Working with Beacon V2: A Comprehensive Guide to Creating, Uploading, and Searching for Variants with Beacons"
zenodo_link: 'https://zenodo.org/records/10658688'
subtopic: ''
questions:
- What does the term "beacon" refer to?
- How can MongoDB be employed to establish a beacon tailored for your institution?
- In what manner can variant data and metadata be readied into a format compatible with beacons?
- What are the steps involved in importing data into a beacon seamlessly?
- How does one perform queries on a beacon to retrieve information about variants?
objectives:
- Gain a comprehensive understanding of beacons
- Master the skill of utilizing MongoDB to construct beacons
- Learn the art of preparing data for beacons, including the transformation of variants and metadata into structures compatible with beacon requirements
- Achieve proficiency in seamlessly importing data into beacons through a step-by-step guide
- Develop the ability to query beacons for variants
requirements:
- type: external
  title: The Unix Shell
  link: "http://swcarpentry.github.io/shell-novice/"
- type: external
  title: Version Control with Git
  link: "https://swcarpentry.github.io/git-novice/"
time_estimation: 2H
key_points:
- Understanding Beacons involves grasping their definition and concept, recognizing their significance in efficient data management, and exploring diverse use cases for these essential data structures.
- Mastering MongoDB for Beacon Creation entails familiarizing yourself with the MongoDB database system, acquiring the skills to construct beacons for advancing institutional data management,
  and exploring tailored MongoDB features for efficient beacon development.
- Effective Data Preparation for Beacons entails mastering techniques for organizing variant data, ensuring proper formatting of metadata to meet beacon requirements,
  and acknowledging the crucial role of well-prepared data in optimizing beacon functionality.
- Achieve seamless data importation into beacons by following a comprehensive step-by-step guide, prioritizing data integrity throughout the process,
  and acquiring skills in effective troubleshooting and error handling.
- Efficiently querying beacons for variants involves mastering diverse querying methods, leveraging MongoDB's capabilities to enhance search precision,
  and emphasizing the necessity for precise and efficient queries within beacon structures.
contributors:
- Khaled196
- kkamieniecka
- poterlowicz-lab
---


Beacon v2 is a data query protocol and API that allows the researcher to seek information about specific genomic variants of biomedical research and clinical applications from the data providers (beacon provider) without 
accessing the original or the whole dataset {% cite Rambla2022 %}. The protocol was developed by the [Global Alliance for Genomics and Health (GA4GH)](https://www.ga4gh.org/) in 2021 as an update for the former Beacon v1.
The second version of Beacon comes with additional features that provide users with more detailed information about the queried variants. Unlike the previous version, 
which only returned a Yes or No response indicating the presence of variants, Beacon v2 presents comprehensive results about the variants being searched. {% cite Rueda2022 %}.

The Beacon v2 comprises two main parts, the framework and the models. The framework specifies how the requests and responses should be formatted, while the models determine the organization of the biological data response 
{% cite Rueda2022 %}.

![Beacon Framework and Model](../../images/Beacon_framework_and_model.png "Schematic representation of Beacon v2 API specification from B2RI Documentation")

MongoDB is a versatile, document-oriented NoSQL database that provides flexibility and scalability for modern applications, allowing developers to efficiently manage and query unstructured data. 

MongoDB databases are created with different security levels. This controls accessing the Beacon through the import and query processes. There are three security levels for Beacon. Public, registered and controlled, in our import and query tools, we have classified them into public and authenticated.


| Security Level | Description                                         |
|----------------|-----------------------------------------------------|
| Public         | Beacon can be accessed by any request               |
| Registered     | Only known users can access the data                |
| Controlled     | Only specifically granted users can access the data |

The Beacon data is divided into two parts. The Metadata and Genomic variants. The metadata is initially saved in different formats like Excel, CSV, etc., while the Genomic variations are kept in VCF files.
So, before creating a Beacon database for them, the Beacon providers must prepare the data in JSON, also known as Beacon-friendly format (BFF). This can be done by extracting the required information from those files following the beacon schemas. 
The beacon providers can use already developed tools such as in the [B2RI tools](https://github.com/EGA-archive/beacon2-ri-tools) or develop their own tools. 

For this tutorial we will show the steps to create a beacon and query it for the [HG00096](https://s3.console.aws.amazon.com/s3/buckets/1000genomes-dragen?region=us-west-2&bucketType=general&prefix=data/dragen-3.5.7b/hg38_altaware_nohla-cnv-anchored/HG00096/&showversions=false) 
biosample structural variants VCF file obtained from the [1000 Genome project](https://s3.console.aws.amazon.com/s3/buckets/1000genomes-dragen?region=us-west-2&bucketType=general&tab=objects) using scripts created by The University of 
Bradford Computational and Data-Driven Science research team.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Create Beacon protocol using MongoDB

The beacon-providing institution creates a Beacon to serve the institution's needs, and it depends if they want to make their Beacon open access or private.
For this training we will create a provate access bracon that with admin user rigts to modefiy the beacon with an additional account for querying the database only. 

We will use docker and docker-compose for this step. If you don't have it installed, please follow this [documentation](https://docs.docker.com/engine/install/) to install it


> <hands-on-title>Create Beacon Database on MongoDB</hands-on-title>
> 1. Create a directory in your local environment and name it as `mongo-init`
> ```bash
> mkdir mongo-init
> ```
> 2. Change your location to the created directory
> ```bash
> cd <mongo-init
> ```
> 3. Use any text editor you are comfortable with to create a new YAML file and name it `docker-compose.yaml`
> ```bash
> nano docker-compose.yaml
> ```
> 4. Copy the text below into the `docker-compose.yaml` file
> >```yaml
> > version: '3.6'
> > services:
> >
> >  mongo-client:
> >    image: mongo:3.6
> >    restart: unless-stopped
> >    volumes:
> >      - ./mongo/db:/data/db
> >      - ./mongo-init:/docker-entrypoint-initdb.d
> >    ports:
> >      - "27017:27017"
> >    environment:
> >      MONGO_INITDB_ROOT_USERNAME: admin
> >      MONGO_INITDB_ROOT_PASSWORD: adminpass
> >
> >  mongo-express:
> >    image: mongo-express
> >    restart: unless-stopped
> >    environment:
> >      - ME_CONFIG_MONGODB_SERVER=mongo-client
> >      - ME_CONFIG_MONGODB_PORT=27017
> >      - ME_CONFIG_BASICAUTH_USERNAME=admin
> >      - ME_CONFIG_BASICAUTH_PASSWORD=adminpass
> >    ports:
> >      - "8081:8081"
> >
> >  mongo-init:
> >    image: mongo:3.6
> >    restart: "no"
> >    depends_on:
> >      - mongo-client
> >    environment:
> >      - MONGO_INITDB_DATABASE=admin
> >      - MONGO_INITDB_ROOT_USERNAME=admin
> >      - MONGO_INITDB_ROOT_PASSWORD=adminpass
> >    volumes:
> >      - ./mongo-init:/docker-entrypoint-initdb.d
> > ```
> 5. Create the path `mongo/db` in your directory using `$mkdir` tool
> ```bash
> mkdir mongo
> mkdir mongo/db
> ```
> You can change the name of that bath, but you have to change that also from the docker-compose.yaml file. 
> We have everything ready for creating the MongoDB server hosted in the docker container. This will create a Beacon database with admin access. Next we will add users. 
> 6. Create a another directory in your mongo-init directory and name it as `mongo-init`
> ```bash
> mkdir mongo-init
> ```
> 7. Change your location to the created directory
> ```bash
> cd <mongo-init
> ```
> 8. Use any text editor you are comfortable with to create a new JS file and name it `create-user.js`
> ```bash
> nano create-user.js
> ```
> 9. Copy the text below into the `docker-compose.yaml` file
> >```js
> > // create_user.js
> >
> > // Connect to the admin database
> > var adminDB = db.getSiblingDB("admin");
> >
> > // Create a new user with read-only access to all databases
> > adminDB.createUser({
> >  user: "query_user",
> >  pwd: "querypassword",
> >  roles: [
> >    { role: "read", db: "admin" },
> >    { role: "read", db: "beacon" }, // Adjust this for your needs
> >    // Add additional read roles as needed
> >  ]
> >});
> > ```
> This will add a user (user name query_user and password querypassword) account with read-only permision to beacon database. This is important to avoid unwanted modefications to the Beacon database. 
> To know more about MongoDB please read the [MongoDB documentation](https://www.mongodb.com/docs/),
> 10. Run the tool `$docker-compose` with the following parameters
> ```bash
> docker-compose up -d
> ```
> 11. Check the created docker containers and test if your docker container is running
> ```bash
> docker ps
> ```
> This will give you a massage similar to this
> ```
> CONTAINER ID   IMAGE           COMMAND                  CREATED       STATUS                         PORTS                                           NAMES
> 96d7886f9cdb   mongo:3.6       "docker-entrypoint.s…"   5 weeks ago   Up 5 weeks                     27017/tcp                                       mongo-init_mongo-init_1
> 1bb520888fbf   mongo:3.6       "docker-entrypoint.s…"   5 weeks ago   Up 5 weeks                     0.0.0.0:27017->27017/tcp, :::27017->27017/tcp   mongo-init_mongo-client_1
> a148fb3db385   mongo-express   "/sbin/tini -- /dock…"   5 weeks ago   Restarting (1) 9 seconds ago                                                   mongo-init_mongo-express_1
> ```
> 12. Test your docker to see if it is running by using `$docker exec` command
> ```bash
> docker exec -it <mongo-client> bash #If the docker image for the Mongo was installed with a different name, Change the name in <mongo-client>
> ```
> This will take you into the docker container. You can exit that by prising `ctrl + d` from your keyboard.
> This will create an empty MongoDB server where we can add the beacon database or any additional databases. 
{: .hands_on}



# Data Preparation
First, start with uploading and preparing the input data to import them itno beacon database. The datasets were optained from [1000genomes-dragen](https://us-west-2.console.aws.amazon.com/s3/buckets/1000genomes-dragen?region=us-west-2&bucketType=general&tab=objects). 

| Name | Format | Data size (MB) |
| ---- | ----- | ------------ |
| HG00096.cnv.vcf | vcf | HG00096.cnv.vcf KB |
| igsr-1000-genomes-30x-on-grch38.tsv | tsv | 405.7 KB |

## Get data

> <hands-on-title>Data upload</hands-on-title>
>
> 1. For this tutorial, make a new history.
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
> 2. Import the data files from
>    [Zenodo](https://zenodo.org/records/10658688):
>
>    ```
>    https://zenodo.org/records/10658688/files/HG00096.cnv.vcf
>    https://zenodo.org/records/10658688/files/igsr-1000-genomes-30x-on-grch38.tsv
>    ```
>   This will download the strucral varient and the metadata file.
> 
>    In some cases the same dataset can be found in the Galaxy shared data library. 
>    Ask the instructor for more details about this.
>
>    The dat aset can also be downloaded a local storage.  
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md format="fastqsanger.gz" %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 4. Give the data meaningful names and tags to facilitate analysis.
>     
>
>    When uploading data from a link, Galaxy names the files after the link address.
>    It might be useful to change or modify the name to something more meaningful.
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
>
> 5. To track the data in the history, it is recommended to tag the datasets by attaching a meaningful tag '#'
>    to them. The tagging will automatically be attached to any file generated
>    from the original tagged dataset.
>   *e.g.*, `#genomicvariants` for structural varients VCF file and
>   *e.g.*, `#phenopacket` metadata tsv file.
>
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
>
{: .hands_on}


## Convert the Genomic Variants VCF file into JSON
We will preprocess the genomic variant data and convert it from VCF format into JSON format.

> <hands-on-title>Convert the Genomic Variants VCF file into JSON</hands-on-title>
> 1. Run {% tool [ CNV VCF2JSON](toolshed.g2.bx.psu.edu/repos/iuc/cnv_vcf2json/cnv_vcf2json/1.0.4+galaxy0) %} on the structural variants VCF file
>       - {% icon param-files %} *"CNV VCF file"*: `HG00096.cnv.vcf` file
>
>
>
> 2. Inspect the *JSON* output produced by the tool
>
>    > <question-title></question-title>
>    >
>    > In the inpected JSON file, what are the key variables and why are they important for genomic analysis?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. **biosampleId**: This variable represents the unique identifier for the biological sample being analyzed.
>    > >    Understanding the biosample allows researchers to trace the origins of the genetic data and contextualize its
>    > >    relevance within a specific study or experiment.
>    > >
>    > > 2. **assemblyId**: The assemblyId denotes the reference genome assembly against which the genetic variants are aligned. 
>    > >    This is crucial for ensuring consistency and accuracy in genomic analyses, as different assemblies may produce varying results.
>    > >
>    > > 3. **variantInternalId**: This identifier specifies the genomic location and type of the variant, providing a standardized format for 
>    > >    referencing genetic variations. It enables researchers to precisely locate and identify specific genetic alterations within the genome.
>    > >
>    > > 4. **variantType**: The variantType indicates the type of genetic variant being described, such as deletion (DEL), insertion (INS), or substitution (SNP). 
>    > >    Understanding the variant type is essential for interpreting its potential impact on gene function and phenotype.
>    > >
>    > > 5. **variantId**: This variable serves as a unique identifier for the variant, often referencing external databases or ontologies. 
>    > >    It facilitates data integration and interoperability across different genomic datasets and analysis platforms.
>    > >
>    > >
>    > > 6. **start** and **end**: These values specify the genomic coordinates of the variant, delineating its precise location within the reference genome. 
>    > >    Knowing the start and end positions is crucial for accurately defining the boundaries of the variant and assessing its potential functional consequences.
>    > >
>    > > 7. **referenceName**: The referenceName indicates the chromosome or genomic contig on which the variant is located. It provides crucial 
>    > >    contextual information for interpreting the genomic coordinates and understanding the genetic context of the variant.
>    > >
>    > >
>    > > 8. **info**: This section contains additional information related to the variant, such as legacy identifiers, copy number counts (cnCount), 
>    > >    and copy number values (cnValue). These supplementary details can provide insights into the structural and functional 
>    > >    implications of the variant within the genome.
>    > >
>    > >
>    > >    Understanding these key variables is essential for conducting effective genomic analyses, as they provide critical information 
>    > >    about the genetic variants under investigation and enable researchers to interpret their biological significance accurately.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}





## Phenopacket Schema

To help us better understand, diagnose, and treat both common and unusual diseases, the Phenopacket Schema is an open standard for exchanging disease and phenotypic data. Building more comprehensive 
disease models is made possible for physicians, biologists, and researchers studying disease and drugs by using Phenopackets, which connect comprehensive phenotype descriptions with patient, disease, 
and genetic data. 

We are using the Biosamles phenopacket, This is a biological material unit from which the substrate molecules (genomic DNA, RNA, proteins, etc.) are extracted for molecular analyses 
(mass spectrometry, array hybridization, sequencing, etc.). Tissue biopsies, single cells from cultures used for single-cell genome sequencing, and protein fractions from gradient centrifugations are a 
few examples. The same Biosample may be referred to by many instances (e.g., technical replicates) or types of research (e.g., both RNA-seq and genomic array experiments).

![v2.0 phenopacket schema](../../images/phenopacket-schema-v2-overview.png "Overview of v2.0 of the schema optained from Phenopacket Schema Documentation")


> <hands-on-title>Create Phenopacket Schema json file from Metadata</hands-on-title>
> 1. Run {% tool [ CNV Phenopacket](toolshed.g2.bx.psu.edu/repos/iuc/cnv_phenopacket/cnv_phenopacket/1.0.2+galaxy0) %} on the structural variants VCF file
>       - {% icon param-files %} *"Metadata file"*: `igsr-1000-genomes-30x-on-grch38.tsv` file
>
>
>
> 2. Inspect the *JSON* output produced by the tool
>
>    > <question-title></question-title>
>    >
>    > In the inpected JSON file, what are the key variables and why are they important for genomic analysis?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. **id**: The id field serves as a unique identifier for the individual described in the Phenopacket. This 
>    > >    identifier allows researchers to track and reference specific individuals across different datasets 
>    > >    and analyses.
>    > >
>    > > 2. **individualId**: This variable represents the individual's unique identifier within the dataset or 
>    > >    study. It helps maintain data integrity and consistency by ensuring accurate identification of the 
>    > >    individual across various analyses and datasets.
>    > >
>    > > 3. **sex**: The sex field specifies the biological sex of the individual, providing important demographic
>    > >    information for genomic analyses. Understanding the individual's sex is crucial for interpreting 
>    > >    genetic variations and assessing their potential relevance to sex-specific traits or diseases.
>    > >
>    > > 4. **description**: This field provides additional information about the individual, such as their 
>    > >    participation in specific genomic projects or studies. It offers context for the individual's 
>    > >    genomic data and helps researchers understand the sources and potential biases of the data.
>    > >
>    > > 5. **procedure**: The procedure field contains information about the procedures or protocols used to 
>    > >    generate the genomic data associated with the individual. This includes codes or identifiers that 
>    > >    reference specific methodologies or experimental protocols, ensuring transparency and 
>    > >    reproducibility in genomic analyses.
>    > >
>    > >
>    > > 6. **files**: This section contains details about the genomic data files associated with the individual. It 
>    > >    includes information such as file identifiers and attributes like the data format (htsFormat) and 
>    > >    genome assembly (genomeAssembly). Understanding these file attributes is essential for 
>    > >    accessing and interpreting the individual's genomic data effectively.
>    > >
>    > > 7. **htsFormat**: The htsFormat field specifies the format of the genomic data files, such as Variant 
>    > >    Call Format (VCF). Knowing the data format is crucial for selecting appropriate analysis tools 
>    > >    and pipelines and ensuring compatibility with downstream analysis workflows.
>    > >
>    > > 8. **genomeAssembly**: This variable denotes the reference genome assembly used for aligning and
>    > >    interpreting the individual's genomic data. It ensures consistency and accuracy in genomic 
>    > >    analyses by specifying the genomic reference against which genetic variants are annotated and interpreted.
>    > >
>    > >
>    > >    By understanding these key variables in the Phenopacket schema, researchers can effectively interpret and analyze the genomic profile 
>    > >    of the individual described, facilitating insights into genetic variations and their potential implications for health and disease.
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}





# Import data into Beacon MongoDB

Now that the data are in the beacon proper format and with creating the Beacon MongoDB server, we are ready to import the data we have into Beacon.

The Beacon import tool is flexible. The tool imports the data from an institutional galaxy or local environment.

> <hands-on-title>Import data into Beacon MongoDB</hands-on-title>
>
> 1. Clone to the usegalaxy-eu, [galaxy-beacon-import](https://github.com/usegalaxy-eu/galaxy-beacon-import) GitHub repo to install the Beacon2-import.py and the beacon2-search.py tools.
> ```bash
> git clone https://github.com/usegalaxy-eu/galaxy-beacon-import.git
> ```
> 2. Change your location to the cloned directory
> ```bash
> cd galaxy-beacon-import #Change the name if the directory is cloned with a different name
> ```
> 3. Install beacon2-import.py dependant tools using *$pip*
> ```bash
> pip3 install -r requirements.txt
> ```
> 4. Run the Beacon2-import.py tool for the genomic variance and phenopacket JSON files with the following parameters.
> ```bash
> python3 beacon2-import.py -i <genomicvarience_JSON_FILE>  -d beacon -c genomecvarience - H < Hostname/IP of the beacon database> -P <Port of the beacon database>
> python3 beacon2-import.py -i <phenopacket_JSON_FILE>  -d beacon -c phenopacket - H < Hostname/IP of the beacon database> -P <Port of the beacon database>
> ```
> The tool also works with a Beacon protocol with authentication. Check the *Advanced Connection to MongoDB* in the tool help section `python beacon2-import.py -h`
{: .hands_on}


# Search Beacon MongoDB

In the last step, we imported the data into the Beacon. Now, we will query the database to look for the samples that match our query. 

We are looking to see if there is a deletion mutation in the gene **located** in `chromosome 1`, which **starts** at `58278107` and **ends** at `58279217`.


> <hands-on-title>Query the Beacon MongoDB</hands-on-title>
>
> 1. Run the Beacon2-search.py tool with the following parameters
> ```bash
> python3 beacon2-search.py range  -d beacon -c genomecvarience -rn 1 -s 2651463 -e 2653075 -v DEL
> ```
> The srarch function will queiry the beacon database and print out the resutls that matches our quiery specifications. In this case it will print something like this. 
> > ```json
> >    {
> >       "biosampleId": "HG00096",
> >       "assemblyId": "GRCh38",
> >       "variantInternalId": "chr1:58278107-58279217:EFO:0030069",
> >       "variantType": "DEL",
> >       "variantId": "EFO:0030069",
> >       "start": 58278107,
> >       "end": 58279217,
> >       "referenceName": "1",
> >       "info": {
> >           "legacyId": "DRAGEN:LOSS:chr1:58278108-58279217",
> >           "cnCount": 0,
> >           "cnValue": 0.0937719
> >       }
> >     }
> > ```
> The tool also works with a Beacon protocol with authentication. Check the *Advanced Connection to MongoDB* in the tool help section `python beacon2-search.py -h`
> 
> > <question-title></question-title>
>    >
>    > What dose variantId "EFO:0030069" means? 
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > EFO:0030069 is a term used to describe the complete genomic deletion. The term was set by the CNV community.
>    > > For more information go to the [CNV annotation formats](https://cnvar.org/resources/CNV-annotation-standards/#cnv-term-use-comparison-in-computational-fileschema-formats)
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}


# Conclusion

Now, you have a general knowledge of Beacon and MongoDB and how to create your own beacon and upload your data. 

You can apply what you learned in this tutorial to create a beacon query for your institution's genomic variant data.

For more information about how to query Beacon databases, please look into [pymongo documantation](https://pymongo.readthedocs.io/en/stable/tutorial.html). Use the documentation steps to query the phenopacket and search for the metadata for our sample `{"id": "HG00096"}`.

We hope you find this tutorial helpful!