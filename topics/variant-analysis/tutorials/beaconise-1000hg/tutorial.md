---
layout: tutorial_hands_on

title: "Working with Beacon V2: A Comprehensive Guide to Creating, Uploading, and Searching for Variants with Beacons"
zenodo_link: ''
subtopic: ''
questions:
- What does the term "beacon" refer to?
- How can MongoDB be employed to establish a beacon tailored for your institution?
- In what manner can variant data and metadata be readied into a format compatible with beacons?
- What are the steps involved in importing data into a beacon seamlessly?
- How does one perform queries on a beacon to retrieve information about variants?
objectives:
- Gain a comprehensive understanding of beacons by unveiling their concept and exploring their significance
- Master the skill of utilizing MongoDB to construct beacons, fostering institutional advancement through effective data management
- Learn the art of preparing data for beacons, including the transformation of variants and metadata into structures compatible with beacon requirements
- Achieve proficiency in seamlessly importing data into beacons through a step-by-step guide, ensuring a smooth and error-free process
- Develop the ability to query beacons for variants, unleashing the formidable capabilities of MongoDB in conducting powerful and precise variant searches
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


# Introduction

Beacon v2 is a data query protocol and API that allows the researcher to seek information about specific genomic variants of biomedical research and clinical applications from the data providers (beacon provider) without 
accessing the original or the whole dataset {% cite Rambla2022 %}. The protocol was developed by the Global Alliance for Genomics and Health (GA4GH) in 2021 as an update for the former Beacon v1 
{% cite GA4GH %}{% cite ELIXIR Beacon %}. The second version of Beacon comes with additional features that provide users with more detailed information about the queried variants. Unlike the previous version, 
which only returned a Yes or No response indicating the presence of variants, Beacon v2 presents comprehensive results about the variants being searched. {% cite Rueda2022 %}.

The Beacon v2 comprises two main parts, namely the framework and the models. The framework specifies how the requests and responses should be formatted, while the models determine the organization of the biological data response 
{cite Rueda2022}.

Beacon MongoDB databases can be created with different security levels for accessing Beacon data through import and query processes. These databases can be classified into three security levels: public, registered, and controlled. 
Alternatively, they can be classified into two security levels in our import and query tools: public and authenticated.

![Beacon Framework and Model](../../images/Beacon_framework_and_model.png "Schematic representation of Beacon v2 API specification")


| Security Level | Description                                         |
|----------------|-----------------------------------------------------|
| Public         | Beacon can be accessed by any request               |
| Registered     | Only known users can access the data                |
| Controlled     | Only specifically granted users can access the data |

The genomic variants data is usually divided into two parts. The Metadata and Genomic variants. The metadata is initially saved in different formats like Excel, CSV, etc., while the Genomic variations are kept in VCF files.
So, before beaconing them, the beacon providers must prepare the data in JSON, also known as Beacon-friendly format (BFF). This can be done by extracting the required information from those files following the beacon schemas. 
The beacon providers can use already developed tools such as in the [B2RI tools](https://github.com/EGA-archive/beacon2-ri-tools) or develop their own tools. 

For this tutorial we will show the steps to beaconing and query the [HG00096](https://s3.console.aws.amazon.com/s3/buckets/1000genomes-dragen?region=us-west-2&bucketType=general&prefix=data/dragen-3.5.7b/hg38_altaware_nohla-cnv-anchored/HG00096/&showversions=false) 
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



# prepare the Data

Before we start working with beacon. First, we need to prepare the data to be in the right format, Beacon protocole desiend to take and search the genomic database and its metadata in BFF format (Json format).
Normaly, The genomic variation are saved in VCF file format and metadata are kept by the resurchers in EXCIL, CSV ,etc. So, before upload them into Beacon database we need to convert them into the proper format
using the availabe Beacon tools in Galaxy. 

we will start with preparing the environmnet by instlling the tools. These scripts we made using LINUX Ubuntu 20.04 system. please modify the scripts for your system if requierd. 

## Setting up the local environment

> <hands-on-title>Preprare the local enviroment</hands-on-title>
>
> 1. Create directory on your local envirment and give it a suitble naming
>
> ```
> mkdir <directory_name>
> ```
>
> 2. Move to the created Directory on the shell
> ```
> cd <directory_name>
> ```
> 3. Create another two directories and name them as inputs, output and scripts
>
> ```
> mkdir inputs
> mkdir outputs
> mkdir scripts
> ```
>
> 4. Crate a conda enviroment and give it a suitable name
>
> ```
> conda create -n <enviroment name>
> ```
> 
> 5. Activate the the enviroment 
> ```
> conda activate <enviroment name>
> ``` 
> 6. add bioconda channel to your channels list
> ```
> conda config --add channels defaults
> conda config --add channels bioconda
> conda config --add channels conda-forge
> conda config --set channel_priority strict
> ``` 
> 7. install the beacon2-ri-tools from conda 
> ```
> conda install beacon2-ri-tools
> ``` 
> with this, we created the directory where we will install the data to, and installed the beacon2-ri-tools wich we will use to preprocess the data.
> 8. Move to the scripts directory and install vcf2json.py, phenopacket.py and galaxy-beacon import tools from Zenodo. 
> ```
> cd scripts
> wget https://zenodo.org/records/10657357/files/phenopacket.py
> wget https://zenodo.org/records/10657357/files/vcf2json.py
> ```
{: .hands_on}



## Convert the Genomic Vareints file into JSON

> <hands-on-title>Convert VCF file into JSON</hands-on-title>
>
> 1. Change the directory to the inputs directory
> ```
> cd path/to/the/inputs
> ```
> 2. Install the genomic varients file from Zenodo using wget tool
> ```
> wget https://zenodo.org/records/10657357/files/HG00096.cnv.vcf
> ``` 
> 3. Run the tool *vcf2json.py* tool
> ```
> python path/to/tools/vcf2json.py -i path/to/inputs/HG00096.cnv.vcf -o path/to/outputs/HG00096.json
> ```
> This will converth the genomic variations VCF file into JSON file
{: .hands_on}



## Preprosses the Metadata Files


> <hands-on-title>JSON Phemopacket Metadata preperation</hands-on-title>
>
> 1. Change the directory to the inputs directory
> ```
> cd path/to/inputs
> ```
> 2. Install the phenopacket metadata from Zenodo using wget tool
> ```
> wget https://zenodo.org/records/10657357/files/igsr-1000-genomes-30x-on-grch38.tsv
> ``` 
> 3. Run the tool *phenopacket.py* with the following parametars
> ```
> python path/to/tools/phenopacket.py -i path/to/inputs/igsr-1000-genomes-30x-on-grch38.tsv -o path/to/outputs/phenopacket.json
> ```
> The tool will extract the informations we need from the tsv file and create a phenopacket json file from them follwing a modefied biosamble schema.
{: .hands_on}




# Create Beacon data discovery protocol using MongoDB

Beacon is created by the Beacon providing instituion to serve and work as the institution prespective and if they want to make their beacon as open access or requers authentication to query the data. 

In this tutorial we will show an example on how to create both of those typs of beacon by creating 2 docker servers for MongoDB one for open access and the other Created from privet docker image.

## Create and open access MongoDB server

> <hands-on-title>Create open access Beacon Database on MongoDB</hands-on-title>
>
> 1. create directory on your local envirment and give it a suitble naming
>
> ```
> mkdir <directory_name>
> ```
>
> 2. Move to the created Directory on the shell
> ```
> cd <directory_name>
> ```
> 3. crate an empty file and name it 'docker-compose.yaml' using any text editor you have.
> ```
> nano docker-compose.yaml
> ``` 
> 4. copy the text below into the 'docker-compose.yaml' file. 
> ```
> version: '3.6'
> services:
> 
>   mongo-client:
>     image: mongo:3.6
>     restart: unless-stopped
>     volumes:
>       - ./beacon/db:/data/db
>     ports:
>       - "27017:27017"
> 
>   mongo-express:
>     image: mongo-express
>     restart: unless-stopped
>     environment:
>       - ME_CONFIG_MONGODB_SERVER=mongo-client
>       - ME_CONFIG_MONGODB_PORT=27017
>       - ME_CONFIG_BASICAUTH_USERNAME=admin
>       - ME_CONFIG_BASICAUTH_PASSWORD=adminpass
>     ports:
>       - "8081:8081"
> ```
> 5. create the path 'beacon/db' in your directory using 'mkdir' tool
> ```
> mkdir beacon
> mkdir beacon/db
> ```
> You can change the name of that bath but you have to change that also from the docker-compose.yaml file. 
> Now we have everything ready for creating the MongoDB server hosted in docker container the only step left is to run docker. 
> 6. Run the tool 'docker-compose' with the following parametars
> ```
> docker-compose up -d
> ```
> 7. Check the created docker containers and test if your dokcer container is runinig
> ```
> docker ps
> ```
> This will give you a massage similaer to this
> ```
> CONTAINER ID   IMAGE           COMMAND                  CREATED      STATUS          PORTS                                           NAMES
> 54aeff806042   mongo:3.6       "docker-entrypoint.s…"   6 weeks ago   Up 6 weeks   0.0.0.0:27017->27017/tcp, :::27017->27017/tcp   beacon_mongo-client_1
> 38a8e4481963   mongo-express   "/sbin/tini -- /dock…"   6 weeks ago   Up 6 weeks   0.0.0.0:8081->8081/tcp, :::8081->8081/tcp       beacon_mongo-express_1
>
> ```
> 8. Test your dokcer to see if it runing by 
> ```
> docker exec -it mongo-client bash
> ```
> this take you into the docker container, you can exit that by prising 'ctrl + d' from your keyboard.  
>
> With this we have created an empty MongoDB server were we can add the beacon database or additional databases to it. 
{: .hands_on}


## Create authentecated beacon MongoDB database

> <hands-on-title>Create authentication requierd Beacon Database on MongoDB</hands-on-title>
>
> 1. create directory on your local envirment and give it a suitble naming
>
> ```
> mkdir <directory_name>
> ```
>
> 2. Move to the created Directory on the shell
> ```
> cd <directory_name>
> ```
> 3. crate an empty file and name it 'docker-compose.yaml' using any text editor you have.
> ```
> nano docker-compose.yaml
> ``` 
> 4. copy the text below into the 'docker-compose.yaml' file. 
> ```
> version: '3.1'
> 
> # networks:
> #   beacon-priv:
> #   idp-priv:
> #   pub:
> 
> services:
> 
>   ###########################################
>   # MongoDB Database
>   ###########################################
> 
>   mongo-client:
>     image: mongo:3.6
>     restart: unless-stopped
>     volumes:
>       - ./beacon/db:/data/db
>     ports:
>       - "27017:27017"
>     environment:
>       MONGO_INITDB_ROOT_USERNAME: root
>       MONGO_INITDB_ROOT_PASSWORD: example
> 
>   mongo-express:
>     image: mongo-express
>     restart: unless-stopped
>     environment:
>       - ME_CONFIG_MONGODB_SERVER=mongo-client
>       - ME_CONFIG_MONGODB_PORT=27017
>       - ME_CONFIG_BASICAUTH_USERNAME=admin
>       - ME_CONFIG_BASICAUTH_PASSWORD=adminpass
>     ports:
>       - "8081:8081"
> ```
> Now we have everything ready for creating the MongoDB server hosted in docker container the only step left is to run docker. 
> 5. run the tool 'docker-compose' with the following parametars
> ```
> docker-compose up -d
> ```
> 6. Check the created docker containers and test if your dokcer container is runinig
> ```
> docker ps
> ```
> This will give you a massage similaer to this
> ```
> 6a6d32af7952   mongo       "/sbin/tini -- /dock…"   8 days ago   Up 8 days       0.0.0.0:27017->27017/tcp, :::27017->27017/tcp   beacon_mongo_1
> 484aebfefcc6   mongo-express   "/sbin/tini -- /dock…"   8 days ago   Up 10 seconds   0.0.0.0:8081->8081/tcp, :::8081->8081/tcp       beacon_mongo-express_1
> ```
>
> There are many ways to create Beacon using MongoDB you can check B2RI MongoDB beacon docker on [github](https://github.com/EGA-archive/beacon2-ri-tools/tree/main) repo for more intormation. 
> In this tutorial we will use the first MongoDB method to create the Beacon protocole. Depend on your requierment you can chose between the two methods to customise your Beacon.   
{: .hands_on}






# Import data into Beacon MongoDB

Now with preparing the data to be in the beacon propper format and with creating the Beacon protocl, we are ready will import the data into Beacon.

The beacon import tool is fexble to import data from institutional galaxy or from local environment,


> <hands-on-title>Import data into Beacon MongoDB</hands-on-title>
>
> 1. Clone to the usegalaxy-eu, [galaxy-beacon-import](https://github.com/usegalaxy-eu/galaxy-beacon-import) GitHub repo to install the Beacon2 import and search tools.
>
> ```
> git clone https://github.com/usegalaxy-eu/galaxy-beacon-import.git
> ```
>
> 2. Move to the created Directory on the shell
> ```
> cd galaxy-beacon-import # change the name if the directory is cloned with a different name
> ```
> 3. Prepare import script
> ```
> pip3 install -r requirements.txt
> ```
> 4. run the Beacon2-import.py tool for the genomic varience and phenopacket json files with the following parametars.
> ```
> python3 beacon2-import.py -i <genomicvarience_JSON_FILE>  -d beacon -c genomecvarience - H < Hostname/IP of the beacon database> -P <Port of the beacon database>
> python3 beacon2-import.py -i <phenopacket_JSON_FILE>  -d beacon -c phenopacket - H < Hostname/IP of the beacon database> -P <Port of the beacon database>
> ```
> Tho tool can work also with a Beacon protocol with authentication, check the *Addvanced Connection to MongoDB* in the tool help section (python beacon2_import.py -h)
{: .hands_on}






# Search Beacon MongoDB

We are done with importing the data into the beacon we created now we will query the database to look for the samples that match our query. 

we are looking to see if there a Delation mutation  in the gene located in chromosome 1 which start at 58278107 and end in 58279217


> <hands-on-title>Import data into Beacon MongoDB</hands-on-title>
>
> 1. run the Beacon2-search.py tool for the genomic varience and phenopacket json files with the following parametars.
> ```
> python3 beacon2-search.py range  -d beacon -c genomecvarience -rn 1 -s 2651463 -e 2653075 -v DEL
> ```
> The srarch function will queiry the beacon database and print out the resutls that matches our quiery specifications in this case it will print something like this. 
> ``` 
>     {
>        "biosampleId": "HG00096",
>        "assemblyId": "GRCh38",
>        "variantInternalId": "chr1:58278107-58279217:EFO:0030069",
>        "variantType": "DEL",
>        "variantId": "EFO:0030069",
>        "start": 58278107,
>        "end": 58279217,
>        "referenceName": "1",
>        "info": {
>            "legacyId": "DRAGEN:LOSS:chr1:58278108-58279217",
>            "cnCount": 0,
>            "cnValue": 0.0937719
>        }
> ```
> Tho tool can work also with a Beacon protocol with authentication, check the *Addvanced Connection to MongoDB* in the tool help section (python beacon2_import.py -h)
{: .hands_on}


# Conclusion

