---
layout: tutorial_hands_on

title: "Querying the University of Bradford GDC Beacon Database for Copy Number Variants (CNVs)"
zenodo_link: 'https://zenodo.org/records/10658688'
subtopic: 'data-management'
priority: 1
questions:
- What does the term "Beacon" refer to in the context of genomic databases?
- How can the Beacon2 CNV tool be utilized to verify the presence of specific Copy Number Variants (CNVs)?
objectives:
- Comprehend the concept of Beacon databases and their role in genomic research.
- Demonstrate proficiency in querying Beacon databases to identify and analyze specific Copy Number Variants (CNVs).
follow_up_training:
- type: internal
  topic_name: variant-analysis
  tutorials:
      - beaconise_1000hg
time_estimation: 30M
level: Introductory
key_points:
- Understanding Beacon Databases involves grasping their purpose and significance in genomic research and data sharing.
- Mastering the Beacon2 CNV Tool includes developing skills to effectively query and analyze specific CNVs within the Beacon2 database.
- Efficient Querying for CNVs focuses on accurate data retrieval and precise query formulation using the Beacon2 tool.
- Hands-On Application provides practical experience in querying Beacons, ensuring data integrity, and troubleshooting queries.
- Recognizing the Importance of CNVs highlights the role of CNVs in genetic variation and their relevance in research using Beacon databases.
contributors:
- khaled196
- poterlowicz-lab
---


The concept of a genomics "Beacon" refers to facilitating connection between genomic data suppliers, developers, and researchers interested in acquiring genetic variation data. The Beacon system was intended to be simple: an API that allows users to query genomic data collections for the presence of specified genetic variations and receive a simple "Yes" or "No" response. The term "Beacon" was chosen to represent the goal of illuminating the hitherto opaque world of genetic data sharing through widespread engagement.The Beacon Project, which became one of the initial Global Alliance for Genomics and Health (GA4GH) Driver Projects, was warmly welcomed by both members of the GA4GH developer community and genomic resource providers {% cite Rambla2022 %}.

By 2016, over 35 institutions worldwide had adopted the Beacon protocol across more than 90 genomic datasets. At this point, ELIXIR, the European bioinformatics infrastructure organization, became interested in moving the community initiative toward a defined definition. Their efforts were focused on improving usability in order to encourage wider adoption of the Beacon system within the biomedical genomics industry.

The Beacon v1.0 protocol was formally accepted as a GA4GH standard in 2018 {% cite Fiume2019 %}. Although this version included some new features, such as limited support for structural variant inquiries, quantitative replies, and interaction with external protocols, it mostly stuck to the basic idea of querying variations and obtaining aggregate results. However, it quickly became clear that the procedure required significant extension, particularly to fulfill the objectives of clinical applications in rare disease and cancer genomics. This understanding prompted the start of the Beacon v2 design process, which attempted to reconfigure the protocol to support a wider range of use cases.

Copy number variations (CNVs) are genomic segments with different numbers of copies. These variations, which include both DNA sequence amplifications and deletions, have been found in every domain of life, including bacteria, archaea, plants, and mammals. CNVs contribute significantly to genomic variety and can play an important role in promoting fast adaptive evolution, as well as the development of inherited and somatic human illnesses such as cancer {% cite Lauer2019 %}.

Subsequently, The classification of CNV gains and losses into more specific categories: copy number gain, low-level copy number gain, high-level copy number gain, focal genome amplification, copy number loss, low-level copy number loss, high-level copy number loss, and complete genomic deletion. This detailed categorization enhances the precision in identifying and describing CNVs. More information can be found on the [ELIXIR hCNVs resource page](https://cnvar.org/resources/CNV-annotation-standards/).

In this tutorial, we will use the Beacon2 Query tool to explore a Beacon built from the [GDC Public Access CNVs Data]((https://portal.gdc.cancer.gov/analysis_page?app=Downloads)) to identify Copy Number Variants (CNVs). We will also perform targeted queries on this dataset using various search criteria to pinpoint specific data points. This tutorial will showcase the tool's versatility and effectiveness in exploring and analyzing genomic data.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Querying the Beacon Database to Retrieve and Display CNV Records for Specific Genomic Locations

A common use case for querying Copy Number Variation (CNV) records is when researchers are interested in identifying CNVs within a targeted chromosome or at a specific genomic location within that chromosome.

To conduct this type of query, you need to specify parameters that precisely define the desired genomic region in the Beacon database.

Those parametars are, "CHROMOSOME", "Start", and "End".


> <hands-on-title>Query the Beacon MongoDB</hands-on-title>
>
> 1. Use {% tool [Beacon2 CNV](toolshed.g2.bx.psu.edu/repos/iuc/beacon2_cnv/beacon2_cnv/2.1.1+galaxy0) %} to query the Beacon **genomicVariations** collection
>   - *"DATABASE HOST"*: `20.108.51.167`
>   - *"DATABASE PORT"*: `27017`
>   - *"DATABASE"*: `beacon`
>   - *"COLLECTION"*: `genomicVariations`
>   - *"CHROMOSOME"*: `1`
>   - *"START"*: `1574102`
>   - *"CHROMOSOME"*: `1674102`
>
> The search function will query the Beacon database to retrieve and display results located on chromosome 1, specifically within the genomic range between positions 1,574,102 and 1,674,102.
> After the query, review the output file to determine how many records match these specific criteria.
>
> > <question-title></question-title>
>    >
>    > What dose variantId "EFO:0030070" means? 
>    >
>    > > ```json
>    > >{'_id': ObjectId('66bc9dcb520b66614b58f024'),
>    > > 'assemblyId': 'GRCh38',
>    > > 'biosampleId': 'ENSG00000196433.13',
>    > > 'definitions': {'Location': {'chromosome': 'X',
>    > >                              'end': 1643081,
>    > >                              'start': 1615059}},
>    > > 'diseaseType': 'adnexal and skin appendage neoplasms',
>    > > 'gene': 'ASMT',
>    > > 'id': 'refvar-66bc9dcb520b66614b58f024',
>    > > 'info': {'cnCount': 3, 'legacyId': 'DUP:chrX:1615059-1643081'},
>    > > 'primarySite': 'breast',
>    > > 'updated': '2024-08-14T10:18:40.152888',
>    > > 'variantInternalId': 'X:1615059-1643081:EFO:0030071',
>    > > 'variantState': {'id': 'EFO:0030071', 'label': 'low-level gain'},
>    > > 'variantType': 'DUP'}
>    > > ```
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > EFO:0030069 is a term used to describe the copy number gain. The term was set by the CNV community.
>    > > For more information go to the [CNV annotation formats](https://cnvar.org/resources/CNV-annotation-standards/#cnv-term-use-comparison-in-computational-fileschema-formats)
>    > >
>    > {: .solution}
>    {: .question}
{: .hands_on}



# Querying the Beacon Database Collection for a specific CNVs. 

Suppose we are searching for specific variants in the database—for instance, records related to a low-level gain in a particular genomic location.

Building on the previous example, we will add the "VARIANT STATE" parameter to the query filter to refine our search and retrieve the exact variants we are interested in.


> <hands-on-title>Query the Beacon MongoDB</hands-on-title>
>
> 1. Use {% tool [Beacon2 CNV](toolshed.g2.bx.psu.edu/repos/iuc/beacon2_cnv/beacon2_cnv/2.1.1+galaxy0) %} to query the Beacon **genomicVariations** collection
>   - *"DATABASE HOST"*: `20.108.51.167`
>   - *"DATABASE PORT"*: `27017`
>   - *"DATABASE"*: `beacon`
>   - *"COLLECTION"*: `genomicVariations`
>   - *"CHROMOSOME"*: `X`
>   - *"START"*: `1574102`
>   - *"END"*: `1674102`
>   - *"*"VARIANT STATE"*"*: `low-level gain`
>
> This will print out the low-level gain CNV records located on chromosome 1 within the genomic range between positions 1,574,102 and 1,674,102.
{: .hands_on}


# Querying the Beacon Database Collection for Specific CNVs in a Targeted Gene.

Suppose we are interested in identifying specific Copy Number Variations (CNVs) within a particular gene in the Beacon database. For example, we might be looking for records related to deletions or duplications within a gene known to be associated with a specific disease or condition.

To achieve this, we will refine our query by focusing on a targeted gene and using the relevant parameters—such as "GENE NAME" and "VARIANT TYPE"—to filter the database records. This approach will help us retrieve CNV records that match our specific criteria, allowing us to analyze variations within the targeted gene effectively.


> <hands-on-title>Query the Beacon MongoDB</hands-on-title>
>
> 1. Use {% tool [Beacon2 CNV](toolshed.g2.bx.psu.edu/repos/iuc/beacon2_cnv/beacon2_cnv/2.1.1+galaxy0) %} to query the Beacon **genomicVariations** collection
>   - *"DATABASE HOST"*: `20.108.51.167`
>   - *"DATABASE PORT"*: `27017`
>   - *"DATABASE"*: `beacon`
>   - *"COLLECTION"*: `genomicVariations`
>   - *"VARIANT STATE"*: `low-level gain`
>   - *"GENE NAME"*: `ASMT`
>
> This query will retrieve and print CNV records associated with low-level gain within the ASMT gene from the Beacon database. By specifying the targeted gene name and variant state, the query will focus on CNVs that are particularly relevant to conditions associated with this gene.
{: .hands_on}


# Querying the Beacon Database Collection for Specific CNVs in a Targeted Gene Based on Primary Site or Disease Type.

In genomic research, it is often essential to focus on specific Copy Number Variations (CNVs) that are not only located within a particular gene but also associated with certain clinical features, such as a primary site of a tumor or a specific disease type. For instance, researchers might be interested in identifying CNVs in the BRCA1 gene that are linked to breast cancer or in the TP53 gene associated with various types of cancers.

To conduct such a detailed analysis, we can query the Beacon database by narrowing down the search to CNVs within a targeted gene, while simultaneously filtering the results based on relevant clinical parameters such as the primary site of the disease or the type of disease. This approach allows us to extract highly specific data that can provide deeper insights into the genetic underpinnings of particular conditions or cancer types, facilitating targeted research and potential clinical applications.


> <hands-on-title>Query the Beacon MongoDB</hands-on-title>
>
> 1. Use {% tool [Beacon2 CNV](toolshed.g2.bx.psu.edu/repos/iuc/beacon2_cnv/beacon2_cnv/2.1.1+galaxy0) %} to query the Beacon **genomicVariations** collection
>   - *"DATABASE HOST"*: `20.108.51.167`
>   - *"DATABASE PORT"*: `27017`
>   - *"DATABASE"*: `beacon`
>   - *"COLLECTION"*: `genomicVariations`
>   - *"VARIANT STATE"*: `low-level loss`
>   - *"GENE NAME"*: `BRCA1`
>   - *"PRIMARY SITE"*: `breast`
>   - *"DISEASE TYPE"*: `adnexal and skin appendage neoplasms` *optional*
>
> This query will retrieve and print CNV records associated with a low-level loss within the BRCA1 gene from the Beacon database. By specifying the gene name, variant state, and filtering based on the primary site as "breast," the query will focus on CNVs that are highly relevant to breast cancer. Additionally, you can refine the search further by including the disease type, such as adnexal and skin appendage neoplasms, to narrow down the results to those most pertinent to your research objectives.
{: .hands_on}



In this tutorial, we used our University of Bradford Beacon server. To query a specific Beacon, update your credentials in the Galaxy user preferences and enter your server's details by adjusting the "DATABASE HOST" and "DATABASE PORT" fields accordingly.

> <comment-title>Use Credentials to Access Specific Beacon</comment-title>
> 1. Make sure you are logged in to Galaxy.
> 2. Go to **User > Preferences** in the top menu bar.
> 3. To add beacon database credentials, click on **Manage Information** and fill in the Beacon2 Account empty fields `db_auth_source`, `db_user` and `db_password`.
> 4. Make the changes and click on the **Save** button at the bottom.
{: .comment}

# Conclusion

Now, you have a general knowledge of Beacons and how to query them.

You can apply what you learned in this tutorial to query Beacon made by your institution or other institutions.

We hope you find this tutorial helpful!