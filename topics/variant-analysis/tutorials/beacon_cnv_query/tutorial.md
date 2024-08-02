---
layout: tutorial_hands_on

title: "Querying a Beacon Database for Copy Number Variants (CNVs)"
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

In this tutorial, we will utilize the Beacon2 Query tool to explore a Beacon built using data from [hg38_altaware_nohla-cnv-anchored](https://us-west-2.console.aws.amazon.com/s3/buckets/1000genomes-dragen?region=us-west-2&bucketType=general&prefix=data/dragen-3.5.7b/hg38_altaware_nohla-cnv-anchored/&showversions=false) to identify Copy Number Variants (CNVs). Additionally, we will perform queries on the same dataset using various search terms to locate specific data points, demonstrating the tool's versatility and effectiveness in genomic data exploration.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Querying the Entire Database Collection for CNVs Without Specifying Genomic Location

We will perform a simple query to identify reads that match the *"VARIANT STATE"*: `copy number gain`


> <hands-on-title>Query the Beacon MongoDB</hands-on-title>
>
> 1. Use {% tool [Beacon2 CNV](toolshed.g2.bx.psu.edu/repos/iuc/beacon2_cnv/beacon2_cnv/2.1.1+galaxy0) %} to query the Beacon **genomicVariations** collection
>   - *"DATABASE HOST"*: `20.108.51.167`
>   - *"DATABASE PORT"*: `27017`
>   - *"DATABASE"*: `beacon`
>   - *"COLLECTION"*: `genomicVariations`
>   - *"*"VARIANT STATE"*"*: `copy number gain`
>
> The search function will query the Beacon database and display the results that match our specified criteria.
>
>
> > <question-title></question-title>
>    >
>    > What dose variantId "EFO:0030070" means? 
>    >
>    > > ```json
>    > >{'_id': ObjectId('66ab66797f2ec2a3f3484ed3'),
>    > > 'assemblyId': 'GRCh38',
>    > > 'biosampleId': 'HG00103',
>    > > 'definitions': {'Location': {'chromosome': '1',
>    > >                              'end': 16728256,
>    > >                              'start': 16716711}},
>    > > 'id': 'refvar-66ab66797f2ec2a3f3484ed3',
>    > > 'info': {'cnCount': 4,
>    > >          'cnValue': 2.17226,
>    > >          'legacyId': 'DRAGEN:GAIN:chr1:16716712-16728256'},
>    > > 'updated': '2024-08-01T10:42:58.106345',
>    > > 'variantInternalId': 'chr1:16716711-16728256:EFO:0030070',
>    > > 'variantState': {'id': 'EFO:0030070', 'label': 'copy number gain'}}
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



# Querying the Entire Database Collection for CNVs at a Specific Genomic Location

We will perform a targeted query to identify reads that match the *"VARIANT STATE"*: `copy number gain` on **Chromosome** `20`starting at position `29313239` and ending at `29340136`


> <hands-on-title>Query the Beacon MongoDB</hands-on-title>
>
> 1. Use {% tool [Beacon2 CNV](toolshed.g2.bx.psu.edu/repos/iuc/beacon2_cnv/beacon2_cnv/2.1.1+galaxy0) %} to query the Beacon **genomicVariations** collection
>   - *"DATABASE HOST"*: `20.108.51.167`
>   - *"DATABASE PORT"*: `27017`
>   - *"DATABASE"*: `beacon`
>   - *"COLLECTION"*: `genomicVariations`
>   - *"CHROMOSOME"*: `20`
>   - *"START"*: `29313239`
>   - *"END"*: `29340136`
>   - *"*"VARIANT STATE"*"*: `copy number gain`
>
> This query will generate a list of all reads that meet the specified parameters, matching the *"VARIANT STATE"* of `copy number gain` on **Chromosome** `20`, within the defined range from position `29313239` to `29340136`.
{: .hands_on}


In this tutorial, we used our testing Beacon server. To query a specific Beacon, update your credentials in the Galaxy user preferences and enter your server's details by adjusting the "DATABASE HOST" and "DATABASE PORT" fields accordingly.

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