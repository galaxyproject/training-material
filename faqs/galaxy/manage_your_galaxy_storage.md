---
title: How do I manage my Galaxy storage?
area: features
layout: faq
box_type: tip
contributors: [dadrasarmin]
---

Now, it is possible to bring your own Storage to Galaxy for computation, storage, and archiving of your results. You can add more storage options to your account by following these steps:

- Click on your Username on top right part of the website and then click on `Preferences`.
- From the middle panel, click on the `Manage Your Galaxy Storage` (previously called `Storage location`).
- Click on the `+ Create` button on top of the page. Here, you get multiple options to connect various storage options to your account.

For all of the possible storage options, you should fill the following fields:

- In the `Name` section, give a name to your storage. This name will be used to choose the storage on Galaxy when you want to select a Storage using `User preferences > Preferred Galaxy Storage`.
- Optionally, you can provide a `Description` for this Storage. This is a note for yourself.


{% include _includes/cyoa-choices.html option1="Onedata Storage" option2="Amazon Web Services S3 Storage" option3="Azure Blob Storage" option4="Google Cloud Storage" option5="Any S3 Compatible Storage" default="Any S3 Compatible Storage" text="Select the Storage you like to add to your Galaxy account." disambiguation="BYOS" %}

<div class="Onedata-Storage" markdown="1">
If you have an account in [Onedata](https://onedata.org/), you can use such an object store as a Storage for your Galaxy datasets; they will be stored in the Onedata space of your choice. The minimal supported Onezone version is 21.02.4. More information on Onedata can be found on [Onedata's website](https://onedata.org/#/home).

There are extensive tutorials for setting up and utilizing of OneData on Galaxy Training Network (GTN). At the moment, we have the following tutorials for Onedata on GTN:
1. [Getting started with Onedata distributed storage]( {% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.md %} )
2. [Onedata user-owned storage]( {% link topics/galaxy-interface/tutorials/onedata-byos/tutorial.md %} )
3. [Setting up a dev Onedata instance]( {% link topics/dev/tutorials/onedata-dev-instance/tutorial.md %} )
4. [Configuring the Onedata connectors (remotes, Object Store, BYOS, BYOD)]( {% link topics/admin/tutorials/onedata-configuration/tutorial.md %} )

In short, you can connect your Galaxy account to an Onedata Storage as follows:

- In the `Onezone domain` field, please fill in the address to your `Onezone` domain. It could be something like "datahub.egi.eu".
- In case you want to disable validation of SSL certificates, you can use `Disable tls certificate validation?` option. However, we strongly recommend you to not use this option unless you know what your are doing.
- Provide name of a space that Galaxy data will be stored on Onedata using `Space Name`. If there is more than one space with the same name, you can explicitly specify which one to select by using the format `<space_name>@<space_id>` (for example `demo@7285220ecc636075ae5759aec7ad65d3cha8f9`).
- If you want to provide a path to store Galaxy data, you can use the `Galaxy root directory` field. If this field is empty, the data will be stored in the space's root directory.
- You should provide an `Access Token` to Galaxy for the Onedata space. Your [access token](https://onedata.org/#/home/documentation/topic/stable/tokens), suitable for REST API access in a Oneprovider service. **Must** allow both read and write data access.
- Click on `Create`.
</div>
<div class="Amazon-Web-Services-S3-Storage" markdown="1">
Amazon's Simple Storage Service (S3) is Amazon's primary cloud storage service. More information on S3 can be found in [Amazon's documentation](https://aws.amazon.com/s3/). You have to create a [bucket](https://docs.aws.amazon.com/AmazonS3/latest/userguide/UsingBucket.html) to use in your AWS web console before using this feature.
- You have to provide an `Access Key ID` to be able to use AWS Storage on Galaxy. A security credential for interacting with AWS services can be created from your [AWS web console](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html). Creating an "Access Key" creates a pair of keys used to identify and authenticate access to your AWS account - the first part of the pair is "Access Key ID" and should be entered here. The second part of your key is the secret part called the "Secret Access Key". Place that in the secure part of this form below.
- Provide the [AWS S3 Bucket](https://docs.aws.amazon.com/AmazonS3/latest/userguide/UsingBucket.html) to store your datasets in the `Bucket` field.
- You should enter the second part of the key you created above, `Access Key ID`, in the `Secret Access Key` section. Read more on access keys on [AWS documentation](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html).
- Click on `Create`.
</div>
<div class="Azure-Blob-Storage" markdown="1">
To setup access to your [Azure Blob Storage](https://learn.microsoft.com/en-us/azure/storage/blobs/storage-blobs-introduction) within the Galaxy, follow the steps:
- Provide the name of your Azure Blob Storage account in the `Container Name` field. More information about container's name could be found on [the Microsoft documentation here](https://learn.microsoft.com/en-us/azure/storage/blobs/storage-blobs-introduction#containers).
- Fill the `Storage Account Name` based on your account. More information is available [on Microsoft website](https://learn.microsoft.com/en-us/azure/storage/common/storage-account-overview).
- Please provide the account access key to your Azur Blob Storage account, using `Account Key` field. This is the documentation on [Managing storage account access keys](https://learn.microsoft.com/en-us/azure/storage/common/storage-account-keys-manage?tabs=azure-portal).
- Click on `Create`.
</div>
<div class="Google-Cloud-Storage" markdown="1">
For the setup you will need to generate [HMAC Keys](https://cloud.google.com/storage/docs/authentication/hmackeys) - these can be linked to your user or a service account. Additionally, you will need to define a [default Google cloud project](https://cloud.google.com/storage/docs/aws-simple-migration#defaultproj) to allow Galaxy to access your Google Cloud Storage via the interfaces described in this FAQs.

- To connect Galaxy to your Google Cloud Storage, you have to generate [HMAC Keys](https://cloud.google.com/storage/docs/authentication/hmackeys). You can use the information after generating the keys to fill the `Access ID` field.
- Use the `Bucket` field to specify the name of [bucket](https://cloud.google.com/storage/docs/buckets) you have created to store your Galaxy data. Documentation for how to create buckets can be found in[ this part of the Google Cloud Storage documentation](https://cloud.google.com/storage/docs/creating-buckets).
- You will receive a `Secret Key` after you generated [HMAC Keys](https://cloud.google.com/storage/docs/authentication/hmackeys). Secret Key should be 40 characters long and look something like the example used the Google documentation - `bGoa+V7g/yqDXvKRqq+JTFn4uQZbPiQJo4pf9RzJ`.
- Click on `Create`.
</div>
<div class="Any-S3-Compatible-Storage" markdown="1">
The APIs used to connect to Amazon's S3 (Simple Storage Service) have become something of an unofficial standard for cloud storage across a variety of vendors and services. Many vendors offer storage APIs compatible with S3. Here, you can configure such service as a Galaxy storage as long as you are able to find the connection details and have the relevant credentials.

- Provide the `Access Key ID`. This is part of your access tokens or access keys that describe the user that is accessing the data. The [Amazon documentation](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html) calls these an "access key ID", the [CloudFlare documentation](https://developers.cloudflare.com/r2/examples/aws/boto3/) describes these as "aws_access_key_id". Internally to Galaxy, we often just call this the "access_key".
- Provide the `Bucket` name. The [bucket](https://docs.aws.amazon.com/AmazonS3/latest/userguide/UsingBucket.html) to store your datasets in. How to setup buckets for your storage will vary from service to service but all S3 compatible storage services should have the concept of a bucket to namespace a grouping of your data together with.
- Using the `S3-Compatible API Endpoint`, you should provide the endpoint URL for your storage service. It is also called "endpoint URL" in some services and the format varies based on the providers. For example, [CloudFlare endpoint URL](https://developers.cloudflare.com/ruleset-engine/rulesets-api/endpoints/) is something like `john.r2.cloudflarestorage.com` and [MinIO](https://min.io/docs/minio/linux/integrations/aws-cli-with-minio.html) endpoint URL is similar to `https://play.min.io:9000`.
- `Secret Access Key` compliment your `Access Key ID` to connect to the S3 compatible storage. The [Amazon documentation](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html) calls these an "secret access key" and the [CloudFlare documentation](https://developers.cloudflare.com/r2/examples/aws/boto3/) describes these as "aws_secret_access_key". Internally to Galaxy, we often just call this the "secret_key".
- Click on `Create`.
</div>

> <tip-title>What can you do after you connected a Storage</tip-title>
>
> You can pick the connected Storage for your analysis as follows:
> 1. Click on your username. Click on `Preferences`.
> 2. Click on `Preferred Galaxy Storage`. Here, you can pick the Storage of your choice. The default option is Galaxy Storage.
>
> Instead of using a default storage location for your account, it is also possible to select it at different levels: per **History**, per **Tool**, and **Workflow**.
> 
> To set a Storage for a specific **History**, you should click on the Galaxy History Storage choice ({% icon galaxy-history-storage-choice %}) icon on the right panel. Then, select the added external storage as the preferred storage location for the **History**. If you execute a **Workflow** in this history, the all results of the workflow will be stored in the external storage (that you selected).
> To verify it, you can click on the Dataset details icon ({% icon details %}) of a job on the right panel and you can see that the user's external storage is used as the "Dataset Storage".
>
> Of course, if instead of a **workflow**, you can run just one **tool** using your connected Storage. To do this, you have to set the Galaxy History Storage choice ({% icon galaxy-history-storage-choice %}) as described above. Then, you can run one (or more) **tool** in this history and the results will be available on your Storage.
{: .tip}