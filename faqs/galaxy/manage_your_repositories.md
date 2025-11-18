---
title: How do I manage my repositories on Galaxy?
area: features
layout: faq
box_type: tip
contributors: [dadrasarmin]
---

Here, we are going to briefly explain how you can Bring-Your-Own-Data to Galaxy or export your dataset, results, or history to 3rd party repositories. In order to add a new repository to your account follow these steps:

- Click on your Username on top right part of the website and then click on `Preferences`.
- From the middle panel, click on the `Manage Your Repositories` (previously called `Manage your remote file sources`).
- Click on the `+ Create` button on top of the page. Here, you get multiple options to connect various repositories to your account.

For all of the possible repositories, you should fill the following fields:

- In the `Name` section, give a name to your repository. This name will be used to choose the repository on Galaxy for importing or exporting datasets.
- Optionally, you can provide a `Description` for this repository. This is a note for yourself.

{% include _includes/cyoa-choices.html option1="Onedata" option2="Amazon Web Services Private Bucket" option3="Amazon Web Services Public Bucket" option4="Azure Blob" option5="Dropbox" option6="eLabFTW" option7="An FTP Server" option8="Export to Google Drive" option9="InvenioRDM" option10="S3 Compatible Storage with Credentials" option11="WebDAV" option12="Zenodo" default="Onedata" text="Select the repository you like to add to your Galaxy account." disambiguation="BYOD" %}

<div class="Onedata" markdown="1">
If you have an [Onedata](https://onedata.org/) account, you can use this repository to import and/or export your data directly from and to Onedata. The minimal supported Onezone version is 21.02.4. More information on Onedata can be found on [Onedata's website](https://onedata.org/#/home).

There are extensive tutorials for setting up and utilizing of OneData on Galaxy Training Network (GTN). At the moment, we have the following tutorials for Onedata on GTN:
1. [Getting started with Onedata distributed storage]( {% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.md %} )
2. [Importing (uploading) data from Onedata]( {% link topics/galaxy-interface/tutorials/onedata-remote-import/tutorial.md %} )
3. [Exporting to Onedata remote]( {% link topics/galaxy-interface/tutorials/onedata-remote-export/tutorial.md %} )
4. [Setting up a dev Onedata instance]( {% link topics/dev/tutorials/onedata-dev-instance/tutorial.md %} )
5. [Configuring the Onedata connectors (remotes, Object Store, BYOS, BYOD)]( {% link topics/admin/tutorials/onedata-configuration/tutorial.md %} )

In short, you can connect your Galaxy account to an Onedata repository as follows:

- In the `Onezone domain` field, please fill in the address to your `Onezone` domain. It could be something like "datahub.egi.eu".
- Using the `Writable?` option you can decide whether to grant access to Galaxy to export (write) to your Onedata or not.
- You should provide an `Access Token` to Galaxy so it can read (import) and write (export) data to your OneData. Read more on [access tokens here](https://onedata.org/#/home/documentation/21.02/user-guide/tokens.html). You can limit the access to read-only data access, unless you wish to export data to your repository (write permissions are needed then).
- In case you want to disable validation of SSL certificates, you can use `Disable tls certificate validation?` option. However, we strongly recommend you to not use this option unless you know what your are doing.
- Click on `Create`.
</div>
<div class="Amazon-Web-Services-Private-Bucket" markdown="1">
To connect an AWS private bucket to your Galaxy account, you need to submit the following information on the form:
- First, read the [Manage access keys for IAM (Identity and Access Management) users](https://docs.aws.amazon.com/IAM/latest/UserGuide/id_credentials_access-keys.html) documentation of AWS. Also, you should be familiar with Buckets ([Buckets overview](https://docs.aws.amazon.com/AmazonS3/latest/userguide/UsingBucket.html)). 
- Please fill in the `Access Key ID` (something like `AKIAIOSFODNN7EXAMPLE`) and `Secret Access Key` (similar to `wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY`) in the corresponding fields on the Galaxy interface.
- Please enter the URL to your Bucket (for example, `https://amzn-s3-demo-bucket.s3.us-west-2.amazonaws.com`) in the `Bucket` section.
- Click on `Create`.
</div>
<div class="Amazon-Web-Services-Public-Bucket" markdown="1">
To connect anonymously to an AWS public bucket using your Galaxy account, you need to enter the Bucket address in the `Bucket` section. For more information about AWS Bucket, please read [AWS documentaion](https://docs.aws.amazon.com/AmazonS3/latest/userguide/UsingBucket.html). Click on `Create`.
</div>
<div class="Azure-Blob" markdown="1">
To setup access to your [Azure Blob Storage](https://learn.microsoft.com/en-us/azure/storage/blobs/storage-blobs-introduction) within the Galaxy, follow the steps:
- Provide the name of your Azure Blob Storage account in the `Container Name` field. More information about container's name could be found on [the Microsoft documentation here](https://learn.microsoft.com/en-us/azure/storage/blobs/storage-blobs-introduction#containers).
- Fill the `Storage Account Name` based on your account. More information is available [on the Microsoft website](https://learn.microsoft.com/en-us/azure/storage/common/storage-account-overview).
- Using the `Hierarchical?` option you can determine whether your storage is hierarchical or not. More information on Data Lake Storage namespaces can be found in the [Azure Blob Storage documentation](https://learn.microsoft.com/en-us/azure/storage/blobs/data-lake-storage-namespace).
- Please provide the account access key to your Azur Blob Storage account, using `Account Key` field. This is the documentation on [Managing storage account access keys](https://learn.microsoft.com/en-us/azure/storage/common/storage-account-keys-manage?tabs=azure-portal).
- If you want to be able to export data to your Azure Blob Storage container, please set `Writable?` option to "Yes".
- Click on `Create`.
</div>
<div class="Dropbox" markdown="1">
- We recommend to first login to your Dropbox account.
- On the Galaxy website, click on the `Create` button of the `Dropbox` section. You will be redirected to the Dropbox website for authentication.
- You have to login there and grant access for the Galaxy.
- Click on `Create`.
</div>
<div class="eLabFTW" markdown="1">
[eLabFTW](https://www.elabftw.net/) is a free and open source electronic lab notebook from [Deltablot](https://www.deltablot.com/). Each lab can either host their own installation or go for Deltablot's hosted solution. Using Galaxy, you can connect to an eLabFTW instance of your choice.
- Provide a URL with the protocol (http or https) and the domain name in the `eLabFTW instance endpoint (e.g. https://demo.elabftw.net)` field.
- If you want to let Galaxy to export data to your eLabFTW, please set the `Allow Galaxy to export data to eLabFTW?` to "Yes" to grant required access to Galaxy. **Keep in mind that your API key must have matching permissions.**
- You should provide an `API Key` to your eLabFTW as well. To do so, navigate to the Settings page on your eLabFTW server and go to the API Keys tab to generate a new key. Choose "Read/Write" permissions to enable both importing and exporting data. "Read Only" API keys still work for importing data to Galaxy, but they will cause Galaxy to error out when exporting data to eLabFTW. You will receive a string (similar to `2-50dd721027f56a2e119b3bdbf64f4b8518b3f82b97e7876d56dad74109c8be73d8919b88097d3c9eb8952`) and you should enter this in the `API Key` field of Galaxy interface.
- Click on `Create`.
</div>
<div class="An-FTP-Server" markdown="1">
You can setup connections to FTP and FTPS servers to import and export files as follows:
- Provide the address to your FTP server using the `FTP Host` field.
- If you want to login with a specific user, provide the username in the `FTP User` field. Leave this blank to connect to the server anonymously (if allowed by the server).
- If you want to export data to this FTP, you should set the `Writable?` option to "Yes".
- Please specify the port that Galaxy should use to connect to your FTP server using the `FTP Port` field.
- In the `FTP Password` field provide the password to connect to the FTP server. Leave this blank to connect to the server anonymously (if allowed by the server).
- Click on `Create`.
</div>
<div class="Export-to-Google-Drive" markdown="1">
- We recommend to login to your Google account first.
- On the Galaxy website, click on `Select` button of `Export to Google Drive`. You will be redirected to the Google.
- Pick the account that you want to connect to Galaxy for import and export. Grant the required permissions.
- You will be back on the Galaxy portal and you can access your Google Drive for import and export (depending on your how you set up your accuont).
- Click on `Create`.
</div>
<div class="InvenioRDM" markdown="1">
[InvenioRDM](https://inveniosoftware.org/) is a research data management platform that allows you to store, share, and publish research data. You can connect to an InvenioRDM instance of your choice by following these steps:
- Please fill the address to your InvenioRDM in the following field: `InvenioRDM instance endpoint` (for example, https://inveniordm.web.cern.ch/). This should include the protocol (http or https).
- Use the `Allow Galaxy to export data to InvenioRDM?` option to give permission to Galaxy to export data to your repository or not.
- Click on `Create`.
</div>
<div class="S3-Compatible-Storage-with-Credentials" markdown="1">
- You should fill `Publication Name` with a name as the "creator" metadata of the records. This could be a person or an organization. You can later modify this. If left blank, an anonymous user will be used as the creator.
- You should also enter your `Personal Access Token`. You can get this information in your InvenioRDM instance. Navigate to Account Settings. Then, go to Applications to generate a new token. This will allow Galaxy to display your draft records and upload files to them.
- Click on `Create`.
</div>
<div class="WebDAV" markdown="1">
Using WebDAV you can connect various services that supports WebDAV protocol such as OwnCloud and NextCloud among others. The configuration of WebDAV is slightly variable from service to service but the general principles apply everywhere.
  
- Provide the server address to this repository in the `Server Domain` field.
- In the `WebDAV server Path`, you have to provide the path on this server to WebDAV.
- In the `Username` field, you should write the username you use to login to this server.
- You can grant write access for this repository using the `Writable?` (set to `Yes`) and therefore make it possible to export datasets, or histories to your connected repository.
- Click on `Create`.

> As an example, if I want to connect my nextCloud repository to my Galaxy account, I should login to my nextCloud server and find the information from `File settings` (bottom left of the page) under the WebDAV section to fill this template. It could be something like: `https://server_address.com/remote.php/dav/files/username_or_text`. Here, the `Server Domain` is `https://server_address.com` and `WebDAV server Path` is `remote.php/dav/files/username_or_text`.

In some cases, you may need to activate some features on your ownCloud or nextCloud to allow this integration. For example, some nextCloud servers require the user to use "App Passwords". This can be done using the `Settings > Security > Devices & sessions > Create new app password`.
</div>
<div class="Zenodo" markdown="1">
[Zenodo](https://zenodo.org/) is an open-access repository for research data, software, publications, and other digital artifacts. It is developed and maintained by [CERN](https://home.cern/) and funded by the European Commission as part of the [OpenAIRE](https://www.openaire.eu/) project. Zenodo provides a free platform for researchers to share and preserve their work, ensuring long-term access and reproducibility. Zenodo is widely used by researchers, institutions, and organizations to share scientific knowledge and comply with open-access mandates from funding agencies.
- Using the `Allow Galaxy to export data to Zenodo?`, you can decide whether you like to give write access to Galaxy or not. Set it to "Yes" if you want to export data from Galaxy to Zenodo, set it to "No" if you only need to import data from Zenodo to Galaxy.
- Provide a name for the "creator" metadata of your records on Zenodo using the `Publication Name` field.  You can always change this value later by editing the records in Zenodo. If left blank, an anonymous user will be used as the creator.
- You have to provide a `Personal Access Token` from your Zenodo account to Galaxy. To do so, you need to log into your account. Then, visit this site: https://zenodo.org/account/settings/applications/. Alternatively, you can click on your username on top right and then click on "Applications". Here, you need to create a "Personal Access Token". This will allow Galaxy to display your draft records and upload files to them. If you enabled the option to export data from Galaxy to Zenodo, make sure to enable the **deposit:write** scope when creating the token.
- Click on `Create`.
</div>

> <tip-title>What can you do after you connected a repository</tip-title>
>
> #### Importing data to your Galaxy account
>
> When you connect a repository to your Galaxy account, you can use it to import data to Galaxy. To do so, you can click on the `Upload` Icon on the left panel. In the poped up window, you can click on `Choose from repository` to select a repository that you have added to your account. Navigate to a file that you want to upload to your Galaxy account, check the box of the file, and click on `Select`. You can determine the format of the file, give it a name, and then click on `Start` to upload the file to your Galaxy account.
> 
> #### Exporting histories, datasets, and results to connected repositories
>
> If you have given Galaxy the permission to write to your repository, you can export your histories, datasets and reulsts in the history to that repository.
>
> ##### Histories
> 
> If you want to export a history, you should click on the History Options icon ({% icon galaxy-history-options %}) on the right panel. Then, you can click on `Export History to File`. Next, you can click on `to repository` on the middle panel. If you click on the `Click to select directory`, there will be a pop up window. Here, you can pick a repository that you have added to your account and when you are in that repository, click on `Select`. You can give a `Name` to your exported history, so you can find it easier in your connected repository. Finally, click on `Export` to write the history to your repository. Similarly, you can use `to RDM repository` or `to Zenodo` instead of the `to repository` option in the middle panel to export your history to connected RDM repositories or Zenodo.
>
> To have more options on exporting your history, you can click on `Show advanced export options` on top of the middle panel. This provides further control over the format and datasets that will be included in your exported history.
>
> ##### Datasets
>
> If you are interested to export a single dataset or results to a connected repository, you can use a tool called {% tool [Export datasets](export_remote) %}.
> 1. Select the desired option from `What would you like to export?`.
> 2. Using the `Directory URI` option, you can `Select` a connected repository. You can also give it a directory name here.
> 3. We recommend to export the metadata with your datasets and results using the `Include metadata files in export?`.
{: .tip}