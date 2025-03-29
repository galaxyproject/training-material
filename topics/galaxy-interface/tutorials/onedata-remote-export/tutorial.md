---
layout: tutorial_hands_on

title: "Exporting to Onedata remote"
subtopic: cleanup
tags:
  - storage
contributors:
  - lopiola
  - bwalkowi

time_estimation: "5m"
level: Introductory

requirements:
 - type: "internal"
   topic_name: galaxy-interface
   tutorials:
     - download-delete-data
     - onedata-getting-started
     - onedata-remote-import

questions:
- How to export Galaxy datasets to Onedata?
- What configuration is needed for Onedata export?
- What permissions are required for Onedata export?
objectives:
- Configure a writable Onedata Remote File Source.
- Learn how to export datasets to Onedata.
- Understand the requirements for Onedata export.
key_points:
- Export can be done through the "Export datasets" tool.
- Remote File Source must be configured as writable to enable export.
- Write access is required in the Onedata access token.
- Export works with both generic and user-configured Remote File Sources.
---


# Prerequisites

1. This tutorial assumes that you have basic knowledge about Onedata and access 
   to a Onedata ecosystem. If needed, follow 
   [this tutorial]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html %})
   first!
2. To use Onedata as a remote file source for data export, you need the
   **domain** of the **Onezone service** and a suitable **access token**. Here
   is the relevant 
   [guide on how to get them]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html#accessing-onedata-services %}). 
3. The Galaxy server must be properly configured by the admins for the Onedata
   remote and/or Onedata BYOD templates to be available. Here is the 
   [corresponding tutorial]({% link topics/admin/tutorials/onedata-configuration/tutorial.html %}).


# Introduction

Onedata Remote File Source can be configured as writable. In such a case, you
can export your Galaxy datasets to a Onedata Space.


# Configuration

Follow the same steps as in the 
[Importing (uploading) data from Onedata]({% link topics/galaxy-interface/tutorials/onedata-remote-import/tutorial.html#configuration %})
tutorial to configure a Onedata Remote File Source.

Make sure that the Remote File Source is writable. In case of you own,
it's as simple as checking the toggle in the configuration. In case of
generic remotes, it's in the hands of admins to mark the Onedata remote
as writable.


# Exporting datasets

Follow these steps:

1. Find the **Export datasets to remote files source** tool in the toolbox.
2. Choose the dataset to be exported and choose your Onedata Remote File Source
   in the **Directory URI** section.
   ![Export dataset](../../images/onedata-remote-export/dataset-export.png)
3. Adjust the configuration if needed.
4. Run the tool.
5. Check your Onedata account to find the exported dataset at the specified path.


# Troubleshooting

In case of errors:

1. Check the section **Use distributed compute resources** in the
   **Manage Information** section of the **Preferences** menu. It's
   possible that the specific Pulsar endpoint is unable to run the
   job. Try using the default settings.
2. Avoid whitespace characters (spaces, tabs) in the **Directory URI**,
   as they are known to cause problems.
3. Take a look at the 
   [troubleshooting]({% link topics/galaxy-interface/tutorials/onedata-getting-started/tutorial.html#troubleshooting %}) 
   guide.


# Related topics

The Onedata Remote file source can also be used for data import — see the 
[tutorial]({% link topics/galaxy-interface/tutorials/onedata-remote-import/tutorial.html %}).
