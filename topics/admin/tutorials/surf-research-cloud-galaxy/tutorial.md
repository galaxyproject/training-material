---
layout: tutorial_hands_on

title: "Galaxy usage on SURF Research Cloud"
zenodo_link: ""
questions:
  - How do I start a Galaxy instance on SURF Research Cloud?
  - How do I link to external storage from my Galaxy workspace?
  - How do I install tools?
objectives:
  - Understand how to work with a Galaxy on-demand instance in SURF Research Cloud
requirements:
  - type: "none"
    title: Access to the SURF Research Cloud
  - type: "none"
    title: SRAM authentication
  - type: "none"
    title: An SSH key connected to your account
  - type: "none"
    title: Have a basic understanding of how Galaxy works
time_estimation: "30m"
level: Introductory
key_points:
  - With SRC you can start your own Galaxy on-demand instance in a secure environment
  - You can log in and use Galaxy via the SRC with your university credentials
contributions:
  authorship:
  - mirelaminkova
  - hexylena
  - dometto
  funding:
  - surf

priority: 10

edam_ontology:
- topic_0605 # Informatics
- topic_3071 # Data Management

abbreviations:
  SRC: SURF Research Cloud
  GDPR: General Data Protection Regulations
  SRAM: SURF Research Access Management
  CO: Collaborative Organisation

follow_up_training:
  - type: "internal"
    topic_name: admin
    tutorials:
      - surf-research-cloud-pulsar

subtopic: cloud
tags:
 - deploying

---

Using Galaxy via the {SRC} allows researchers to start Galaxy instances on-demand and analyze their data in a secure environment following the {GDPR}. The instance provides secure authentication, where users must have a SURF Research account prior to this tutorial, have set the {SRAM} authentication method, and connect an SSH key to their accounts. In case you are not familiar with {SRC} and need help in setting up your accounts, please follow the instructions on the [SURF Knowledge Base](https://servicedesk.surf.nl/wiki/display/WIKI/SURF+Research+Cloud)

In this training the aim is to focus on:

- Understanding how to start Galaxy on {SRC}
- Attaching external volumes for storage

Galaxy instances can be started and stopped on demand, depending on personal cases and requirements. Inside the SRC members should have access to all publicly available catalog items. If you are not able to create a catalog item, please [contact SURF servicedesk](mailto:servicedesk@surf.nl).

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Prerequisites

This tutorial assumes you are member of a {CO} in {SRAM} that has access to {SRC} and a wallet with budget in SRC with enough sources to create Galaxy and Pulsar catalog items. (For more information please refer to the [SURF Knowledge Base](https://servicedesk.surf.nl/wiki/display/WIKI/Budgets%2C+wallets%2C+contracts).

You must have previous experience working with data inside Galaxy.

# Starting a Galaxy instance inside SRC step-by-step

> <hands-on-title>Access the SRC</hands-on-title>
> 1. The first and most important step is to have access to the SURF Research Cloud.
> 2. You will need to login to the [portal](https://portal.live.surfresearchcloud.nl).
>    ![login screen for SURF research cloud](./images/1.png)
>
> 3. Assuming you are a Dutch researcher, you should expect to see your educational institution listed on the authentication provider screen.
>    ![login with screen, showing an authentication provider being chosen.](./images/2.png)
> 4. Log in with your institute credentials and follow the SRAM authentication on your phone. Inside the environment, you can see if there are any active virtual machines, the options to create a workspace, new storage or request new wallet for your projects.
>    ![SRC dashboard, options like create new workspace, storage, or request a new wallet are available. one workspace is listed as running: galaxy test tim.](./images/3.png)
{: .hands_on}

With this, you're ready to get started launching resources within {SRC}.

## Create an external storage

The purpose of the storage is to attach it to your workspace, where you will save your data. This way, if the machine is deleted, your data will not be lost. We strongly advise anyone using the SRC to create storages and attach them to their machines!


> <hands-on-title>Create a Storage Volume</hands-on-title>
> 1. Locate the "Create new storage" option
>    ![top half of the SRC dashboard showing options like create new workspace, storage, or request a new wallet.](./images/4.png)
> 2. Click on "Create new".
> 3. Select the collaborative organisation you want the storage to be part of.
> 4. Select the desired cloud provider and the size of the volume you want to use.
>    > <tip-title>Unsure of the size required?</tip-title>
>    > If you are unsure about the size of the storage you need, contact your administrator. Consider also consulting [Galaxy from an Administrator's Point of View]({% link topics/admin/tutorials/introduction/slides.html %}#10)
>    {: .tip}
>
>    ![create your storage dashboard showing a selection of providers like surf and aws and azure, and then a selection of storage sizes below ranging from 5GB to 1.5TB](./images/5.png)
>
> 6. Name your storage however you like and press on submit.
>
>    ![form to provide name and preview what will be created.](./images/6.png)
>
> 7. You will be redirected to the main page, where you will be able to track the status of the storage creation.
>
>    ![bottom half of home page showing the Storage tab, and a resource named storage-galaxy in state creating.](./images/7.png)
{: .hands_on}

Once the storage is deployed, you can attach it to any of your machines! 🎉

Be aware, that in order to attach an external storage to an existing machine, you would need to pause the workspace first! (If you are not sure how to do that please refer to [External storage volumes](https://servicedesk.surf.nl/wiki/display/WIKI/External+storage+volumes))

# Create a new workspace

Depending on user’s rights, you can create a new workspace for your collaboration.

> <tip-title>A Workspace?</tip-title>
> Workspaces may also be known as *Virtual Machines/VMs* or *Instances* in other cloud providers.
{: .tip}

> <hands-on-title>Create a Workspace</hands-on-title>
> 1. On the right side under the "Workspaces" click on the "Add" button.
>
>    ![dashboard of workspaces showing a single running workspace titled galaxy test tim, and on the right a button to add a new workspace.](./images/8.png)
>
> 2. Then, you will be redirected to a new page, where you must first choose the collaborative organisation in which you want to create your workspace (in case you are a member of multiple organisations). Once chosen, you a new page will be loaded, where you will have access to all available catalog items.
>
>    ![Create your workspace screen showing 9 catalog items on the first page with pages more. Visible are Galaxy Pulsar node amongst some docker and ubuntu options.](./images/9.png)
>
> 3. Use the magnifying glass on the right side of the panel and search for Galaxy. Two catalog items will appear - the Galaxy instance designed for SURF and a Galaxy Pulsar node that can be connected to the instance. Select "Galaxy at SURF".
>
>    ![similar to previous image, but the items described above](./images/10.png)
>
> 4. SURF Research Cloud allows researchers to host their catalog items on different cloud providers. The Galaxy catalog item is currently supported only on the HPC Cloud and for Ubuntu 22.04. On this page, users can select the number of cores and RAM that they want on their machine. More sizes can be added in the future, and on request. Choose wisely and in case you are not certain contact your administrator!
>
>    ![workspace creation screen, choose the cloud provider is locked to SURF HPC, flavour is locked to 22.04, and the size options are available from 1GB/8CPU to 60C/750GB Hi-mem](./images/11.png)
>
> 5. Select the storage you have created earlier so it is attached to the new workspace.
>
>    ![create your workspace screen, storage only lists one option named storage-galaxy](./images/12.png)
>
> 6. Lastly, before the workspace is deployed, you need to choose for how long the machine will run.
>
>    > <tip-title>Expiration Date</tip-title>
>    > The standard life-time of the VM is 5 days. If you need it for longer, this option can be changed once the machine is running.
>    > Note, that once the machine is expired and deleted it cannot be restored! Plan accordingly and migrate your data in time to prevent data loss!
>    >
>    > This is an incredibly useful feature as it saves you from forgetting to destroy a VM. Especially for GPU nodes it can help you ensure that they disappear after your computation is complete.
>    {:.tip}
>
> 7. In this form, you can also select how to name your workspace, add a description (if you want) and specify the hostname. Then, scroll down, customise your Galaxy name, add your storage name and submit the workspace creation. 
>    ![almost there, some final details reads the heading. an expiration date is set and galaxy at surf is provided in the name. the hostname is set to galaxyatsurf and the description is blank](./images/17.png)
>    Change `my_storage` to the storage you attached earlier. In this case, `galaxy-storage`. Bonus: You can also rename the Galaxy to your liking! 
>    ![customise your Galaxy name and attach your storage. change the name of my_data to the storage you attached earlier, in this case, galaxy-storage](./images/18.png)
{: .hands_on}

With that, your workspace will start being created.

![A galaxy instance named galaxy demo is in status creating.](./images/14.png "Creating the workspace can take some time. You will see on the main page the status of your workspace.")

> <tip-title>Also launching Pulsar?</tip-title>
> While you are waiting, you can create your Galaxy Pulsar node simultaneously to save time.
{: .tip}

Once the workspace is running, press on the down arrow to check any information such as the ID, the owner of the workspace, the collaboration it is part of, when it is created and when it expires and the wallet it is using. Additionally, you can pause the workspace when you are not using it. Pressing on the storage option will display the storage that has been attached to the workspace.

![the details panel of the galaxy demo instance is shown including fields like an id, owner, date created and expiration date, url., etc](./images/15.png)

To access the Galaxy instance, press on the "Access" button. This will redirect you to a new webpage. The {SRAM} login will be displayed again, and after authenticating, you will be redirected to the Galaxy at SURF Research Cloud page!

![A galaxy instance homepage is shown, with a brand indicating that it was deployed on surf research cloud featuring a cute surfer emoji.](./images/16.png)

Congratulations! You have started a Galaxy instance inside the SURF Research Cloud! 🏄‍♀️🎉

If you have administrator rights, you can download tools that can be used by the members of the collaborative organisation. To obtain administrator rights you need to be part of the `src_co_admin` group inside the collaborative organisation. In case you are not sure how to download tools as a Galaxy administrator, please refer to [Installing Tools into Galaxy](https://galaxyproject.org/admin/tools/add-tool-from-toolshed-tutorial/).

## Going Further with Pulsar?

Another important step when starting a local Galaxy instance on {SRC} is connecting it to a Galaxy Pulsar Node, which you can do [in this Pulsar SRC tutorial]({% link topics/admin/tutorials/surf-research-cloud-pulsar/tutorial.md %}). The Pulsar is an useful component that is used to execute your tasks on remote servers/clusters. In this case, the Pulsar is started as a new workspace, allowing users to also run interactive tools.

In SRC, this can also enable you to save costs. By choosing a smaller Galaxy instance, and attaching a large Pulsar node when you need to do computations, you can be sure to optimise your budget expenditures.
