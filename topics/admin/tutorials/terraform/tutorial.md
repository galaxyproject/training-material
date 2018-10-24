---
layout: tutorial_hands_on

title: "Terraform"
zenodo_link: ""
questions:
  - What is Terraform?
  - In which situations is it good/bad?
  - How to use it for managing your VM cluster
objectives:
  - Learn Terraform basics
  - Launch a VM with Terraform
  - Launch and tear down a cluster with Terraform
time_estimation: "60m"
key_points:
  - Terraform lets you develop implement infrastructure-as-code within your organisation
  - It can drastically simplify management of large numbers of VMs
contributors:
  - erasche
---

# Overview
{:.no_toc}

In this tutorial we have briefly cover what Terraform is and how you can leverage it for your needs. This will not make you an expert on Terraform but will give you the tools you need in order to maintain your cloud infrastructure as code.

This will be a very practical training with emphasis on looking at examples from modules and becoming self sufficient. This tutorial uses OpenStack.

> ### Agenda
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Cattle vs. Pets

For those of you with more experience as System Administrators, you may have heard the term "cattle vs. pets."

Pets
:    These servers are managed by hand, with admins SSHing in to make changes directly on the machine, and no log of this. If a "pet" dies, everyone is sad.

Cattle
:    These servers are managed in bulk, you no longer own a beloved pet but hundreds of identical ones. You stop caring if an individual server lives or dies as you have many identical ones ready to take its place, and you can easily replace them.


# What is Terraform?

Terraform is a tool you can use to "sync" your infrastructure (VMs in Clouds, DNS records, etc.) with code you have written in configuration files. This is known as **infrastructure as code**.
Knowing that your infrastructure is exactly what you expect it to be can simplify your operations significantly. You can have confidence that if anything changes, any images crash or are accidentally deleted, that you can immediately re-build your infrastructure.

UseGalaxy.eu used Terraform to migrate easily between clouds, when our old cloud shut down and our new cloud was launched. We simply updated the credentials for connecting to the cloud, and ran Terraform. It recognised that none of our existing resources were there and recreated all of them. This made life incredibly easy for us, we knew that our infrastructure would be exactly correct relative to what it looked like in the previous cloud.

## What can be managed

Terraform supports [many providers](https://www.terraform.io/docs/providers/index.html), allowing you to easily manage resources no matter where they are located. Primarily this consists of resources like virtual machines and DNS records.

If you have VMs you are managing, whether for individual projects or for large scale infrastructure like UseGalaxy.eu, your work or research can be simpler and more reproducible if you are doing it with automation like Terraform.

## What should be managed

Support for certain databases and other software is available, but whether or not this is a good idea depends heavily on your workflow. UseGalaxy.eu found that managing databases and database users was not a good fit with our workflow. We launch VMs with terraform and then provision them with ansible. To then have a second step where we re-connect and provision the database with terraform was an awkward workflow, so we mostly let Ansible manage these resources.

Some groups use Ansible or Bash scripts in order to launch VMs. This can be a better workflow for certain cases. Launching VMs that manage themselves or are automatically scaled by some external process is not a good fit. Terraform expects the resource to be there the next time, with the same ID, with the same values, and will alert you if it isn't.

# Managing a Single VM

We will start small, by managing a single VM in our cloud account. Make sure you have your OpenStack credentials available.

## Keypair

Terraform reads all files with the extension `.tf` in your current directory. Resources can be in a single file, or organised across several different files. We had the best experience by separating out logical groups of resources:

- DNS entries in a single file
- Security groups in separate `.tf` files
- Keypairs in a single file
- Servers/Instances in individual files.

> ### {% icon hands_on %} Hands-on: Setting up Terraform
>
> 1. [Install Terraform](https://www.terraform.io/downloads.html)
>
> 2. Create a new directory and go into it.
>
> 3. Create a `providers.tf` file with the following contents:
>
>    ```hcl
>    provider "openstack" {}
>    ```
>
> 4. Run `terraform init`
>
>    > ### {% icon question %} Question
>    >
>    > What did the output look like?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > If this is not the first time you've run it, you will see Terraform download any missing providers:
>    > >
>    > > ```
>    > > Initializing provider plugins...
>    > > - Checking for available provider plugins on https://releases.hashicorp.com...
>    > > - Downloading plugin for provider "openstack" (1.10.0)...
>    > >
>    > > The following providers do not have any version constraints in configuration,
>    > > so the latest version was installed.
>    > >
>    > > To prevent automatic upgrades to new major versions that may contain breaking
>    > > changes, it is recommended to add version = "..." constraints to the
>    > > corresponding provider blocks in configuration, with the constraint strings
>    > > suggested below.
>    > >
>    > > * provider.openstack: version = "~> 1.10"
>    > >
>    > > Terraform has been successfully initialized!
>    > >
>    > > You may now begin working with Terraform. Try running "terraform plan" to see
>    > > any changes that are required for your infrastructure. All Terraform commands
>    > > should now work.
>    > >
>    > > If you ever set or change modules or backend configuration for Terraform,
>    > > rerun this command to reinitialize your working directory. If you forget, other
>    > > commands will detect it and remind you to do so if necessary.
>    > > ```
>    > >
>    > {: .solution }
>    {: .question}
{:.hands_on}

This is the minimum to get started with Terraform, defining which providers we will use and any parameters for them. We have two options now:

1. We can load the OpenStack credentials as environment variables
2. We can code the OpenStack credentials in the `providers.tf` file

We recommend the first option, as often Terraform plans are [made publicly available](https://github.com/usegalaxy-eu/infrastructure) in order to collaborate on them, and this prevents accidentally committing your credentials to the git repository where.

We will start by managing your SSH keypair for the cloud as this is an easy thing to add.

> ### {% icon tip %} Your Operating System
>
> If you are on Windows and do not know your public key, please skip to the "Adding an Instance" section, as you probably do not have the tools installed to do this, and we cannot document how to do it. Instead you can simply reference the key by name later. We will point out this location.
>
{: .tip}

> ### {% icon hands_on %} Hands-on: Keypairs
>
> 1. Find the public key you will use for connecting to the your new VM. It is usually known as `id_rsa.pub`
>
>    > ### {% icon tip %} No public key
>    > If you can find the private key file (possibly a `cloud.pem` file you downloaded earlier from OpenStack), then you can find the public key by running the command:
>    >
>    > ```shell
>    > $ ssh-keygen -y -f /path/to/key.pem
>    > ```
>    {: .tip}
>
> 2. Create a new file, `main.tf` with the following structure. In `public_key`, write the complete value of your public key that you found in the first step.
>
>    ```hcl
>    resource "openstack_compute_keypair_v2" "my-cloud-key" {
>      name       = "my-key"
>      public_key = "ssh-rsa AAAAB3Nz..."
>    }
>    ```
>
> 3. Run `terraform plan`
>
>    > ### {% icon question %} Question
>    >
>    > What does the output look like? Did it execute successfully?
>    >
>    > > ### {% icon solution %} Solution
>    > >
>    > > If you have not sourced your OpenStack credentials file, you will see something like the following:
>    > >
>    > > ```
>    > > Error: Error running plan: 1 error(s) occurred:
>    > >
>    > > * provider.openstack: One of 'auth_url' or 'cloud' must be specified
>    > > ```
>    > >
>    > > You should source your openstack credentials first, and then re-run `terraform plan`:
>    > >
>    > > ```
>    > > An execution plan has been generated and is shown below.
>    > > Resource actions are indicated with the following symbols:
>    > >   + create
>    > >
>    > > Terraform will perform the following actions:
>    > >
>    > >   + openstack_compute_keypair_v2.my-cloud-key
>    > >       id:          <computed>
>    > >       fingerprint: <computed>
>    > >       name:        "my-key"
>    > >       private_key: <computed>
>    > >       public_key:  "ssh-rsa AAAAB...."
>    > >       region:      <computed>
>    > >
>    > >
>    > > Plan: 1 to add, 0 to change, 0 to destroy.
>    > > ```
>    > >
>    > > If you see this, then everything ran successfully
>    > >
>    > {: .solution }
>    {: .question}
>
{: .hands_on }

Let's look at the output in detail:

```
An execution plan has been generated and is shown below.
Resource actions are indicated with the following symbols:
  + create

Terraform will perform the following actions:

  + openstack_compute_keypair_v2.my-cloud-key
      id:          <computed>
      fingerprint: <computed>
      name:        "my-key"
      private_key: <computed>
      public_key:  "ssh-rsa AAAAB..."
      region:      <computed>


Plan: 1 to add, 0 to change, 0 to destroy.

------------------------------------------------------------------------

Note: You didn't specify an "-out" parameter to save this plan, so Terraform
can't guarantee that exactly these actions will be performed if
"terraform apply" is subsequently run.
```

Terraform informs us first about the different symbols used. Here it tells us that it will `+ create` a resource. Sometimes it will `- delete` or `~ update in-place`. Next it describes in detail what will be done and why. For creating a resource it does not give us much information, we will see more types of changes later.

Lastly it informs us that we did not save our plan. Terraform can maintain a concept of what the remote resource's state looks like between your `terraform plan` and `terraform apply` steps. This is a more advanced feature and will not be covered today.

Now that we've defined a keypair (or cannot, but have an existing one in OpenStack)

> ### {% icon hands_on %} Applying our plan
>
> Now that you've reviewed the plan and everything looks good, we're ready to apply our changes.
>
> 1. Run `terraform apply`. Because we did not re-use our plan from the previous step, terraform will re-run the plan step first, and request your confirmation:
>
>    ```
>    An execution plan has been generated and is shown below.
>    Resource actions are indicated with the following symbols:
>      + create
>
>    Terraform will perform the following actions:
>
>      + openstack_compute_keypair_v2.my-cloud-key
>          id:          <computed>
>          fingerprint: <computed>
>          name:        "my-key"
>          private_key: <computed>
>          public_key:  "" => "ssh-rsa AAAAB3..."
>          region:      <computed>
>
>
>    Plan: 1 to add, 0 to change, 0 to destroy.
>
>    Do you want to perform these actions?
>      Terraform will perform the actions described above.
>      Only 'yes' will be accepted to approve.
>
>      Enter a value:
>    ```
>
> 2. Confirm that everything looks good and enter the value 'yes', and hit enter.
>
>    ```
>    openstack_compute_keypair_v2.my-cloud-key: Creating...
>      fingerprint: "" => "<computed>"
>      name:        "" => "my-key"
>      private_key: "" => "<computed>"
>      public_key:  "" => "ssh-rsa AAAAB3..."
>      region:      "" => "<computed>"
>    openstack_compute_keypair_v2.my-cloud-key: Creation complete after 3s (ID: my-key)
>
>    Apply complete! Resources: 1 added, 0 changed, 0 destroyed.
>    ```
{: .hands_on}

We should now have a keypair in our cloud!


## Adding an Instance

We've now have:

- Terraform installed
- The OpenStack plugin intialized
- A keypair in OpenStack

> ### {% icon tip %} Tip: Doing this training outside of a training event
> If you are doing this training outside of an event, then you will likely need to do some additional steps specific to your Cloud:
>
> 1. Identify a small image flavor that you can use for your testing.
>
> 2. Upload this `raw` format image available from [UseGalaxy.eu](https://usegalaxy.eu/static/vgcn/denbi-centos7-j10-2e08aa4bfa33-master.raw) to your OpenStack
>
> You will need to do both of these things before you can proceed.
{: .tip}

You are now ready to launch an instance!

> ### {% icon hands_on %} Hands-on: Launching an Instance
>
> 1. Open your `main.tf` in your text editor and add the following.
>
>    > ### {% icon warning %} Warning: Correct image/flavor/network/security_group names
>    > The documentation below notes some specific values for the `image_name`, `flavor_name`, `security_groups`, and `network` properties. These *may not be correct* for your training, instead your instructor will provide these values to you.
>    {: .warning-box}
>
>    ```hcl
>    resource "openstack_compute_instance_v2" "test" {
>      name            = "test-vm"
>      image_name      = "denbi-centos7-j10-2e08aa4bfa33-master"
>      flavor_name     = "m1.tiny"
>      key_pair        = "${openstack_compute_keypair_v2.cloud2.name}"
>      security_groups = ["default"]
>
>      network {
>        name = "public"
>      }
>    }
>    ```
>
>    There are some important things to note, we reproduce the above configuration but with comments:
>
>    "openstack_compute_instance_v2" is the 'type' of a resource, 'test' is the name of this specific resource.
>    ```
>    resource "openstack_compute_instance_v2" "test" {
>    ```
>
>    This will become the server name
>    ```
>      name            = "test-vm"
>    ```
>
>    The image which will be launched
>    ```
>      image_name      = "denbi-centos7-j10-2e08aa4bfa33-master"
>    ```
>
>    Which flavor
>    ```
>      flavor_name     = "m1.tiny"
>    ```
>
>    This uses Terraform's knowledge of the
>    `openstack_compute_keypair_v2.my-cloud-key` resource which you
>    previously described. Then it accesses the `.name` attribute.
>    This allow your to update the name at any time, and still have your
>    resource definitions be correct.
>    ```
>      key_pair        = "${openstack_compute_keypair_v2.my-cloud-key.name}"
>    ```
>
>    The security group to apply. This is a comma separated list, but
>    `default` should work for most OpenStack clouds.
>    ```
>      security_groups = ["default"]
>    ```
>
>    A network that is accessible to you.
>    ```
>      network {
>        name = "public"
>      }
>    ```
>
> 2. Run `terraform apply`. Running the `plan` step is not necessary, it is just useful to see what changes will be applied without starting the apply process.
>
>    Note that first terraform "refreshes" the state of the remote resources that it already manages, before checking what changes need to be made.
>
>    ```
>    $ terraform apply
>    openstack_compute_keypair_v2.my-cloud-key: Refreshing state... (ID: my-key)
>
>    An execution plan has been generated and is shown below.
>    Resource actions are indicated with the following symbols:
>      + create
>
>    Terraform will perform the following actions:
>
>      + openstack_compute_instance_v2.test
>          id:                         <computed>
>          access_ip_v4:               <computed>
>          access_ip_v6:               <computed>
>          all_metadata.%:             <computed>
>          availability_zone:          <computed>
>          flavor_id:                  <computed>
>          flavor_name:                "m1.tiny"
>          force_delete:               "false"
>          image_id:                   <computed>
>          image_name:                 "denbi-centos7-j10-2e08aa4bfa33-master"
>          key_pair:                   "my-key"
>          name:                       "test-vm"
>          network.#:                  "1"
>          network.0.access_network:   "false"
>          network.0.fixed_ip_v4:      <computed>
>          network.0.fixed_ip_v6:      <computed>
>          network.0.floating_ip:      <computed>
>          network.0.mac:              <computed>
>          network.0.name:             "public"
>          network.0.port:             <computed>
>          network.0.uuid:             <computed>
>          power_state:                "active"
>          region:                     <computed>
>          security_groups.#:          "1"
>          security_groups.3814588639: "default"
>          stop_before_destroy:        "false"
>
>
>    Plan: 1 to add, 0 to change, 0 to destroy.
>
>    Do you want to perform these actions?
>      Terraform will perform the actions described above.
>      Only 'yes' will be accepted to approve.
>
>      Enter a value:
>    ```
>
> 3. If everything looks good, enter 'yes'.
>
>    ```
>    openstack_compute_instance_v2.test: Creating...
>      access_ip_v4:               "" => "<computed>"
>      access_ip_v6:               "" => "<computed>"
>      all_metadata.%:             "" => "<computed>"
>      availability_zone:          "" => "<computed>"
>      flavor_id:                  "" => "<computed>"
>      flavor_name:                "" => "m1.tiny"
>      force_delete:               "" => "false"
>      image_id:                   "" => "<computed>"
>      image_name:                 "" => "denbi-centos7-j10-2e08aa4bfa33-master"
>      key_pair:                   "" => "my-key"
>      name:                       "" => "test-vm"
>      network.#:                  "" => "1"
>      network.0.access_network:   "" => "false"
>      network.0.fixed_ip_v4:      "" => "<computed>"
>      network.0.fixed_ip_v6:      "" => "<computed>"
>      network.0.floating_ip:      "" => "<computed>"
>      network.0.mac:              "" => "<computed>"
>      network.0.name:             "" => "public"
>      network.0.port:             "" => "<computed>"
>      network.0.uuid:             "" => "<computed>"
>      power_state:                "" => "active"
>      region:                     "" => "<computed>"
>      security_groups.#:          "" => "1"
>      security_groups.3814588639: "" => "default"
>      stop_before_destroy:        "" => "false"
>    openstack_compute_instance_v2.test: Still creating... (10s elapsed)
>    openstack_compute_instance_v2.test: Still creating... (20s elapsed)
>    openstack_compute_instance_v2.test: Creation complete after 30s (ID: 7a2ed5ba-0801-49b5-bf1b-9bf9cec733fa)
>
>    Apply complete! Resources: 1 added, 0 changed, 0 destroyed.
>    ```
{: .hands_on}

You now have a running instance! We do not know the IP address so we cannot login yet. You can obtain that from the OpenStack dashboard, or via the `terraform show` command


> ### {% icon hands_on %} `terraform show` and logging in
>
> 1. Run `terraform show`
>
>    The values in the following output will not match yours.
>
>    ```
>    openstack_compute_instance_v2.test:
>      id = 7a2ed5ba-0801-49b5-bf1b-9bf9cec733fa
>      access_ip_v4 = 192.52.32.231
>      access_ip_v6 =
>      all_metadata.% = 0
>      availability_zone = nova
>      flavor_id = a9a58bcf-29a3-47f4-b4c3-e71e40127833
>      flavor_name = m1.tiny
>      force_delete = false
>      image_id = 922ae172-275b-4e6a-bb9d-c7ace52fc8d4
>      image_name = denbi-centos7-j10-2e08aa4bfa33-master
>      key_pair = my-key
>      name = test-vm
>      network.# = 1
>      network.0.access_network = false
>      network.0.fixed_ip_v4 = 192.52.32.231
>      network.0.fixed_ip_v6 =
>      network.0.floating_ip =
>      network.0.mac = fa:16:3e:e6:cd:6c
>      network.0.name = public
>      network.0.port =
>      network.0.uuid = 6a68b40c-356f-4d4c-95e9-e5c72855fb35
>      power_state = active
>      region = Freiburg
>      security_groups.# = 1
>      security_groups.3814588639 = default
>      stop_before_destroy = false
>    openstack_compute_keypair_v2.my-cloud-key:
>      id = my-key
>      fingerprint = 4a:72:2c:7e:d2:b7:f2:c4:5a:9e:bb:ca:c7:8f:d7:a3
>      name = my-key
>      private_key =
>      public_key = ssh-rsa AAAAB3...
>      region = Freiburg
>    ```
>
> 2. Here we can obtain the IP address. In the above output, it is `192.52.32.231`
>
> 3. Login with `ssh`: `ssh centos@<your-ip-address>`. You should see something like the following:
>
>    ```
>    The authenticity of host '192.52.32.231 (192.52.32.231)' can't be established.
>    ECDSA key fingerprint is SHA256:oEHgotonwT3a47K0kLdsNVd/QqctBjDK7F829J0VnwQ.
>    Are you sure you want to continue connecting (yes/no)? yes
>    Warning: Permanently added '192.52.32.231' (ECDSA) to the list of known hosts.
>          _        _   _ ____ _____
>         | |      | \ | |  _ \_   _| ██      ██
>       __| | ___  |  \| | |_) || |     ██████
>      / _' |/ _ \ | . ' |  _ < | |     ██████
>     | (_| |  __/_| |\  | |_) || |_  ██  ██  ██
>      \__,_|\___(_)_| \_|____/_____| ██  ██  ██
>    ====================================================
>    Hostname..........: test-vm.novalocal
>    Release...........: CentOS Linux release 7.5.1804 (Core)
>    Build.............: de.NBI Generic denbi-centos7-j10-2e08aa4bfa33-master
>    Build Date........: 2018-10-24T14:28:38Z
>    Current user......: centos
>    Load..............: 0.00 0.07 0.05 2/148 1644
>    Uptime............: up 5 minutes
>    ====================================================
>    [centos@test-vm ~]$
>    ```
>
{: .hands_on}

With this, you're done. You've launched a single VM with Terraform. The image we are using comes with:

- docker
- singularity
- conda
- CVMFS connected to Galaxyproject's biological databases

This is the standard image that UseGalaxy.eu uses for all of its compute nodes.

# Managing a Cluster of VMs

When you are done playing around with the individual VM, we're ready to explore launching an entire cluster.

- We will launch a VM to act as our HTCondor central manager
- We will launch a VM to act as an NFS server for our "cluster"
- We will launch two execution nodes to process jobs

We use a very similar setup for our production cloud cluster that powers UseGalaxy.eu's compute and training infrastructure.

## Configuration

We will start by setting up the new nodes and launching them as a test:

> ### {% icon hands_on %} Hands-on: Adding configuration for other nodes
>
> 1. Update the name of the instance you already have, `test` to `name = "central-manager"`
>
>    ```
>    resource "openstack_compute_instance_v2" "test" {
>      name            = "central-manager"
>      ...
>    }
>    ```
>
> 2. Add the NFS server. It is mostly the same configuration you've done before, just with a different `name` and resource name.
>
>    ```
>    resource "openstack_compute_instance_v2" "nfs" {
>      name            = "nfs-server"
>      image_name      = "denbi-centos7-j10-2e08aa4bfa33-master"
>      flavor_name     = "m1.tiny"
>      key_pair        = "${openstack_compute_keypair_v2.my-cloud-key.name}"
>      security_groups = ["default"]
>
>      network {
>        name = "public"
>      }
>    }
>    ```
>
> 3. Lastly, we'll add the execution nodes. Here we will launch two servers at once by using the `count` parameter. We can add `count` in the main portion of a Terraform resource definition and it will create that many copies of that same resource. In the case where we wish to be able to distinguish between them, we can use `${count.index}` in the name or another field to distinguish them.
>
>    ```ini
>    resource "openstack_compute_instance_v2" "exec" {
>      name            = "exec-${count.index}"
>      image_name      = "denbi-centos7-j10-2e08aa4bfa33-master"
>      flavor_name     = "m1.tiny"
>      key_pair        = "${openstack_compute_keypair_v2.my-cloud-key.name}"
>      security_groups = ["default"]
>      count           = 2
>
>      network {
>        name = "public"
>      }
>    }
>    ```
>
> 4. Run `terraform apply`
>
>    You will see something like the following:
>
>    ```
>    openstack_compute_keypair_v2.my-cloud-key: Refreshing state... (ID: my-key)
>    openstack_compute_instance_v2.test: Refreshing state... (ID: 7a2ed5ba-0801-49b5-bf1b-9bf9cec733fa)
>
>    An execution plan has been generated and is shown below.
>    Resource actions are indicated with the following symbols:
>      + create
>      ~ update in-place
>
>    Terraform will perform the following actions:
>
>      + openstack_compute_instance_v2.exec[0]
>          id:                         <computed>
>          access_ip_v4:               <computed>
>          access_ip_v6:               <computed>
>          all_metadata.%:             <computed>
>          availability_zone:          <computed>
>          ...
>
>      ~ openstack_compute_instance_v2.test
>          name:                       "test-vm" => "central-manager"
>
>
>    Plan: 3 to add, 1 to change, 0 to destroy.
>
>    Do you want to perform these actions?
>      Terraform will perform the actions described above.
>      Only 'yes' will be accepted to approve.
>
>      Enter a value:
>    ```
>
>    This time, Terraform is able to modify the `openstack_compute_instance_v2.test` resource in place, it can just update the name of the VM.
>
> 5. Enter `yes` and terraform will make the requested changes
>
{: .hands_on}

## Cloud-init

blah blah blah

## Login

## HTCondor

## Submitting a Test Job

# Tearing Everything Down

## Destroy time provisioners
