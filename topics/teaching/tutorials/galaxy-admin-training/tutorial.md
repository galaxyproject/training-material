---
layout: tutorial_hands_on
redirect_from:
  - /topics/instructors/tutorials/galaxy-admin-training/tutorial

title: Galaxy Admin Training
subtopic: prepare
time_estimation: 60m
questions:
  - How do I organise a Galaxy Admin Training (GAT)
  - What do I need to set up?
  - What should I know during the training?
objectives:
  - Interact with the UseGalaxy.eu admins to arrange for infrastructure
  - Run a great training!
key_points:
  - Infrastructure is available for running GATs for free from UseGalaxy.eu
  - This can be very convenient and easy to use
  - EU provides appropriate DNS entries so you can run trainings with ITs.
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - terraform
contributors:
  - hexylena
---

# Introduction


![GAT logo is the GTN logo over a space background and text reading galaxy admin traiing.](../../images/gat.png)

Setting up and running a Galaxy Admin Training is not a very complicated process thanks to a significant amount of work that has been put into making it easy and quick.

This tutorial has multiple audiences who are all adressed within this one tutorial. We encourage everyone involved in hosting a GAT event to read this, so you are aware of all of the moving parts which are required for such an event.

> <agenda-title></agenda-title>
>
> In this tutorial, we will see:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Planning Your Training

First consider the requirements for your training:

- How many days will the training be?
- What topics do you want to cover?
- How many helpers will you need?
- Will the training be online or in person? Teaching online is slower, so you may not be able to cover as many topic as in person; you may also need more helpers.
- Do you have your own infrastructure, or do you need to request infrastructure?

We recommend [checking out an example schedule from the GAT repository](https://github.com/galaxyproject/admin-training/tree/2020-barcelona) and modifying it based on your needs.

## Requesting Infrastructure

If you do not have a cloud available that can, at minimum:

- setup VMs which are accessible from where your students are
- expose ports 20, 80, 443 (optional), 8080 (optional but nice)

Then consider [contacting UseGalaxy.eu](mailto:galaxy@informatik.uni-freiburg.de?subject=Hosting%20a%20GAT) and let us know you'd like to host a GAT event. Tell us the number of students and dates of the event.

# Provisioning the Infrastructure

Skip to the section on setting up below that is appropriate for you

## [EU admins] Setting up VMs

> <hands-on-title>EU: Setting up VMs</hands-on-title>
>
> 1. Open a PR editing the [count variable](https://github.com/usegalaxy-eu/infrastructure/blob/master/gat.tf#L1) and increasing it to the desired number of VMs.
>
> 2. Merge
>
> 3. Once deployed, clone/pull the repository locally and run `./bin/process-training-output.sh`
>
>    ```
>    $ ./bin/process-training-output.sh
>    ubuntu  gat-0.training.galaxyproject.eu ...
>    ubuntu  gat-1.training.galaxyproject.eu ...
>    ```
>
> 4. Place this in the [GAT Machines](http://gxy.io/gatmachines) spreadsheet in the correct tab.
>
{: .hands_on}

## [Non-EU admins] Setting up the VMs

> <hands-on-title>Non-EU: Setting up VMs</hands-on-title>
>
> 1. If you have your own cloud, setup VMs per student. We strongly recommend:
>
>    - Ubuntu 18.04+
>    - 8Gb RAM
>    - 2 VCPUs
>
>    If you're using terraform, you can take inspiration from [UseGalaxy.eu's terraform configuration](https://github.com/usegalaxy-eu/infrastructure/blob/master/gat.tf)
>
> 2. Place the list of usernames and IPs and passwords in the [GAT Machines](http://gxy.io/gatmachines) spreadsheet in the correct tab.
>
{: .hands_on}

## [Everyone] Bootstrapping the VMs

Once your VMs are running, great! Now you'll need to bootstrap the instances and prepare them for student use.

The GAT team maintains some infrastructure to handle the bootstrapping in [the GAT repository](https://github.com/galaxyproject/admin-training/tree/2020-bcc/bootstrap-instances)

We use a more complicated hosts file in the project, as we have multiple pools of VMs across global regions. This is seen in the `workshop_eu` and `workshop_oz` groups, and the corresponding commands in the Makefile. However, this distinction of regions is not necessary. The most simple hosts file looks like:

```
# Your Machines
[workshop_instances]
192.0.2.1
192.0.2.2
192.0.2.3
...

# Some variables for those machines
[workshop_instances:vars]
ansible_host_key_checking = false
ansible_user = ubuntu
ansible_become = true
ansible_ssh_private_key_file = ~/admintraining.key
set_password = true      # Generate a random password
```

> <hands-on-title>Everyone: Setting up VMs</hands-on-title>
>
> 1. `git clone https://github.com/galaxyproject/admin-training/` and change into that repo
>
> 2. `cd bootstrap-instances`
>
> 3. In the `hosts` file, edit the `[workshop_eu]` section to list all of your IPs or DNS entries.
>
>    You can specify `ansible_user=something` if it uses a different username, and `ansible_password=password` for each machine. E.g.
>
>    A range of DNS entries (`gat-0` to `gat-39`)
>
>    ```
>    [workshop_]
>    gat-[0:39].training.galaxyproject.eu
>    ```
>
>    Some IPs
>
>    ```
>    [workshop_eu]
>    192.0.2.1
>    192.0.2.2
>    192.0.2.3
>    192.0.2.4
>    ```
>
>    Some IPs with additional information
>
>    ```
>    [workshop_eu]
>    192.0.2.1 ansible_password=2121
>    192.0.2.2 ansible_password=1212
>    192.0.2.3 ansible_password=4321
>    192.0.2.4 ansible_password=1235
>
>    [workshop_instances:vars]
>    ansible_user = myuser
>    ```
>
>    Remove the `workshop_oz` from the top of the file under `[workshop_instances:children]`
>
> 4. Setup a python virtualenv and activate it.
>
> 5. Install the requirements (listed in `requirements.txt`)
>
> 6. Run `make all`
>
{: .hands_on}

The last command runs the `playbook.yml`, which in turn does a large number of things:

- Bootstraps python on the machine if it isn't available already (Ubuntu stopped including it for some reason.)
- (Generates if needed) and copies an SSH key to all of the machines to make login easier (this is stored in `id_rsa` and `id_rsa.pub` in the same directory as the makefile.)
- This key is set as the ubuntu user's key, and added to their authorized_keys, to permit them to run ansible easily on that machine (if they don't configure the local connection correctly.)
	- (When pulsar is in use) the pulsar machines are provisioned identically to the ones where Galaxy is setup, so the students can login passwordlessly to their pulsar machine.
- Updates packages (slow)
- Installs basic dependencies (emacs, vim, nano, git, etc.)
- Adds the `gat-cli` script to `/usr/bin/gat`
- Optionally sets a password to the machine
	- If you set `set_password=true` in the hosts file, you can set a password on machines.
	- This is useful when your machines only have an SSH key on them, and no password set for students.
- Reboots the machines

# Testing

Once your VMs are ready, you should test ALL of the lessons you intend to teach. Most of the training should be fine, but sometimes changes in e.g. Galaxy versions or the availability of newer versions of some ansible modules will indicate you should test the training as you plan to teach it, and update the training materials where relevant.

# Starting Your Training

We recommend providing a website [similar to our GitHub repository](https://github.com/galaxyproject/admin-training/tree/2020-barcelona) (or using our repository! Ask us for a branch.) with at minimum the following links:

- **Q&A** pointing to a Google document where students can ask questions.
	- In our experience this is an excellent format for discussion: it allows students to ask as many questions as they have, and responses can be given in real-time right below each question.
	- Additionally images and complex formatting is easy
	- Lastly, it's *anonymous* which many students prefer.
- **Chat** pointing to `https://gxy.io/gatchat` (or your preferred channel)
- **VM List** pointing to `https://gxy.io/gatmachines` (or your own spreadsheet)
- The slides and tutorials for your training.

# During Your Training

During planning for BCC2020 we found that monitoring student progress would be extremely difficult, it was not a situation we had encountered before as this was the first remote GAT during the pandemic. So we developed a small utility, the `gat-cli`, which assists in monitoring students' progress.

```
$ gat
Galaxy Admin Training (gat) tool:

  gat status-ansible     [Admin] Check status of ansible training
  gat status-galaxy      [Admin] Check status of ansible-galaxy training
  gat status-cvmfs       [Admin] Check status of cvmfs training
  gat status-pulsar      [Admin] Check status of pulsar training
```

Each of these commands run a couple checks on their local machine. For example the `status-galaxy` command checks:

- Check that the `postgresql` service is running
- Does a db named `galaxy` exist
- Does http://localhost:8080/api/version respond with *some* content
- Check that the `galaxy` service is running
- Check that the `nginx` service is running
- Does https://localhost:443/api/version respond with *some* content

However, the `gat` command is only available on the VMs, so we've written an ansible task (wrapped in the Makefile) which SSHs into every machine, and runs the gat command for you:

```
$ make check-galaxy-eu
gat-19    postgres ✘ ✘ galaxy(http) ✘ SysD-gxy ✘ SysD-nginx ✘ galaxy(ssl) ✘
gat-22    postgres ✔ ✘ galaxy(http) ✘ SysD-gxy ✘ SysD-nginx ✘ galaxy(ssl) ✘
gat-2     postgres ✔ ✔ galaxy(http) ✔ SysD-gxy ✘ SysD-nginx ✘ galaxy(ssl) ✘
gat-15    postgres ✔ ✘ galaxy(http) ✘ SysD-gxy ✘ SysD-nginx ✘ galaxy(ssl) ✘
gat-24    postgres ✔ ✘ galaxy(http) ✘ SysD-gxy ✘ SysD-nginx ✘ galaxy(ssl) ✘
gat-6     postgres ✘ ✘ galaxy(http) ✘ SysD-gxy ✘ SysD-nginx ✘ galaxy(ssl) ✘
```

You can see a checkmark reported for every step completed by students, giving you a nice overview of how many students have completed each step, and if you're ready to move on. Additionally you know precisely which students you should reach out to check in with, if they aren't progressing. The `-eu` again refers to the specific pool of machines, and if you're using a different hosts file, with different group names, you may need to edit the Makefile accordingly.

There are other `make check-` commands for each of the `gat` status commands. Run `make` to list all of them.

# After Your Training

Once your training is concluded, go through the questions students have asked in the Google Doc, and consider contributing them back to the training materials with Tips and Question boxes covering these student questions. Here is [an example pull request](https://github.com/galaxyproject/training-material/pull/1922) when we did this after BCC2020.

# Conclusion

