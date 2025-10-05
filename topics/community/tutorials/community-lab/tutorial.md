---
layout: tutorial_hands_on
title: Creation of the labs in the different Galaxy instances for your community
level: Introductory
redirect_from:
- /topics/dev/tutorials/community-lab/tutorial
subtopic: sig

questions:
- How can I create a lab for a Galaxy community?
objectives:
- Create labs for Galaxy communities
time_estimation: 1H
key_points:
- The Community lab is a centralised webpage built on the Galaxy framework that enables communities to use specific tools, workflows and tutorials on different Galaxy servers.
tags:
- Community
- SIG
- CoDex
contributions:
  authorship:
    - bebatut
    - paulzierep
    - scorreard
    - BirdmanRidesAgain

---

The **Community lab**, a centralised webpage that enables communities to rapidly aggregate, curate, integrate, display, and launch relevant tools, workflows, and training on different Galaxy servers. This user-friendly interface, built on the Galaxy framework, provides community members with data analysis capacity without requiring programming expertise. Users can run individual tools or create complex workflows, with full provenance tracking to ensure reproducibility, designed specifically for the community research (Nasr et al., 2024).
  *For example, [the microgalaxy lab (Europe)](https://microbiology.usegalaxy.eu).*


The aim is this tutorial is to create the files necessary to display the labs in each Galaxy instance.

You can also use the [Galaxy Labs engine](https://labs.usegalaxy.org.au).

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Create the lab specific files (automatic)

The labs are composed of various files, some that are specific to your community and some that are common between all the labs.
To automatically create the necessary files from a set of templates, a script was generated ("`sources/bin/populate_labs.sh`)

This script will create the following structure and files:
communities/<your community>/lab/
    CONTRIBUTORS
    README.md
    base.yml
    conclusion.html
    intro.html
    usegalaxy.eu.yml
    usegalaxy.fr.yml
    usegalaxy.org.yml

communities/<your community>/lab/sections_templates/
        1_data_import_and_preparation.yml
        2_tools.yml
        3_workflows.yml
        4_tutorials.yml
        5_support_and_help.yml
        6_community.yml

{% include _includes/cyoa-choices.html option1="GitHub_Actions" option2="Command-line" default="GitHub_Actions"
       text="The files can be created in two wyas, using GitHub actions (on the GitHub interface) or using the command line to launch a script." %}

<div class="GitHub_Actions" markdown="1">

> <hands-on-title>Generate the files using Github Actions</hands-on-title>
>
> 1. Go to the [Galaxy Codex repo](https://github.com/galaxyproject/galaxy_codex)
> 2. Go to the file `.github/workflows/populate_labs.yaml`
> 3. On the right, click on the pen to "dit this file"
> 4. Duplicate section called "Populate biodiversity Lab"
>
>    For example:
>    ```
>    - name: Populate biodiversity Lab
>      run: | 
>        bash sources/bin/populate_labs.sh
>      env: 
>        COMMUNITY: biodiversity
>    ```
> 5. Change the name and community name for your community
> 6. Commit changes the changes to a new branch that you name "new lab for <Community-name>"
> 7. On the top panel of the Galaxy codex repository, click on "Actions"
> 8. On the left panel, click on "Populate labs"
> 9. In the table, click on "Run workflow" and select the branch you just created "new lab for <Community-name>"
> 10. Click on "Run workflow"
>
> This will create a Pull Request with the new files in the appropriate community folder. To find the Pull Request, on the top panel of the Galaxy codex repository, click on "Pull requests".
{: .hands_on}

</div>

<div class="Command-line" markdown="1">

> <hands-on-title>Generate the files using the command line</hands-on-title>
>
> 1. Fork and clone the [Galaxy Codex repo](https://github.com/galaxyproject/galaxy_codex)
> 2. Set the `COMMUNITY` variable : `export $COMMUNITY=<your community>`
> 3. Lauch `sources/bin/populate_labs.sh`
>
> This will create a bunch of the new files in the appropriate community folder (communities/<your community>/lab/).
{: .hands_on}

</div>

# Modify the generated files to personalize your community lab

Depending on how you generated the files (GitHub actions or command-line), you should be able to find your community specific lab files.
If you generated the files using the GitHub Actions, this will be in the automatically generated pull request.
If you generated the files using the command line, this will be in the local folder that were generated.

Files to update : 
- "communities/<your community>/lab/CONTRIBUTORS" --> Add the handles of everyone who contributed in the lab
- "communities/<your community>/lab/README.md" --> Change all >your-community> by your community name
- "communities/<your community>/lab/base.yml" --> Change all >your-community> by your community name
- "communities/<your community>/lab/intro.html" --> Include a description of your community.

Files that do not require a manual update : 
- "communities/<your community>/lab/conclusion.html" --> No change required.
- "communities/<your community>/lab/sections/*" --> No change required.

The files in the section folders contain the code for each table visible in the lab.
You can check different labs for inspiration, such as the [microgalaxy lab](https://microbiology.usegalaxy.eu).

# Create community specific sections to personalize your community lab

If you want additional sections, for example, the "Microbial isolates" and "Microbiome" sections in the [microgalaxy lab](https://microbiology.usegalaxy.eu/).

> <hands-on-title>Add community specific sections</hands-on-title>
> 1. Find an appropriate template. You can use other files in your section folder, or browse other labs, such as [microgalaxy lab code](https://github.com/galaxyproject/galaxy_codex/tree/main/communities/microgalaxy/lab), and copy the raw code as a template.
> 2. Create the file in your section folder. Name it with a digit (numerical order) and a descriptive name (i.e. "7_microbial_isolate").
> 3. Copy the raw code in this newly created file.
> 4. Update the code to display what you want.
> 5. Save the file.
> 6. Open `communities/<your community>/lab/base.yml`
> 7. Add the previously created file in the sections (see [microgalaxy base file](https://github.com/galaxyproject/galaxy_codex/blob/main/communities/biodiversity/lab/base.yml) for example)
> 8. Save `communities/<your community>/lab/base.yml`
> 9. Commit the changes, create the pull request (if not done previously)
{: .hands_on}

The Pull Request will be reviewed. Make sure to respond to any feedback.

# Include the labs in the different instances

For the lab to be accessible from the different instances, you need to add files in each instance independently.

For the French instance (https://usegalaxy.fr), all the steps are described in [Issue 297](https://gitlab.com/ifb-elixirfr/usegalaxy-fr/infrastructure/-/issues/297).
You can use the [merge request done for the biodiversity lab](https://gitlab.com/ifb-elixirfr/usegalaxy-fr/infrastructure/-/merge_requests/1302) as a reference.

For the European instance (https://usegalaxy.eu), this tutorial will be updated later.
You can use the [pull request done for the biodiversity lab](https://github.com/usegalaxy-eu/infrastructure-playbook/pull/1629) as a reference.

For the American instance (https://usegalaxy.org), this tutorial will be updated later.
You can use the [pull request done for the biodiversity lab](https://github.com/galaxyproject/usegalaxy-playbook/pull/427) as a reference.

For the Australian instance (https://usegalaxy.org.au), this tutorial will be updated later.
You can use the [pull request done for the biodiversity lab](https://github.com/usegalaxy-au/infrastructure/issues/2786#issuecomment-3330995427) as a reference.

# Conclusion

You now have all the files necessary to display your lab in several instances and they should soon be available to your users.
