---
layout: tutorial_hands_on

title: "Contributing with GitHub via its interface"
questions:
  - "How can I contribute to an open-source project with GitHub and its interface"
objectives:
  - "Edit a file via GitHub interface"
  - "Create a pull request"
  - "Update a pull request"
time_estimation: "20M"
key_points:
  - "You can't add or push directly to the `master` branch, so you need to create a pull request"
  - "1 pull request = 1 branch"
  - "The pull request is the foundation of the collaborative development of the training material"
contributors:
  - bebatut
---

# Introduction
{:.no_toc}

All the training material which you find on [{{ site.url }}{{ site.baseurl }}]({{ site.url }}{{ site.baseurl }}) is stored on a [GitHub](https://github.com) repository ([{{ site.github_repository }}]({{ site.github_repository }})), a code hosting platform for version control and collaboration. GitHub interface is quite intuitive and simplifies the contributions from anyone.

> ### Agenda
>
> In this tutorial, you will learn how to use GitHub interface to contribute to the training material:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# GitHub

The GitHub repository for the training material is: [{{ site.github_repository }}]({{ site.github_repository }}).

> ### {% icon hands_on %} Hands-on: Familiarization with GitHub
>
> 1. Go on the GitHub repository: [{{ site.github_repository }}]({{ site.github_repository }})
> 
>    ![GitHub interface](../../images/github_interface.png "Interface of the GitHub repository of the training material")
>
> 2. Click on `CONTRIBUTORS.yaml` file
>
>    ![Click on CONTRIBUTOR](../../images/github_contributor.png)
>
> 3. View the file 
> 
>    You should see something like:
>
>    ![Click on CONTRIBUTOR](../../images/github_contributor_interface.png "CONTRIBUTOR file")
>
{: .hands_on}

This `CONTRIBUTORS.yaml` file is where we collect the information (name, email, etc) about the different contributors to display them on our [Hall of Fame]({{ site.baseurl }}/hall-of-fame). You will add your information there. But first you need to sign in to GitHub to be able to change this file.

> ### {% icon hands_on %} Hands-on: Sign in to GitHub
>
> 1. Create a GitHub account (if you do not have one already)
> 
>    ![GitHub signup](../../images/github_sign_up.png)
>
> 2. Sign in (once you have a GitHub account)
>
{: .hands_on}

# Edit a file

You can now modify the `CONTRIBUTORS.yaml` file to add your information there

> ### {% icon hands_on %} Hands-on: Edit a file
>
> 1. Open the `CONTRIBUTORS.yaml` file on GitHub
> 2. Click on {% icon hands_on %} icon (top right of the file)
>
>    ![GitHub edit](../../images/github_edit.png)
>
>    A new page will open: 
>
>    ![GitHub contributor edit page](../../images/github_contributor_edit.png "CONTRIBUTOR file in edit mode")
>
> 4. Modify the `CONTRIBUTORS.yaml` by adding yourself
>
>    You should use your GitHub username and add it followed by `:` (the `:` is important) at the correct position given the alphabetical order.
>
> 5. Scroll down to the bottom of the file
> 6. Fill the **Propose file change** form 
>
>    It can also be named **Commit changes** for the ones with write accesses to the repository
>
>    1. Fill the box "Update CONTRIBUTORS.yaml" with "Add < GitHub username > as contributor" (replace < GitHub username > by your GitHub username)
>    
>       > ### {% icon comment %} Commit messages
>       > This a commit message: a description explaining why a particular change was made. Theses messages capture the history of the changes, so other contributors can understand what have been done and why
>       {: .comment} 
>
>    2. Leave "Add an optional extended description..." empty
>
>    ![Propose file change form](../../images/github_propose_change_form.png)
>
> 7. Click on **Propose file change**
>
>    ![Submit propose file change form](../../images/github_propose_change_form_submit.png)
>
{: .hands_on}

Without realizing it, GitHub let you create your first branch (named here `patch-1`) and add a changement on this branch.

> ### {% icon comment %} Branching
> Branching is the way to work on different versions of a repository at one time. By default your repository has one branch named `master` which is considered to be the definitive branch. When you create a branch off the master branch, you're making a copy, or snapshot, of `master` as it was at that point in time.
>
> By changing a file in this branch, it will diverge from the `master` branch. It will contain data that is only on this new branch.
{: .comment}

# Open a Pull Request

Then the addition of your information in the `CONTRIBUTORS.yaml` file is currently only on your branch `patch-1`. Not on the `master` branch and so not only on the [Hall of Fame]({{ site.baseurl }}/hall-of-fame). You can't add or push directly to the `master` branch, so you need to create what we call a pull request.

The GitHub interface guides you through this process: after clicking on **Propose file change**, a new page opens:

![Pull request form](../../images/github_pr_interface.png "Pull request form")

> ### {% icon hands_on %} Hands-on: Edit a file
>
> 1. Open and read the [CONTRIBUTING.md]({{ site.github_repository }}/blob/master/CONTRIBUTING.md) file
> 1. Come back to the pull request
> 2. Fill in the pull request description
>
>    ![Pull request description](../../images/github_pr_form.png)
>
>    1. Add a title for the Pull Request
>    1. Add a message explaining the changes you made
>
>       This message is a good way to introduce yourself and to explain the message you made. Be kind and descriptive. It helps the reviewers to understand why you did your changes and if it should be intergrated into the `master` branch (and then website).
>
>       > ### {% icon comment %} Pull request messages
>       > GitHub uses [Markdown](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet), a simple Markup language, to render the Pull request messages. You can then add bold test, lists, images, etc.
>       {: .comment} 
>
> 1. Scroll down and check the changes you made
>
>    - In green with `+`: what you added
>    - In red with `-`: what you deleted
>
> 1. Click on **Create pull request**
{: .hands_on}

You have created your first pull request!

Your pull request will be reviewed. There are two possible outcomes:

1. Your pull request is accepted. Congratulations! Your changes will be merged into the master branch of the original repository. The website will be re-built and you will be in the [Hall of Fame]({{ site.baseurl }}/hall-of-fame)
2. Your pull request needs modifications: the reviewers will ask for some changes, possibly because the automatic tests are failing.

# Update a Pull Request

One of the reviewers of your pull request asked you to add your name after your GitHub username in the `CONTRIBUTORS.yaml` file.

> ### {% icon hands_on %} Hands-on: Update a Pull Request
> 
> 1. Go to the list of pull request tab on GitHub
> 2. Click on your pull request
> 
>     You can see here the comments the reviewers left you
>
> 3. Click on **Files changed** tab and see the changes you made
>
>    ![Pull request files changed tab](../../images/github_pr_file_changed.png)
>
> 4. Click on {% icon hands_on %} icon
> 5. Add your name below your GitHub username
>
>    It should look like:
>
>    ```
>    bebatut:
>         name: Bérénice Batut
>    ```
> 
> 5. Navigate to the bottom of the file
> 6. Fill the **Commit changes** form, similarly to the **Propose file change** form before
> 7. Make sure the **Commit directly to the `patch-1` branch** is selected
>
>    ![Commit directly to the `patch-1` branch](../../images/github_pr_commit.png)
> 
> 8. Click on **Commit changes**
>
>    ![Submit propose file change form](../../images/github_pr_commit_changes.png)
> 
>    The pull request should be automatically updated
>
> 9. Check that the new changes are added to the pull request on GitHub
> 
{: .hands_on}

# Conclusion
{:.no_toc}

With this tutorial, you learn how to use GitHub to change a file, create a pull request and then contribute to the training material. What you have learned here can be applied to any file.

> ### {% icon tip %} More about GitHub
> Via the GitHub interface, you can also go further: create file, branch directly, etc.
> To learn that, we recommend you to read the [GitHub guide](https://guides.github.com/activities/hello-world/)
{: .tip}
