---
layout: tutorial_hands_on

title: "Contributing with GitHub via command-line"
questions:
  - "How can I contribute to an open-source project with GitHub?"
  - "What is the GitHub flow?"
objectives:
  - "Fork a repository on GitHub"
  - "Clone a remote repository locally"
  - "Create a branch"
  - "Commit changes"
  - "Push changes to a remote repository"
  - "Create a pull request"
  - "Update a pull request"
time_estimation: "30m"
subtopic: contribute
key_points:
  - "Work on your fork"
  - "1 branch = 1 feature"
  - "1 pull request = 1 branch"
  - "Be careful of the branch you're working on"
  - "Check the automated tests"
contributors:
  - cmonjeau
  - bebatut
---

# Introduction


Most of the content is written in [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/) with some metadata (or variables) found in [YAML](http://yaml.org/) files. Everything is stored on a [GitHub](https://github.com) repository: [{{ site.github_repository }}]({{ site.github_repository }}).

The process of development of new content is open and transparent, using git and following the [GitHub flow](https://guides.github.com/introduction/flow/):

![Open source development process](../../images/oss_development_process.png "Open source development process")

1. Create a fork
1. Clone your fork of this repository to create a local copy on your computer
1. Create a new branch in your local copy for each significant change
2. Commit the changes in that branch
1. Push that branch to your fork on GitHub
1. Submit a pull request from that branch to the [original repository]({{ site.github_repository }})
1. If you receive feedback, make changes in your local clone and push them to your branch on GitHub: the pull request will update automatically
1. Pull requests will be merged by the training team members after at least one other person has reviewed the Pull request and approved it.

> <agenda-title></agenda-title>
>
> In this tutorial, you will learn how to contribute to the GitHub repository:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Create a fork of this repository on GitHub

A fork is a copy of a repository. [Forking](https://help.github.com/articles/fork-a-repo/) a repository allows you to freely experiment with changes without affecting the original project:

![Explanation of the forking](../../images/PR_schema_01.png "Creation of a fork"){: width="900px"}

Forking a repository is a simple two-step process:

> <hands-on-title>Fork the repository</hands-on-title>
>
> 1. Go on the GitHub repository: [{{ site.github_repository }}]({{ site.github_repository }})
> 2. Click on **Fork** (top-right corner of the page)
>
>    ![Fork](../../images/PR_fork.jpg)
>
> > <comment-title></comment-title>
> > When you click the **Fork** button GitHub will show you a list with your user account and any groups where you can create a fork of this repository (starting with `@`). If you see yourself only below "You have existing forks of this repository:", it means you already have a fork and you have nothing to do.
> {: .comment}
>
{: .hands_on}

# Clone the GitHub repository on your computer

To modify the content of the repository, you need a copy of it on your computer. This step of importing a git repository is called "cloning":

![Explanation of the cloning](../../images/PR_schema_02.png "Cloning a repository"){: width="900px"}

> <hands-on-title>Clone the GitHub repository</hands-on-title>
>
> 1. Get the URL of your fork. You can find this on the main page of your fork under the green button:
>    1. Click on **Code** (right)
>
>       ![Get fork URL](../../images/clone_code_button.png)
>
>    2. Copy the URL
>
>       It should be something like `https://github.com/< Your GitHub Username >/training-material.git`
>
> 1. Open a terminal
> 2. Navigate with `cd` to the folder in which you will clone the repository
> 3. Clone the repository with the command:
>
>    ```
>    $ git clone https://github.com/< Your GitHub Username >/training-material.git
>    ```
>
> 4. Navigate to the repository
>
>    ```
>    $ cd training-material
>    ```
{: .hands_on}

> <comment-title></comment-title>
> If you already have a local copy of the GitHub repository, you need to update it before doing any changes. To learn how to do that, please follow the [last section](#stay-up-to-date).
{: .comment}

# Create a new branch

You have now your repository locally and you want to modify it. For this example tutorial, you will add yourself as contributor of the project to appear on the [Hall of Fame]({% link hall-of-fame.md %}).

In GitHub flow, there is a concept: one new feature or change = one branch.

When you're working on a project, you're going to have a bunch of different features or ideas in progress at any given time – some of which are ready to go, and others which are not. Branching exists to help you manage this workflow. You should develop different features on different branches to help keep the changes you make simple and easy to review.

![Divergence of a branch compared to main](../../images/PR_branches_01.png "When you create a branch in your project, you're creating an environment where you can try out new ideas. Changes you make on a branch don't affect the main branch")

Here for this tutorial, you will create a branch called "my_new_branch" in which you will modify the `CONTRIBUTORS.yaml` file, the file used to generate the [Hall of Fame]({% link hall-of-fame.md %}).

> <hands-on-title>Create a branch</hands-on-title>
>
> 1. List the existing branch
>
>    ```
>    $ git branch
>      * main
>    ```
>
>    The branch on which you are is shown with the `*`
>
> 2. Create a new branch
>
>    ```
>    $ git checkout -b my_new_branch
>    Switched to a new branch 'my_new_branch'
>    ```
>
> 3. List the existing branch to check that the branch has been created and you are now on it
{: .hands_on}

This branch is added to your local copy:

![Creation of a branch on the local copy of the repository](../../images/PR_schema_03.png "Creation of a branch on the local copy of the repository")

# Make your changes on this branch

You have created your first branch! Now you want to make the change in the `CONTRIBUTING.yaml` file. By changing a file in this branch, it will diverge from the `main` branch. It will contain data that is only on this new branch:

![Divergence of the branch compared to main](../../images/PR_branches_02.png "The changes on your branch will not be on the main branch")

> <hands-on-title>Make changes in a branch</hands-on-title>
>
> 1. Open with your favorite text editor the `CONTRIBUTORS.yaml` file that is on your computer
> 2. Add yourself in the `CONTRIBUTORS.yaml` file
>
>    You should use your GitHub username and add it followed by `:` at the correct position given the alphabetical order
>
> 2. Save the file
> 3. Check the changes you made
>
>    ```
>    $ git status
>    On branch my_new_branch
>    Changes not staged for commit:
>      (use "git add/rm <file>..." to update what will be committed)
>      (use "git checkout -- <file>..." to discard changes in working directory)
>
>        modified:   CONTRIBUTORS.yaml
>
>    no changes added to commit (use "git add" and/or "git commit -a")
>    ```
>
> 4. Add the file and commit the changes
>
>    ```
>    $ git add CONTRIBUTORS.yaml
>    $ git commit -m "Add ..."
>    ```
>
> 5. Check that there are no more changes to commit with `git status`
>
{: .hands_on}

# Push your branch on your GitHub repository

The changes you made on your branch are only on the local copy of the repository. To propagate them online, you need to push them on your fork on GitHub:

> <hands-on-title>Push the changes</hands-on-title>
>
> 1. Push the changes to the GitHub repository
>
>    ```
>    $ git push origin my_new_branch
>    ```
>
>    > <details-title>"push origin"</details-title>
>    > When you `git push`ed, you specified `origin`. Git repositories can know that forks exist in multiple places. When you clone one, it creates a "remote" (a remote repository) which it names `origin`, set to the URL that you used when you cloned. By having multiple remotes, you can manage more complex workflows.
>    {: .details}
>
> 2. Go to your GitHub repository
> 3. Change to the "my_new_branch" branch:
>    1. Click on **Branch: main** (left)
>
>       ![Selecting branch on GitHub](../../images/PR_branch_github.png)
>
>    2. Select the branch "my_new_branch"
>
> 4. Check that your name is in the `CONTRIBUTORS.yaml` file
>
{: .hands_on}

![Pushing changes to the fork from the local copy](../../images/PR_schema_04.png "Pushing changes from the local copy to the fork on GitHub")

# Open a pull request

You pushed your changes to GitHub, but currently they are only on your fork. You want to have these changes in the main GitHub repository in order to appear on our [Hall of Fame]({% link hall-of-fame.md %}) online. You can't add or push directly the main GitHub repository, so you need to create what we call a pull request:

![Pull request](../../images/PR_schema_05.png "Pull Requests provide a way to notify project maintainers about the changes you'd like them to consider")

> <hands-on-title>Create a pull request</hands-on-title>
>
> 2. Go to your GitHub repository
> 1. Click on **Compare & pull request**
>
>    !["Opening a pull request"](../../images/PR_button.png)
>
> 3. Check that the selected branch are correct: **main** on the left and your branch name on the right
>
>    ![Branches in PR](../../images/PR_branch_check.png)
>
> 2. Fill in the pull request description
>
>    ![PR description](../../images/PR_description.png)
>
>    1. Add a title for the Pull Request
>    2. Add a message explaining the changes you made (Be kind <i class="fa fa-smile-o" aria-hidden="true"></i>)
>    3. Click on **Create pull request**
> 3. Go to **Pull requests** to check if it is there
{: .hands_on}

Once the pull is open, it will be reviewed. There are two possible outcomes:

1. Your pull request is accepted. Congratulations! Your changes will be merged into the main branch of the original repository. The website will be re-built and you will be in the [Hall of Fame]({% link hall-of-fame.md %})
2. Your pull request needs modifications: the reviewers will ask for some changes, possibly because the automatic tests are failing.

# Make the requested changes

One of the reviewers of your pull request asked you to add your name after your GitHub username in the `CONTRIBUTORS.yaml` file

> <hands-on-title>Make further changes</hands-on-title>
>
> 1. Make the requested changes in the `CONTRIBUTORS.yaml` file
>
>    It should look like
>
>    ```yaml
>    bebatut:
>         name: Bérénice Batut
>    ```
>
> 2. Check the changes that you made
>
>    ```
>    $ git status
>    On branch my_new_branch
>    Changes not staged for commit:
>      (use "git add/rm <file>..." to update what will be committed)
>      (use "git checkout -- <file>..." to discard changes in working directory)
>
>        modified:   CONTRIBUTORS.yaml
>
>    no changes added to commit (use "git add" and/or "git commit -a")
>    ```
>
> 3. Add the file and commit the changes
>
>    ```
>    $ git add CONTRIBUTORS.yaml
>    $ git commit -m "Add ..."
>    ```
>
> 4. Check that there are no more changes to commit with `git status`
> 5. Push the new changes to GitHub
>
>    ```
>    $ git push origin my_new_branch
>    ```
>
>    The pull request should be automatically updated
>
> 6. Check that the new changes are added to the pull request on GitHub
>
{: .hands_on}

# Check the automatic tests

When a pull request is opened, some automated tests are automatically launched on [Travis](http://travis-ci.com/) to be sure that the changes do not break the website, the URL are valid, etc.

On the bottom of your pull request, you can see the status of the tests:

- Yellow (with circle)

    ![Running](../../images/PR_test_yellow.png)

    The tests are still running

- Red (with cross)

    ![Failed tests](../../images/PR_test_red.png)

    When it is red, you can investigate why by clicking on **Details**. You will be redirected on [Travis](http://travis-ci.com/) where you can see the logs of the tests. Get in touch with us on [Gitter]({{ site.gitter_url }}) if you need help to understand the issue.

- Green (with tick)

    ![Passed tests](../../images/PR_test_green.png)

    The tests passed. Good job!

    Even it is green, we recommend to check the result of the tests, as some of tests are allowed to fail (to avoid too much noise).

# Stay up to date

You now want to work on a new tutorial or make some other new changes. However since you get a local copy, some changes have happened to the original GitHub repository. You need then to update your local copy of the repository before changing anything.

> <hands-on-title>Update the local copy</hands-on-title>
>
> 1. Move to the `main` branch
>
>    ```
>    $ git checkout main
>    ```
>
> 2. Add a reference to the original GitHub repository
>
>    ```
>    $ git remote add upstream https://github.com/galaxyproject/training-material.git
>    ```
>
>    > <comment-title>Error "remote upstream already exists"</comment-title>
>    > If you have done step 2 before and try to `remote add` again, git will tell you that a "remote upstream already exists". In this case you can safely continue to step 4.
>    {: .comment}
>
> 4. Update the local copy of the repository by "pulling" in the content of the original GitHub repository
>
>    ```
>    $ git pull upstream main
>    ```
>
{: .hands_on}

You can now restart the GitHub flow to propose new changes: start by [creating a new branch](#create-a-new-branch).


# Close the Pull Request

Great! You now know how to make pull request on GitHub, and how to make changes after a review.
Reviewers can now approve and merge your pull request.

Because this was just a practice pull request, let's close it again.


> <hands-on-title>Close the Pull Request</hands-on-title>
>
> Once you have run through all these steps, please close the pull request again.
>
> 1. Go to the [list of pull request tab on GitHub](https://github.com/galaxyproject/training-material/pulls)
> 2. Click on your pull request
> 3. Scroll to the bottom of the page
> 3. Click on "Close pull request" button
>
> Whenever you add your first real contribution, you can add yourself to the `CONTRIBUTORS.yaml` file in that PR.
>
{: .hands_on}


# Conclusion


With this tutorial, you have learned some basics `git` commands and principles:

![Summary of the links between GitHub, fork and local repository](../../images/PR_global_schema.png "Summary of the links between GitHub, fork and local repository")

You also learned the GitHub flow and its cycle:

1. Create a new branch in your local copy
2. Commit the changes in that branch
1. Push that branch to your fork on GitHub
1. Submit a pull request from that branch to the main repository
1. Wait for feedbacks and make requested changes
1. Update your local copy
1. Restart the cycle

You can now contribute and help us to improve our tutorials!

> <details-title>Git in Depth</details-title>
> This tutorial was a quick introduction to explain the basics of contributing to the training material. We recommend that everyone follow a more detailed git tutorials:
> - [Software Carpentry tutorial](http://swcarpentry.github.io/git-novice/) which explains a bit more in detail some git commands and the concept of remote, local and stagging
> - [Learn Git Branching](https://learngitbranching.js.org/) by GitHub to test the different git commands
>
> You should also download, print and keep always with you the [Git Cheat Sheet](https://services.github.com/on-demand/downloads/github-git-cheat-sheet.pdf)
{: .details}
