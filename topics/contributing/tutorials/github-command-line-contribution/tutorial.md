---
layout: tutorial_hands_on
topic_name: contributing
tutorial_name: github-command-line-contribution
---

# Introduction
{:.no_toc}

Most of the content is written in [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/) with some metadata (or variables) found in [YAML](http://yaml.org/) files. Everything is stored on a [GitHub](https://github.com) repository: [{{ site.github_repository }}]({{ site.github_repository }}).

The process of development of new content is open and transparent, using git and following the [GitHub flow](https://guides.github.com/introduction/flow/):

![Open source development process](../../images/oss_development_process.png "Open source development process")

1. Create a fork
1. Clone your fork of this repository to create a local copy on your computer 
1. Create a new branch in your local copy for each significant change
2. Commit the changes in that branch
1. Push that branch to your fork on GitHub
1. Submit a pull request from that branch to the [master repository]({{ site.github_repository }})
1. If you receive feedback, make changes in your local clone and push them to your branch on GitHub: the pull request will update automatically
1. Pull requests will be merged by the training team members after at least one other person has reviewed the Pull request and approved it.

> ### Agenda
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

> ### {% icon hands_on %} Hands-on: Fork the repository
>
> 1. Go on the GitHub repository: [{{ site.github_url }}]({{ site.github_url }})
> 2. Click on **Fork** (top-right corner of the page)
>   
>    ![Fork](../../images/PR_fork.jpg)
>
{: .hands_on}

# Clone the GitHub repository on your computer

To modify the content of the repository, you need to it on your computer. This step of importing a git repository is called "cloning":

![Explanation of the cloning](../../images/PR_schema_02.png "Cloning a repository"){: width="900px"}

> ### {% icon hands_on %} Hands-on: Clone the GitHub repository
>
> 1. Get the URL of your fork
>    1. Click on **Clone or download** (right)
>
>       ![Get fork URL](../../images/PR_github_url.png)
>
>    2. Copy the URL
>
>       It should be something like `https://github.com/< GitHub username >/training-material.git`
>
> 1. Open a terminal
> 2. Navigate with `cd` to the folder in which cloning the repository
> 3. Clone the repository
>
>    ```
>    $ git clone https://github.com/< GitHub username >/training-material.git
>    ```
>
> 4. Navigate to the repository
>
>    ```
>    $ cd training-material
>    ```
{: .hands_on}

# Create a new branch

You have now your repository locally and you want to modify it. For example for this tutorial, you want to add yourself as contributor of the project to appear on the [Hall of Fame]({{ site.baseurl }}/hall-of-fame).

In the GitHub flow, there is a concept: one new feature/change = one branch.

When you're working on a project, you're going to have a bunch of different features or ideas in progress at any given time – some of which are ready to go, and others which are not. Branching exists to help you manage this workflow.

![Divergence of a branch compared to master](../../images/PR_branches_01.png "When you create a branch in your project, you're creating an environment where you can try out new ideas. Changes you make on a branch don't affect the master branch")

Here for this tutorial, you will create a branch called "my_new_branch" in which you will modify the `CONTRIBUTORS.yaml` file, the file used to generate the [Hall of Fame]({{ site.baseurl }}/hall-of-fame).

> ### {% icon hands_on %} Hands-on: Create a branch
>
> 1. List the existing branch
>
>    ```
>    $ git branch 
>      * master
>    ```
>
> 2. Create a new branch
>
>    ```
>    $ git checkout -b my_new_branch
>    Switched to a new branch 'my_new_branch'
>    ```
{: .hands_on}

This branch is added to your local copy:

![Creation of a branch on the local copy of the repository](../../images/PR_schema_03.png "Creation of a branch on the local copy of the repository")

# Make your changes on this branch

You created a branch. Now you want to add the change in the `CONTRIBUTING.yaml` file. By doing that, you will diverge from the `master` branch:

![Divergence of the branch compared to master](../../images/PR_branches_02.png "The changes on your branch will not be on the master branch")

> ### {% icon hands_on %} Hands-on: Make changes in a branch
>
> 1. Modify the `CONTRIBUTORS.yaml` by adding yourself
>
>    You should use your GitHub username and add it followed by `:` at the correct position given the alphabetical order
>
> 2. Check the changes you made
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
> 2. Commit the changes
>
>    ```
>    $ git add CONTRIBUTORS.yaml
>    $ git commit -m "Add ..."
>    ```
>
> 3. Check there is no changes to add again with `git status`
> 
{: .hands_on}

# Push your branch on your GitHub repository

The changes you made on your branch are only on the local copy of the repository. To propagate them online, you need to push them on your fork on GitHub:

> ### {% icon hands_on %} Hands-on: Push the changes
>
> 1. Push the changes to the GitHub repository
>
>    ```
>    $ git push origin my_new_branch
>    ```
>
> 2. Go on your GitHub repository
> 3. Move to the "my_new_branch" branch:
>    1. Click on **Branch: master** (left)
>
>       ![Selecting branch on GitHub](../../images/PR_branch_github.png)
>
>    2. Select the branch "my_new_branch"
>
> 4. Check you are in the `CONTRIBUTORS.yaml` file
{: .hands_on}

![Pushing changes to the fork from the local copy](../../images/PR_schema_04.png "Pushing changes from the local copy to the fork on GitHub")

# Open a pull request

You pushed your changes to GitHub, but on your fork, You want now to have them in the main GitHub repository to appear on our [Hall of Fame]({{ site.baseurl }}/hall-of-fame) online. You can't add or push directly the main GitHub repository. You need to create what we call a pull request:

![Pull request](../../images/PR_schema_05.png "Pull Requests provide a way to notify project maintainers about the changes you'd like them to consider")

> ### {% icon hands_on %} Hands-on: Create a pull request
>
> 2. Go on your GitHub repository
> 1. Click on **Compare & pull request**
>
>    !["Opening a pull request"](../../images/PR_button.png)
>
> 3. Check that the selected branch are correct: **master** on the left and your branch name on the right
>
>    ![Branches in PR](../../images/PR_branch_check.png)
> 
> 2. Fill the pull request description
>    
>    ![PR description](../../images/PR_description.png)
>
>    1. Add a title for the Pull Request
>    2. Add a message explaining the changes you made (Be kind :smile:)
>    3. Click on **Create pull request**
> 3. Go to **Pull requests** to check if it is there
{: .hands_on}

Once the pull is open, it will be reviewed. 2 cases are then possible:

- Your pull request is accepted. Congratulations! Your changes will be merged into the master branch of the original repository. The website will be re-built and you will be in the [Hall of Fame]({{ site.baseurl }}/hall-of-fame)
- Your pull request need modifications: the reviewers ask for some changes or the automatic tests are failing.

# Make the requested changes

One of the reviewers of your pull request asked you to add your name after your GitHub username in the `CONTRIBUTORS.yaml` file

> ### {% icon hands_on %} Hands-on: Make further changes
>
> 2. Do the requested changes in the `CONTRIBUTORS.yaml` file 
>
>    It should look like
>    
>    ```
>    bebatut:
>         name: Bérénice Batut
>    ```
>
> 2. Check the changes you made
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
> 2. Commit the changes
>
>    ```
>    $ git add CONTRIBUTORS.yaml
>    $ git commit -m "Add ..."
>    ```
>
> 3. Check there is no changes to add again with `git status`
> 4. Push the new changes to GitHub
>
>    ```
>    $ git push origin my_new_branch
>    ```
>
>    The pull request should be automatically updated
>
> 5. Check the new changes has been added to the pull request on GitHub
>
{: .hands_on}

# Check the automatic tests

Once pull Request is open, some automated tests are automatically launched on [Travis](http://travis-ci.org/) to be sure that the changes do not break the website, the URL are valid, etc.

On the bottom of your pull request, you can see the status of the tests:

- Yellow (with circle)

    ![Running](../../images/PR_test_yellow.png)

    The tests are still running

- Red (with cross)
    
    ![Failed tests](../../images/PR_test_red.png)

    When it is red, you can investigate why by clicking on **Details**. You will be redirected on [Travis](http://travis-ci.org/) where you can see the logs of the tests. Get in touch with us on [Gitter]({{ site.gitter_url }}) if you need help to understand the issue.

- Green (with tick)

    ![Passed tests](../../images/PR_test_green.png)

    The tests passed. Good point!

    Even it is green, we recommend to check the result of the tests, as some of tests are allowed to fail (to avoid too much noise).

# Conclusion
{:.no_toc}

With this tutorial, you learn the basics to contribute using GitHub:

![Summary of the links between GitHub, fork and local repository](../../images/PR_global_schema.png "Summary of the links between GitHub, fork and local repository")

We can now contribute and help us to improve our tutorials

This tutorial is a quick introduction to explain the basics to contribute to the training material. We recommend anyone to follow a git tutorial, e.g. the one of [Software Carpentry](), and test yourself using the [GitHub test]().
