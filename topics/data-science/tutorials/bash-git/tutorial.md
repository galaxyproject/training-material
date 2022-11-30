---
layout: tutorial_hands_on
title: "Version Control with Git"
time_estimation: "65M"
questions:
- "What is version control and why should I use it?"
- "How do I get set up to use Git?"
- "Where does Git store information?"
- "How do I record changes in Git?"
- "How do I check the status of my version control repository?"
- "How do I record notes about what changes I made and why?"
- "How can I identify old versions of files?"
- "How do I review my changes?"
- "How can I recover old versions of files?"
objectives:
- "Understand the benefits of an automated version control system."
- "Understand the basics of how automated version control systems work."
- "Configure `git` the first time it is used on a computer."
- "Understand the meaning of the `--global` configuration flag."
- "Create a local Git repository."
- "Describe the purpose of the `.git` directory."
- "Go through the modify-add-commit cycle for one or more files."
- "Explain where information is stored at each stage of that cycle."
- "Distinguish between descriptive and non-descriptive commit messages."
- "Explain what the HEAD of a repository is and how to use it."
- "Identify and use Git commit numbers."
- "Compare various versions of tracked files."
- "Restore old versions of files."
key_points:
- "Version control is like an unlimited 'undo'."
- "Version control also allows many people to work in parallel."
- "Use `git config` with the `--global` option to configure a user name, email address, editor, and other preferences once per machine."
- "`git init` initializes a repository."
- "Git stores all of its repository data in the `.git` directory."
- "`git status` shows the status of a repository."
- "Files can be stored in a project's working directory (which users see), the staging area (where the next commit is being built up) and the local repository (where commits are permanently recorded)."
- "`git add` puts files in the staging area."
- "`git commit` saves the staged content as a new commit in the local repository."
- "Write a commit message that accurately describes your changes."
- "`git diff` displays differences between commits."
- "`git checkout` recovers old versions of files."
contributions:
  authorship:
    - Sofokli5
  editing:
    - fpsom
    - shiltemann
    - hexylena
  funding:
    - carpentries
    - erasmusplus

---

Version control is a way of tracking the change history of a project, and `git` is one of the most popular systems for doing that! This tutorial will guide you through the basics of using git for version control.

> <agenda-title></agenda-title>
>
> In this tutorial, you will learn how to create a git repo, and begin working with it.
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Basics

We'll start by exploring how version control can be used
to keep track of what one person did and when.
Even if you aren't collaborating with other people,
automated version control is much better than this situation:

![Cartoon titled 'final'.doc, showing a grad student and their advisor going through multiple revisions. The first named final.doc, then final_rev.2.doc, final_rev.6.comments.doc, a long filename with the revision number 18, until a final filename, revision 22, with special characters indicating frustration where the file name includes the text 'why did I come to grad school'.](../../images/bash-git/phd101212s.png "Piled Higher and Deeper by Jorge Cham, http://www.phdcomics.com")

We've all been in this situation before: it seems unnecessary to have
multiple nearly-identical versions of the same document. Some word
processors let us deal with this a little better, such as Microsoft
Word's
[Track Changes](https://support.office.com/en-us/article/Track-changes-in-Word-197ba630-0f5f-4a8e-9a77-3712475e806a),
Google Docs' [version history](https://support.google.com/docs/answer/190843?hl=en), or
LibreOffice's [Recording and Displaying Changes](https://help.libreoffice.org/Common/Recording_and_Displaying_Changes).

Version control systems start with a base version of the document and
then record changes you make each step of the way. You can
think of it as a recording of your progress: you can rewind to start at the base
document and play back each change you made, eventually arriving at your
more recent version.

![Changes Are Saved Sequentially, graphic shows three documents with text being added in each new revision.](../../images/bash-git/play-changes.svg)

Once you think of changes as separate from the document itself, you
can then think about "playing back" different sets of changes on the base document, ultimately
resulting in different versions of that document. For example, two users can make independent
sets of changes on the same document.

![Different Versions Can be Saved, showing a document splitting into two, with different changes.](../../images/bash-git/versions.svg)

Unless multiple users make changes to the same section of the document - a conflict - you can
incorporate two sets of changes into the same base document.

![Multiple Versions Can be Merged, shows two documents with different changes merging into a final document with both changes.](../../images/bash-git/merge.svg)

A version control system is a tool that keeps track of these changes for us,
effectively creating different versions of our files. It allows us to decide
which changes will be made to the next version (each record of these changes is
called a [commit](https://git-scm.com/docs/gitglossary#def_commit), and keeps useful metadata
about them. The complete history of commits for a particular project and their
metadata make up a [repository](https://git-scm.com/docs/gitglossary#def_repository).
Repositories can be kept in sync across different computers, facilitating
collaboration among different people.

> <tip-title>The Long History of Version Control Systems</tip-title>
>
> Automated version control systems are nothing new.
> Tools like [RCS](https://en.wikipedia.org/wiki/Revision_Control_System), [CVS](https://en.wikipedia.org/wiki/Concurrent_Versions_System), or [Subversion](https://en.wikipedia.org/wiki/Apache_Subversion) have been around since the early 1980s and are used by
> many large companies.
> However, many of these are now considered legacy systems (i.e., outdated) due to various
> limitations in their capabilities.
> More modern systems, such as Git and [Mercurial](https://swcarpentry.github.io/hg-novice/),
> are *distributed*, meaning that they do not need a centralized server to host the repository.
> These modern systems also include powerful merging tools that make it possible for
> multiple authors to work on
> the same files concurrently.
{: .tip}

> <question-title>Paper Writing</question-title>
>
> *   Imagine you drafted an excellent paragraph for a paper you are writing, but later ruin
>     it. How would you retrieve the *excellent* version of your conclusion? Is it even possible?
>
> *   Imagine you have 5 co-authors. How would you manage the changes and comments
>     they make to your paper?  If you use LibreOffice Writer or Microsoft Word, what happens if
>     you accept changes made using the `Track Changes` option? Do you have a
>     history of those changes?
>
> > <solution-title></solution-title>
> >
> > *   Recovering the excellent version is only possible if you created a copy
> >     of the old version of the paper. The danger of losing good versions
> >     often leads to the problematic workflow illustrated in the PhD Comics
> >     cartoon at the top of this page.
> >
> > *   Collaborative writing with traditional word processors is cumbersome.
> >     Either every collaborator has to work on a document sequentially
> >     (slowing down the process of writing), or you have to send out a
> >     version to all collaborators and manually merge their comments into
> >     your document. The 'track changes' or 'record changes' option can
> >     highlight changes for you and simplifies merging, but as soon as you
> >     accept changes you will lose their history. You will then no longer
> >     know who suggested that change, why it was suggested, or when it was
> >     merged into the rest of the document. Even online word processors like
> >     Google Docs or Microsoft Office Online do not fully resolve these
> >     problems.
> {: .solution}
{: .question }


Before diving in the tutorial, we need to open {% tool [RStudio](interactive_tool_rstudio) %}. If you do not know how or never interacted with RStudio, please follow the [dedicated tutorial]({% link topics/galaxy-interface/tutorials/rstudio/tutorial.md %}).

{% snippet faqs/galaxy/interactive_tools_rstudio_launch.md %}

> <hands-on-title>Installing git</hands-on-title>
> The R Console and other interactive tools like RStudio are great for prototyping code and exploring data, but sooner or later we will want to use our program in a pipeline or run it in a shell script to process thousands of data files. This is one of those cases and, in order to do that, we will use the terminal provided by the RStudio itself.
> We go to "Tools" and pick the "Shell..." option and we are good to go. Our workspace is the left, terminal window that just opened.
>
> Fortunately, [miniconda](https://docs.conda.io/en/latest/miniconda.html) is already installed. Miniconda is a package manager that simplifies the installation processes. We can and will use it to install every essential package for our tutorial. However, it is of critical importance that we do that in an new environment within our existing base and install our packages in said environment.
>
> > <code-in-title>Environment and Packages</code-in-title>
> > ```bash
> > $ conda create -n name_of_your_env nano git
> > $ conda activate name_of_your_env
> > ```
> {: .code_in}
>
>
> | Software | Version | Manual | Available for | Description |
> | -------- | ------------ | ------ | ------------- | ----------- |
> | [git](https://git-scm.com/) | 2.35.3 | [git Manual](https://git-scm.com/book/en/v2) | Linux, MacOS | Git is a free and open source distributed version control system designed to handle everything from small to very large projects with speed and efficiency. |
> | [GNU nano](https://www.nano-editor.org/) | 2.9.8 | [Nano manual](https://www.nano-editor.org/dist/latest/nano.html) | Linux, MacOS | GNU nano is a small and friendly text editor. |
>
{: .hands_on}

# Setting up git

When we use Git on a new computer for the first time,
we need to configure a few things. Below are a few examples
of configurations we will set as we get started with Git:

*   our name and email address,
*   what our preferred text editor is,
*   and that we want to use these settings globally (i.e. for every project).

On a command line, Git commands are written as `git verb options`,
where `verb` is what we actually want to do and `options` is additional optional information which may be needed for the `verb`. So here is how
Sherlock sets up his new laptop:

> <code-in-title>Setting up with bash</code-in-title>
> ```bash
> $ git config --global user.name "Sherlock Holmes"
> $ git config --global user.email "sherlock@baker.street"
> ```
{: .code-in}

Please use your own name and email address instead of Sherlocks's. This user name and email will be associated with your subsequent Git activity,
which means that any changes pushed to
[GitHub](https://github.com/),
[BitBucket](https://bitbucket.org/),
[GitLab](https://gitlab.com/) or
another Git host server
after this lesson will include this information.

For this lesson, we will be interacting with [GitHub](https://github.com/) and so the email address used should be the same as the one used when setting up your GitHub account. If you are concerned about privacy, please review [GitHub's instructions for keeping your email address private][git-privacy].

> <tip-title>Keeping your email private</tip-title>
>
> If you choose to use a private email address with GitHub, then use that same email address for the `user.email` value, e.g. `username@users.noreply.github.com` replacing `username` with your GitHub one.
{: .tip}


> <tip-title>Line Endings</tip-title>
>
> As with other keys, when you hit <kbd>Enter</kbd> or <kbd></kbd> or on Macs, <kbd>Return</kbd>, on your keyboard,
> your computer encodes this input as a character.
> Different operating systems use different character(s) to represent the end of a line.
> (You may also hear these referred to as newlines or line breaks.)
> Because Git uses these characters to compare files,
> it may cause unexpected issues when editing a file on different machines.
> Though it is beyond the scope of this lesson, you can read more about this issue
> [in the Pro Git book](https://www.git-scm.com/book/en/v2/Customizing-Git-Git-Configuration#_core_autocrlf).
>
> You can change the way Git recognizes and encodes line endings
> using the `core.autocrlf` command to `git config`.
> The following settings are recommended:
>
> > <code-in-title>Settings with bash</code-in-title>
> > On macOS and Linux:
> >
> >```bash
> >$ git config --global core.autocrlf input
> >```
> >
> >And on Windows:
> >
> >```bash
> >$ git config --global core.autocrlf true
> >```
> {: .code-in}
{: .tip}

Sherlock also has to set his favorite text editor, following this table:

| Editor             | Configuration command                            |
|:-------------------|:-------------------------------------------------|
| Atom | `$ git config --global core.editor "atom --wait"`|
| nano               | `$ git config --global core.editor "nano -w"`    |
| BBEdit (Mac, with command line tools) | `$ git config --global core.editor "bbedit -w"`    |
| Sublime Text (Mac) | `$ git config --global core.editor "/Applications/Sublime\ Text.app/Contents/SharedSupport/bin/subl -n -w"` |
| Sublime Text (Win, 32-bit install) | `$ git config --global core.editor "'c:/program files (x86)/sublime text 3/sublime_text.exe' -w"` |
| Sublime Text (Win, 64-bit install) | `$ git config --global core.editor "'c:/program files/sublime text 3/sublime_text.exe' -w"` |
| Notepad (Win)    | `$ git config --global core.editor "c:/Windows/System32/notepad.exe"`|
| Notepad++ (Win, 32-bit install)    | `$ git config --global core.editor "'c:/program files (x86)/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"`|
| Notepad++ (Win, 64-bit install)    | `$ git config --global core.editor "'c:/program files/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"`|
| Kate (Linux)       | `$ git config --global core.editor "kate"`       |
| Gedit (Linux)      | `$ git config --global core.editor "gedit --wait --new-window"`   |
| Scratch (Linux)       | `$ git config --global core.editor "scratch-text-editor"`  |
| Emacs              | `$ git config --global core.editor "emacs"`   |
| Vim                | `$ git config --global core.editor "vim"`   |
| VS Code                | `$ git config --global core.editor "code --wait"`   |

It is possible to reconfigure the text editor for Git whenever you want to change it.

> <tip-title>Exiting Vim</tip-title>
> Note that Vim is the default editor for many programs. If you haven't used Vim before and wish to exit a session without saving
> your changes, press <kbd>Esc</kbd> then type `:q!` and hit <kbd>Enter</kbd> or <kbd></kbd> or on Macs, <kbd>Return</kbd>.
> If you want to save your changes and quit, press <kbd>Esc</kbd> then type `:wq` and hit <kbd>Enter</kbd> or <kbd></kbd> or on Macs, <kbd>Return</kbd>.
{: .tip}

Git (2.28+) allows configuration of the name of the branch created when you
initialize any new repository.  Sherlock decides to use that feature to set it to `main` so
it matches the cloud service he will eventually use.

> <code-in-title>Configure the name of the created branch</code-in-title>
> ```bash
> $ git config --global init.defaultBranch main
> ```
{: .code-in}

> <tip-title>Default Git branch naming</tip-title>
>
> Source file changes are associated with a "branch."
> For new learners in this lesson, it's enough to know that branches exist, and this lesson uses one branch.
> By default, Git will create a branch called `master`
> when you create a new repository with `git init` (as explained in the next Episode). This term evokes
> the racist practice of human slavery and the
> [software development community](https://github.com/github/renaming)  has moved to adopt
> more inclusive language.
>
> In 2020, most Git code hosting services transitioned to using `main` as the default
> branch. As an example, any new repository that is opened in GitHub and GitLab default
> to `main`.  However, Git has not yet made the same change.  As a result, local repositories
> must be manually configured have the same main branch name as most cloud services.
>
> For versions of Git prior to 2.28, the change can be made on an individual repository level.  The
> command for this is in the next episode.  Note that if this value is unset in your local Git
> configuration, the `init.defaultBranch` value defaults to `master`.
>
{: .tip}

The five commands we just ran above only need to be run once: the flag `--global` tells Git
to use the settings for every project, in your user account, on this computer.

You can check your settings at any time:

> <code-in-title>Checking your settings</code-in-title>
> ```bash
> $ git config --list
> ```
{: .code-in}

You can change your configuration as many times as you want: use the
same commands to choose another editor or update your email address.

> <tip-title>Proxy</tip-title>
>
> In some networks you need to use a
> [proxy](https://en.wikipedia.org/wiki/Proxy_server). If this is the case, you
> may also need to tell Git about the proxy:
>
> > <code-in-title>Git and proxy</code-in-title>
> > ```bash
> > $ git config --global http.proxy proxy-url
> > $ git config --global https.proxy proxy-url
> > ```
> {: .code-in}
>
> > <code-in-title>To disable the proxy, use</code-in-title>
> >
> > ```bash
> > $ git config --global --unset http.proxy
> > $ git config --global --unset https.proxy
> > ```
> {: .code-in}
{: .tip}

> <tip-title>Git Help and Manual</tip-title>
>
> Always remember that if you forget the subcommands or options of a `git` command, you can access the relevant list of options typing `git <command> -h` or access the corresponding Git manual by typing
> `git <command> --help`, e.g.:
>
> > <code-in-title>Access help </code-in-title>
> >```bash
> >$ git config -h
> >$ git config --help
> >```
> {: .code-in}
>
> While viewing the manual, remember the `:` is a prompt waiting for commands and you can press <kbd>Q</kbd> to exit the manual.
>
> More generally, you can get the list of available `git` commands and further resources of the Git manual typing:
>
> > <code-in-title>Access available commands</code-in-title>
> >```bash
> >$ git help
> >```
>
> {: .code-in}
>
{: .tip}

[git-privacy]: https://help.github.com/articles/keeping-your-email-address-private/

# Create a Repo

Once Git is configured,
we can start using it.

We will continue with the story of Sherlock who is investigating a crime and is collecting information about suspects.

![Cartoon of sherlock examining the git logo with his magnifying glass.](../../images/bash-git/sherlock_git.png)


First, let's create a directory for our work and then move into that directory:

> <code-in-title>Create a workspace</code-in-title>
> ```bash
> $ mkdir suspects
> $ cd suspects
> ```
{: .code-in}

Then we tell Git to make `suspects` a [repository](https://git-scm.com/docs/gitglossary#def_repository)
-- a place where Git can store versions of our files:

> <code-in-title>Turn our workspace into directory</code-in-title>
> ```bash
> $ git init
> ```
{: .code-in}

It is important to note that `git init` will create a repository that
includes subdirectories and their files---there is no need to create
separate repositories nested within the `suspects` repository, whether
subdirectories are present from the beginning or added later. Also, note
that the creation of the `suspects` directory and its initialization as a
repository are completely separate processes.

If we use `ls` to show the directory's contents,
it appears that nothing has changed:

> <code-in-title>Show directory content</code-in-title>
> ```bash
> $ ls
> ```
{: .code-in}

But if we add the `-a` flag to show everything,
we can see that Git has created a hidden directory within `suspects` called `.git`:

> <code-in-title>Show everything in our directory</code-in-title>
> ```bash
> $ ls -a
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> .	..	.git
> ```
{: .code-out}

Git uses this special subdirectory to store all the information about the project,
including all files and sub-directories located within the project's directory.
If we ever delete the `.git` subdirectory, we will lose the project's history.

Next, we will change the default branch to be called `main`.
This might be the default branch depending on your settings and version
of git.

> <code-in-title>Rename the branch</code-in-title>
> ```bash
> $ git checkout -b main
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> Switched to a new branch 'main'
> ```
{: .code-out}


We can check that everything is set up correctly
by asking Git to tell us the status of our project:

> <code-in-title>Check</code-in-title>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
>
> No commits yet
>
> nothing to commit (create/copy files and use "git add" to track)
> ```
{: .code-out}

If you are using a different version of `git`, the exact
wording of the output might be slightly different.

## Places to Create Git Repositories

 Along with tracking information about suspects (the project we have already created), Sherlock would also like to track information specific clues. So, Sherlock creates a `clues` project inside his `suspects` project with the following sequence of commands:

> <code-in-title>Create a project within a project</code-in-title>
> ```bash
> $ cd ~/Desktop   # return to Desktop directory
> $ cd suspects     # go into suspects directory, which is already a Git repository
> $ ls -a          # ensure the .git subdirectory is still present in the suspects directory
> $ mkdir clues    # make a subdirectory suspects/clues
> $ cd clues       # go into clues subdirectory
> $ git init       # make the clues subdirectory a Git repository
> $ ls -a          # ensure the .git subdirectory is present indicating we have created a new Git repository
> ```
{: .code-in}

> <question-title>Tracking in a subdirectory</question-title>
> Is the `git init` command, run inside the `clues` subdirectory, required for tracking files stored in the `clues` subdirectory?
>
> > <solution-title></solution-title>
> >
> > No. Sherlock does not need to make the `clues` subdirectory a Git repository
> > because the `suspects` repository will track all files, sub-directories, and
> > subdirectory files under the `suspects` directory.  Thus, in order to track
> > all information about clues, Sherlock only needed to add the `clues` subdirectory
> > to the `suspects` directory.
> >
> {: .solution}

{: .question}

> <tip-title>"Nested" repositories</tip-title>
> Additionally, Git repositories can interfere with each other if they are "nested":
> the outer repository will try to version-control
> the inner repository. Therefore, it's best to create each new Git
> repository in a separate directory. To be sure that there is no conflicting
> repository in the directory, check the output of `git status`. If it looks
> like the following, you are good to go to create a new repository as shown
> above:
>
> > <code-in-title>Before creating a new repository</code-in-title>
> > ```bash
> > $ git status
> > ```
> {: .code-in}
> > <code-out-title></code-out-title>
> >```
> >fatal: Not a git repository (or any of the parent directories): .git
> >```
> {: .code-out}
>
{: .tip}

## Correcting `git init` Mistakes
Dr. Watson explains to Sherlock how a nested repository is redundant and may cause confusion down the road. Sherlock would like to remove the nested repository. How can Sherlock undo his last `git init` in the `clues` subdirectory?

### USE `rm` WITH CAUTION!

Removing files from a Git repository needs to be done with caution. But we have not learned yet how to tell Git to track a particular file; we will learn this in the next episode. Files that are not tracked by Git can easily be removed like any other "ordinary" files with
> <code-in-title>Remove files </code-in-title>
> ```bash
> $ rm filename
> ```
{: .code-in}

Similarly a directory can be removed using `rm -r dirname` or `rm -rf dirname`. If the files or folder being removed in this fashion are tracked by Git, then their removal becomes another change that we will need to track, as we will see later.

Git keeps all of its files in the `.git` directory. To recover from this little mistake, Sherlock can just remove the `.git` folder in the clues subdirectory by running the following command from inside the `suspects` directory:

> <code-in-title>Remove the `.git` folder</code-in-title>
> ```bash
> $ rm -rf clues/.git
> ```
{: .code-in}

But be careful! Running this command in the wrong directory will remove the entire Git history of a project you might want to keep. Therefore, always check your current directory using the command `pwd`.

# Tracking Changes

First let's make sure we're still in the right directory.
You should be in the `suspects` directory.

> <code-in-title>Check directory</code-in-title>
> ```bash
> $ cd ~/suspects
> ```
{: .code-in}

Let's create a file called `colonel.txt` that contains some notes about the colonel Smith's connection to the case.
We'll use `nano` to edit the file;
you can use whatever editor you like.
In particular, this does not have to be the `core.editor` you set globally earlier. But remember, the bash command to create or edit a new file will depend on the editor you choose (it might not be `nano`). For a refresher on text editors, check out ["Which Editor?"](https://swcarpentry.github.io/shell-novice/03-create/) in [The Unix Shell](https://swcarpentry.github.io/shell-novice/) lesson.

> <code-in-title>Edit file with nano in bash</code-in-title>
> ```bash
> $ nano colonel.txt
> ```
{: .code-in}

Type the text below into the `colonel.txt` file:

```
No alibi for the night of murder.
```


Let's first verify that the file was properly created by running the list command (`ls`):

> <code-in-title>Check new file</code-in-title>
> ```bash
> $ ls
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> colonel.txt
> ```
{: .code-out}

`colonel.txt` contains a single line, which we can see by running:

> <code-in-title>Inspect the new file</code-in-title>
> ```bash
> $ cat colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> No alibi for the night of murder.
> ```
{: .code-out}

If we check the status of our project again,
Git tells us that it's noticed the new file:

> <code-in-title>Check</code-in-title>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
>
> No commits yet
>
> Untracked files:
>   (use "git add <file>..." to include in what will be committed)
>
> colonel.txt
>
> nothing added to commit but untracked files present (use "git add" to track)
> ```
{: .code-out}

The "untracked files" message means that there's a file in the directory that Git isn't keeping track of. We can tell Git to track a file using `git add`:

> <code-in-title>Track file</code-in-title>
> ```bash
> $ git add colonel.txt
> ```
{: .code-in}

and then check that the right thing happened:

> <code-in-title>Check</code-in-title>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
>
> No commits yet
>
> Changes to be committed:
>   (use "git rm --cached <file>..." to unstage)
>
>	new file:   colonel.txt
>
> ```
{: .code-out}

Git now knows that it's supposed to keep track of `colonel.txt`, but it hasn't recorded these changes as a commit yet. To get it to do that, we need to run one more command:

> <code-in-title>Record changes as a commit with a descriptive title</code-in-title>
> ```bash
> $ git commit -m "Start notes for colonel Smith as a suspect"
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
>[main (root-commit) f22b25e] Start notes on colonel Smith as a suspect
> 1 file changed, 1 insertion(+)
> create mode 100644 colonel.txt
> ```
{: .code-out}

When we run `git commit`,
Git takes everything we have told it to save by using `git add`
and stores a copy permanently inside the special `.git` directory.
This permanent copy is called a [commit](https://git-scm.com/docs/gitglossary#def_commit)
(or [revision](https://git-scm.com/docs/gitglossary#def_revision)) and its short identifier is `f22b25e`. Your commit may have another identifier.

We use the `-m` flag (for "message")
to record a short, descriptive, and specific comment that will help us remember later on what we did and why.
If we just run `git commit` without the `-m` option,
Git will launch `nano` (or whatever other editor we configured as `core.editor`)
so that we can write a longer message.

[Good commit messages][commit-messages] start with a brief <50 characters) statement about the
changes made in the commit. Generally, the message should complete the sentence "If applied, this commit will" <commit message here>.
If you want to go into more detail, add a blank line between the summary line and your additional notes. Use this additional space to explain why you made changes and/or what their impact will be.

> <code-in-title>If we run `git status` now:</code-in-title>
>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
> nothing to commit, working directory clean
> ```
{: .code-out}

it tells us everything is up to date.
If we want to know what we've done recently,
we can ask Git to show us the project's history using `git log`:

> <code-in-title>Access project's history</code-in-title>
> ```bash
> $ git log
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
> Author: Sherlock Holmes <sherlock@baker.street>
> Date:   Thu Aug 22 09:51:46 2013 -0400
>
>   Start notes on colonel as a suspect
> ```
{: .code-out}

`git log` lists all commits  made to a repository in reverse chronological order.
The listing for each commit includes
the commit's full identifier
(which starts with the same characters as
the short identifier printed by the `git commit` command earlier),
the commit's author,
when it was created,
and the log message Git was given when the commit was created.

> <tip-title>Where Are My Changes?</tip-title>
>
> If we run `ls` at this point, we will still see just one file called `colonel.txt`.
> That's because Git saves information about files' history
> in the special `.git` directory mentioned earlier
> so that our filesystem doesn't become cluttered
> (and so that we can't accidentally edit or delete an old version).
{: .tip}

Now suppose Sherlock adds more information to the file.
(Again, we'll edit with `nano` and then `cat` the file to show its contents;
you may use a different editor, and don't need to `cat`.)

> <code-in-title>Edit with nano</code-in-title>
> ```bash
> $ nano colonel.txt
> $ cat colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> No alibi for the night of murder.
> No clear motive. Seems high unlikely.
> ```
{: .code-out}

When we run `git status` now,
it tells us that a file, it already knows about, has been modified:

> <code-in-title>Find about current status</code-in-title>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
> Changes not staged for commit:
>  (use "git add <file>..." to update what will be committed)
>  (use "git checkout -- <file>..." to discard changes in working directory)
>
>	modified:   colonel.txt
>
> no changes added to commit (use "git add" and/or "git commit -a")
> ```
{: .code-out}

The last line is the key phrase:
"no changes added to commit".
We have changed this file,
but we haven't told Git we will want to save those changes
(which we do with `git add`)
nor have we saved them (which we do with `git commit`).
So let's do that now. It is good practice to always review
our changes before saving them. We do this using `git diff`.
This shows us the differences between the current state
of the file and the most recently saved version:

> <code-in-title>Review the changes </code-in-title>
> ```bash
> $ git diff
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> diff --git a/colonel.txt b/colonel.txt
> index df0654a..315bf3a 100644
> --- a/colonel.txt
> +++ b/colonel.txt
> @@ -1 +1,2 @@
>  No alibi for the night of murder.
> +No clear motive. Seems high unlikely.
> ```
{: .code-out}

The output is cryptic because
it is actually a series of commands for tools like editors and `patch`
telling them how to reconstruct one file given the other.
If we break it down into pieces:

1.  The first line tells us that Git is producing output similar to the Unix `diff` command
    comparing the old and new versions of the file.
2.  The second line tells exactly which versions of the file
    Git is comparing;
    `df0654a` and `315bf3a` are unique computer-generated labels for those versions.
3.  The third and fourth lines once again show the name of the file being changed.
4.  The remaining lines are the most interesting, they show us the actual differences
    and the lines on which they occur.
    In particular,
    the `+` marker in the first column shows where we added a line.

After reviewing our change, it's time to commit it:

> <code-in-title>Commit after checking changes</code-in-title>
> ```bash
> $ git commit -m "Add concerns about existence of motive for colonel"
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
> Changes not staged for commit:
>   (use "git add <file>..." to update what will be committed)
>   (use "git checkout -- <file>..." to discard changes in working directory)
>
> 	modified:   colonel.txt
>
> no changes added to commit (use "git add" and/or "git commit -a")
> ```
{: .code-out}

Whoops:
Git won't commit because we didn't use `git add` first.
Let's fix that:

> <code-in-title>First `add` then `commit`</code-in-title>
> ```bash
> $ git add colonel.txt
> $ git commit -m "Add concerns about existence of motive for colonel"
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> [main 34961b1] Add concerns about existence of motive for colonel
>  1 file changed, 1 insertion(+)
> ```
{: .code-out}

Git insists that we add files to the set we want to commit
before actually committing anything. This allows us to commit our
changes in stages and capture changes in logical portions rather than
only large batches.
For example,
suppose we're adding a few citations to relevant research to our thesis.
We might want to commit those additions,
and the corresponding bibliography entries,
but *not* commit some of our work drafting the conclusion
(which we haven't finished yet).

To allow for this,
Git has a special *staging area*
where it keeps track of things that have been added to
the current [changeset](https://git-scm.com/docs/gitglossary#def_changeset)
but not yet committed.

> <tip-title>Staging Area</tip-title>
>
> If you think of Git as taking snapshots of changes over the life of a project,
> `git add` specifies *what* will go in a snapshot
> (putting things in the staging area),
> and `git commit` then *actually takes* the snapshot, and
> makes a permanent record of it (as a commit).
> If you don't have anything staged when you type `git commit`,
> Git will prompt you to use `git commit -a` or `git commit --all`,
> which is kind of like gathering *everyone* to take a group photo!
> However, it's almost always better to
> explicitly add things to the staging area, because you might
> commit changes you forgot you made. (Going back to the group photo simile,
> you might get an extra with incomplete makeup walking on
> the stage for the picture because you used `-a`!)
> Try to stage things manually,
> or you might find yourself searching for "git undo commit" more
> than you would like!
{: .tip}

![The Git Staging Area cartoon, a document is shown going into the staging area via "git add", and then into the repository via "git commit".](../../images/bash-git/git-staging-area.svg)

Let's watch as our changes to a file move from our editor
to the staging area
and into long-term storage.
First,
we'll add another line to the file:

> <code-in-title>Add and review new line with `nano` and `cat`</code-in-title>
> ```bash
> $ nano colonel.txt
> $ cat colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> No alibi for the night of murder.
> No clear motive. Seems high unlikely.
> Fingerprints on victims glasses.
> ```
{: .code-out}

> <code-in-title>Review changes</code-in-title>
> ```bash
> $ git diff
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> diff --git a/colonel.txt b/colonel.txt
> index 315bf3a..b36abfd 100644
> --- a/colonel.txt
> +++ b/colonel.txt
> @@ -1,2 +1,3 @@
>  No alibi for the night of murder.
>  No clear motive. Seems high unlikely.
> +Fingerprints on victims glasses.
> ```
{: .code-out}

So far, so good:
we've added one line to the end of the file
(shown with a `+` in the first column).
Now let's put that change in the staging area
and see what `git diff` reports:

> <code-in-title>See changes</code-in-title>
> ```bash
> $ git add colonel.txt
> $ git diff
> ```
{: .code-in}

But, there is no output:
as far as Git can tell,
there's no difference between what it's been asked to save permanently
and what's currently in the directory.
However,
if we do this:

> <code-in-title>See what is in the staging area</code-in-title>
> ```bash
> $ git diff --staged
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> diff --git a/colonel.txt b/colonel.txt
> index 315bf3a..b36abfd 100644
> --- a/colonel.txt
> +++ b/colonel.txt
> @@ -1,2 +1,3 @@
>  No alibi for the night of murder.
>  No clear motive. Seems high unlikely.
> +Fingerprints on victims glasses.
> ```
{: .code-out}

it shows us the difference between
the last committed change
and what's in the staging area.
Let's save our changes:

> <code-in-title>Save changes</code-in-title>
> ```bash
> $ git commit -m "Make notes about colonel's fingerprints"
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> [main 005937f] Make notes about colonel's fingerprints
>  1 file changed, 1 insertion(+)
> ```
{: .code-out}

> <code-in-title>Check new status</code-in-title>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
> nothing to commit, working directory clean
> ```
{: .code-out}

> <code-in-title>and look at the history of what we've done so far: </code-in-title>
> ```bash
> $ git log
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> commit 005937fbe2a98fb83f0ade869025dc2636b4dad5 (HEAD -> main)
> Author: Holmes Sherlock <Holmes@tran.sylvan.ia>
> Date:   Thu Aug 22 10:14:07 2013 -0400
>
>     Make notes about colonel's fingerprints
>
> commit 34961b159c27df3b475cfe4415d94a6d1fcd064d
> Author: Sherlock Holmes <sherlock@baker.street>
> Date:   Thu Aug 22 10:07:21 2013 -0400
>
>     Add concerns about existence of motive for colonel
>
> commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
> Author: Sherlock Holmes <sherlock@baker.street>
> Date:   Thu Aug 22 09:51:46 2013 -0400
>
>     Start notes on colonel as a suspect
> ```
{: .code-out}

> <tip-title>Word-based diffing</tip-title>
>
> Sometimes, e.g. in the case of the text documents a line-wise
> diff is too coarse. That is where the `--color-words` option of
> `git diff` comes in very useful as it highlights the changed
> words using colors.
{: .tip}

> <tip-title>Paging the Log</tip-title>
>
> When the output of `git log` is too long to fit in your screen,
> `git` uses a program to split it into pages of the size of your screen.
> When this "pager" is called, you will notice that the last line in your
> screen is a `:`, instead of your usual prompt.
>
> *   To get out of the pager, press <kbd>Q</kbd>.
> *   To move to the next page, press <kbd>Spacebar</kbd>.
> *   To search for `some_word` in all pages,
>     press <kbd></kbd>
>     and type `some_word`.
>     Navigate through matches pressing <kbd>N</kbd>.
{: .tip}

> <tip-title>Limit Log Size</tip-title>
>
> To avoid having `git log` cover your entire terminal screen, you can limit the
> number of commits that Git lists by using `-N`, where `N` is the number of
> commits that you want to view. For example, if you only want information from
> the last commit you can use:
>
> > <code-in-title>See specific number of commits</code-in-title>
> > ```bash
> > $ git log -1
> > ```
> {: .code-in}
>
> > <code-out-title></code-out-title>
> > ```
> > commit 005937fbe2a98fb83f0ade869025dc2636b4dad5 (HEAD -> main)
> > Author: Sherlock Holmes <sherlock@baker.street>
> > Date:   Thu Aug 22 10:14:07 2013 -0400
> >
> >    Make notes about colonel's fingerprints
> > ```
> {: .code-out}
>
> You can also reduce the quantity of information using the
> `--oneline` option:
>
> > <code-in-title>See only basic information</code-in-title>
> > ```
> > $ git log --oneline
> > ```
> {: .code-in}
>
> > <code-out-title></code-out-title>
> > ```
> > 005937f (HEAD -> main) Make notes about colonel's fingerprints
> > 34961b1 Add concerns about existence of motive for colonel
> > f22b25e Start notes on colonel as a base
> > ```
> {: .code-out}
>
> You can also combine the `--oneline` option with others. One useful
> combination adds `--graph` to display the commit history as a text-based
> graph and to indicate which commits are associated with the
> current `HEAD`, the current branch `main`, or
> [other Git references][git-references]:
>
> > <code-in-title>Extra options</code-in-title>
> > ```bash
> > $ git log --oneline --graph
> > ```
> {: .code-in}
> > <code-out-title></code-out-title>
> > ```
> > * 005937f (HEAD -> main) Make notes about colonel's fingerprints
> > * 34961b1 Add concerns about existence of motive for colonel
> > * f22b25e Start notes on colonel as a base
> > ```
> {: .code-out}
{: .tip}

> <tip-title>Directories </tip-title>
>
> Two important facts you should know about directories in Git.
>
> 1. Git does not track directories on their own, only files within them.
>    Try it for yourself:
>
>    > <code-in-title>Test directory tracking  </code-in-title>
>    > ```bash
>    > $ mkdir mysteries
>    > $ git status
>    > $ git add mysteries
>    > $ git status
>    > ```
>    {: .code-in}
>
>    > Note, our newly created empty directory `mysteries` does not appear in
>    the list of untracked files even if we explicitly add it (_via_ `git add`) to our
>    repository. This is the reason why you will sometimes see `.gitkeep` files
>    in otherwise empty directories. Unlike `.gitignore`, these files are not special
>    and their sole purpose is to populate a directory so that Git adds it to
>    the repository. In fact, you can name such files anything you like.
>
> 2. If you create a directory in your Git repository and populate it with files,
>    you can add all files in the directory at once by:
>
>     ```
>     git add <directory-with-files>
>     ```
>
>    Try it for yourself:
>    > <code-in-title>Add multiple files</code-in-title>
>    > ```bash
>    > $ touch mysteries/belgravia mysteries/reichenbachfall
>    > $ git status
>    > $ git add mysteries
>    > $ git status
>    > ```
>    {: .code-in}
>
>    Before moving on, we will commit these changes.
>
>    > <code-in-title>Commit the new changes</code-in-title>
>    > ```
>    > $ git commit -m "Add some initial thoughts on mysteries"
>    > ```
>    {: .code-in}
>
{: .tip}

To recap, when we want to add changes to our repository,
we first need to add the changed files to the staging area
(`git add`) and then commit the staged changes to the
repository (`git commit`):

![The Git Commit Workflow](../../images/bash-git/git-committing.svg)

# Let's put us to the test

> <question-title>Choosing a Commit Message</question-title>
>
> Which of the following commit messages would be most appropriate for the
> last commit made to `colonel.txt`?
>
> 1. "Changes"
> 2. "Added line 'Fingerprints on victims glasses.' to colonel.txt"
> 3. "Make notes about colonel's fingerprints"
>
> > <solution-title></solution-title>
> > Answer 1 is not descriptive enough, and the purpose of the commit is unclear;
> > and answer 2 is redundant to using "git diff" to see what changed in this commit;
> > but answer 3 is good: short, descriptive, and imperative.
> {: .solution}
{: .question }

> <question-title>Committing Changes to Git</question-title>
>
> Which command(s) below would save the changes of `myfile.txt`
> to my local Git repository?
>
> 1. ```
>    $ git commit -m "my recent changes"
>    ```
>
> 2. ```
>    $ git init myfile.txt
>    $ git commit -m "my recent changes"
>    ```
> 3. ```
>    $ git add myfile.txt
>    $ git commit -m "my recent changes"
>    ```
> 4. ```
>    $ git commit -m myfile.txt "my recent changes"
>    ```
>
> > <solution-title></solution-title>
> >
> > 1. Would only create a commit if files have already been staged.
> > 2. Would try to create a new repository.
> > 3. Is correct: first add the file to the staging area, then commit.
> > 4. Would try to commit a file "my recent changes" with the message myfile.txt.
> {: .solution}
{: .question }

> <question-title>Committing Multiple Files</question-title>
>
> The staging area can hold changes from any number of files
> that you want to commit as a single snapshot.
>
> 1. Add some text to `colonel.txt` noting your suspicion on judge Brown.
> 2. Create a new file `judge.txt` with your initial thoughts
> about the judge as a suspect.
> 3. Add changes from both files to the staging area,
> and commit those changes.
>
> > <solution-title></solution-title>
> >
> > The output below from `cat colonel.txt` reflects only content added during
> > this exercise. Your output may vary.
> >
> > First we make our changes to the `colonel.txt` and `judge.txt` files:
> > > <code-in-title>Edit `colonel.txt`</code-in-title>
> > > ```bash
> > > $ nano colonel.txt
> > > $ cat colonel.txt
> > > ```
> > {: .code-in}
> > > <code-out-title></code-out-title>
> > > ```
> > > Maybe judge Brown should also be considerable as a suspect.
> > > ```
> > {: .code-out}
> > > <code-in-title>Create and edit `judge.txt`</code-in-title>
> > > ```bash
> > > $ nano judge.txt
> > > $ cat judge.txt
> > > ```
> > {: .code-in}
> > > <code-out-title></code-out-title>
> > > ```
> > > Judge seems like a nice guy, but has a shady past.
> > > ```
> > {: .code-out}
> >
> > Now you can add both files to the staging area. We can do that in one line:
> > > <code-in-title>Add both files</code-in-title>
> > > ```bash
> > > $ git add colonel.txt judge.txt
> > > ```
> > {: .code-in}
> >
> > > <code-in-title>Or with multiple commands: </code-in-title>
> > > ```bash
> > > $ git add colonel.txt
> > > $ git add judge.txt
> > > ```
> > {: .code-in}
> >
> > Now the files are ready to commit. You can check that using `git status`. If you are ready to commit use:
> > > <code-in-title>Commit changes</code-in-title>
> > > ```bash
> > > $ git commit -m "Write plans to start a base on judge"
> > > ```
> > {: .code-in}
> > > <code-out-title></code-out-title>
> > > ```
> > > [main cc127c2]
> > >  Write plans to start a base on judge
> > >  2 files changed, 2 insertions(+)
> > >  create mode 100644 judge.txt
> > > ```
> > {: .code-out}
> {: .solution}
{: .question }

> <question-title>`bio` Repository</question-title>
>
> 1. Create a new Git repository on your computer called `bio`.
> 2. Write a three-line biography for yourself in a file called `me.txt`,
> commit your changes
> 3. Modify one line, add a fourth line
> 4. Display the differences
> between its updated state and its original state.
>
> > <solution-title></solution-title>
> >
> > If needed, move out of the `suspects` folder:
> > > <code-in-title>Change directory</code-in-title>
> > > ```bash
> > > $ cd ..
> > > ```
> > {: .code-in}
> >
> > Create a new folder called `bio` and 'move' into it:
> > > ### Create folder
> > > ```bash
> > > $ mkdir bio
> > > $ cd bio
> > > ```
> > {: .code-in}
> >
> > > <code-in-title>Initialise git:</code-in-title>
> > > ```bash
> > > $ git init
> > > ```
> > {: .code-in}
> >
> > Create your biography file `me.txt` using `nano` or another text editor.
> > Once in place, add and commit it to the repository:
> >
> > > <code-in-title>Create file and edit it</code-in-title>
> > > ```bash
> > > $ git add me.txt
> > > $ git commit -m "Add biography file"
> > > ```
> > {: .code-in}
> >
> > Modify the file as described (modify one line, add a fourth line).
> > To display the differences
> > between its updated state and its original state, use `git diff`:
> >
> > > <code-in-title>Display the differences</code-in-title>
> > > ```bash
> > > $ git diff me.txt
> > > ```
> > {: .code-in}
> {: .solution}
{: .question }

[commit-messages]: https://chris.beams.io/posts/git-commit/
[git-references]: https://git-scm.com/book/en/v2/Git-Internals-Git-References

# History Exploring

As we saw in the previous episode, we can refer to commits by their
identifiers.  You can refer to the _most recent commit_ of the working
directory by using the identifier `HEAD`.

We've been adding one line at a time to `colonel.txt`, so it's easy to track our
progress by looking, so let's do that using our `HEAD`s.  Before we start,
let's make a change to `colonel.txt`, adding yet another line.

> <code-in-title>Edit `colonel.txt`</code-in-title>
> ```bash
> $ nano colonel.txt
> $ cat colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> No alibi for the night of murder.
> No clear motive. Seems high unlikely.
> Fingerprints on victims glasses.
> Maybe judge Brown should also be considerable as a suspect.
> ```
{: .code-out}

Now, let's see what we get.

> <code-in-title>Display the changes</code-in-title>
> ```bash
> $ git diff HEAD colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> diff --git a/colonel.txt b/colonel.txt
> index b36abfd..0848c8d 100644
> --- a/colonel.txt
> +++ b/colonel.txt
> @@ -1,3 +1,4 @@
>  No alibi for the night of murder.
>  No clear motive. Seems high unlikely.
>  Fingerprints on victims glasses.
> +Maybe judge Brown should also be considerable as a suspect.
> ```
{: .code-out}

which is the same as what you would get if you leave out `HEAD` (try it).  The
real goodness in all this is when you can refer to previous commits.  We do
that by adding `~1`
(where "~" is "tilde", pronounced [**til**-d*uh*])
to refer to the commit one before `HEAD`.

> <code-in-title>Refer to previous commits</code-in-title>
> ```bash
> $ git diff HEAD~1 colonel.txt
> ```
{: .code-in}

If we want to see the differences between older commits we can use `git diff`
again, but with the notation `HEAD~1`, `HEAD~2`, and so on, to refer to them:

> <code-in-title>Refer to previous commits</code-in-title>
> ```bash
> $ git diff HEAD~3 colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> diff --git a/colonel.txt b/colonel.txt
> index df0654a..b36abfd 100644
> --- a/colonel.txt
> +++ b/colonel.txt
> @@ -1 +1,4 @@
>  No alibi for the night of murder.
> +No clear motive. Seems high unlikely.
> +Fingerprints on victims glasses.
> +Maybe judge Brown should also be considerable as a suspect.
> ```
{: .code-out}

We could also use `git show` which shows us what changes we made at an older commit as
well as the commit message, rather than the _differences_ between a commit and our
working directory that we see by using `git diff`.

> <code-in-title>`git show` and `git diff` differences </code-in-title>
> ```bash
> $ git show HEAD~3 colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
> Author: Sherlock Holmes <sherlock@baker.street>
> Date:   Thu Aug 22 09:51:46 2013 -0400
>
>     Start notes on colonel as a base
>
> diff --git a/colonel.txt b/colonel.txt
> new file mode 100644
> index 0000000..df0654a
> --- /dev/null
> +++ b/colonel.txt
> @@ -0,0 +1 @@
> +No alibi for the night of murder.
> ```
{: .code-out}

In this way,
we can build up a chain of commits.
The most recent end of the chain is referred to as `HEAD`;
we can refer to previous commits using the `~` notation,
so `HEAD~1`
means "the previous commit",
while `HEAD~123` goes back 123 commits from where we are now.

We can also refer to commits using
those long strings of digits and letters
that `git log` displays.
These are unique IDs for the changes,
and "unique" really does mean unique:
every change to any set of files on any computer
has a unique 40-character identifier.
Our first commit was given the ID
`f22b25e3233b4645dabd0d81e651fe074bd8e73b`,
so let's try this:

> <code-in-title>Display specific commit</code-in-title>
> ```bash
> $ git diff f22b25e3233b4645dabd0d81e651fe074bd8e73b colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> diff --git a/colonel.txt b/colonel.txt
> index df0654a..93a3e13 100644
> --- a/colonel.txt
> +++ b/colonel.txt
> @@ -1 +1,4 @@
>  No alibi for the night of murder.
> +No clear motive. Seems high unlikely.
> +Fingerprints on victims glasses.
> +Maybe judge Brown should also be considerable as a suspect.
> ```
{: .code-out}

That's the right answer,
but typing out random 40-character strings is annoying,
so Git lets us use just the first few characters (typically seven for normal size projects):

> <code-in-title>Shorter alternative</code-in-title>
> ```bash
> $ git diff f22b25e colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> diff --git a/colonel.txt b/colonel.txt
> index df0654a..93a3e13 100644
> --- a/colonel.txt
> +++ b/colonel.txt
> @@ -1 +1,4 @@
>  No alibi for the night of murder.
> +No clear motive. Seems high unlikely.
> +Fingerprints on victims glasses.
> +Maybe judge Brown should also be considerable as a suspect.
> ```
{: .code-out}

All right! So
we can save changes to files and see what we've changed. Now, how
can we restore older versions of things?
Let's suppose we change our mind about the last update to
`colonel.txt` (the "ill-considered change").

`git status` now tells us that the file has been changed,
but those changes haven't been staged:

> <code-in-title>Check</code-in-title>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
> Changes not staged for commit:
>   (use "git add <file>..." to update what will be committed)
>   (use "git checkout -- <file>..." to discard changes in working directory)
>
>     modified:   colonel.txt
>
> no changes added to commit (use "git add" and/or "git commit -a")
> ```
{: .code-out}

We can put things back the way they were
by using `git checkout`:

> <code-in-title>Restore with `git checkout`</code-in-title>
> ```bash
> $ git checkout HEAD colonel.txt
> $ cat colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> No alibi for the night of murder.
> No clear motive. Seems high unlikely.
> Fingerprints on victims glasses.
> ```
{: .code-out}

As you might guess from its name,
`git checkout` checks out (i.e., restores) an old version of a file.
In this case,
we're telling Git that we want to recover the version of the file recorded in `HEAD`,
which is the last saved commit.
If we want to go back even further,
we can use a commit identifier instead:

> <code-in-title>Restrore specific commit</code-in-title>
> ```bash
> $ git checkout f22b25e colonel.txt
> ```
{: .code-in}

> <code-in-title>Check</code-in-title>
> ```bash
> $ cat colonel.txt
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> No alibi for the night of murder.
> ```
{: .code-out}

> <code-in-title>Check</code-in-title>
> ```bash
> $ git status
> ```
{: .code-in}

> <code-out-title></code-out-title>
> ```
> On branch main
> Changes to be committed:
>   (use "git reset HEAD <file>..." to unstage)
>
>     modified:   colonel.txt
>
> ```
{: .code-out}

Notice that the changes are currently in the staging area.
Again, we can put things back the way they were
by using `git checkout`:

> <code-in-title>Restore </code-in-title>
> ```bash
> $ git checkout HEAD colonel.txt
> ```
{: .code-in}

> <warning-title>Don't Lose Your HEAD</warning-title>
>
> Above we used
>
> ```
> $ git checkout f22b25e colonel.txt
> ```
>
> to revert `colonel.txt` to its state after the commit `f22b25e`. But be careful!
> The command `checkout` has other important functionalities and Git will misunderstand
> your intentions if you are not accurate with the typing. For example,
> if you forget `colonel.txt` in the previous command.
>
> > <code-in-title>Error recipe</code-in-title>
> > ```bash
> > $ git checkout f22b25e
> > ```
> {: .code-in}
>
> > <code-out-title></code-out-title>
> > ```
> > Note: checking out 'f22b25e'.
> >
> > You are in 'detached HEAD' state. You can look around, make experimental
> > changes and commit them, and you can discard any commits you make in this
> > state without impacting any branches by performing another checkout.
> >
> > If you want to create a new branch to retain commits you create, you may
> > do so (now or later) by using -b with the checkout command again. Example:
> >
> > git checkout -b <new-branch-name>
> >
> > HEAD is now at f22b25e Start notes on colonel as a base
> > ```
> {: .code-out}
>
> The "detached HEAD" is like "look, but don't touch" here,
> so you shouldn't make any changes in this state.
> After investigating your repo's past state, reattach your `HEAD` with `git checkout main`.
{: .warning}

It's important to remember that
we must use the commit number that identifies the state of the repository
*before* the change we're trying to undo.
A common mistake is to use the number of
the commit in which we made the change we're trying to discard.
In the example below, we want to retrieve the state from before the most
recent commit (`HEAD~1`), which is commit `f22b25e`:

![Git Checkout](../../images/bash-git/git-checkout.svg)

So, to put it all together,
here's how Git works in cartoon form:

![https://figshare.com/articles/How_Git_works_a_cartoon/1328266](../../images/bash-git/git_staging.svg)

> <tip-title>Simplifying the Common Case</tip-title>
>
> If you read the output of `git status` carefully,
> you'll see that it includes this hint:
>
> ```
> (use "git checkout -- <file>..." to discard changes in working directory)
> ```
>
> As it says,
> `git checkout` without a version identifier restores files to the state saved in `HEAD`.
> The double dash `--` is needed to separate the names of the files being recovered
> from the command itself:
> without it,
> Git would try to use the name of the file as the commit identifier.
{: .tip}

The fact that files can be reverted one by one
tends to change the way people organize their work.
If everything is in one large document,
it's hard (but not impossible) to undo changes to the introduction
without also undoing changes made later to the conclusion.
If the introduction and conclusion are stored in separate files,
on the other hand,
moving backward and forward in time becomes much easier.

> <question-title>Recovering Older Versions of a File</question-title>
>
> Jennifer has made changes to the Python script that she has been working on for weeks, and the
> modifications she made this morning "broke" the script and it no longer runs. She has spent
> ~ 1hr trying to fix it, with no luck...
>
> Luckily, she has been keeping track of her project's versions using Git! Which commands below will
> let her recover the last committed version of her Python script called
> `data_cruncher.py`?
>
> 1. `$ git checkout HEAD`
>
> 2. `$ git checkout HEAD data_cruncher.py`
>
> 3. `$ git checkout HEAD~1 data_cruncher.py`
>
> 4. `$ git checkout <unique ID of last commit> data_cruncher.py`
>
> 5. Both 2 and 4
>
>
> > <solution-title></solution-title>
> >
> > The answer is (5)-Both 2 and 4.
> >
> > The `checkout` command restores files from the repository, overwriting the files in your working
> > directory. Answers 2 and 4 both restore the *latest* version *in the repository* of the file
> > `data_cruncher.py`. Answer 2 uses `HEAD` to indicate the *latest*, whereas answer 4 uses the
> > unique ID of the last commit, which is what `HEAD` means.
> >
> > Answer 3 gets the version of `data_cruncher.py` from the commit *before* `HEAD`, which is NOT
> > what we wanted.
> >
> > Answer 1 can be dangerous! Without a filename, `git checkout` will restore **all files**
> > in the current directory (and all directories below it) to their state at the commit specified.
> > This command will restore `data_cruncher.py` to the latest commit version, but it will also
> > restore *any other files that are changed* to that version, erasing any changes you may
> > have made to those files!
> > As discussed above, you are left in a *detached* `HEAD` state, and you don't want to be there.
> {: .solution}
{: .question }

> <question-title>Reverting a Commit</question-title>
>
> Jennifer is collaborating with colleagues on her Python script.  She
> realizes her last commit to the project's repository contained an error, and
> wants to undo it.  Jennifer wants to undo correctly so everyone in the project's
> repository gets the correct change. The command `git revert [erroneous commit ID]` will create a
> new commit that reverses the erroneous commit.
>
> The command `git revert` is
> different from `git checkout [commit ID]` because `git checkout` returns the
> files not yet committed within the local repository to a previous state, whereas `git revert`
> reverses changes committed to the local and project repositories.
>
> Below are the right steps and explanations for Jennifer to use `git revert`,
> what is the missing command?
> 1. `________ # Look at the git history of the project to find the commit ID`
>
> 2. Copy the ID (the first few characters of the ID, e.g. 0b1d055).
>
> 3. `git revert [commit ID]`
>
> 4. Type in the new commit message.
>
> 5. Save and close
>
>
> > <solution-title></solution-title>
> >
> > The command `git log` lists project history with commit IDs.
> >
> > The command `git show HEAD` shows changes made at the latest commit, and lists
> > the commit ID; however, Jennifer should double-check it is the correct commit, and no one
> > else has committed changes to the repository.
> {: .solution}
{: .question }

> <question-title>Understanding Workflow and History</question-title>
>
> What is the output of the last command in
>
> ```
> $ cd suspects
> $ echo "judge has unresolved military issues" > judge.txt
> $ git add judge.txt
> $ echo "judge has enemies in the city" >> judge.txt
> $ git commit -m "Comment on judge as a suspect"
> $ git checkout HEAD judge.txt
> $ cat judge.txt #this will print the contents of judge.txt to the screen
> ```
>
> 1. ```
>    judge has enemies in the city
>    ```
> 2. ```
>    judge has unresolved military issues
>    ```
> 3. ```
>    judge has unresolved military issues
>    judge has enemies in the city
>    ```
> 4. ```
>    Error because you have changed judge.txt without committing the changes
>    ```
>
> > <solution-title></solution-title>
> >
> > The answer is 2.
> >
> > The command `git add judge.txt` places the current version of `judge.txt` into the staging area.
> > The changes to the file from the second `echo` command are only applied to the working copy,
> > not the version in the staging area.
> >
> > So, when `git commit -m "Comment on judge as an unsuitable base"` is executed,
> > the version of `judge.txt` committed to the repository is the one from the staging area and
> > has only one line.
> >
> >  At this time, the working copy still has the second line (and
> >  `git status` will show that the file is modified). However, `git checkout HEAD judge.txt`
> >  replaces the working copy with the most recently committed version of `judge.txt`.
> >
> >  So, `cat judge.txt` will output
> >  ```
> >  judge has unresolved military issues.
> >  ```
> {: .solution}
{: .question }

> <question-title>Checking Understanding of `git diff</question-title>
>
> Consider this command: `git diff HEAD~9 colonel.txt`. What do you predict this command
> will do if you execute it? What happens when you do execute it? Why?
>
> Try another command, `git diff [ID] colonel.txt`, where [ID] is replaced with
> the unique identifier for your most recent commit. What do you think will happen,
> and what does happen?
{: .question }

> <question-title>Getting Rid of Staged Changes</question-title>
>
> `git checkout` can be used to restore a previous commit when unstaged changes have
> been made, but will it also work for changes that have been staged but not committed?
> Make a change to `colonel.txt`, add that change, and use `git checkout` to see if
> you can remove your change.
{: .question }

> <question-title>Explore and Summarize Histories</question-title>
>
> Exploring history is an important part of Git, and often it is a challenge to find
> the right commit ID, especially if the commit is from several months ago.
>
> Imagine the `suspects` project has more than 50 files.
> You would like to find a commit that modifies some specific text in `colonel.txt`.
> When you type `git log`, a very long list appeared.
> How can you narrow down the search?
>
> > <solution-title></solution-title>
> > Recall that the `git diff` command allows us to explore one specific file,
> > e.g., `git diff colonel.txt`. We can apply a similar idea here.
> >
> > ```
> > $ git log colonel.txt
> > ```
> >
> >
> > Unfortunately some of these commit messages are very ambiguous, e.g., `update files`.
> > How can you search through these files?
> >
> > Both `git diff` and `git log` are very useful and they summarize a different part of the history
> > for you.
> > Is it possible to combine both? Let's try the following:
> >
> > ```
> > $ git log --patch colonel.txt
> > ```
> >
> >You should get a long list of output, and you should be able to see both commit messages and
> >the difference between each commit.
> >
> >Question: What does the following command do?
> >
> > ```
> > $ git log --patch HEAD~9 *.txt
> > ```
> {: .solution}
{: .question }
