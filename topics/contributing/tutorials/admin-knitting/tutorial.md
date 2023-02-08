---
layout: tutorial_hands_on

title: "Updating diffs in admin training"
questions:
  - "How does it work?"
objectives:
  - "Update some text in an earlier commit"
time_estimation: "5m"
key_points:
  - "The knitting script works OK?"
  - "Learn to rebase now, if you haven't yet"
  - "git am is not fun"
  - "If you have issues, ping me."
contributors:
  - hexylena
---

# Updating Diffs


The admin training was recently converted completely to use diffs. This was
done so we could build a git repository, and time travel to any point within
admin training, in support of students and testing.

Here's some common editing scenarios to help you out.

> <agenda-title></agenda-title>
>
> In this tutorial, we will:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## How it Works

Tutorials are written with diffs that look like this:

    ```diff
    --- /dev/null
    +++ b/ansible.cfg
    @@ -0,0 +1,4 @@
    +[defaults]
    +interpreter_python = /usr/bin/python3
    +inventory = hosts
    +retry_files_enabled = false
    ```
    {: data-commit="Add ansible.cfg"}

The `data-commit` will be turned into the commit message, so, you can write something nice here. Otherwise it's a pretty standard diff. For tutorials which are written in this way, we can extract all of the diffs.

### Frogging

In knitting you sometimes need to rip out the stitches, a process sometimes called frogging.

> If they are not secured, the loops of a knitted course will come undone when their yarn is pulled; this is known as ripping out, unravelling knitting, or humorously, frogging (because you 'rip it', this sounds like a frog croaking: 'rib-bit'). [via Wikipedia](https://en.wikipedia.org/wiki/Knitting)
{: .quote}

```console
$ python bin/knit-frog.py topics/admin/tutorials/singularity/tutorial.md /tmp/03-singularity
/tmp/03-singularity-commit-0000-add-golang-and-singulary-ansible-roles.patch
/tmp/03-singularity-commit-0001-configure-golang-and-singularity.patch
/tmp/03-singularity-commit-0002-add-the-roles-to-the-playbook.patch
/tmp/03-singularity-commit-0003-configure-the-container-and-dependency-resolvers.patch
/tmp/03-singularity-commit-0004-configure-the-dependency-resolvers.patch
/tmp/03-singularity-commit-0005-configure-the-container-resolver.patch
/tmp/03-singularity-commit-0006-update-the-job-conf-xml-with-singularity-destination.patch
```

This rips out the diffs and writes them into patch files.

### Knitting

With a series of patch files, we have a similar script to knit them back together.

```console
$ python bin/knit.py topics/admin/tutorials/singularity/tutorial.md --patches /tmp/03-singularity*
```

This takes the patch files and lines them up in the tutorial. It is **not intelligent**. You should have the same number of patch files as you do diffs with `data-commit`. Up to you to manage that.

## Automation

But mostly you don't need to care about the `knit.py` or `knit-frog.py`, we have `knit-automated.sh` which takes care of calling those properly. You can call it with `export` to write the patches into `/tmp/git-gat/`

```
$ ./bin/knit-automated.sh export
/tmp/git-gat/0-ansible-galaxy-commit-0000-add-requirements.patch
/tmp/git-gat/0-ansible-galaxy-commit-0001-add-ansible-cfg.patch
/tmp/git-gat/0-ansible-galaxy-commit-0002-add-hosts.patch
/tmp/git-gat/0-ansible-galaxy-commit-0003-add-initial-group-variables-file.patch
/tmp/git-gat/0-ansible-galaxy-commit-0004-add-initial-galaxy-playbook.patch
/tmp/git-gat/0-ansible-galaxy-commit-0005-add-pip--miniconda--galaxy-to-playbook.patch
/tmp/git-gat/0-ansible-galaxy-commit-0006-configure-miniconda-and-galaxy.patch
[...]
/tmp/git-gat/6-pulsar-commit-0000-add-requirements.patch
/tmp/git-gat/6-pulsar-commit-0001-configure-rabbitmq.patch
/tmp/git-gat/6-pulsar-commit-0002-add-role.patch
/tmp/git-gat/6-pulsar-commit-0003-add-pulsar-group-variables.patch
/tmp/git-gat/6-pulsar-commit-0004-add-pulsar-host.patch
/tmp/git-gat/6-pulsar-commit-0005-add-pulsar-playbook.patch
/tmp/git-gat/6-pulsar-commit-0006-add-pulsar-plugin.patch
/tmp/git-gat/6-pulsar-commit-0007-add-pulsar-destination.patch
/tmp/git-gat/6-pulsar-commit-0008-send-bwa-and-bwa-mem-to-pulsar.patch
```

And there you can build a repository from them!

```
[hxr@mk:/tmp/git-gat]$ git init
[hxr@mk:/tmp/git-gat]$ git am -3 *
Applying: admin/ansible-galaxy/0000: Add requirements
applying to an empty history
Applying: admin/ansible-galaxy/0001: Add ansible.cfg
Applying: admin/ansible-galaxy/0002: Add hosts
Applying: admin/ansible-galaxy/0003: Add initial group variables file
Applying: admin/ansible-galaxy/0004: Add initial galaxy playbook
Applying: admin/ansible-galaxy/0005: Add pip, miniconda, galaxy to playbook
[...]
Applying: admin/pulsar/0004: Add pulsar host
Applying: admin/pulsar/0005: Add pulsar playbook
Applying: admin/pulsar/0006: Add pulsar plugin
Applying: admin/pulsar/0007: Add pulsar destination
Applying: admin/pulsar/0008: Send bwa and bwa-mem to pulsar
```

Et voila!

```console
$ rm *.patch
$ tree
.
├── ansible.cfg
├── files
│   └── galaxy
│       ├── config
│       │   ├── dependency_resolvers_conf.xml
│       │   └── tool_destinations.yml
│       ├── dynamic_job_rules
│       │   ├── map_resources.py
│       │   └── my_rules.py
│       └── tools
│           └── testing.xml
├── galaxy.yml
├── group_vars
│   ├── all.yml
│   ├── galaxyservers.yml
│   └── pulsarservers.yml
├── hosts
├── pulsar.yml
├── requirements.yml
└── templates
    ├── galaxy
    │   └── config
    │       ├── container_resolvers_conf.xml.j2
    │       ├── job_conf.xml.j2
    │       └── job_resource_params_conf.xml.j2
    └── nginx
        ├── galaxy.j2
        └── redirect-ssl.j2

10 directories, 18 files
```

When you're done editing, you can run the `./bin/knit-automated.sh import`, it will

- extract the history as patches (using `git format-patch`)
- knit them into the tutorials

> <warning-title>Danger: Don't mess with the commit names at this step</warning-title>
> They're used for figuring out which commit a patch file belongs to. Sorry!
{: .warning}


# Git-Gat Knick-Knacks: Some Common Scenarios

## I want to update the text of a specific commit

> <question-title>Example problem</question-title>
> I want to edit the commit titled "Add production facing vars" to add an addtional variable in admin/ansible-galaxy.
{: .question}

> <hands-on-title>Make the change you want to see in the GAT</hands-on-title>
> 1. Start by exporting
>
>    ```
>    $ bash bin/knit-automated.sh export
>    ```
>
> 2. Build the repository, and apply the patches to get the current state.
>
>    ```
>    [hxr@mk:/tmp/git-gat]$ git init
>    [hxr@mk:/tmp/git-gat]$ git am -3 *.patch
>    Applying: admin/ansible-galaxy/0000: Add requirements
>    applying to an empty history
>    Applying: admin/ansible-galaxy/0001: Add ansible.cfg
>    Applying: admin/ansible-galaxy/0002: Add hosts
>    Applying: admin/ansible-galaxy/0003: Add initial group variables file
>    Applying: admin/ansible-galaxy/0004: Add initial galaxy playbook
>    Applying: admin/ansible-galaxy/0005: Add pip, miniconda, galaxy to playbook
>    Applying: admin/ansible-galaxy/0006: Configure miniconda and galaxy
>    Applying: admin/ansible-galaxy/0007: Configure galaxy config
>    [...]
>    [hxr@mk:/tmp/git-gat]$ rm *.patch
>    ```
>
> 3. Find the commit:
>
>    ```
>    $ git log --oneline | grep "Add production facing vars"
>    2f5b416 admin/ansible-galaxy/0022: Add production facing vars
>    ```
> 4. Edit rebase one upstream of that commit
>
>    ```
>    $ git rebase -i 2f5b416~
>    ```
>
>    Mark it for editing
>
>    ```
>    edit 2f5b416 admin/ansible-galaxy/0022: Add production facing vars
>    pick a352bec admin/ansible-galaxy/0023: Add nginx x-accel-redir and g-i-g webhook config to nginx
>    pick 59f9ff0 admin/singularity/0000: Add golang and singulary ansible roles
>    pick f38ddbf admin/singularity/0001: Configure golang and singularity
>    ```
>
>    It should stop there
>
>    ```
>    Stopped at 2f5b416...  admin/ansible-galaxy/0022: Add production facing vars
>    You can amend the commit now, with
>
>      git commit --amend
>
>    Once you are satisfied with your changes, run
>
>      git rebase --continue
>    ```
>
> 5. Open your file, make your change. Here I'm showing the change I made with `git diff`
>
>    ```
>    $ git diff | cat
>    diff --git a/group_vars/galaxyservers.yml b/group_vars/galaxyservers.yml
>    index 15fbf1e..8fe4343 100644
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -53,6 +53,8 @@ galaxy_config:
>         expose_user_name: true
>         expose_dataset_path: true
>         expose_potentially_sensitive_job_metrics: true
>    +    # NFS workarounds
>    +    retry_job_output_collection: 3
>         # Debugging
>         cleanup_job: onsuccess
>         allow_user_impersonation: true
>    ```
>
>
> 6. The status probably looks like this, the file that was changed in that commit is marked as modified
>
>    ```
>    $ git status
>    interactive rebase in progress; onto 7284036
>    Last command done (1 command done):
>       edit 2f5b416 admin/ansible-galaxy/0022: Add production facing vars
>    Next commands to do (46 remaining commands):
>       pick 3c6a72e admin/ansible-galaxy/0023: Add nginx x-accel-redir and g-i-g webhook config to nginx
>       pick a893955 admin/singularity/0000: Add golang and singulary ansible roles
>      (use "git rebase --edit-todo" to view and edit)
>    You are currently editing a commit while rebasing branch 'main' on '7284036'.
>      (use "git commit --amend" to amend the current commit)
>      (use "git rebase --continue" once you are satisfied with your changes)
>
>    Changes to be committed:
>      (use "git restore --staged <file>..." to unstage)
>            modified:   group_vars/galaxyservers.yml
>    ```
>
> 7. Ammend the commit, and continue the rebase
>
>    ```
>    $ git commit --amend
>    [detached HEAD ed290de] admin/ansible-galaxy/0022: Add production facing vars
>     Author: The Galaxy Training Network <galaxytrainingnetwork@gmail.com>
>     Date: Mon Feb 15 14:06:56 2021 +0100
>     1 file changed, 22 insertions(+)
>    $ git rebase --continue
>    Successfully rebased and updated refs/heads/main.
>    ```
>
> 8. Import the changes back in the GTN repository
>
>    ```
>    ./bin/knit-automated.sh import
>    ```
>
{: .hands_on}


The above shows commands with outputs, the rest will be more abbreviated

## Upgrading an existing tutorial

> <question-title>Example problem</question-title>
> I have a new tutorial in the schedule and I want to upgrade it to use this feature
{: .question}

> <hands-on-title>Make the change you want to see in the GAT</hands-on-title>
>
> 1. Start by exporting
>
>    ```
>    $ bash bin/knit-automated.sh export
>    ```
>
> 2. Build the repository... **up to the point where your tuto is**
>
>    ```
>    [hxr@mk:/tmp/git-gat]$ git init
>    [hxr@mk:/tmp/git-gat]$ git am -3 0* 1* 2*
>    Applying: admin/ansible-galaxy/0000: Add requirements
>    applying to an empty history
>    Applying: admin/ansible-galaxy/0001: Add ansible.cfg
>    Applying: admin/ansible-galaxy/0002: Add hosts
>    Applying: admin/ansible-galaxy/0003: Add initial group variables file
>    Applying: admin/ansible-galaxy/0004: Add initial galaxy playbook
>    Applying: admin/ansible-galaxy/0005: Add pip, miniconda, galaxy to playbook
>    Applying: admin/ansible-galaxy/0006: Configure miniconda and galaxy
>    Applying: admin/ansible-galaxy/0007: Configure galaxy config
>    [...]
>    ```
>
> 3. Go diff-by-diff in your tutorial, making the change, and then making a new commit on top for each diff. For each of these, add the `{: data-commit="Something usful"}` below the commit.
>
>    Do they have a `{ % endraw % }`? then move that above the backticks.
>
> 4. Format them as patches
>
>    ```
>    git format-patch <the last commit before you started committing>
>    ```
>
> 5. Knit them together
>
>    ```
>    python bin/knit.py topics/admin/tutorials/<yourtuto>/tutorial.md --patches /tmp/<the-patches-you-exported>*
>    ```
>
> 6. That's part one done! Now we need to finish the rest of the patches.
>
> 7. Update `./bin/knit-automated.sh` to include your tutorial in the right order.
>
> 8. Go through an export-import cycle
>
>    ```
>    $ rm -rf /tmp/git-gat/
>    $ ./bin/knit-automated.sh export
>    $ cd /tmp/git-gat/
>    $ git am -3 *.patch # Resolve any conflicts
>    $ cd -
>    $ ./bin/knit-automated.sh import
>    ```
>
> 9. Done!
{: .hands_on}

## I want to add a new commit

> <hands-on-title>Make the change you want to see in the GAT</hands-on-title>
>
> 1. Add your commit where it belongs in the tutorial
>
>    <pre>
>    {% raw %}
>    { raw }
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -22,3 +22,7 @@
>       version: 0.0.2
>     - src: galaxyproject.slurm
>       version: 0.1.3
>    +- name: usegalaxy_eu.rabbitmq
>    +  version: 0.1.0
>    +- src: galaxyproject.pulsar
>    +  version: 1.0.7
>    { endraw }
>    ```
>    {: data-commit="Add requirements"}
>    {% endraw %}
>    </pre>
>
>    Use your imagination to pretend the raw/endraw tags are correct with the missing percent signs.
>
> 2. Now in the above diff we've written some --- and +++ and more metadata we don't really have. So: just make it up, the commit will fail to apply, and we'll go fix it. You can really write *anything* here you want, and the diff will fail which is what we want.
>
> 3. Now, go apply the commits
>
>    ```
>    $ ./bin/knit-automated.sh export
>    $ cd /tmp/git-gat/
>    $ git am -3 *.patch # Resolve any conflicts
>    ```
>
> 4. This will fail. Here's how to fix it:
>
>    1. Run `git am --show-current-patch` to see the above diff you've inserted (or the otherwise incorrect diff that's a downstream result of your new change).
>    2. Make the correct change to the file on disk.
>    3. `git add <path>`
>    3. `git am --continue`
>
> 5. You'll probably have to do this quite a few more times, diffs have a habit of affecting downstream diffs (e.g. changing a version number or being included in the diff context.)
>
> 6. Once you've finished all of them, you'll just `./bin/knit-automated.sh import` and commit the resulting changes!
{: .hands_on}

The alternative to the above approach is

1. Export
2. git rebase and make the corresponding change in the history
3. Let the rebase finish
4. Insert a 'dummy' diff (raw/open code+diff/endraw/close code/data-commit tag) without any content and then run `import` which should insert it in the correct place.

# Errors that can occur and how to fix them


## sha1 information is lacking or useless

This is because there's a mismatch in the commits you're trying to apply and what the merge commit is seeing. Usually you won't notice. **fix by rebasing, and roundtripping the commits.**

```
Applying: admin/pulsar/0012: Send bwa and bwa-mem to pulsar
Tagging step-17
Committing 18-gxadmin-commit-0000-add-requirement.patch
18-gxadmin-commit-0002-add-the-gxadmin-role.patch
Applying: admin/gxadmin/0000: Add requirement
error: sha1 information is lacking or useless (requirements.yml).
error: could not build fake ancestor
Patch failed at 0001 admin/gxadmin/0000: Add requirement
hint: Use 'git am --show-current-patch' to see the failed patch
When you have resolved this problem, run "git am --continue".
If you prefer to skip this patch, run "git am --skip" instead.
To restore the original branch and stop patching, run "git am --abort".
```

Here you grumble at how picky it is.

```
$ git am --show-current-patch
From: The Galaxy Training Network <galaxytrainingnetwork@gmail.com>
Date: Mon, 15 Feb 2021 14:06:56 +0100
Subject: admin/gxadmin/0000: Add requirement


--- a/requirements.yml
+++ b/requirements.yml
@@ -26,3 +26,5 @@
   version: 0.1.0
 - src: galaxyproject.pulsar
   version: 1.0.6
+- src: galaxyproject.gxadmin
+  version: 0.0.8
--
2.25.1
```

That looks right, doesn't it? But if we grep we'll see

```
$ ag -Q '1.0.6' topics/admin
topics/admin/tutorials/gxadmin/tutorial.md
55:>       version: 1.0.6

topics/admin/tutorials/monitoring/tutorial.md
80:>       version: 1.0.6

topics/admin/tutorials/pulsar/tutorial.md
149:>    +  version: 1.0.6
```

again right, but is it? No. Github is testing the merge commit. Rebase your commits, and now you'll see:


```
$ ag -Q '1.0.6' topics/admin
topics/admin/tutorials/gxadmin/tutorial.md
55:>       version: 1.0.6

topics/admin/tutorials/monitoring/tutorial.md
80:>       version: 1.0.6

$ ag -Q '1.0.7' topics/admin
topics/admin/tutorials/pulsar/tutorial.md
149:>    +  version: 1.0.7
```

Aha, two different version. This is what breaks the hashes (hence the error message.)

Roundtripping the commits and fixing up everywhere there was a single character difference in the diff context (UGH) will fix it.
