---
title: How to read a Diff
area: diffs
box_type: tip
layout: faq
contributors: [hexylena]
---

If you haven't worked with diffs before, this can be something quite new or different.

If we have two files, let's say a grocery list, in two files. We'll call them 'a' and 'b'.


> > <code-in-title>Old</code-in-title>
> > ```
> > $ cat old
> > 🍎
> > 🍐
> > 🍊
> > 🍋
> > 🍒
> > 🥑
> > ```
> {: .code-in}
>
> > <code-out-title>New</code-out-title>
> > ```
> > $ cat new
> > 🍎
> > 🍐
> > 🍊
> > 🍋
> > 🍍
> > 🥑
> > ```
> {: .code-out}
{: .code-2col}

We can see that they have some different entries. We've removed 🍒 because they're awful, and replaced them with an 🍍

Diff lets us compare these files

```bash
$ diff old new
5c5
< 🍒
---
> 🍍
```

Here we see that 🍒 is only in a, and 🍍 is only in b. But otherwise the files are identical.

There are a couple different formats to diffs, one is the 'unified diff'

```bash
$ diff -U2 old new
--- old	2022-02-16 14:06:19.697132568 +0100
+++ new	2022-02-16 14:06:36.340962616 +0100
@@ -3,4 +3,4 @@
 🍊
 🍋
-🍒
+🍍
 🥑
```

This is basically what you see in the training materials which gives you a lot of context about the changes:

- `--- old` is the 'old' file in our view
- `+++ new` is the 'new' file
- @@ these lines tell us where the change occurs and how many lines are added or removed.
- Lines starting with a - are removed from our 'new' file
- Lines with a + have been added.

So when you go to apply these diffs to your files in the training:

1. Ignore the header
2. Remove lines starting with - from your file
3. Add lines starting with + to your file

The other lines (🍊/🍋 and 🥑) above just provide "context", they help you know where a change belongs in a file, but **should not be edited** when you're making the above change. Given the above diff, you would find a line with a 🍒, and replace it with a 🍍

#### Added & Removed Lines

Removals are very easy to spot, we just have removed lines

```bash
--- old	2022-02-16 14:06:19.697132568 +0100
+++ new	2022-02-16 14:10:14.370722802 +0100
@@ -4,3 +4,2 @@
 🍋
 🍒
-🥑
```

And additions likewise are very easy, just add a new line, between the other lines in your file.

```bash
--- old	2022-02-16 14:06:19.697132568 +0100
+++ new	2022-02-16 14:11:11.422135393 +0100
@@ -1,3 +1,4 @@
 🍎
+🍍
 🍐
 🍊
```

#### Completely new files

Completely new files look a bit different, there the "old" file is `/dev/null`, the empty file in a Linux machine.

```bash
$ diff -U2 /dev/null old
--- /dev/null	2022-02-15 11:47:16.100000270 +0100
+++ old	2022-02-16 14:06:19.697132568 +0100
@@ -0,0 +1,6 @@
+🍎
+🍐
+🍊
+🍋
+🍒
+🥑
```

And removed files are similar, except with the new file being /dev/null

```bash
--- old	2022-02-16 14:06:19.697132568 +0100
+++ /dev/null	2022-02-15 11:47:16.100000270 +0100
@@ -1,6 +0,0 @@
-🍎
-🍐
-🍊
-🍋
-🍒
-🥑
```
