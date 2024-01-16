---
title: My Rscript tool generates a strange R error on STDOUT about an invalid operation on a closure called 'args' ?
box_type: question
layout: faq
contributors: [fubar2]
---

Did your code declare the `args vector` with something like `args = commandArgs(trailingOnly=TRUE)` before it tried to access args[1] ?
See the plotter tool for a sample

