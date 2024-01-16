---
title: What is Gene Ontology (GO)?
box_type: question
layout: faq
contributors: [combesf]
---

A very commonly used way of specifying these sets is to gather genes/proteins that share
the same Gene Ontology (GO) term, as specified by the Gene Ontology Consortium.

The GO project provides an ontology that describes gene products and their relations
in three non-overlapping domains of molecular biology, namely “Molecular Function”,
“Biological Process”, and “Cellular Component”. Genes/proteins are annotated by one
or several GO terms, each composed of a label, a definition and a unique identifier.
GO terms are organized within a classification scheme that supports relationships, and
formalized by a hierarchical structure that forms a directed acyclic graph (DAG). In
such a graph is used the notions of child and parent, where a child inherits from one
or multiple parents, child class having a more specific annotation than parent class
(e.g. “glucose metabolic process” inherits from “hexose metabolic” parent term which
itself inherits from “monosaccharide metabolic process” etc.). In this graph, each node corresponds to a GO term composed of genes/proteins sharing the same annotation, while
directed edges between nodes represents their relation (e.g. ‘is a’, ‘part of’) and
their roles in the hierarchy (i.e. parent and child).

[Further reading](http://www.geneontology.org/page/introduction-go-resource)

