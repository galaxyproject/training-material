---
title: "GTN ADR: Why Jekyll and not another Static Site Generator (SSG)"
area: adr
box_type: tip
layout: faq
contributors: [hexylena, shiltemann, bebatut]
---

[FAQ: What is an ADR?]({% link topics/contributing/faqs/gtn-adr.md %})

## Context and Problem Statement

We needed a static site generator for the GTN, one had to be chosen. We chose Jekyll because of it's good integration with GitHub and GitHub Pages. Over time our requirements have changed but we still need one SSG.

## Decision Drivers

* Must be easy for contributors to setup and use
* Needs to be relatively performant (full rebuilds may not take more than 2 minutes.)
* Must allow us to develop custom plugins

## Considered Options

* Jekyll
* Hugo
* A javascript option
* Another SSG.

## Decision Outcome

Chosen option: "Jekyll", because of the amount of time and effort we have sunk into it over the years has made it a good platform for us, despite limitations.

Over time we have invested heavily into Jekyll, any choice to switch *must* take that into consideration. Consider the following output of `scc _plugins bin/`


Language   | Files  | Lines | Blanks | Comments | Code  | Complexity
--------   | -----  | ----- | ------ | -------- | ----  | ----------
YAML       | 117    | 9830  | 71     | 33       | 9726  | 0
Ruby       | 90     | 14471 | 1795   | 2617     | 10059 | 1163
JSON       | 48     | 3075  | 0      | 0        | 3075  | 0
Python     | 24     | 3693  | 284    | 272      | 3137  | 310
Shell      | 21     | 1529  | 175    | 262      | 1092  | 84
JavaScript | 5      | 299   | 38     | 19       | 242   | 48
Markdown   | 4      | 76    | 19     | 0        | 57    | 0
Dockerfile | 2      | 60    | 15     | 1        | 44    | 14
Plain      | Text   | 2     | 18     | 0        | 0     | 18         | 0
BASH       | 1      | 51    | 8      | 4        | 39    | 1
CSS        | 1      | 3     | 0      | 0        | 3     | 0
Docker     | ignore | 1     | 1      | 0        | 0     | 1          | 0
gitignore  | 1      | 123   | 0      | 0        | 123   | 0
Total      | 317    | 33229 | 2405   | 3208     | 27616 | 1620

- Estimated Cost to Develop (organic) $880,671
- Estimated Schedule Effort (organic) 13.11 months
- Estimated People Required (organic) 5.97
- Processed 1081253 bytes, 1.081 megabytes (SI)

This is a *lot* of code that would need to be rewritten if another language was ever chosen.

The YAML comprises our Kwalify Schemas. There is a good argument for moving to JSON Schema instead. The Ruby however is the bulk of the code that would need to be rewritten. It does a significant number of complex things:

- collecting and collating files off disk / in Jekyll's Page model into "Learning Materials", very large objects with hundreds of properties that are used to render each and every template.
- Generating hundreds of pages with a multitude of calculated properties. These would all need to be hand translated.

Additionally any layouts would need to be rewritten from our existing Liquid templates. Note that this is not the full set of templates.


Language | Files | Lines | Blanks | Comments | Code | Complexity
-------- | ----- | ----- | ------ | -------- | ---- | ----------
HTML     | 69    | 5937  | 830    | 96       | 5011 | 0
Markdown | 4     | 125   | 1      | 0        | 124  | 0
Total    | 73    | 6062  | 831    | 96       | 5135 | 0

- Estimated Cost to Develop (organic) $150,543
- Estimated Schedule Effort (organic) 6.70 months
- Estimated People Required (organic) 2.00


### Consequences

* Good, because it works well for us and has scaled sufficiently to an incredible number of output pages (~7k html/22k files in a full GTN production deployment.) with acceptable build times (<5 minutes in prod, most of the action execution is taken up by contacting other servers, dependencies, and uploading the results.)
* Good, because it has a well supported ecosystem of plugins we can leverage for common tasks
* Good, because we can easily write our own plugins for many tasks.
* Bad, because we it remains difficult to install
* Bad, because people must know Ruby and very few people do (but it isn't that hard to learn!)

## Pros and Cons of the Options

### Hugo

* Good, because it would be a single binary, easier to install
* Bad, because plugins do not exist, it does not have a way to hook the internals and work with them which we use extensively.
* Bad, because what plugins do exist, only exist as 'shortcodes' that are written in Go templates which are not as powerful as Ruby.

### A JavaScript option

* Good, because we could re-use code from other places
* Bad, because the average lifetime of a JavaScript SSG is maybe one year.
* Bad, because they are also quite slow on average (Hub compile times are on the order of 10 minutes.)
