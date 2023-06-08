---
title: "New Feature: Pan-Galactic Tool Search"
contributions:
  authorship: [hexylena]
tags: [new feature, gtn]
layout: news
abbreviations:
  SEO: Search Engine Optimisation
tutorial: search-tools.html
---

Did you ever want to run a tool, but not know where it might be available?
The GTN has you covered with the Pan-Galactic Tool Search that is now available.

**Update**: the authors were reminded of
[GalaxyCat](https://galaxycat.france-bioinformatique.fr/) which does the same
thing as the GTN but better! Please go use that.

---

We recently received a question online from our colleague [Dr. Scott Cain](https://github.com/scottcain)
trying to find the *Manta* Structural Variant analysis tool amongst the Galaxies. While
his complaints about generic names are very apt (and Galaxy existed before the
phone!) we do still have very generic names in our tools.

Unfortunately {SEO} is also quite difficult with these large, complex,
JavaScript based web applications. So, the GTN has added a Pan Galactic Tool
Search!

<iframe src="https://genomic.social/@scottcain/110498601155255378/embed" width="100%" height="250">
Scott: Arg! Please don't use common words to name your software projects! First (and a long time ago) you have the fine folks @galaxyproject using that word, then you have the SV caller Manta, making it impossible to answer with a simple google whether manta is available at a public Galaxy hub.
</iframe>

## How it Works

The GTN uses the list of [Public Galaxy
Servers](https://galaxyproject.org/use/) available from the Hub to construct
the dropdowns you see on Tutorials indicating which servers support a given
tutorial.

Every time the GTN gets deployed (at minimum once a day), we collect a list of
tools that are available from each of those servers. This metadata is available
directly from the GTN (`/api/psl.json`).

When we render a tutorial, we check the list of tools used in that tutorial,
both those annotated by authors, as well as whichever ones are mentioned in the
associated workflows, and list the intersection of what's available.

Since we had this metadata already available, it was trivial to generate a tool
search which simply searches through this file and lists the relevant servers.

## Implementation

As the GTN is a static site, so is our search. It takes a second to load the
tool search interface (~7MB uncompressed), but after that you can search more
or less instantaneously through all 15k versions of the 6k known tools we see
across the Galaxy.
