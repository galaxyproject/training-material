---
title: "New Feature: Easy Abbreviation"
tags: [new feature]
contributors: [hexylena, rikeshi, simonbray]
tutorial: "topics/dev/tutorials/bioblend-dev/tutorial.html"
layout: news
---

Thanks to the great tutorial developed by first time contributor {% include _includes/contributor-badge.html id="rikeshi" %} and edited by {% include _includes/contributor-badge.html id="simonbray" %}, we noticed that they defined a number of abbreviations and re-used those throughout their tutorial.

As the GTN is intended to be easy for contributors and easy for learners, we wanted to make use of the [`<abbr>`](https://developer.mozilla.org/en-US/docs/Web/HTML/Element/abbr) tag which allows you to define commonly re-used terms in your HTML. However this is a bit clumsy to write every time, so we've implemented a tag and some metadata which permits easily defining and referencing abbreviations throughout your text.

In your tutorial metadata you can add an abbreviations section like:

```yaml
---
title: My awesome tutorial
...
abbreviations:
  API: Application Programming Interface
  JSON: JavaScript Object Notation
---
```

And in your text you can use braces to refer to the term

> > ### {% icon code-in %} Input: Markdown
> > <code>
> > The `/jobs` &lbrace;API&rbrace; will return &lbrace;JSON&rbrace;. When we call the &lbrace;API&rbrace; we'll get back this result &lbrace;JSON&rbrace;.
> > </code>
> {: .code-in}
>
> > ### {% icon code-out %} Output
> >
> > The `/jobs` Application Programming Interface (API) will return JavaScript Object Notation (JSON). When we call the <abbr title="Application Programming Interface">API</abbr> we'll get back this result <abbr title="JavaScript Object Notation">JSON</abbr>.
> >
> {: .code-out}
{: .code-2col}

These will even generate an automatic Glossary at the end. Check out the use of this new feature in the BioBlend Dev Tutorial!
