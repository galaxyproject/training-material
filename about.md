---
layout: page
title: About the GTN
hide_title: true
---

<div id='hero'>
  <p id="hero-headline">
  The GTN is a collection of free, FAIR, open-source, reusable e-Learning materials for life sciences and beyond
  </p>
  <p id="hero-subtitle">
  Our mission: provide a platform for accessible, interactive training materials for everyone.
  </p>
</div>

<h2>Our Pillars</h2>

<div class="row" id="button-row">
<div class="col-md-2">
<div>Free & Open</div>
<a href="#free--open-source">
<img src="{% link assets/images/undraw_open_source.svg %}" alt="cartoon of MIT and CC license">
</a>
</div>
<div class="col-md-2">
<div>Automatically FAIR</div>
<a href="#findable-accessible-interoperable-reusable">
<img src="{% link assets/images/undraw_automation.svg %}" alt="cartoon of a factory">
</a>
</div>
<div class="col-md-2">
<div>Community Driven</div>
<a href="#community">
<img src="{% link assets/images/undraw_community.svg %}" alt="cartoon of several humans chatting">
</a>
</div>
<div class="col-md-2">
<div>Accessible</div>
<a href="#accessibility">
<img src="{% link assets/images/undraw_accessible.svg %}" alt="cartoon of a blind person crossing a street">
</a>
</div>
<div class="col-md-2">
<div>Easy to Contribute</div>
<a href="#easy-contributions">
<img src="{% link assets/images/undraw_developer.svg %}" alt="cartoon of a person at a laptop working">
</a>
</div>
<div class="col-md-2">
<div>Classroom or Self-study</div>
<a href="#e-learning-materials">
<img src="{% link assets/images/undraw_online_learning.svg %}" alt="cartoon of a person studying at a laptop">
</a>
</div>
</div>

## A global view

Our materials have grown significantly over the years. We have a global audience, with users from over 200 countries. Our materials are used by thousands of learners every month. See [more detailed statistics]({% link stats/index.md %}) on our dedicated page.

<div markdown=0>
{% snippet faqs/gtn/gtn_stats.md %}
</div>

## Findable, Accessible, Interoperable, Reusable

<div class="row">
<div class="col-md-4">
<figure>
  <img alt="A screenshot of a radar plot visibly from fair-checker's site showing a perfect score." src="{% link news/images/FAIR-Checker-us.png %}">
  <figcaption>FAIR Checker reports a 100% FAIRness score for GTN materials, better than some competing e-learning material infrastructures.</figcaption>
</figure>
</div>
<div class="col-md-8" markdown=1>

The GTN infrastructure has been developed in accordance with the FAIR (Findable, Accessible, Interoperable, Reusable) principles for training materials ({% cite Garcia2020 %}). Following these principles enables trainers and trainees to find, reuse, adapt, and improve the available tutorials.

By contributing to the GTN, you can be certain that materials are 100% FAIR and compliant with the 10 Simple Rules for FAIR training materials. Every learning resource contributed to the GTN **automatically** receives:

- Automatic Persistent URLs
- Automatic TeSS registration
- Automatic BioSchemas compliant LD+JSON metadata
- Automatic monthly archival copies on [the GTN archive](https://training.galaxyproject.org/archive/)
- Page view statistics via [Plausible](https://plausible.galaxyproject.eu/training.galaxyproject.org)

We have an open and welcoming collaborative peer-review and curation process as users contribute materials to the GTN.

The GTN is integrated with the [ELIXIR Training e-Support System (TeSS)](https://tess.elixir-europe.org/), which is a platform for discovering and sharing learning materials and events. Every tutorial, slide deck, event, and learning pathway in the GTN is automatically registered with TeSS, making it easy for users to find and access the materials.

</div>
</div>

## Accessibility

The GTN strives to meet WCAG 2.0 AA compliance, and WCAG 2.0 AAA where possible. We regularly test our materials for their screen-reader accessibility and have implemented <a href="{% link accessibility.md %}">a range of features</a> to ensure the website is accessible to everyone, regardless of disability.

## e-Learning Materials

<div class="row">
<div class="col-md-8" markdown=1>

Materials in the GTN are designed to be interactive, engaging, and accessible. Many are additionally designed to reproduce a published analysis from literature, helping document and showcase reproducibility.

Tutorials are designed around rich metadata, every tutorial including extensive metadata that help learners and educators decide if the materials are useful to them:

- Audience
- Difficulty Level
- SMART Learning Objectives
- Key Points
- Pre-requisite materials
- License (always open source!)
- Duration (for teaching)
- and more

Every tutorial is highly structured with both expository information, clearly demarcated hands-on exercises, and then formative assessments to help learners check their understanding.

GTN e-Learning materials are complemented with excellent technical support such as input datasets available on Zenodo, as well as workflows which match the tutorial's steps. Some tutorials additionally feature recordings from previous times this material was taught. Galaxy tutorials also indicate which public servers can be used for the tutorial, identifying servers which support every tool used in a given learning material.

</div>
<div class="col-md-4">
<figure>
  <img alt="Screenshot of a GTN learning material showing metadata and authors for a tutorial" src="{% link assets/images/tutorial-overview.png %}">
  <figcaption>GTN e-Learning materials feature robust metadata and a clear, easy to follow structure.</figcaption>
</figure>
</div>
</div>

## Proven Impact

Please read further about our community driven approach and the infrastructure we have built to support both educators and learners in {% cite Batut2018 %} and {% cite hiltemann2023galaxy %}.

## Community

Our community has been key to the GTN's success, providing materials, feedback, and support to learners. We have a growing number of contributors, who have helped to create our diverse set of materials.

### Easy Contributions

We aim to make contribution easy: all contributions to the GTN can be done in a web browser, without the installation of any additional software to decrease burdens of contributors. Additionally we take many steps to make it easier for contributors to contribute. Based on their feedback we support adding events, news, FAQs, and more via <a href="{{ site.baseurl }}/news/2024/07/17/google-forms.html">Google forms</a>.

### Stay in Touch

You can follow us on <a rel="me" href="https://mstdn.science/@gtn"><span style="fill: var(--hyperlink);">{{ "assets/images/mastodon.svg" | load_svg }}</span> Mastodon</a> or <a rel="me" href="https://bsky.app/profile/galaxytraining.bsky.social"><span style="fill: var(--hyperlink);">{{ "assets/images/bluesky-logo.svg" | load_svg }}</span> Bluesky</a>, or get in touch with us via <a rel="me" href="https://matrix.to/#/%23Galaxy-Training-Network_Lobby%3Agitter.im">Matrix Chat</a>.

{% include _includes/contributor-quilt.html %} 

We are always looking for new contributors, so if you have a tutorial you would like to share, please get in touch!

{% include _includes/cta.html title="How You Can Help" subtitle="Join the GTN today, sharing your learning and expertise!" url="/training-material/topics/contributing/index.html" call="Get Involved" %}


## History

<div class="row">
<div class="col-md-8" markdown=1>

The GTN started in July of 2015 as a collection of just a handful of <i>Galaxy-specific</i> training materials. A few slides, and a few hands-on tutorials were available covering topics like RNA-Seq, ChIP-Seq, and other NGS topics.

The GTN is a community-driven project, with contributions from people around the world.

Over time, it has grown to include materials for a wide range of topics, in bioinformatics and beyond! We now cover a range of scientific topics and methodologies.

</div>

<div class="col-md-4">
<figure>
  <img alt="screenshot of the gtn on the left and galaxy on the right, both are very old looking" src="{% link shared/images/interactive_training-old.png %}">
  <figcaption>The GTN started out in July 2015, and has changed significantly since then</figcaption>
</figure>
</div>
</div>


## Free & Open Source

Both the GTN Framework, and all materials in the GTN are openly licensed and free to use. The GTN code and libraries is available under [MIT](https://github.com/galaxyproject/training-material/blob/main/LICENSE.md), and most tutorials are available under a Creative Commons license [CC-BY-4.0](https://spdx.org/licenses/CC-BY-4.0")


<style>
#hero {
    border: 2px solid var(--border-light);
    margin: 1rem;
    padding: 1rem 3rem;
    text-align: center;
    border-radius: 25px;
}
#hero-headline {
    font-size: 2.5rem;
    font-weight: bold;
}
#hero-subtitle {
    font-size: 1.5rem;
}
.main-content {
    font-size: 1.25;
}
figure {
margin: 0;
}
.tutorial_tag {
  margin-inline-end: 0.2rem;
}

#button-row div.col-md-2 div {
    font-size: 1.5rem;
}
#button-row a {
    display: block;
    border: none;
}
#button-row a:hover {
    border: 1px solid var(--hyperlink);
}
</style>

## References

{% bibliography --cited %}
