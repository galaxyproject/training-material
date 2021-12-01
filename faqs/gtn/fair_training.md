---
title: How does the GTN ensure our training materials are FAIR?
area: contributors
layout: faq
---

This infrastructure has been developed in accordance with the FAIR (Findable, Accessible, Interoperable, Reusable) principles for training materials {% cite Garcia2020 %}. Following these principles enables trainers and trainees to find, reuse, adapt, and improve the available tutorials.

10 Simple Rules                                                            | Implementation in GTN framework
---------------                                                            | -------------------------------
Plan to share your training materials online                               | Online [training material portfolio](https://training.galaxyproject.org/), managed via [a public GitHub repository](https://github.com/galaxyproject/training-material)
Improve findability of your training materials by properly describing them | Rich metadata associated with each tutorial that are visible and accessible via schema.org on each tutorial webpage.
Give your training materials a unique identity                             | URL persistency with redirection in case of renaming of tutorials. Data used for tutorials stored on Zenodo and associated with a Digital Object Identifiers (DOI)
Register your training materials online                                    | Tutorials automatically registered on TeSS, the [ELIXIR's Training e-Support System](https://tess.elixir-europe.org/)
If appropriate, define access rules for your training materials            | Online and free to use without registration
Use an interoperable format for your training materials                    | Content of the tutorials and slides written in Markdown. Metadata associated with tutorials stored in YAML, and workflows in JSON. All of this metadata is available [from the GTN's API](https://training.galaxyproject.org/training-material/api/)
Make your training materials (re-)usable for trainers                      | Online. Rich metadata associated with each tutorial: title, contributor details, license, description, learning outcomes, audience, requirements, tags/keywords, duration, date of last revision. Strong technical support for each tutorial: workflow, data on Zenodo and also available as data libraries on UseGalaxy.\*, tools installable via the Galaxy Tool Shed, list of possible Galaxy instances with the needed tools.
Make your training materials (re-)usable for trainees                      | Online and easy to follow hands-on tutorials. Rich metadata with "Specific, Measurable, Attainable, Realistic and Time bound" (SMART) learning outcomes following Bloom's taxonomy. Requirements and follow-up tutorials to build learning path. List of Galaxy instances offering needed tools, data on Zenodo and also available as data libraries on UseGalaxy.\*. Support chat embedded in tutorial pages.
Make your training materials contribution friendly and citable             | Open and collaborative infrastructure with contribution guidelines, a CONTRIBUTING file and a chat. Details to cite tutorials and give credit to contributors available at the end of each tutorial.
Keep your training materials up-to-date                                    | Open, collaborative and transparent peer-review and curation process. Short time between updates.
