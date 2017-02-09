Topic name
==========

*Short introduction about the topic, why it is important to learn it (the questions behind this topic), what are the objectives, ...*

# Slides

# Tutorials

## 16S Microbial Analysis with Mothur: MiSeq SOP

**Trainer/admin instructions:**

Instructions for configuring your Galaxy instance for this tutorial:

- Tools:
  - mothur suite (MTS, owner: iuc)
  - xy_plot (MTS, owner: devteam)
  - collapse_collections (MTS, owner nml)
  - newick_display (MTS, owner dcorreia)
  - krona_text (MTS, owner saskia_hiltemann)

- Galaxy Config
  - Requires at least Galaxy version 16.04 (for Phinch external display application)
  - For Phinch viewing, make sure your webserver is configured to allow requests from the Phinch server. For example in nginx add following to your config:

    ```
    location /galaxy {
      add_header Access-Control-Allow-Origin "http://www.bx.psu.edu";
      ..
    }
    ```
  - For viewing of the html outputs in Galaxy (e.g. Krona, tools that output SVGs wrapped in HTML), change the following line in your galaxy.ini file  `sanitize_all_html = False`.
  - For viewing SVGs directly withing galaxy, set `serve_xss_vulnerable_mimetypes = True` in `galaxy.ini`

## Input datasets

**Metagenomics with Mothur: MiSeq SOP** The input datasets for the tutorials are available on
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.165147.svg)](https://doi.org/10.5281/zenodo.165147)


## Galaxy instance

For these tutorials, you can use the dedicated Docker image:

```
docker run -d -p 8080:80 <name-of-docker-image>
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080 ).

# References

## Books

## Papers

## Websites

# Contributors

This material is maintained by:

- Saskia Hiltemann
- Maintainer 2

For any question related to this topic and the content, you can contact them.

The following individuals have contributed to this training material:

- Saskia Hiltemann
- Stefan Boers
