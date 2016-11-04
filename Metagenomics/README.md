Topic name
==========

*Short introduction about the topic, why it is important to learn it (the questions behind this topic), what are the objectives, ...*

# Slides

Information on where the slide decks can be found, generally `http://bgruening.github.io/training-material/<topic_name>/slides/index.html` and what is inside the slide decks.

# Tutorials

## Metagenomics with Mothur: MiSeq SOP

**Trainer/admin instructions:**

The following contains instructions on configuring your Galaxy instance for this tutorial:

- Tools:
  - mothur suite (MTS, owner: iuc)
  - regex find and replace (MTS, owner: jjohnson)
  - xy_plot (MTS, owner: devteam)
  - collapse_collections (MTS, owner nml)

- Data
  - Get from zenodo here: TODO
  - Optional: put in a data library named *Galaxy training: Metagenomics with Mothur - MiSeq SOP* to prevent mass uploading during a class (training manual instruct users to check for presence of this data library before uploading data themselves). Create a subfolder named *Input Data* and put the contents of `input_data.zip` there. Also create another folder named *Reference Data* and place contents of reference_data.zip there.

- Galaxy Config
  - Requires at least Galaxy version 16.04 (for Phinch external display application)
  - For Phinch viewing, make sure your webserver is configured to allow requests from the Phinch server. For example in nginx add following to your config:

    ```
    location /galaxy {
      add_header Access-Control-Allow-Origin "http://www.bx.psu.edu";
      ..
    }
    ```

## Input datasets

*Zenodo DOI badge, e.g. [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.60520.svg)](http://dx.doi.org/10.5281/zenodo.60520)*

The input datasets for the tutorials are available on
[Zenodo with a dedicated DOI](http://dx.doi.org/10.5281/zenodo.60520).

## Galaxy instance

For these tutorials, you can use the dedicated Docker image:

```
docker run -d -p 8080:80 <name-of-docker-image>
```

It will launch a flavored Galaxy instance available on
[http://localhost:8080](http://localhost:8080).

# References

## Books

**Name et al:** [*Book name*](link/to/the/book)

    One line to describe why the book is interesting for this topic

**Name et al:** [*Book name*](link/to/the/book)

    One line to describe why the book is interesting for this topic

## Papers

**Name et al:** [*Article name*](link/to/the/article)

    One line to describe why the article is interesting for this topic

**Name et al:** [*Article name*](link/to/the/article)

    One line to describe why the article is interesting for this topic

## Websites

# Contributors

This material is maintained by:

- Maintainer 1
- Maintainer 2

For any question related to this topic and the content, you can contact them.

The following individuals have contributed to this training material:

- Name 1
- Name 2
