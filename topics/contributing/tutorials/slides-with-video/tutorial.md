---
layout: tutorial_hands_on

title: "Adding auto-generated video to your slides"
questions:
  - "How can we add auto-generated video?"
  - "How does it work?"
  - "What do I need to do to make it optimal for viewers?"
objectives:
  - "Adding a video to a set of slides"
time_estimation: "20m"
subtopic: extras
key_points:
  - "Thanks to the GTN, videos are easy to add"
  - "Be mindful of your captions. Short sentences are good!"
contributors:
  - hexylena
---

# Video Lectures
{:.no_toc}

Based on the work by Delphine Larivière and James Taylor with their [COVID-19 Lectures](https://github.com/galaxyproject/video-lectures/) we have implemented a similar feature in the Galaxy Training Network.

> ### Agenda
>
> In this tutorial, we will:
>
> 1. TOC
> {:toc}
>
{: .agenda}

## How it Works

We wrote a short script which does the following:

*Locally and in production*:

- Extracts a 'script' from the slides. We extract every presenter comment in the slidedeck, and turn this into a text file.
- Every line of this text file is then narrated by [Amazon Polly](https://aws.amazon.com/polly/) (if you have money) or [MozillaTTS](https://github.com/synesthesiam/docker-mozillatts) (free).
- The slide deck is converted to a PDF, and then each slide is extracted as a PNG.
- Captions are extracted from the audio components.
- The narration is stitched together into an mp3
- The images are stitched together into an mp4 file
- The video, audio, and captions are muxed together into a final mp4 file

*In production*

- We use Amazon Polly, paid for by the Galaxyproject
- The result is uploaded to an S3 bucket

# Enabling Video

We have attempted to simplify this process as much as possible, but making good slides which work well is up to you.

## Writing Good Captions

Every slide must have some narration in the presenter notes. It does not make sense for students to see a slide without commentary. For each slide, you'll need to write presenter notes in full, but short sentences.

### Sentence Structure

Use simple and uncomplex sentences whenever possible. Break up ideas into easy to digest bits. Students will be listening to this spoken and possibly reading the captions.

*2021-05-01* There used to be a limit of ~120 characters per sentence, but this is no longer an issue. We now break up sentences which are too long in the captions and show them over multiple timepoints. So if you need to write a really long sentence, you can, but we still advise to simplify sentences where possible.

### Captions per Slide

Every slide must have some speaker notes in this system, **NO exceptions**.

### Punctuation

Sentences should end with punctuation like `.` or `?` or even `!` if you're feeling excited.

### Abbreviations

These are generally fine as-is. (e.g. `e.g.`/`i.e.` is fine as-is, `RNA` is fine, etc.) Make sure abbreviations are all caps though.

> **Good**
> This role deploys CVMFS.
{: .code-in}

### "Weird" Names

In the captions you will want to teach the GTN how to pronounce these words by editing `bin/ari-map.yml` to provide your definition.

E.g.

Word       | Pronunciation
---------- | ---
SQLAlchemy | SQL alchemy
FastQC     | fast QC
nginx      | engine X
gxadmin    | GX admin
/etc       | / E T C

The same applies to the many terms we read differently from how they are written, e.g. 'src' vs 'source'. Most of us would pronounce it like the latter, even though it isn't spelt that way. Our speaking robot doesn't know what we mean, so we need to spell it out properly.

So we write the definition in the [`bin/ari-map.yml`](https://github.com/galaxyproject/training-material/blob/master/bin/ari-map.yml) file.

### Other Considerations

(*Written 2020-12-16, things may have changed since.*)

Be sure to check the pronunciation of the slides. There are known issues with [heteronyms](https://en.wikipedia.org/wiki/Heteronym_(linguistics)), words spelt the same but having different pronunciation and meaning. Consider "read" for a classic example, or ["analyses"](https://en.wiktionary.org/wiki/analyses#English) for one that comes up often in the GTN. "She analyses data" and "Multiple analyses" are pronounced quite differently based on their usage in sentences. See the [wiktionary](https://en.wiktionary.org/wiki/analyses#English) page for more information, or the [list of English heteronyms](https://en.wiktionary.org/wiki/Category:English_heteronyms) you might want to be aware of.

This becomes an issue for AWS Polly and Mozilla's TTS which both don't have sufficient context sometimes to choose between the two pronunciations. You'll find that "many analyses" is pronounced correctly while "multiple analyses" isn't.

Oftentimes the services don't understand part of speech, so by adding adjectives to analyses, you confuse the engine in to thinking it should be the third person singular pronunciation. This is probably because it only has one or two words of context ahead of the word to be pronounced.

## Enable the Video

Lastly, we need to tell the GTN framework we would like videos to be generated.

> ### {% icon hands_on %} Hands-on: Enable video
>
> 1. Edit the `slides.html` for your tutorial
> 2. Add `video: true` to the top
{: .hands_on}

That's it! With this, videos can be automatically generated.


# How it works: In Detail

1. We take our markdown slides, e.g. [`topics/introduction/tutorials/galaxy-intro-short/slides.html`](https://github.com/galaxyproject/training-material/blob/main/topics/introduction/tutorials/galaxy-intro-short/slides.html)
2. In order for them to be processed, slides must have an annotation saying `video: true` in the header metadata, and then 'speaker notes' (everything after the ??? before the ---)
3. This is turned into our 'plain text slides' which just renders the markdown a bit more nicely ([example](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-short/slides-plain.html))
4. Then we run ari.sh which does the following:

	- `make video` is run which runs [`bin/ari-make.sh`](https://github.com/galaxyproject/training-material/blob/main/bin/ari-make.sh)
	- This builds PDFs for any slides which have changed
	- And runs `./bin/ari.sh` with the PDF, the original Slides, and where the mp4 should be saved.
		- In [`./bin/ari.sh`](https://github.com/galaxyproject/training-material/blob/main/bin/ari.sh)
		- It [extracts metadata](https://github.com/galaxyproject/training-material/blob/main/bin/ari.sh#L38) from the tutorial (title, authors, etc.)
		- It [builds a 'script'](https://github.com/galaxyproject/training-material/blob/main/bin/ari.sh#L51), a json document with blocks for every line of the speaker notes that were in the slides.
		- [Those get converted into mp3 files](https://github.com/galaxyproject/training-material/blob/main/bin/ari.sh#L55) by AWS Polly (or MozillaTTS), one per slide.
		- The PDFs get turned into [a series of PNG images](https://github.com/galaxyproject/training-material/blob/main/bin/ari.sh#L60)
		- We take the timings of the mp3 files together with the json 'script' to [write out webvtt / srt subtitles](https://github.com/galaxyproject/training-material/blob/main/bin/ari.sh#L69) which get embedded into the video, and supplied next to it.
		- [editly is used](https://github.com/mifi/editly) to knit together the PNGs + mp3s with appropriate delay


All of this is run on cron by [`.github/workflows/video.yml`](https://github.com/galaxyproject/training-material/blob/main/.github/workflows/video.yml) which handles building all of these videos and then later uploading them to s3.

Many of the scripts internally are prefixed with `ari`, we named our internal version after [github.com/jhudsl/ari/](https://github.com/jhudsl/ari/) which inspired it, but we wanted a version that would be more closely tied to the GTN and integrate with our infrastructure nicely, so we ended up writing our own.

# Conclusion
{:.no_toc}
