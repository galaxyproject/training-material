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

- Extracts a 'script' from the slides. We extract every presenter comment in the slidedeck, and turn this into a text file.
- Every line of this text file is then narrated by [Amazon Polly](https://aws.amazon.com/polly/)

  - *NB*: this currently means it is non-free, and requires AWS credentials in order to build the videos locally, which is not an option for everyone. We do not have plans for integrating an open source voice engine but we would welcome it!

- The slide deck is converted to a PDF, and then each slide is extracted as a PNG.
- Captions are extracted from the audio components.
- The narration is stitched together into an mp3
- The images are stitched together into an mp4 file
- The video, audio, and captions are muxed together into a final mp4 file
- This is uploaded to an S3 bucket

# Enabling Video

We have attempted to simplify this process as much as possible, but making good slides which work well is up to you.

## Writing Good Captions

Every slide must have some narration in the presenter notes. It does not make sense for students to see a slide without commentary. For each slide, you'll need to write presenter notes in full, but short sentences.

### Sentence Structure

Use short and uncomplex sentences whenever possible. Break up ideas into easy to digest bits. Students will be listening to this spoken and possibly reading the captions.

The captioning process is completely automated, but it means that for very long sentences, we do not currently break them up into multiple captions. So please keep your sentences under ~80 characters where possible.

> > **Good**
> > - Configuration management manages the configuration of machines.
> > - It specifies what software should be installed, and how it should be configured.
> {: .code-in}
>
> > **Bad**
> > Configuration management manages the configuration of machines; it specifies what software should be installed, and how it should be configured
> {: .code-out}
{: .code-2col}


### Punctuation

Every sentence must end with punctuation like `.` or `?` or even `!` if you're feeling excited. We merge all lines of the presenter notes before having them spoken by Polly. If you forget the punctuation it will be one long run-on sentence.

### Abbreviations and "Weird" Names

In the captions you will want to spell out abbreviations

> > **Good**
> > This role deploys C V M F S.
> {: .code-in}
>
> > **Bad**
> > This role deploys CVMFS.
> {: .code-out}
{: .code-2col}

For many terms we read them differently than they're written, e.g. 'src' vs 'source'. Most of us would pronounce it like the latter, even though it isn't spelt that way. Our speaking robot doesn't know what we mean, so we need to spell it out properly.

> > **Good**
> > Copy copies a file from src on localhost...
> {: .code-in}
>
> > **Bad**
> > Copy copies a file from source on localhost...
> {: .code-out}
{: .code-2col}


## Enable the Video

Lastly, we need to tell the GTN framework we would like videos to be generated.

> ### {% icon hands_on %} Hands-on: Enable video
>
> 1. Edit the `slides.html` for your tutorial
> 2. Add `video: true` to the top
{: .hands_on}

That's it! With this, videos can be automatically generated.


# Conclusion
{:.no_toc}
