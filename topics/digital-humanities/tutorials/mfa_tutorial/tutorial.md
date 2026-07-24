---
layout: tutorial_hands_on
title: Montreal Forced Aligner
level: Introductory
zenodo_link: https://zenodo.org/records/21371207
questions:
  - How can you use Montreal Forced Aligner to align an audio and a TextGrid file in Galaxy?
objectives:
- Use the Montreal Forced Aligner to break down speech into its smallest sounds and align them with their corresponding transcription.
time_estimation: 1H
contributions:
  authorship:
    - wallajess    
  reviewing:
    - Sch-Da
  testing:
    - Sch-Da
---

<!-- This is a comment. -->

This tutorial explains how to use the Montreal Forced Aligner (MFA) {% cite mcauliffe17_interspeech %} on the Galaxy platform. 

The aligner breaks down speech into its smallest possible sounds (phones) and aligns the indivdual 
sounds to their corresponding orthographic transcription. For example, the word 'fox' is broken down
into four sounds: F AA K S or /f ɑː k s/ (depending on whether you use ARPAbet or IPA for the transcription). The phonetic transcription of 
each sound is given a precise time boundary to match it to the corresponding sound in the audio file. 
Each sound can then be analyzed individually without the need for lengthy manual segmentation. For example, researchers can then 
measure and compare each pronunciation of a certain sound in very little time.

To achieve this, MFA has acoustic models and pronunciation dictionaries for numerous languages, 
which can be selected in Galaxy. 

The tool requires an audio or video file and an orthographic transcription as a TextGrid file. 
Audio can be transcribed using e.g. Elan {% cite wittenburg-etal-2006-elan %} or Praat {% cite Praat %}. Speech should be segmented into breath groups (the speech in between two breaths) and the transcript exported as a TextGrid. The audio file and TextGrid should then be compressed into one .zip file to be uploaded into Galaxy.


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Preprocessing

Audio or video files must first be segmented and transcribed using e.g. Elan {% cite wittenburg-etal-2006-elan %} or Praat
{% cite Praat %} and the transcript exported as a TextGrid (alternatively as a .lab/.txt file but these must be pasted in as a single line). The TextGrid should contain a single 
tier with short intervals (a breath gorup). Each interval should correspond to an utterance in the audio file. The tier name can be a speaker ID (though this is not necessary).


## Required data formats
The filenames of the audio/video file and the transcript must be identical except 
for the extension (e.g. .wav or .TextGrid). It is helpful to have the speaker ID 
prefix to the filename (e.g. JW_interview.wav, JW_interview.TextGrid), then you 
can instruct the aligner to use the corresponding number of characters (2 in the previous example) as the speaker information. This speaker ID will be retained in the output, allowing you to later identify and analyze speech patterns by speaker.

These should then be compressed into one .zip file to be uploaded into Galaxy.

## Data Validation Checklist

Before uploading, verify your files meet these requirements:

- {% icon param-check %} Audio file and TextGrid have identical base filenames (e.g., `JW_interview.wav` and `JW_interview.TextGrid`)
- {% icon param-check %} Audio format is WAV or MP3 (WAV is best for quality)
- {% icon param-check %} TextGrid contains a single tier with short intervals, each representing one breath group
- {% icon param-check %} Transcription text matches the actual spoken content (typos will cause alignment errors)
- {% icon param-check %} Speaker ID prefix is consistent, if present (e.g. all files start with 2-3 character speaker code)
- {% icon param-check %} Files are compressed into a single .zip archive

The following example file matches the above criteria and will be used in the next steps.

## Uploading the data
> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial
> 2. Import your files or if you want to follow the tutorial use the files from [Zenodo]({{ page.zenodo_link }}) or from
>    the shared data library (`GTN - Material` -> `{{ page.topic_name }}`
>     -> `{{ page.title }}`):
>
>    ```
>    https://zenodo.org/records/21371207/files/the_fox_and_the_grapes.zip?download
>    
>    ```
>    
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets to a more meaningful name (optional):
>    - Click on the dataset produced by MFA Find OOVs
>    - Click the {% icon galaxy-pencil %} (edit) icon
>    - Choose a meaningful name, for example, `the_fox_and_the_grapes.zip`
>    
> 4. Change the datatype of the dataset to mfa_corpus_model.zip
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

# Align the files

Before we align the files, we need to first check whether all of the words are 
in our dictionary. MFA offers the Find OOVs (out of vocabulary words) tool to do this.


## Sub-step with **MFA Find OOVs**

> <hands-on-title> Find out-of-vocabulary words </hands-on-title>
>
> 1. {% tool [MFA Find OOVs](mfa_find_oovs) %} with the following parameters:
>    - {% icon param-file %} *"Corpus Archive"*: `the_fox_and_the_grapes.zip` (Uploaded dataset)
>    - *"Dictionary Source"*: `Use a built-in MFA Dictionary`
>        - *"Select Dictionary"*: `English Us (ARPA)`
>    - *"Speaker Characters"*: `0`
>
> 2. Rename the output file for clarity:
>    - Click on the dataset you uploaded
>    - Click the {% icon galaxy-pencil %} (edit) icon
>    - Change the name to `OOVs_english_us` or similar
>    - This will help you identify this as the out-of-vocabulary words file
>
> **Expected output:** A .zip file containing a text file (`oovs_found_english_us_arpa.txt`) listing all words from your transcription that are not in the English US ARPA dictionary. If no OOVs are found, the file will be empty (if that is the case, you can skip the next three steps and go straight to MFA Align).
>
{: .hands_on}

## Sub-step with **Unzip**

In order to generate the  pronunciations of the OOVs, we need the single file oovs_found_english_us_arpa.txt. We therefore need to unzip the output of Find OOVs so we can work with the file. The missing pronunciations can be generated automatically in the next step (MFA G2P) or manually using the editor.

> <hands-on-title> Unzipping your file </hands-on-title>
>
> 1. {% tool [Unzip a file](toolshed.g2.bx.psu.edu/repos/imgteam/unzip/unzip/6.0+galaxy5) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `OOVs_english_us` (output of **MFA Find OOVs** {% icon tool %})
>    - *"What to extract"*: `Single file`
>        - *"File path"*: `oovs_found_english_us_arpa.txt`
>
> 2. Rename the output file for clarity:
>    - Click on the dataset produced by MFA Find OOVs
>    - Click the {% icon galaxy-pencil %} (edit) icon
>    - Change the name to `unknown_words` or similar
>    - This will help you identify this as the out-of-vocabulary words file
>
>    > <comment-title> Check the OOVs </comment-title>
>    >
>    > Double-check the OOVs found file and make sure that none of the OOVs are typos. To do so, click on the {% icon galaxy-eye %} icon next to the dataset name to look at the file. If you do find a typo, you must correct this in the transcription, then upload the .TextGrid again and rerun the previous steps.
>   {: .comment}
>
{: .hands_on}


## Sub-step with **MFA G2P**

Now we can generate pronunciations for each of the out of vocabulary words using MFA's built in grapheme-to-phoneme (G2P) tool. Select the same model you used for finding OOVs.

> <hands-on-title> Generate pronunciations with G2P </hands-on-title>
>
> 1. {% tool [MFA G2P](mfa_g2p) %} with the following parameters:
>    - *"Input Mode"*: `Word List (Text File)`
>        - {% icon param-file %} *"Word List"*: `unknown_words` (output of **Unzip** {% icon tool %})
>    - *"G2P Model Source"*: `Use a built-in MFA Model`
>        - *"Select G2P Model"*: `English Us (ARPA)`
>
>    > <comment-title> Check the pronunciations manually </comment-title>
>    >
>    > It is important to check the pronunciations manually. MFA can also offer multiple possible pronunciations for each word. Correct any incorrect pronunciations and delete any additional pronunciations that do not apply. To do so, you can use the editor which can be found in Visualization on the side bar. There you can select the generated dictionary and correct the pronunciation.
>    {: .comment}
>
>    > <tip-title> Phonetic notation </tip-title>
>    >
>    > This tutorial uses ARPAbet phonetic notation for US English (e.g., AE, IH, K). Other models use IPA, which is the international standard for most other languages and varieties.
>    {: .tip}
>
> 2. Rename the output file for clarity:
>    - Click on the `output_dictionary` dataset produced by MFA G2P
>    - Click the {% icon galaxy-pencil %} (edit) icon
>    - Change the name to `generated_pronunciations_oovs` or `generated_dictionary`
>    - This clearly identifies this as the file containing generated pronunciations for OOV words
>
> **Expected output:** A dictionary file in MFA format (one entry per line: word followed by phonetic transcription). File size depends on the number of OOVs.
>
{: .hands_on}


## Sub-step with **MFA Merge**

Now we will add the new words to the built-in dictionary so that they are recognized when we run the aligner.

> <hands-on-title> Add new words to your dictionary </hands-on-title>
>
> 1. {% tool [MFA Merge](mfa_merge) %} with the following parameters:
>    - *"Dictionary Source"*: `Use a built-in MFA Dictionary`
>        - *"Select Dictionary"*: `English Us (ARPA)`
>    - In *"Dictionary to add"*:
>        - {% icon param-repeat %} *"Insert Dictionary to add"*
>            - {% icon param-file %} *"MFA dictionary"*: `generated_pronunciations_oovs` (output of **MFA G2P** {% icon tool %})
>
>
>    > <comment-title> If corrections were uploaded </comment-title>
>    >
>    > If you uploaded a file with corrected pronunciations, you will have to select the uploaded file with the corrected dictionary entries as your dictionary to add.
>    {: .comment}
>
> 2. Rename the output file for clarity:
>    - Click on the `output` dataset produced by MFA Merge
>    - Click the {% icon galaxy-pencil %} (edit) icon
>    - Change the name to `updated_dictionary_english_us` or `merged_dictionary`
>    - This clearly identifies this as the merged dictionary containing both built-in entries and newly generated pronunciations
>
> **Expected output:** A merged dictionary file (single file, not a .zip) containing all entries from the built-in English US ARPA dictionary plus your newly generated pronunciations for OOV words. This is the complete dictionary your audio will be aligned against.
>
{: .hands_on}

## Sub-step with **MFA Align**

Now we're ready to align the pronunciations with the audio. Since we only have one spekaer and no spekaer ID, we leave the number of spekaer characters at 0. If you have a spekaer ID at the start of your file name, then select the number of characters corresponding to the speaker ID.

> <hands-on-title> Align Audio and Pronounciation </hands-on-title>
>
> 1. {% tool [MFA Align](mfa_align) %} with the following parameters:
>    - {% icon param-file %} *"Corpus Archive"*: `the_fox_and_the_grapes.zip` (Input dataset)
>    - *"Dictionary Source"*: `Upload from History`
>        - {% icon param-file %} *"Dictionary File"*: `updated_dictionary_english_us` (output of **MFA Merge** {% icon tool %})
>    - *"Acoustic Model Source"*: `Use a built-in MFA Model`
>        - *"Select Acoustic Model"*: `English Us (ARPA)`
>
>
>    > <comment-title> Check the output in Praat </comment-title>
>    >
>    > Always validate the alignment output:
>    >
>    > 1. Download the .TextGrid file from the output archive (see steps below)
>    > 2. Open the original audio file + aligned TextGrid together in Praat
>    > 3. Check: Do the boundaries align with the actual speech transitions?
>    > 4. Listen for obvious misalignments: does silence have phoneme labels? Are vowels properly segmented?
>    > 5. Minor misalignment is normal, but if you have more serious or widespread errors, check your transcription or OOVs
>    >
>    > **Praat playback tip:** Right-click on a TextGrid interval and select "Play" to hear just that phoneme segment.
>    {: .comment}
>
> 2. Rename the output file for clarity:
>    - Click on the output archive produced by MFA Align
>    - Click the {% icon galaxy-pencil %} (edit) icon
>    - Change the name to `aligned_textgrids` or `mfa_output_aligned`
>    - This clearly identifies this as the final aligned output containing the TextGrid files with time boundaries for each phoneme
>
> **Expected output:** A .zip archive containing TextGrid files (one per input file) with three tiers: words, phones (individual sounds), and word confidence scores. Total runtime depends on audio length.
>
{: .hands_on}

## Downloading and Exporting Results

> <hands-on-title> Extract and download your aligned files </hands-on-title>
>
> 1. Download the aligned archive:
>    - Click the download icon ({% icon download %}) next to the `aligned_textgrids` dataset
>    - Save the .zip file to your computer
>
> 2. Extract the files on your computer
>   
> 3. Open in Praat for analysis:
>    - In Praat, open the original audio file (File → Open → Read from file...)
>    - Then open the corresponding TextGrid (Objects window → Read → Read from file...)
>    - Select both and click "View & Edit" to visualize the alignment
>
>
{: .hands_on}


# Next Steps and Use Cases

Once you have aligned TextGrids, you can:

- **Measure vowel formants:** Extract and analyze vowel characteristics from aligned vowel intervals. To do so, you can use new-fave.
- **Compare pronunciations:** Analyze variation across speakers or time
- **Assess speech rate:** Measure duration of phonemes or words to study articulation speed
- **Study consonant voicing:** Use the precise time boundaries to measure voice onset time (VOT)


# Conclusion

In this tutorial, we uploaded a sample audio file and transcription. 
We checked for any out-of-vocabulary words (OOVs) missing from our dictionary and generated the pronunciations using G2P. 
After checking and, if necessary, correcting these pronunciations, we merged the generated pronunciations with the built-in dictionary and ran the aligner on our files, producing aligned TextGrid files with precise phoneme-level time boundaries.

You now have a fully aligned and labeled corpus ready for phonetic analysis.
