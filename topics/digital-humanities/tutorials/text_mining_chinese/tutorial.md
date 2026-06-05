---
layout: tutorial_hands_on
redirect_from:
  - /topics/statistics/tutorials/text_mining_chinese/tutorial

title: Text-Mining Differences in Chinese Newspaper Articles
level: Intermediate
zenodo_link: https://doi.org/10.5281/zenodo.14899614
questions:
  - How can I automatically compare two Chinese newspaper articles?
  - What characters were censored in a Chinese newspaper published in Hong Kong in the 1930s?
objectives:
  - Learn to clean and compare two texts.
  - Extract specific information from texts.
  - Visualise your results in a word cloud.
requirements:
  - type: "internal"
    topic_name: digital-humanities
    tutorials:
      - introduction_to_dh
time_estimation: 1H
key_points:
  - The *diff* tool allows comparing two similar texts automatically
  - The word cloud shows redactions from the texts at a glance
tags:
- text mining
contributions:
  authorship:
    - Sch-Da
  funding:
    - deKCD
    - mwk
    - deNBI
answer_histories:
  - label: "UseGalaxy.eu"
    history: https://usegalaxy.eu/u/schnda/h/example-answer-history-text-mining-differences-in-chinese-newspaper-articles-3
    date: 2026-05-06
---


The British Hong Kong Government censored Chinese newspapers before their publication in the colony in the 1930s ({% cite Ng2022 %}).
Replacement characters like × visibly marked those redactions, making them obvious even to those who did not read any Chinese.

![Example of a Chinese article with symbol × marking censored characters]({% link topics/digital-humanities/tutorials/text_mining_chinese/images/Example_x_censored_Chinese_article.svg %} "Example of a Chinese article with symbol × marking censored characters")

The schematic example in Figure 1, adapted from {% cite Schneider2024 %}, shows what such a censored Chinese article looked like. It is read from right to left and top to bottom. The two more prominent lines on the right are the article title, and the following text is the article's main body. It contains the character × several times, indicating various instances where it was censored.

Despite this obvious form of censorship, no research is looking into what Chinese characters the × replaced. My dissertation ({% cite Schneider2024 %}), which informs this workflow, started at this point. Through extensive research, I found several articles censored in the Hong Kong edition of Da gong bao (大公報) as uncensored versions. Those mostly came from Chinese editions printed in mainland China, where different censorship regulations applied. Those articles from China were not censored and openly showed the characters redacted in the Hong Kong versions. An example of a censored article could be {% cite Anon..16.10.1938_5598 %} and of the uncensored version {% cite TKP.16.10.1938_18864 %}.

This tutorial uses text mining to compare two newspaper articles, a censored and an uncensored text to answer what characters were censored in a Chinese newspaper published in Hong Kong in the 1930s.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Upload the texts

The machine-readable versions of the Chinese newspaper articles I originally used in my dissertation come from a proprietary database and can not be shared here. Instead, I generated a dummy article with a similar setup in GPT and manually adapted the censorship symbols based on my research. As a result, those dummy articles differ in style and punctuation, as the articles in my research. Therefore, the input files are two texts in traditional Chinese. The first is censored, containing ×, and the second one is uncensored and does not contain replacement symbols. Both texts slightly differ in their layout, which will be unified later.

> <hands-on-title> Upload the Texts </hands-on-title>
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the files from [Zenodo](https://zenodo.org/records/14899614)
>
>
>    ```
>    https://zenodo.org/records/14899614/files/Example_Chinese_newspaper_censored.txt
>    https://zenodo.org/records/14899614/files/Example_Chinese_newspaper_uncensored.txt
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
> 3. Rename the datasets, if necessary
> 4. Check that the datatype is `txt`.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
{: .hands_on}

> <question-title>What do the uploaded texts look like?</question-title>
>
> 1. Name at least one visual difference you notice between the two texts. (No need to understand their content here.)
>
> > <solution-title></solution-title>
> >
> > 1. Visual differences you notice might be:
> > - The texts have different headlines
> > - The uncensored text contains more paragraphs
> > - The censored text contains the symbol × several times
> > - The censored text contains additional symbols at the end
> > - The texts use slightly different punctuation (impressive catch!)
> {: .solution}
>
{: .question}


# Clean both texts

We pre-process and clean both texts to make the comparison easier and more apparent. Regular Expressions (Regex) help unify the text and find and replace certain text parts. Here, the Regex also helps to restructure the text.


> <details-title> More on Regular Expressions (optional) </details-title>
>
> Regular Expressions are a powerful tool to modify text automatically based on word patterns. To read more about Regex's specifics, see its [documentation](https://docs.python.org/3/library/re.html).
> You can check what your input matches in [Regex 101](https://regex101.com). Enter the Regular Expression you want to try in the top field. Insert a sample text below to see what content your Regex catches and to find out how to adapt it.
>
{: .details}

We will use Regular Expressions in a tool called "Replace text" with four different sub-steps. Those will vary if you upload different texts. Apply this step first to the censored and then to the uncensored text to get two cleaned texts.

> <hands-on-title> Cleaning the Text with Regular Expressions </hands-on-title>
>
> 1. {% tool [Replace Text in entire line](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `Example_Chinese_newspaper_censored.txt` (censored text you uploaded)
>    - In *"Replacement"*:
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `\r`
>            - *"Additional sed commands before replacement"*: `:a;N;$!ba;`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `\n`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `\s`
>        - {% icon param-repeat %} *"Insert Replacement"*
>            - *"Find pattern"*: `(.)`
>            - *"Replace with:"*: `\1\n`
>    - Click *"Run Tool"*.
>    
>    Now repeat those steps with the uncensored text.
>
>    > <comment-title>Explaining the above Regular Expressions</comment-title>
>    > Regular expressions can not only find particular words, as you might be familiar with from regular text editors.
>    > It is more powerful and can find particular patterns, for example, only capitalised words or all numbers.
>    > In this step, we mostly delete unnecessary placeholders.
>    > The first pattern we want to find is `\r`. It catches a specific form of invisible line breaks that would create unwanted gaps in the comparison later.
>    > We delete those by leaving the optional "Replace with" field blank.
>    > The additional sed commands before replacement `:a;N;$!ba;` catch all blank spaces with this tool.
>    > It is necessary only once to ensure that particular end-of-line characters are removed consistently.
>    > Similarly, `\n` marks linebreaks. We also delete those by leaving the optional "Replace with" field blank.
>    > The next expression we search for is `\s`. Those are spaces as you see them between words on your computer. We delete those.
>    > As a result, there are no gaps in our text anymore.
>    > In the last step, we want to choose each character with `(.)` and reformat it. We want to have one character per line. Therefore, we replace all characters with `\1\n`. `\1` means the responding character, and `\n` adds a linebreak after each.
 >    > The result is a clean and reformatted text.
>    {: .comment}
>
>    {% snippet faqs/galaxy/analysis_regular_expressions.md %}
>
{: .hands_on}

Remember to apply those steps to both the censored and the uncensored text. Rename both new texts in a meaningful way. Choose, for example, "Cleaned Censored Text" and "Cleaned Uncensored Text".

> <question-title>Take a look at the texts</question-title>
>
> 1. What do the texts look like now? How have they changed?
>
> > <solution-title></solution-title>
> >
> > 1. Before, the texts had many sentences per line. Now, both texts show one character per line and have more lines. The layout is entirely different.
> >
> {: .solution}
>
{: .question}



# Compare the censored and uncensored text

We can now compare the two cleaned texts. This will visualise the differences between the two texts and mark them by colour. Make sure to upload the cleaned, censored text with the replacement characters like ‘×’ first. As text two, upload the cleaned, uncensored text without the replacement characters.  This version (HTML version) creates an HTML file, which colour codes differences as additions (green) or extractions (red) when comparing the texts.

## Create a _diff_ file for researchers

> <hands-on-title> Comparing the texts using <em>diff</em> tool </hands-on-title>
>
> 1. {% tool [diff](toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.10+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"First input file"*: `Cleaned Censored Text` (output of **Replace Text** {% icon tool %})
>    - {% icon param-file %} *"Second input file"*: `Cleaned Uncensored Text` (output of **Replace Text** {% icon tool %})
>    - *"Choose a report format"*: `Generates an HTML report to visualize the differences`
>    - Click *"Run Tool"*
{: .hands_on}

The result is two files: an HTML file and a raw output as a by-product.

> <question-title> Take a look at the HTML report</question-title>
>
> 1. What can you see in lines 19 / 21?
>
> > <solution-title></solution-title>
> >
> > 1. Line 19  on the left shows a red comma on the left side of the table. It corresponds with line 21 on the right, which contains a colon in green. This means the punctuation in the file differs. The censored version contains a comma, while the uncensored one includes a colon.
> >
> {: .solution}
>
{: .question}

The HTML file could look like this:

![Screenshot of the diff tool comparing the censored and uncensored text]({% link topics/digital-humanities/tutorials/text_mining_chinese/images/Diff_WF_HTML.jpg %} "Example of the HTML file comparing the censored and uncensored text")

It shows what passages differ in the two texts. Red parts indicate deletions, and green areas indicate additions.
If passages are identical in both texts, they remain uncoloured.
This output is very convenient for researchers, as it shows differences quickly. 
If all you want to do is compare the two texts - you have achieved your goal and can stop here. Congratulations! 

If you want to extract the censored passages in the next steps, this HTML output is not helpful for further technical processing with Galaxy. For this, we run this _diff_ tool a second time with slightly changed parameters. The output is the basis for our further analysis steps.


## Create a _diff_ file for further processing
This step runs the text comparison line by line again to create another raw file that the computer can work with.
It is less intuitive to understand at first glance. Again, select the cleaned, censored text containing `×` first and the uncensored text second.
But this time, choose another report format.

> <hands-on-title> Run another <em>diff</em> tool </hands-on-title>
>
> 1. {% tool [diff](toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.10+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"First input file"*: `Cleaned Censored Text` (output of **Replace Text** {% icon tool %})
>    - {% icon param-file %} *"Second input file"*: `Cleaned Uncensored Text` (output of **Replace Text** {% icon tool %})
>    - *"Choose a report format"*: `Text file, side-by-side (-y)`
>    - Click *"Run Tool"*
>
>    > <comment-title> The output format </comment-title>
>    > As you can see above, the output of this file is a text file, compared to the HTML output in the last step.
>    {: .comment}
>
>2. Rename your file `diff raw`.
{: .hands_on}


> <question-title>Look at the diff raw file</question-title>
>
> 1. How does this file differ from the HTML file in the last step?
>
> > <solution-title>Answer</solution-title>
> >
> > 1. This output is a text file. Changes between the two texts are not coloured.
> > Instead, they are marked by additional symbols like '>' or '<'.
> >
> {: .solution}
>
{: .question}

Great, we can use this file for our next steps.

# Select only censored lines
Now we want to extract only specific lines. To determine what content was redacted in the first text, we filter the last steps' raw output for lines containing the censorship symbol ×. 

> <hands-on-title> Select text </hands-on-title>
>
> 1. {% tool [Select lines that match an expression](Grep1) %} with the following parameters:
>    - {% icon param-file %} *"Select lines"*: `diff_raw` (output of **diff** {% icon tool %})
>    - *"that are"*: `Matching`
>    - *"the pattern"*: `×`
>    - Click *"Run Tool"*
>
>2. Rename your file `Selected lines`.
>
>
{: .hands_on}


> <question-title>What output do you get?</question-title>
>
> 1. How many lines contain ×?
>
> > <solution-title></solution-title>
> >
> > 1. 13
> >
> {: .solution}
>
{: .question}

Check the datatype of this file, it should be `txt`.
The tools we use in the next steps need tabular files instead. 
Therefore, we will change the datatype of this file to `tabular` to enable the next steps:

> <hands-on-title>Change Datatype</hands-on-title>
>
> 1. **Set the datatype** of `Selected lines` {% icon galaxy-pencil %} to `tabular`.
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
{: .hands_on}

When you have changed the format successfully, you should see × in Column 1 and the Chinese characters in Column 9.


# Ensure consistent file format

After filtering for censored lines, we add a sub-step to ensure smooth computation. The previous setup could cause an error if the characters filtered in the last step were erased. Then the extracted file would lack the last column, causing an error. This is invisible to the researchers in the file. The compute step covers this potential error and ensures all necessary columns exist.

> <hands-on-title> Compute to ensure all columns exist </hands-on-title>
>
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `Selected lines` (output of **Select** {% icon tool %})
>    - *"Input has a header line with column names?"*: `No`
>        - In *"Expressions"*:
>                - *"Add expression"*: `c9`
>                - *"Mode of the operation"*: `Replace`
>                    - *"Use new column to replace column number"*: `9`
>    - In *"Error handling"*:
>        - *"If an expression cannot be computed for a row"*: `Produce an empty column value for the row`
>    - Click *"Run Tool"*
>
> 2. Rename the result to `Censored lines`.
>
{: .hands_on}



# Summarise your findings

This step summarises how often each character appeared in the table previously.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Datamash](toolshed.g2.bx.psu.edu/repos/iuc/datamash_ops/datamash_ops/1.9+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `Censored lines` (output of **Compute** {% icon tool %})
>    - *"Group by fields"*: `9`
>    - *"Sort input"*: `Yes`
>    - In *"Operation to perform on each group"*:
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"Type"*: `count`
>            - *"On column"*: `column 9`
>   - Click *"Run Tool"*
>
> 2. Rename the output file `Quantified Results`.
>
{: .hands_on}


> <question-title></question-title>
>
> 1. How many lines does the file have now?
>
> > <solution-title></solution-title>
> >
> > 1. Three lines showing different characters: 寇, 敵 and 日.
> {: .solution}
>
{: .question}


# Sort your findings

Particularly if you get a long list in the last step, sorting the results from the most to least frequent characters is necessary.
If you are only interested in the quantitative results, this can be your final output.

> <hands-on-title> Sort </hands-on-title>
>
> 1. {% tool [Sort](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sort_header_tool/9.5+galaxy3) %} with the following parameters:
>    - {% icon param-file %} *"Sort Query"*: `Quantified Results` (output of **Datamash** {% icon tool %})
>    - *"Sort on column"*: `Column 2`
>    - *"in"*: `Descending order`
>    - *"using sort flavor"*: `Fast numeric sort`
>    - Click *"Run Tool"*
>
> 2. Rename the output file `Sorted Results`. 
>
>    > <comment-title> How to sort? </comment-title>
>    >
>    > Select column `c2` because it contains the character frequency.
>    {: .comment}
>
{: .hands_on}

> <question-title> Check your results </question-title>
>
> 1. How often did the most frequently censored character appear?
>
> > <solution-title></solution-title>
> >
> > 1. The character 敵, meaning enemy, was censored 10 times in the first text.
> >
> {: .solution}
>
{: .question}

Why would the British Hong Kong Government consistently censor this character? Jump to the conclusion to find out!

# Cut out the censored characters only

If you want to visualise your results, this step gets you there. We use the file `Censored lines` that we created in the "Ensure consistent file format" step.
Unlike the last step, this tabular file contains all 13 lines in 9 columns.
From this file, we cut out all uncensored characters. The result is a single column with multiple rows of Chinese characters.
It allows scaling words by frequency in the word cloud in the next step. As a result, characters that appear more often appear bigger, making the results evident at first sight.

> <hands-on-title> Select the censored characters </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c9`
>    - {% icon param-file %} *"From"*: `Censored lines` (output of **Compute** {% icon tool %})
>    - Click *"Run Tool"*
>
> 2. Rename the output file `Uncensored characters`.
>
>    > <details-title> What do I select here? (optional) </details-title>
>    >
>    > `c9` means column 9. It contains the uncensored characters from text two and is, therefore, cut out in this step.
>    {: .details}
>
{: .hands_on}

# Generate a word cloud

The last step is to visualise the results within a word cloud. It shows which characters were censored in the first text. The bigger the word, the more often it appeared in the text.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Generate a word cloud](toolshed.g2.bx.psu.edu/repos/bgruening/wordcloud/wordcloud/1.9.6+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `Uncensored characters` (output of **Cut** {% icon tool %})
>    - *"Do you want to select a special font?": `Select from a list of fonts`: `Noto Sans Traditional Chinese`
>    - *"Smallest font size to show"*: `8`
>    - *"Color option"*: `Colormap`
>    - *"Ratio of times to try horizontal fitting as opposed to vertical"*: `1.0`
>    - *"Scaling of words by frequency (0 - 1)"*: `0.9`
>    - Click *"Run Tool"*
>
> 2. Rename the output file `Wordcloud of uncensored characters`. 
>
>    > <details-title> Optimise your word cloud (optional) </details-title>
>    >
>    > You can choose different colours to suit your needs. The higher the "Ratio of times to try horizontal fitting as opposed to vertical" is towards "1", the more likely the character or word will appear horizontally.
>    > "Scaling of words by frequency (0 - 1)" allows you to scale the words according to their frequency. The smaller this number, the more equal-sized the characters in your word cloud will be, no matter their number.
>    {: .details}
>
{: .hands_on}

Your word cloud should look similar to this:

![Screenshot of the above Workflow in Galaxy]({% link topics/digital-humanities/tutorials/text_mining_chinese/images/Wordcloud_censored_characters.png %})


# Conclusion

This tutorial used text mining to extract censored characters from a Chinese newspaper.

![Screenshot of the above Workflow in Galaxy]({% link topics/digital-humanities/tutorials/text_mining_chinese/images/Workflow_Screenshot.jpg %} "Screenshot of the above workflow showing each of the steps explained above")

The uploaded dummy texts differed in several ways. They used slightly different punctuation, and some sentences and characters differed. The most obvious difference is that the second text was published uncensored in China, whereas the original, published in Hong Kong, contained censorship symbols. This allowed us to extract what characters were censored in the text from Hong Kong.

Within this workflow, we first unified the layout of both texts, showing one character per line for an easier comparison with the *diff* tool. The tool marked the differences between the two texts in colour. Afterwards, we extracted only lines censored with ×. The extraction of the results ran in two strands: One was counting and sorting the results. This will answer what characters the British Hong Kong Government censored in Chinese newspapers in the 1930s: Based on the (simplified) dummy texts, the characters were 敵 (enemy), 寇 (brave) and 日 (Japan). The character for *enemy* dominates and was censored five times more often than the character for *brave*.

What do those findings tell us? The British Hong Kong Government avoided publishing newspapers with a strong stand against Japan. Why? The British colony of Hong Kong, with a large Chinese population, is located very close to the Chinese mainland. Especially after the Japanese army invaded China in the summer of 1937, the British had to walk a tightrope. They tried to support the Chinese efforts without offending Japan. As a British outpost, Hong Kong had little military power and would not withstand a Japanese attack for long. Therefore, the British tried to appease the Japanese Government and avoid an attack. Calling them brave or enemies openly would have been dangerous. The one redaction of 日 (Japan) is very uncommon. This shows that the censorship practices were adaptable and not always unified. Censoring Hong Kong's newspapers to avoid anti-Japanese content is therefore a practical example of how the British Government's appeasement policies were implemented locally. This newspaper comparison is consistent with the findings in archival sources that I also researched for my dissertation ({% cite Schneider2024 %}) and lays the censored characters open for the first time.
