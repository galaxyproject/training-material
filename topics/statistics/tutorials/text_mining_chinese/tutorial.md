---
layout: tutorial_hands_on

title: Text-Mining Differences in Chinese Newspaper Articles
level: Introductory
zenodo_link: https://doi.org/10.5281/zenodo.14899614
questions:
  - How can I automatically compare two Chinese newspaper articles?
  - What characters were censored in a Chinese newspaper published in Hong Kong in the 1930s?
objectives:
  - Learn to clean and compare two texts.
  - Extract specific information out of texts.
  - Visualise your results in a word cloud.
time_estimation: 1H
key_points:
  - The *diff* tool allows comparing two similar texts automatically
  - The word cloud shows redactions from the texts at one glance
tags:
- Humanities
- Text_mining
contributions:
  authorship:
  - Sch-Da
---


The British Hong Kong Government censored Chinese newspapers before their publication in the colony in the 1930s ({% cite Ng2022 %}).
Replacement characters like × visibly marked those redactions, making them visible even to those who did not read any Chinese.

![Example of a Chinese article with symbol × marking censored characters]({% link topics/statistics/tutorials/text_mining_chinese/images/Example_x_censored_Chinese_article.svg %} "Example of a Chinese article with symbol × marking censored characters")

The schematic example, adapted from {% cite Schneider2024 %}, shows what such a censored Chinese article looked like. It is read from right to left and top to bottom. The two more prominent lines on the right are the article title, and the following text is the article's main body. It contains the character × several times, indicating various instances where it was censored.

Despite this obvious form of censorship, no research is looking into what Chinese characters the × replaced. My dissertation ({% cite Schneider2024 %}), which informs this workflow, started at this point. Through extensive research, I found several articles censored in the Hong Kong edition of Da gong bao (大公報) as uncensored versions. Those mostly came from Chinese editions printed in mainland China, where different censorship regulations applied. Those articles from China were not censored and openly showed the characters redacted in the Hong Kong versions. An example of a censored article could be {% cite Anon..16.10.1938_5598 %} and of the uncensored version {% cite TKP.16.10.1938_18864 %}.

The tutorial uses text mining to compare censored and uncensored text and to answer the following question: What characters were censored in a Chinese newspaper published in Hong Kong in the 1930s?


> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}


# Uploading the texts

The machine-readable versions of the Chinese newspaper articles I originally used in my dissertation come from a proprietary database and can not be shared here. Instead, I generated a dummy article with a similar setup in GPT and manually adapted the censorship symbols based on my research. The articles differ in style and punctuation, as is consistent with the articles in my research. Therefore, the input files are two texts in traditional Chinese. The first is censored, containing ×, and the second one is uncensored and does not contain replacement symbols. Both texts slightly differ in their layout, which will be unified later.

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
> 3. Rename the datasets
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


# Pre-processing

We pre-process and clean both texts to make the comparison easier and more apparent. This step will unify their layout. It uses Regular Expressions (Regex) to find and replace certain text parts. Here, the Regex also helps to restructure the text.


> <details-title> More on Regular Expressions (optional) </details-title>
>
> Regular Expressions are a powerful tool to modify text automatically based on word patterns. To read more about Regex's specifics, see its [documentation](https://docs.python.org/3/library/re.html).
> You can check what your input matches in [Regex 101](https://regex101.com). Enter the Regular Expression you want to try in the top field. Insert a sample text below to see what content your Regex catches and to find out how to adapt it.
>
{: .details}


## Clean up both texts

We will use Regular Expressions in a tool called "Replace text". It contains four different sub-steps. Those will vary if you upload different texts. Apply this step first to the censored and then to the uncensored text to get two cleaned ones.

> <hands-on-title> Cleaning the Text with Regular Expressions </hands-on-title>
>
> 1. {% tool [Replace Text](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"File to process"*: `output` (Input dataset)
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
>
>    > <comment-title>Explaining the above Regular Expressions</comment-title>
>    > Regular expressions can not only find particular words, as you might be familiar with from regular text editors.
>    > It is more powerful and can find particular patterns, for example, only capitalised words or all numbers.
>    > In this step, we mostly delete unnecessary placeholders.
>    > The first pattern we want to find is `\r`. It catches a specific form of invisible linebreaks that would create unwanted gaps in the comparison later.
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
> > 1. Before, the texts showed many sentences in one line. Now, both texts show one character per line and have more lines. The layout is entirely different.
> >
> {: .solution}
>
{: .question}



## Comparing the censored and uncensored text

We can now compare the two cleaned texts. This will visualise the differences between the two texts and mark them by colour. Make sure to upload the cleaned censored text with the replacement characters like ‘×’ first. As text two, upload the cleaned uncensored text without the replacement characters.  This version (HTML version) creates an HTML file, which colour codes differences as additions (green) or extractions (red) when comparing the texts.

### Create a _diff_ file for researchers

> <hands-on-title> Comparing the texts using <em>diff</em> tool </hands-on-title>
>
> 1. {% tool [diff](toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.10+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"First input file"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - {% icon param-file %} *"Second input file"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - *"Choose a report format"*: `Generates an HTML report to visualize the differences`
>
{: .hands_on}

> <question-title> Take a look at the HTML file</question-title>
>
> 1. What can you see in line 20 / 23?
>
> > <solution-title></solution-title>
> >
> > 1. Line 20  on the left shows a red comma on the left side of the table. It corresponds with line 23 on the right, which contains a colon in green. This means the punctuation in the file differs. The censored version contains a comma, while the uncensored one includes a colon.
> >
> {: .solution}
>
{: .question}

The HTML file could look like this:

![Screenshot of the diff tool comparing the censored and uncensored text]({% link topics/statistics/tutorials/text_mining_chinese/images/Diff_WF_HTML.jpg %} "Example of the HTML file comparing the censored and uncensored text")

It shows what passages differ in the two texts. Red parts show deletions and green-coloured areas are additions.
This output is very convenient for researchers, as it shows differences quickly. However, it is not helpful for further processing with Galaxy. For this, we run this tool a second time with slightly changed parameters. The output is the basis for our further analysis.


### Create a _diff_ file for further processing
This step runs the text comparison line by line again to create a raw file that the computer can work with.
It is less intuitive to understand at first glance. Again, clean the censored text with the replacement characters like ‘×’ first and the uncensored text second.

> <hands-on-title> Run another <em>diff</em> tool </hands-on-title>
>
> 1. {% tool [diff](toolshed.g2.bx.psu.edu/repos/bgruening/diff/diff/3.10+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"First input file"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - {% icon param-file %} *"Second input file"*: `outfile` (output of **Replace Text** {% icon tool %})
>    - *"Choose a report format"*: `Text file, side-by-side (-y)`
>
>    > <comment-title> The output format </comment-title>
>    > As you can see above, the output of this file is a text file, compared to the HTML in the last step.
>    {: .comment}
>
{: .hands_on}


> <question-title>Look at this step's output</question-title>
>
> 1. How does this file differ from the HTML file in the last step?
>
> > <solution-title>Answer</solution-title>
> >
> > 1. This output is a raw file and shows several columns. Changes between the two texts are not coloured.
> > Instead, they are marked by various symbols in column 8.
> >
> {: .solution}
>
{: .question}


# Select only censored lines
In the next step, we want to extract specific lines only. To determine what content was redacted in the first text, we filter the last step's raw output for lines containing the censorship symbol ×.


> <hands-on-title> Filter text </hands-on-title>
>
> 1. {% tool [Filter](Filter1) %} with the following parameters:
>    - {% icon param-file %} *"Filter"*: `diff_file` (output of **diff** {% icon tool %})
>    - *"With following condition"*: `ord(c1) == 215`
>
>
>
>    > <details-title> How to select the correct characters here (optional) </details-title>
>    >
>    > The condition "ord(c1) == 215" means that column c1's lines, which contain the censored text, are selected if they match ×. The symbol × is unspecific. Therefore, the Unicode identifier of the character (215) is used for clarity in this condition.
>    {: .details}
>
>
>    > <comment-title> Filter for other characters </comment-title>
>    > Add another Unicode here if you want to select a different character, for example, '□' or '△'.
>    > For example, you can get the respective code on [Character Code Finder](https://www.mauvecloud.net/charsets/CharCodeFinder.html).
>    > Copy the character you want to filter for in the "input" window and select "Decimal Character Codes" as an output. If you do this for symbol ×, you get 215.
>    {: .comment}
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


# Ensure consistent file format

After filtering for the censored lines, we insert a sub-step to ensure smooth computing. The previous setup could cause an error if the characters filtered in the last step were erased. Then, the extracted file would miss the last column, which would cause an error. This is invisible to the researchers in the file. The compute step covers this potential error and ensures all necessary columns exist.

> <hands-on-title> Compute to ensure all columns exist </hands-on-title>
>
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `out_file1` (output of **Filter** {% icon tool %})
>    - *"Input has a header line with column names?"*: `No`
>        - In *"Expressions"*:
>            - {% icon param-repeat %} *"Insert Expressions"*
>                - *"Add expression"*: `c9`
>                - *"Mode of the operation"*: `Replace`
>                    - *"Use new column to replace column number"*: `9`
>    - In *"Error handling"*:
>        - *"If an expression cannot be computed for a row"*: `Produce an empty column value for the row`
>
{: .hands_on}


# Summarise your findings

This step sums up how often each character appeared in the table before.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Datamash](toolshed.g2.bx.psu.edu/repos/iuc/datamash_ops/datamash_ops/1.8+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input tabular dataset"*: `out_file1` (output of **Compute** {% icon tool %})
>    - *"Group by fields"*: `9`
>    - *"Sort input"*: `Yes`
>    - In *"Operation to perform on each group"*:
>        - {% icon param-repeat %} *"Insert Operation to perform on each group"*
>            - *"On column"*: `c9`
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
> 1. {% tool [Sort](sort1) %} with the following parameters:
>    - {% icon param-file %} *"Sort Dataset"*: `out_file` (output of **Datamash** {% icon tool %})
>    - *"on column"*: `c2`
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

Why would the British Hong Kong Government consistently censor this character? Jump to the conclusion to find out.

# Cut out the censored characters only

If you want to visualise your results, this step gets you there. We select only the uncensored characters from text two. The result is only one column with different rows of Chinese characters.
It allows scaling words by frequency in the word cloud in the next step. As a result, characters that appear more often appear bigger, making the results evident at first sight.

> <hands-on-title> Select the censored characters </hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} with the following parameters:
>    - *"Cut columns"*: `c9`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Compute** {% icon tool %})
>
>    > <details-title> What do I select here? (optional) </details-title>
>    >
>    > `c9` means column 9. It contains the uncensored characters from text two and is, therefore, cut out in this step.
>    {: .details}
>
{: .hands_on}

# Generate a word cloud

The last step is to visualise the results within a word cloud. It shows, which characters were censored in the first text. The bigger the word, the more often it appeared in the text.

> <hands-on-title> Task description </hands-on-title>
>
> 1. {% tool [Generate a word cloud](toolshed.g2.bx.psu.edu/repos/bgruening/wordcloud/wordcloud/1.9.4+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `out_file1` (output of **Cut** {% icon tool %})
>    - *"Do you want to select a special font?": `Select from a list of fonts`: `Noto Sans Traditional Chinese`
>    - *"Smallest font size to use"*: `8`
>    - *"Color option"*: `Color`
>    - *"Ratio of times to try horizontal fitting as opposed to vertical"*: `1.0`
>    - *"Scaling of words by frequency (0 - 1)"*: `0.9`
>
>    > <details-title> Optimise your word cloud (optional) </details-title>
>    >
>    > You can choose different colours to suit your needs. The higher the "Ratio of times to try horizontal fitting as opposed to vertical" is towards "1", the more likely the character or word will appear horizontally.
>    > "Scaling of words by frequency (0 - 1)" allows you to scale the words according to their amount. The smaller this number, the more equal-sized the characters in your word cloud will be, no matter their amount.
>    {: .details}
>
{: .hands_on}

Your word cloud should look similar to this:

![Screenshot of the above Workflow in Galaxy]({% link topics/statistics/tutorials/text_mining_chinese/images/Wordcloud_censored_characters.png %})


# Conclusion

This tutorial used text mining to extract censored characters from a Chinese newspaper.

![Screenshot of the above Workflow in Galaxy]({% link topics/statistics/tutorials/text_mining_chinese/images/Workflow_Screenshot.jpg %} "Screenshot of the above workflow showing each of the steps explained above")

The uploaded dummy texts contained several differences. They used slightly different punctuation, and some sentences and characters differed. The most obvious difference is that the second text was published uncensored in China, while the original text, published in Hong Kong, contained censorship symbols. This allowed us to extract what characters were censored in the text from Hong Kong.

Within this workflow, we first unified the layout of both texts, showing one character per line for an easier comparison with *diff* tool. The tool marked the differences between both texts in colour. Afterwards, we extracted only lines censored with ×. The extraction of the results ran in two strands: One was counting and sorting the results. This will answer what characters the British Hong Kong Government censored in their Chinese newspapers in the 1930s: Based on the (simplified) dummy texts, the characters were 敵 (enemy), 寇 (brave) and 日 (Japan). The character for *enemy* dominates and was censored five times more often than the character for *brave*.

What do those findings tell us? The British Hong Kong Government avoided publishing newspapers with a strong stand against Japan. Why? Because the British colony Hong Kong, with a large Chinese population, is located very close to the Chinese mainland. Especially after the Japanese army invaded China in the summer of 1937, the British had to walk a tightrope. They tried to support the Chinese efforts without offending Japan. As a British outpost, Hong Kong had little military power and would not withstand a Japanese attack for long. Therefore, the British tried to appease the Japanese Government and avoid an attack. Calling them brave or enemy openly would have been dangerous. The one redaction of 日 (Japan) is very uncommon. This shows that the censorship practices were adaptable and not always unified. Censoring Hong Kong's newspapers to avoid anti-Japanese content is, therefore, a practical example of how appeasement policies from the British Government were implemented locally. This newspaper comparison is consistent with the findings in archival sources that I also researched for my dissertation ({% cite Schneider2024 %}) and lays the censored characters open for the first time.
