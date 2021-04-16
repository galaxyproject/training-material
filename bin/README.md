# Generate a draft for a tutorial paper

## How to setup the environmment?

1. Get `format-tutorial-for-paper` branch on the training material
2. Install Latex, bibtex
    - Ubuntu: `sudo apt install texlive-latex-extra`
3. Run `make create-pandoc-env`

## How to generate a 1st draft for a tutorial paper?

1. Run `make format-for-article TUTO:=topics/<topic_name>/tutorials/<tutorial_name> JOURNAL:=<f1000research|plos>` by replacing
    - `<topic_name>` by the name of the topic folder
    - `<tutorial_name>` by the name of the tutorial folder
    - `<f1000research|plos>` by the name of the journal to submit to: F1000 Research or PlOS

    It will create a `article` folder in the tutorial folder with several files inside

2. Check the `article_3.pdf`
3. Make any changes to `article_3.tex`

    A TODO list for formatting is written on the top of the file

4. Recompile the latex with `pdflatex` (inside the article folder)

## How does that work?

The command `make format-for-article`:

1. Run a 1st Python script (`bin/format_tutorial_for_article.py`) to create article folder with `article_1.md` formatted for pandoc

    To test it, you need to

    1. Activate usual GTN conda environment (`conda activate galaxy_training_material`)
    2. Run it from the root of the repository

2. Create the tex files from `article_1.md`
    
    1. Run Pandoc using Plos template to generate 1st latex file (`article_2.tex`) using the command

        ``` 
        $ pandoc -f markdown -t latex --template ${ROOT}/${ARTICLE_TEMPLATE}/${JOURNAL}/template.tex -o article_2.tex article_1.md
        ```

    2. Copy the template files (from `share/article-template` folder) to `article` folder in the tutorial folder 
    2. Run pdflatex and bibtex on `article_2` files to 

    To test these steps, you need to 

    1. Activate a dedicated conda environment (`conda activate gtn_paper`)
    2. Move into the `article` folder in the tutorial folder 

3. Run a 2nd Python script (`bin/format_article_latex.py`) to create, from `article_2.tex`, `article_3.tex` with some extra formatting and copy of references

    To test it, you need to

    1. Activate usual GTN conda environment (`conda activate galaxy_training_material`)

5. Run pdflatex again to generate `article_3.pdf`

