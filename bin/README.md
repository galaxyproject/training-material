## To generate a 1st draft for a tutorial paper

1. Get `format-tutorial-for-paper` branch on the training material
2. Install Latex, bibtex
    - Ubuntu: `sudo apt install texlive-latex-extra`
3. Run `make create-pandoc-env`
4. Run `make format-for-article TUTO:=topics/<topic_name>/tutorials/<tutorial_name> JOURNAL:=<f1000research|plos>` by replacing
    - `<topic_name>` by the name of the topic folder
    - `<tutorial_name>` by the name of the tutorial folder
    - `<f1000research|plos>` by the name of the journal to submit to: F1000 Research or PlOS

    It will create a `article` folder in the tutorial folder with several files inside:

    1. Run a 1st Python script to create article folder with `article_1.md` formatted for pandoc
    2. Run Pandoc using Plos template to generate 1st latex file (`article_2.tex`)
    3. Run pdflatex and bibtex
    4. Run a 2nd Python script to create `article_3.tex` with some extra formatting and copy of references
    5. Run pdflatex again to generate `article_3.pdf`

5. Check the `article_3.pdf`
6. Make any changes to `article_3.tex`

    A TODO list for formatting is written on the top of the file

7. Recompile the latex with `pdflatex` (inside the article folder)

 