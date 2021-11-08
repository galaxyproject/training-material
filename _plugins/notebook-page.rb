require 'json'
require 'mkmf'
require 'fileutils'
require 'yaml'
require "kramdown"

module MakeMakefile::Logging
  @logfile = File::NULL
end

module Jekyll
  class JupyterNotebookGenerator < Generator
    safe true

    def convertToNotebook(tutorial, language)
      notebook_json = `notedown --match=#{language} #{tutorial}`
      JSON.parse(notebook_json)
    end

    def extractAndFixMetadataCell(notebook)
      # Strip out %%bash since we'll use the bash kernel
      metadata = nil
      contributors = nil

      File.open('CONTRIBUTORS.yaml', 'r') do |f2|
        contributors = YAML.load(f2.read)
      end

      notebook['cells'].map.with_index {|cell, i|
        if i == 0
          metadata = YAML.load(cell['source'].join(''))
          offset = cell['source'].slice(1..-1).index("---\n")

          by_line = metadata['contributors'].map{|c|
            "[#{contributors.fetch(c, {"name" => c}).fetch('name', c)}](https://training.galaxyproject.org/hall-of-fame/#{c}/)"
          }.join(", ")

          meta_header = [
            "<div style=\"border: 2px solid #8A9AD0; margin: 1em 0.2em; padding: 0.5em;\">\n\n",
            "# #{metadata['title']}\n",
            "\n",
            "by #{by_line}\n",
            "\n",
            "**Objectives**\n",
            "\n",
          ] + metadata['questions'].map{|q| "- #{q}\n"} + [
            "\n",
            "**Objectives**\n",
            "\n",
          ] + metadata['objectives'].map{|q| "- #{q}\n"} + [
            "\n",
            "**Time Estimation: #{metadata['time_estimation']}**\n",
            "\n",
            "</div>\n",
          ]
          cell['source'] = meta_header + cell['source'].slice(offset + 2..-1)
        end
        cell
      }
      return notebook, metadata
    end

    def fixRNotebook(notebook)
      # Set the bash kernel
      notebook['metadata'] = {
        "kernelspec" => {
         "display_name" => "R",
         "language" => "R",
         "name" => "r"
        },
        "language_info" => {
         "codemirror_mode" => "r",
         "file_extension" => ".r",
         "mimetype" => "text/x-r-source",
         "name" => "R",
         "pygments_lexer" => "r",
         "version" => "4.1.0"
        }
      }
      # Strip out %%R since we'll use the bash kernel
      notebook['cells'].map{|cell|
        if cell.fetch('cell_type') == 'code'
          if cell['source'][0] == "%%R\n"
            cell['source'] = cell['source'].slice(1..-1)
          end
        end
        cell
      }
      notebook
    end

    def fixBashNotebook(notebook)
      # Set the bash kernel
      notebook['metadata'] = {
        "kernelspec" =>  {
          "display_name" =>  "Bash",
          "language" =>  "bash",
          "name" =>  "bash"
        },
        "language_info" =>  {
          "codemirror_mode" =>  "shell",
          "file_extension" =>  ".sh",
          "mimetype" =>  "text/x-sh",
          "name" =>  "bash"
        }
      }
      # Strip out %%bash since we'll use the bash kernel
      notebook['cells'].map{|cell|
        if cell.fetch('cell_type') == 'code'
          if cell['source'][0] == "%%bash\n"
            cell['source'] = cell['source'].slice(1..-1)
          end
        end
        cell
      }
      notebook
    end

    def fixSqlNotebook(notebook)
      # Add in a %%sql at the top of each cell
      notebook['cells'].map{|cell|
        if cell.fetch('cell_type') == 'code'
          if cell['source'].join('').index('load_ext').nil?
            cell['source'] = ["%%sql\n"] + cell['source']
          end
        end
        cell
      }
      notebook
    end

    def markdownify(site, text)
      site.find_converter_instance(
        Jekyll::Converters::Markdown
      ).convert(text.to_s)
    end

    def renderMarkdownCells(site, notebook, metadata, page)
      colors = {
        "overview" => "#8A9AD0",
        "agenda" => "#86D486; display: none",
        "keypoints" => "#FFA1A1",
        "tip" => "#FFE19E",
        "warning" => "#de8875",
        "comment" => "#ffecc1",
        "hands_on" => "#dfe5f9",
        "question" => "#8A9AD0",
        "solution" => "#B8C3EA; color: white",
        "details" => "#ddd",
        "feedback" => "#86D486",
        "code-in" => "#86D486",
        "code-out" => "#fb99d0",
      }

      seen_abbreviations = Hash.new

      notebook['cells'].map{|cell|
        if cell.fetch('cell_type') == 'markdown'

          # The source is initially a list of strings, we'll merge it together
          # to make it easier to work with.
          source = cell['source'].join("")

          # Here we replace individual `s with codeblocks, they screw up
          # rendering otherwise by going through rouge
          source = source.gsub(/ `([^`]*)`([^`])/, ' <code>\1</code>\2')
            .gsub(/([^`])`([^`]*)` /, '\1<code>\2</code> ')

          # Replace all the broken icons that can't render, because we don't
          # have access to the full render pipeline.
          cell['source'] = markdownify(site, source)
            .gsub(/{% icon tip %}/, '💡')
            .gsub(/{% icon code-in %}/, '⌨️')
            .gsub(/{% icon code-out %}/, '🖥')
            .gsub(/{% icon question %}/, '❓')
            .gsub(/{% icon solution %}/, '👁')
            .gsub(/{% icon warning %}/, '⚠️')
            .gsub(/{% icon comment %}/, '💬')
            .gsub(/{% icon hands_on %}/, '✏️')

          if metadata.key?('abbreviations')
            metadata['abbreviations'].each{|abbr, defn|
              cell['source'].gsub(/\{#{abbr}\}/) {
                if seen_abbreviations.key?(abbr) then
                  firstdef = false
                else
                  firstdef = true
                  seen_abbreviations[abbr] = true
                end

                if firstdef then
                  "#{defn} (#{abbr})"
                else
                  "<abbr title=\"#{defn}\">#{abbr}</abbr>"
                end
              }
            }
          end

          # Here we give a GTN-ish styling that doesn't try to be too faithful,
          # so we aren't spending time keeping up with changes to GTN css,
          # we're making it 'our own' a bit.

          colors.each{ |key, val|
            cell['source'].gsub!(/<blockquote class="#{key}">/, "<blockquote class=\"#{key}\" style=\"border: 2px solid #{val}; margin: 1em 0.2em\">")
          }

          # Images are referenced in the GTN through relative URLs which is
          # fab, but in a notebook this doesn't make sense as it will live
          # outside of the GTN. We need real URLs.
          #
          # But we'll only do this if it's in "deploy" mode, in order to have
          # the images also work locally.
          if site.config['url'] == "https://training.galaxyproject.org"
            cell['source'].gsub!(/<img src=\"\.\./, '<img src="' + site.config['url'] + site.config['baseurl'] + page.url.split('/')[0..-2].join('/') + '/..')
          end

          # Strip out the highlighting as it is bad on some platforms.
          cell['source'].gsub!(/<pre class="highlight">/, '<pre style="color: inherit; background: white">')
          cell['source'].gsub!(/<div class="highlight">/, '<div>')

          # add a 'hint' to the solution boxes which have blanked out text.
          cell['source'].gsub!(/(<h3 id="-icon-solution--solution">)/, '<div style="color: #555; font-size: 95%;">Hint: Select the text with your mouse to see the answer</div>\1')

          # There is some weirdness in the processing of $s in Jupyter. After a
          # certain number of them, it will give up, and just render everything
          # like with a '<pre>'. We remove this to prevent that result.
          cell['source'].gsub!(/^\s*</, '<')
          # Additionally leading spaces are sometimes interpreted as <pre>s and
          # end up causing paragraphs to be rendered as code. So we wipe out
          # all leading space.
          cell['source'].gsub!(/\$/, '&#36;')
        end
        cell
      }
      notebook
    end


    def generate(site)
      if find_executable('notedown').nil?
        puts "[GTN/Notebooks] We could not find the notedown executable, so, notebooks will not be rendered."
        return
      end

      # For every tutorial with the 'notebook' key in the page data
      site.pages.select{|page| page.data['layout'] == 'tutorial_hands_on' and page.data.has_key?('notebook')}.each do |page|
        # We get the path to the tutorial source
        dir = File.dirname(File.join('.', page.url))
        fn = File.join('.', page.url).sub(/html$/, 'md')
        notebook_language = page.data['notebook'].fetch('language', 'python')

        puts "[GTN/Notebooks] Rendering #{notebook_language} #{fn}"

        # Here we read use `notedown` to convert the tutorial to a Hash
        # representing the notebook
        notebook = convertToNotebook(fn, notebook_language)

        # This extracts the metadata yaml header and does manual formatting of
        # the header data to make for a nicer notebook.
        notebook, metadata = extractAndFixMetadataCell(notebook)

        # Apply language specific conventions
        if notebook_language == 'bash'
          notebook = fixBashNotebook(notebook)
        elsif notebook_language == 'sql'
          notebook = fixSqlNotebook(notebook)
        elsif notebook_language == 'r'
          notebook = fixRNotebook(notebook)
        end


        # Here we loop over the markdown cells and render them to HTML. This
        # allows us to get rid of classes like {: .tip} that would be left in
        # the output by Jupyter's markdown renderer, and additionally do any
        # custom CSS which only seems to work when inline on a cell, i.e. we
        # can't setup a style block, so we really need to render the markdown
        # to html.
        notebook = renderMarkdownCells(site, notebook, metadata, page)

        # Here we add a close to the notebook
        notebook['cells'] = notebook['cells'] + [{
          "cell_type" => "markdown",
          "id" => "final-ending-cell",
          "metadata" => {},
          "source" => [
            "# Key Points\n\n",
          ] + metadata['key_points'].map{|k| "- #{k}\n"} + [
            "\n# Congratulations on successfully completing this tutorial!\n\n",
              "Please [fill out the feedback on the GTN website](https://training.galaxyproject.org/training-material#{page.url}#feedback) and check there for further resources!\n"
          ]
        }]

        # Create the JSON file and inject the data
        page2 = PageWithoutAFile.new(site, "", dir, "tutorial.md.ipynb")
        page2.content = JSON.pretty_generate(notebook)
        page2.data["layout"] = nil
        site.pages << page2
      end
    end
  end
end
