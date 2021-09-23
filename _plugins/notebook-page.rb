require 'json'
require "kramdown"

module Jekyll
  class JupyterNotebookGenerator < Generator
    safe true

    def convertToNotebook(tutorial, language)
      notebook_json = `notedown --match=#{language} #{tutorial}`
      JSON.parse(notebook_json)
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

    def markdownify(site, text)
      site.find_converter_instance(
        Jekyll::Converters::Markdown
      ).convert(text.to_s)
    end

    def renderMarkdownCells(site, notebook)
      notebook['cells'].map{|cell|
        if cell.fetch('cell_type') == 'markdown'
          cell['source'] = markdownify(site, cell['source'].join(""))
        end
        cell
      }
      notebook
    end


    def generate(site)
      # layout: tutorial_slides
      # layout: base_slides

      site.pages.select{|page| page.data['layout'] == 'tutorial_hands_on' and page.data.has_key?('notebook')}.each do |page|
        dir = File.dirname(File.join('.', page.url))
        fn = File.join('.', page.url).sub(/html$/, 'md')
        notebook_language = page.data['notebook'].fetch('language', 'python')

        notebook = convertToNotebook(fn, notebook_language)
        # Apply language specific conventions
        if notebook_language == 'bash'
          notebook = fixBashNotebook(notebook)
        end

        notebook = renderMarkdownCells(site, notebook)

        # Create the JSON file and inject the data
        f = File.new("_site/training-material/#{dir}/tutorial.md.ipynb", "w+")
        f.puts(JSON.generate(notebook))

      end
    end
  end
end
