require 'fileutils'
require 'yaml'
require "kramdown"
require 'json'

module Jekyll
  class RmarkdownGenerator < Generator
    safe true

    def group_doc_by_first_char(data)
      out = []
      first_char = nil
      val = []
      data = data.split("\n")

      # Here we collapse running groups of `>` into single blocks.
      data.each{|line|
        if first_char.nil?
          first_char = line[0]
          val = [line]
        else
          if line[0] == first_char
            val.append(line)
          elsif line[0..1] == '{:' && first_char == '>'
            val.append(line)
          else
            # flush
            out.append(val)
            if line.size > 0
              first_char = line[0]
            else
              first_char = ""
            end
            val = [line]
          end
        end
      }
      # final flush
      out.append(val)

      out.select!{|v|
        !(v[0][0] == '>' && v[-1][0..1] == '{:' && v[-1].match(/.agenda/))
      }
      out.map!{|v|
        if v[0][0] == '>' && v[-1][0..1] == '{:'
          cls = v[-1][2..-2].strip
          res = [":::{#{cls}}"]
          res += v[0..-2].map{|c| c.sub(/^>\s*/, '')}
          res += [":::"]
          res
        else
          v
        end
      }

      out.flatten(1).join("\n")
    end

    def generate(site)

      # For every tutorial with the 'notebook' key in the page data
      site.pages.select{|page| page.data['layout'] == 'tutorial_hands_on'}
        .select{|page| page.data.has_key?('notebook')}
        .select{|page| page.data['notebook']['language'].downcase == "r" }
        .each do |page|

        # We get the path to the tutorial source
        dir = File.dirname(File.join('.', page.url))
        fn = File.join('.', page.url).sub(/html$/, 'Rmd')

        by_line = page.data['contributors'].map{|c|
          "[#{site.data['contributors'].fetch(c, {"name" => c}).fetch('name', c)}](https://training.galaxyproject.org/hall-of-fame/#{c}/)"
        }.join(", ")

        # Replace top level `>` blocks with fenced `:::`
        content = group_doc_by_first_char(page.content)

        # Re-run a second time to catch singly-nested Q&A?
        content = group_doc_by_first_char(content)
            .gsub(/{% icon tip %}/, '💡')
            .gsub(/{% icon code-in %}/, '⌨️')
            .gsub(/{% icon code-out %}/, '🖥')
            .gsub(/{% icon question %}/, '❓')
            .gsub(/{% icon solution %}/, '👁')
            .gsub(/{% icon warning %}/, '⚠️')
            .gsub(/{% icon comment %}/, '💬')
            .gsub(/{% icon hands_on %}/, '✏️')

        content = content + %Q(\n\n# References\n\n<div id="refs"></div>\n)

        # https://raw.githubusercontent.com/rstudio/cheatsheets/master/rmarkdown-2.0.pdf
        # https://bookdown.org/yihui/rmarkdown/

        puts "[GTN/Notebooks/R] Rendering RMarkdown #{fn}"
        fnparts = fn.split('/')
        rmddata = {
          'title' => page.data['title'],
          'author' => "#{by_line}, #{page.data.fetch('license', 'CC-BY')} licensed content from the [Galaxy Training Network](https://training.galaxyproject.org/)",
          'bibliography' => "#{fnparts[2]}-#{fnparts[4]}.bib",
          'output' => {
            'html_notebook' => {
              'toc' => true,
              'toc_depth' => 2,
              'css' => "gtn.css",
              'toc_float' => {
                "collapsed" => false,
                "smooth_scroll" => false,
              },
              'theme' => {'bootswatch' => 'journal'}
            },
            'word_document' => {
              'toc' => true,
              'toc_depth' => 2,
              'latex_engine' => 'xelatex',
            },
            'pdf_document' => {
              'toc' => true,
              'toc_depth' => 2,
              'latex_engine' => 'xelatex',
            },
          },
          'date' => begin page.last_modified.to_s rescue Time.new.to_s end,
          'link-citations' => true,
          'anchor_sections' => true,
          'code_download' => true,
        }
        rmddata['output']['html_document'] = JSON::parse(JSON::generate(rmddata['output']['html_notebook']))

        final_content = [
          "# Introduction\n",
          content.gsub(/```r/, "```{r}"),
          "# Key Points\n",
        ] + page.data['key_points'].map{|k| "- #{k}"} + [
            "\n# Congratulations on successfully completing this tutorial!\n",
              "Please [fill out the feedback on the GTN website](https://training.galaxyproject.org/training-material#{page.url}#feedback) and check there for further resources!\n"
        ]

        page2 = PageWithoutAFile.new(site, "", dir, "tutorial.Rmd")
        page2.content = rmddata.to_yaml(:line_width => rmddata['author'].size + 10) + "---\n" + final_content.join("\n")
        page2.data["layout"] = nil
        page2.data["citation_target"] = 'R'
        site.pages << page2
      end
    end
  end
end
