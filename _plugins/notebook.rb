require 'digest'
require 'json'
require 'fileutils'
require 'yaml'
require 'base64'


module GTNNotebooks
  COLORS = {
    'overview' => '#8A9AD0',
    'agenda' => '#86D486',
    'keypoints' => '#FFA1A1',
    'tip' => '#FFE19E',
    'warning' => '#de8875',
    'comment' => '#ffecc1',
    'hands_on' => '#dfe5f9',
    'question' => '#8A9AD0',
    'solution' => '#B8C3EA',
    'details' => '#ddd',
    'feedback' => '#86D486',
    'code-in' => '#86D486',
    'code-out' => '#fb99d0',
  }
  COLORS_EXTRA = {
    'agenda' => 'display: none',
  }

  ICONS = {
    'tip' => '💡',
    'code-in' => '⌨️',
    'code-out' => '🖥',
    'question' => '❓',
    'solution' => '👁',
    'warning' => '⚠️',
    'comment' => '💬',
    'feedback' => '⁉️',
    'details' => '💬',
    'hands_on' => '✏️',
  }

  ICONS_FA = {
    "far fa-keyboard" => "code-in",
    "fas fa-laptop-code" => "code-out",
    "far fa-comment-dots" => "comment",
    "fas fa-info-circle" => "details",
    "far fa-comments" => "feedback",
    "fas fa-pencil-alt" => "hands_on",
    "fas fa-pencil-alt" => "hands_on",
    "far fa-question-circle" => "question",
    "far fa-eye" => "solution",
    "far fa-lightbulb" => "tip",
    "fas fa-exclamation-triangle" => "warning",
  }

  def self.generate_css
    COLORS.map{|key, val|
      ".#{key} { padding: 0 1em; margin: 1em 0.2em; border: 2px solid #{val} }"
    }.join("\n")
  end

  def self.convert_notebook_markdown(content, language)
    out = []
    inside_block = false
    val = []
    data = content.split("\n")
    data.each.with_index{|line, i|
      if line == "```#{language}"
        if inside_block
          puts data[i-2..i+2]
          raise "[GTN/Notebook] L#{i} Error! we're already in a block:"
        end
        # End the previous block
        out.push([val, inside_block])
        val = []

        inside_block = true
      elsif inside_block && line == '```'
        # End of code block
        out.push([val, inside_block])
        val = []
        inside_block = false
      else
        val.push(line)
      end
    }
    # final flush
    if ! val.nil?
      out.push([val, inside_block])
    end

    notebook = {
     "metadata" => {},
     "nbformat" => 4,
     "nbformat_minor" => 5,
    }

    notebook['cells'] = out.map.with_index{|data, index|
      res = {
        'id' => "cell-#{index}",
        'source' => data[0].map{|x| x.rstrip + "\n"}
      }
      # Strip the trailing newline in the last cell.
      if res['source'].length > 0
        res['source'][-1] = res['source'][-1].rstrip
      end

      # Remove any remaining language tagged code blocks, e.g. in
      # tip/solution/etc boxes. These do not render well.
      res['source'] = res['source'].map{|x| x.gsub(/```#{language}/, '```')}

      if data[1]
        res = res.update({
          'cell_type' => 'code',
          'execution_count' => nil,
          'outputs' => [],
          'metadata' => {
            'attributes' => {
              'classes' => [
                language
              ],
              'id' => '',
            }
          }
        })
      else
        res['cell_type'] = 'markdown'
      end
      res
    }
    notebook
  end

  def self.group_doc_by_first_char(data)
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
          val.push(line)
        elsif line[0..1] == '{:' && first_char == '>'
          val.push(line)
        else
          # flush
          out.push(val)
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
    out.push(val)

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

  def self.construct_byline(metadata)
    if metadata.has_key?('contributors')
      folks = metadata['contributors']
    else
      folks = metadata['contributions'].map{|k, v| v }.flatten
    end

    contributors = nil
    File.open('CONTRIBUTORS.yaml', 'r') do |f2|
      contributors = YAML.load(f2.read)
    end

    folks.map{|c|
      "[#{contributors.fetch(c, {"name" => c}).fetch('name', c)}](https://training.galaxyproject.org/hall-of-fame/#{c}/)"
    }.join(", ")
  end

  def self.add_metadata_cell(notebook, metadata)
    by_line = self.construct_byline(metadata)

    meta_header = [
      "<div style=\"border: 2px solid #8A9AD0; margin: 1em 0.2em; padding: 0.5em;\">\n\n",
      "# #{metadata['title']}\n",
      "\n",
      "by #{by_line}\n",
      "\n",
      "#{metadata.fetch('license', 'CC-BY')} licensed content from the [Galaxy Training Network](https://training.galaxyproject.org/)\n",
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
    metadata_cell = {
      'id' => 'metadata',
      'cell_type' => 'markdown',
      'source' => meta_header
    }
    notebook['cells'].unshift(metadata_cell)
    return notebook
  end

  def self.fixRNotebook(notebook)
    # Set the bash kernel
    notebook['metadata'] = {
      'kernelspec' => {
       'display_name' => 'R',
       'language' => 'R',
       'name' => 'r'
      },
      'language_info' => {
       'codemirror_mode' => 'r',
       'file_extension' => '.r',
       'mimetype' => 'text/x-r-source',
       'name' => 'R',
       'pygments_lexer' => 'r',
       'version' => '4.1.0'
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

  def self.fixBashNotebook(notebook)
    # Set the bash kernel
    notebook['metadata'] = {
      'kernelspec' =>  {
        'display_name' =>  'Bash',
        'language' =>  'bash',
        'name' =>  'bash'
      },
      'language_info' =>  {
        'codemirror_mode' =>  'shell',
        'file_extension' =>  '.sh',
        'mimetype' =>  'text/x-sh',
        'name' =>  'bash'
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

  def self.fixSqlNotebook(notebook)
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

  def self.markdownify(site, text)
    begin
      site.find_converter_instance(
        Jekyll::Converters::Markdown
      ).convert(text.to_s)
    rescue
      require 'kramdown'
      Kramdown::Document.new(text).to_html
    end
  end

  def self.notebook_filter(data, language=nil)
    data['layout'] == 'tutorial_hands_on' \
      and data.has_key?('notebook') \
      and (language.nil? or data['notebook']['language'].downcase == language)
  end

  def self.render_rmarkdown(page_data, page_content, page_url, page_last_modified, fn)
    by_line = self.construct_byline(page_data)

    # Replace top level `>` blocks with fenced `:::`
    content = group_doc_by_first_char(page_content)

    # Re-run a second time to catch singly-nested Q&A?
    content = group_doc_by_first_char(content)

    ICONS.each{ |key, val|
      content.gsub!(/{% icon #{key} %}/, val)
    }
    ICONS_FA.each{ |key, val|
      content.gsub!(/<i class="#{key}" aria-hidden="true"><\/i>/, ICONS[val])
    }

    content = content + %Q(\n\n# References\n\n<div id="refs"></div>\n)

    # https://raw.githubusercontent.com/rstudio/cheatsheets/master/rmarkdown-2.0.pdf
    # https://bookdown.org/yihui/rmarkdown/

    fnparts = fn.split('/')
    rmddata = {
      'title' => page_data['title'],
      'author' => "#{by_line}, #{page_data.fetch('license', 'CC-BY')} licensed content from the [Galaxy Training Network](https://training.galaxyproject.org/)",
      'bibliography' => "#{fnparts[2]}-#{fnparts[4]}.bib",
      'output' => {
        'html_notebook' => {
          'toc' => true,
          'toc_depth' => 2,
          'css' => 'gtn.css',
          'toc_float' => {
            'collapsed' => false,
            'smooth_scroll' => false,
          },
          # 'theme' => {'bootswatch' => 'journal'}
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
      'date' => page_last_modified.to_s,
      'link-citations' => true,
      'anchor_sections' => true,
      'code_download' => true,
    }
    rmddata['output']['html_document'] = JSON::parse(JSON::generate(rmddata['output']['html_notebook']))

    final_content = [
      "# Introduction\n",
      content.gsub(/```r/, "```{r}"),
      "# Key Points\n",
    ] + page_data['key_points'].map{|k| "- #{k}"} + [
        "\n# Congratulations on successfully completing this tutorial!\n",
          "Please [fill out the feedback on the GTN website](https://training.galaxyproject.org/training-material#{page_url}#feedback) and check there for further resources!\n"
    ]

    rmddata.to_yaml(:line_width => rmddata['author'].size + 10) + "---\n" + final_content.join("\n")
  end

  def self.render_jupyter_notebook(data, content, url, last_modified, notebook_language, site, dir)
    # Here we read use internal methods to convert the tutorial to a Hash
    # representing the notebook
    notebook = self.convert_notebook_markdown(content, notebook_language)

    # This extracts the metadata yaml header and does manual formatting of
    # the header data to make for a nicer notebook.
    notebook = self.add_metadata_cell(notebook, data)

    # Apply language specific conventions
    if notebook_language == 'bash'
      notebook = self.fixBashNotebook(notebook)
    elsif notebook_language == 'sql'
      notebook = self.fixSqlNotebook(notebook)
    elsif notebook_language == 'r'
      notebook = self.fixRNotebook(notebook)
    end

    # Here we loop over the markdown cells and render them to HTML. This
    # allows us to get rid of classes like {: .tip} that would be left in
    # the output by Jupyter's markdown renderer, and additionally do any
    # custom CSS which only seems to work when inline on a cell, i.e. we
    # can't setup a style block, so we really need to render the markdown
    # to html.
    notebook = self.renderMarkdownCells(site, notebook, data, url, dir)

    # Here we add a close to the notebook
    notebook['cells'] = notebook['cells'] + [{
      "cell_type" => "markdown",
      "id" => "final-ending-cell",
      "metadata" => {"editable" => false, "collapsed" => false},
      "source" => [
        "# Key Points\n\n",
      ] + data['key_points'].map{|k| "- #{k}\n"} + [
        "\n# Congratulations on successfully completing this tutorial!\n\n",
          "Please [fill out the feedback on the GTN website](https://training.galaxyproject.org/training-material#{url}#feedback) and check there for further resources!\n"
      ]
    }]
    notebook
  end

  def self.renderMarkdownCells(site, notebook, metadata, page_url, dir)
    seen_abbreviations = Hash.new
    notebook['cells'].map{|cell|
      if cell.fetch('cell_type') == 'markdown'

        # The source is initially a list of strings, we'll merge it together
        # to make it easier to work with.
        source = cell['source'].join("").strip

        # Here we replace individual `s with codeblocks, they screw up
        # rendering otherwise by going through rouge
        source = source.gsub(/ `([^`]*)`([^`])/, ' <code>\1</code>\2')
          .gsub(/([^`])`([^`]*)` /, '\1<code>\2</code> ')

        # Strip out includes
        source = source.gsub(/{% include .* %}/, '')

        # Replace all the broken icons that can't render, because we don't
        # have access to the full render pipeline.
        cell['source'] = self.markdownify(site, source)

        ICONS.each{ |key, val|
          cell['source'].gsub!(/{% icon #{key} %}/, val)
        }

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

        COLORS.each{ |key, val|
          if COLORS_EXTRA.has_key? key
            val = val + ";" + COLORS_EXTRA[key]
          end

          cell['source'].gsub!(/<blockquote class="#{key}">/, "<blockquote class=\"#{key}\" style=\"border: 2px solid #{val}; margin: 1em 0.2em\">")
        }

        # Images are referenced in the GTN through relative URLs which is
        # fab, but in a notebook this doesn't make sense as it will live
        # outside of the GTN. We need real URLs.
        #
        # So either we'll embed the images directly via base64 encoding (cool,
        # love it) or we'll link to the production images and folks can live
        # without their images for a bit until it's merged.

        if cell['source'].match(/<img src="\.\./)
          cell['source'].gsub!(/<img src="(\.\.[^"]*)/) { |img|
            path = img[10..-1]
            image_path = File.join(dir, path)

            if img[-3..-1].downcase == 'png'
              #puts "[GTN/Notebook/Images] Embedding png: #{img}"
              data = Base64.encode64(File.open(image_path, "rb").read)
              %Q(<img src="data:image/png;base64,#{data}")
            elsif img[-3..-1].downcase == 'jpg' or img[-4..-1].downcase == 'jpeg'
              #puts "[GTN/Notebook/Images] Embedding jpg: #{img}"
              data = Base64.encode64(File.open(image_path, "rb").read)
              %Q(<img src="data:image/jpeg;base64,#{data}")
            elsif img[-3..-1].downcase == 'svg'
              #puts "[GTN/Notebook/Images] Embedding svg: #{img}"
              data = Base64.encode64(File.open(image_path, "rb").read)
              %Q(<img src="data:image/svg+xml;base64,#{data}")
            else
              #puts "[GTN/Notebook/Images] Fallback for #{img}"
              # Falling back to non-embedded images
              '<img src="https://training.galaxyproject.org/training-material/' + page_url.split('/')[0..-2].join('/') + '/..'
            end
          }
        end

        # Strip out the highlighting as it is bad on some platforms.
        cell['source'].gsub!(/<pre class="highlight">/, '<pre style="color: inherit; background: transparent">')
        cell['source'].gsub!(/<div class="highlight">/, '<div>')
        cell['source'].gsub!(/<code>/, '<code style="color: inherit">')

        # There is some weirdness in the processing of $s in Jupyter. After a
        # certain number of them, it will give up, and just render everything
        # like with a '<pre>'. We remove this to prevent that result.
        cell['source'].gsub!(/^\s*</, '<')
        # Additionally leading spaces are sometimes interpreted as <pre>s and
        # end up causing paragraphs to be rendered as code. So we wipe out
        # all leading space.
        # 'editable' is actually CoCalc specific but oh well.
        cell['metadata'] = {"editable" => false, "collapsed" => false}
        cell['source'].gsub!(/\$/, '&#36;')
      end
      cell
    }
    notebook
  end



end
