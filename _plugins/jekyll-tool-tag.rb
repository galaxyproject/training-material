# frozen_string_literal: true

module Jekyll
  ##
  # Tags are useful in liquid to process data and access internal functions. Ruby functions are *always* faster than liquid templates, so, when possible, consider writing a custom Ruby filter or tag.
  #
  # 
  # - Jekyll::Tags::BibTag - {% bibliography %}
  # - Jekyll::Tags::CiteTag - {% cite hiltemann2023galaxy, %}
  # - Jekyll::Tags::CiteUrlTag - {% cite_url Batut2018 %}
  # - Jekyll::Tags::ColorPickerTag - {% color_picker #ff0000 %}
  # - Jekyll::Tags::CustomLinkTag - {% link file.md %}
  # - Jekyll::Tags::DumpSearchDataTag - {% dump_search_view testing %}
  # - Jekyll::Tags::FileExistsTag - {% file_exists path.md %}
  # - Jekyll::Tags::IconTag - {% icon email %}
  # - Jekyll::Tags::IconTagVar - {% icon var1 %}
  # - Jekyll::Tags::SnippetIncludeTag - {% snippet %}
  # - Jekyll::Tags::ToolTag - {% tool [My Tool](Grouping1) %}
  # - Jekyll::Tags::WorkflowTag - unused?
  module Tags

    # The tool tag which allows us to do fancy tool links
    class ToolTag < Liquid::Tag
      def initialize(tag_name, text, tokens)
        super
        @text = text.strip
      end

      ##
      # {% tool %} - the Tool rendering tag
      #
      # The first part should be any text, the last should be the tool_id/tool_version
      #
      # Examples:
      #
      #   {% tool Group on Column %}
      #   {% tool [Group on Column](Grouping1) %}
      #   {% tool [Group on Column](Grouping1/1.0.0) %}
      #   {% tool [Group on Column](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72+galaxy1) %}
      def render(context)
        format = /\[(.*)\]\((.*)\)/


        m = @text.match(format)

        if m
          # check if a variable was provided for the tool id
          tool = context[m[2].tr('{}', '')] || m[2]
          version = tool.split('/').last

          if tool.count('/').zero?
            "<span class=\"tool\" data-tool=\"#{tool}\" title=\"#{m[1]} tool\" aria-role=\"button\">" \
              '<i class="fas fa-wrench" aria-hidden="true"></i> ' \
              "<strong>#{m[1]}</strong>" \
              '</span>'
          else
            "<span class=\"tool\" data-tool=\"#{tool}\" title=\"#{m[1]} tool\" aria-role=\"button\">" \
              '<i class="fas fa-wrench" aria-hidden="true"></i> ' \
              "<strong>#{m[1]}</strong> " \
              '(' \
              '<i class="fas fa-cubes" aria-hidden="true"></i> ' \
              "Galaxy version #{version}" \
              ')' \
              '</span>'
          end
        else
          %(<span><strong>#{@text}</strong> <i class="fas fa-wrench" aria-hidden="true"></i></span>)
        end
      end
    end

    ##
    # The (unused?) workflow rendering tag.
    class WorkflowTag < Liquid::Tag
      def initialize(tag_name, text, tokens)
        super
        @text = text.strip
      end

      # Call the tag.
      # It MUST be this format:
      # {% workflow [Main Workflow](topics/x/tutorials/y/material/workflows/main.ga) %}
      def render(_context)
        format = /\[(?<title>.*)\]\((?<url>.*)\)/
        m = @text.match(format)
        # puts "Found #{@text} => #{m[:title]}, #{m[:url]}"

        "<span class=\"workflow\" data-workflow=\"#{m[:url]}\"><strong>#{m[:title]}</strong> " \
          '<i class="fas fa-share-alt" aria-hidden="true"></i></span>'
      end
    end
  end
end

Liquid::Template.register_tag('tool', Jekyll::Tags::ToolTag)
Liquid::Template.register_tag('workflow', Jekyll::Tags::WorkflowTag)
