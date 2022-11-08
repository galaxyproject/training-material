module Jekyll
  class ToolTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
      format = /\[(.*)\]\((.*)\)/

      # The first part should be any text, the last should be the tool_id/tool_version
      # {% tool Group on Column %}
      # {% tool [Group on Column](Grouping1) %}
      # {% tool [Group on Column](Grouping1/1.0.0) %}
      # {% tool [Group on Column](toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.72+galaxy1) %}

      m = @text.match(format)

      if m
        # check if a variable was provided for the tool id
        if context[m[2].tr('{}','')]
          tool = context[m[2].tr('{}','')]
        else
          tool = m[2]
        end

        %Q(<span class="tool" data-tool="#{tool}" title="Tested with #{tool}"><strong>#{m[1]}</strong> <i class="fas fa-wrench" aria-hidden="true"></i><i aria-hidden="true" class="fas fa-cog"></i><span class="visually-hidden">Tool: #{tool}</span></span> )
      else
        %Q(<span><strong>#{@text}</strong> <i class="fas fa-wrench" aria-hidden="true"></i></span> )
      end

    end
  end
end

Liquid::Template.register_tag('tool', Jekyll::ToolTag)
