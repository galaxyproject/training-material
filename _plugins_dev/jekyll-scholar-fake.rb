module Jekyll
  class CiteTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
      begin
        citation_text = context.registers[:site].config['cached_citeproc'].render(:citation, id: @text)
        res = %Q(<span class="citation">#{citation_text}</span>)
      rescue
        puts "[GTN/scholar] Could not render #{@text}"
        res = %Q(<span>ERROR INVALID CITATION #{@text}</span>)
      end

      res
    end

    private
  end

  class BibTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
        %Q([bibliography])
    end

    private
  end

end

Liquid::Template.register_tag('cite', Jekyll::CiteTag)
Liquid::Template.register_tag('bibliography', Jekyll::BibTag)
