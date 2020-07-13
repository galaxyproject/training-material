module Jekyll
  class CiteTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
        %Q(<span class="text-muted">[citation hidden; run 'make serve-full' to show]</span>)
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
