module Jekyll
  class IconTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
      cfg = get_config(context)
      icon = cfg[@text] || ''

      if icon.empty?
        raise SyntaxError.new(
          "No icon defined for: '#{@text}'. " +
          "Please define it in `_config.yml` (under `icon-tag:`)."
        )
      end

      %Q(<i class="fa #{icon}" aria-hidden="true"></i>)
    end

    private

    def get_config(context)
      context.registers[:site].config['icon-tag']
    end
  end
end

Liquid::Template.register_tag('icon', Jekyll::IconTag)
