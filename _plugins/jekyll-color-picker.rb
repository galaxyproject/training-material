module Jekyll
  class ColorPickerTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(context)
      "<span style='background-color:#{@text};border-radius:3px;border:1px solid #000;width:12px;height:12px;margin-right:5px;'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>"
    end

  end
end

Liquid::Template.register_tag('color_picker', Jekyll::ColorPickerTag)
