# frozen_string_literal: true

module Jekyll
  # Convert a color into a cute little box
  class ColorPickerTag < Liquid::Tag
    def initialize(tag_name, text, tokens) # :nodoc:
      super
      @text = text.strip
    end

    ##
    # This function renders the color box
    # Params:
    # +context+:: The context of the page
    #
    # Example:
    #  {% color_picker #ff0000 %}
    def render(_context)
      "<span style='background-color:#{@text};border-radius:3px;border:1px solid #000;width:12px;" \
        "height:12px;margin-right:5px;'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span>"
    end
  end
end

Liquid::Template.register_tag('color_picker', Jekyll::ColorPickerTag)
