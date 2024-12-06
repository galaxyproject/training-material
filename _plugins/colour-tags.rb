# frozen_string_literal: true

require 'digest'
require './_plugins/gtn/hsluv'

# Our automatic colour tag generator
# It makes tags colourful in a reproducible way
module ColourTag
  ##
  # This function generates the CSS for a colour tag
  # Params
  # +contents+:: The contents of the tag
  #
  # Returns
  # +String+:: The CSS for the tag
  #
  # Example
  #  ColourTag.colour_tag("test") => "--color-primary: #f799ff; --color-darker: #f571ff; --color-dimmed: #f686ff;"
  def self.colour_tag(contents)
    d = (Digest::SHA256.hexdigest contents).to_i(16)

    hue = ((d >> 4) % 360).abs
    lightnessOffset = 75
    lightness = lightnessOffset + (hash & 0xf)

    # randomly make yellow tags bright
    lightness += (100 - lightness) * 0.75 if (hue > 70) && (hue < 96) && ((d & 0x100) == 0x100)

    primary = Hsluv.hsluv_to_hex(hue, 100, lightness)
    darker = Hsluv.hsluv_to_hex(hue, 100, lightness * 0.9)
    dimmed = Hsluv.hsluv_to_hex(hue, 100, lightness * 0.95)

    "--color-primary: #{primary}; --color-darker: #{darker}; --color-dimmed: #{dimmed};"
  end
end

module Jekyll # :nodoc:
  # The jekyll implementation of the colour tag
  module ImplColourTag
    def cache
      @@cache ||= Jekyll::Cache.new('ColorTags')
    end

    def colour_tag(contents)
      cache.getset(contents) do
        ColourTag.colour_tag(contents)
      end
    end
  end
end

Liquid::Template.register_filter(Jekyll::ImplColourTag)
