require 'digest'
require './_plugins/gtn/hsluv'

module ColourTag
  def self.colour_tag(contents)
    d = (Digest::SHA256.hexdigest contents).to_i(16)

    hue = ((d >> 4) % 360).abs
    lightnessOffset = 75
    lightness = lightnessOffset + (hash & 0xf)

    # randomly make yellow tags bright
    if hue > 70 and hue < 96 and ((d & 0x100) === 0x100)
      lightness += (100 - lightness) * 0.75
    end

    primary = Hsluv.hsluv_to_hex(hue, 100, lightness)
    darker = Hsluv.hsluv_to_hex(hue, 100, lightness * 0.9)
    dimmed = Hsluv.hsluv_to_hex(hue, 100, lightness * 0.95)

    return "--color-primary: #{primary}; --color-darker: #{darker}; --color-dimmed: #{dimmed};"
  end
end



module Jekyll
  module ImplColourTag
    def cache
      @@cache ||= Jekyll::Cache.new("ColorTags")
    end

    def colour_tag(contents)
      cache.getset(contents) do
        ColourTag.colour_tag(contents)
      end
    end
  end
end

Liquid::Template.register_filter(Jekyll::ImplColourTag)
