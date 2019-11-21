module Jekyll
  module ColourTag
    def colour_tag(contents)
      d = (Digest::SHA256.hexdigest contents).to_i(16)

      hue = ((d >> 4) % 360).abs
      saturation = 0.8
      lightness = 85
      bgColor = "hsl(#{hue}, #{saturation * 100}%, #{lightness}%)"
      brColor = "hsl(#{hue}, #{saturation * 100}%, #{lightness - 40}%)"

      r, g, b = hsl2rgb(hue, saturation, lightness / 100.0)
      fgColor = contrasting_colour(r, g, b)

      return "background-color: #{bgColor}; color: #{fgColor}; border: 1px solid #{brColor}";
    end

    def contrasting_colour(r, g, b)
      # Implement W3C contrasting color algorithm
      # http://www.w3.org/TR/AERT#color-contrast
      # Assumes r, g, b are in the set [0, 1]
      o  = (r * 255 * 299 + g * 255 * 587 + b * 255 * 114) / 1000;
      if o > 125 then
        "#333"
      else
        "#ccc"
      end
    end

    def hue2rgb(p, q, t)
      if t < 0 then
        t = t + 1
      end

      if t > 1 then
        t = t - 1
      end

      if t < 1/6 then
        return p + (q - p) * 6 * t
      elsif t < 1/2 then
        return q
      elsif t < 2/3 then
        return p + (q - p) * (2 / 3 - t) * 6
      end

      return p
    end

    def hsl2rgb(h, s, l)
      # Converts an HSL color value to RGB. Conversion formula
      # adapted from http://en.wikipedia.org/wiki/HSL_color_space.
      # Assumes h, s, and l are contained in the set [0, 1] and
      # returns r, g, and b in the set [0, 1].
      r = 0
      g = 0
      b = 0

      if s == 0 then
        r = l
        g = l
        b = l

        return r, g, b
      end

      if l < 0.5 then
        q = l * (1 + s)
      else
        q = l + s - l * s
      end
      p = 2 * l - q

      r = hue2rgb(p, q, h + 1 / 3)
      g = hue2rgb(p, q, h)
      b = hue2rgb(p, q, h - 1 / 3)

     return r, g, b
    end

  end
end

Liquid::Template.register_filter(Jekyll::ColourTag)
