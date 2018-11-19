module Jekyll
  module NoTranslate
    def no_translate(text, words)
        escaped_words = words.map { |word| Regexp.escape(word) }
        text = text.gsub(/(#{escaped_words.join("|")})/i, "<span class=\"notranslate\">\\1</span>")
        text
    end
  end
end

Liquid::Template.register_filter(Jekyll::NoTranslate)
