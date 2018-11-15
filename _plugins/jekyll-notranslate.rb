module Jekyll
  module NoTranslate
    def no_translate(text, words)
        for word in words do
            word = word.downcase
            text = text.gsub(word, "<span class=\"notranslate\">#{word}</span>")
            word = word.capitalize
            text = text.gsub(word, "<span class=\"notranslate\">#{word}</span>")
            word = word.upcase
            text = text.gsub(word, "<span class=\"notranslate\">#{word}</span>")
        end
        text
    end
  end
end

Liquid::Template.register_filter(Jekyll::NoTranslate)
