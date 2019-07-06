module Jekyll
  module NoTranslate
    # `text` is in HTML
    def no_translate(text, words)
      escaped_words = words.map { |word| Regexp.escape(word) }
      text = text.gsub(
        /(\b#{escaped_words.join("|")}\b)(?![^<]*?>)/i,
        "<span class=\"notranslate\">\\1</span>"
      )
      text
    end
  end
end

if $0 == __FILE__
  require 'test/unit'

  class NoTranslateTest < Test::Unit::TestCase
    include Jekyll::NoTranslate

    def test_no_translate_marks_a_word_to_not_translate
      text = 'hello galaxy'
      words = ['galaxy']
      expected = 'hello <span class="notranslate">galaxy</span>'

      assert_equal(expected, no_translate(text, words))
    end

    def test_no_translate_is_case_insensitive
      text = 'hello GalAxY'
      words = ['galaxy']
      expected = 'hello <span class="notranslate">GalAxY</span>'

      assert_equal(expected, no_translate(text, words))
    end

    def test_no_translate_marks_multiple_words_to_not_translate
      text = 'hello galaxy, BED, galaxy'
      words = ['galaxy', 'BED']
      expected = [
        'hello <span class="notranslate">galaxy</span>',
        '<span class="notranslate">BED</span>',
        '<span class="notranslate">galaxy</span>'
      ].join(', ')

      assert_equal(expected, no_translate(text, words))
    end

    def test_no_translate_does_not_change_urls
      text = '<a href="https://galaxyproject.org">hello</a>'
      words = ['galaxy']
      expected = text

      assert_equal(expected, no_translate(text, words))
    end

    def test_no_translate_does_not_change_title_attributes
      text = '<div title="this is from galaxy">hello</div>'
      words = ['galaxy']
      expected = text

      assert_equal(expected, no_translate(text, words))
    end

    def test_no_translate_does_not_change_parts_of_a_word
      text = 'it is galaxyproject'
      words = ['galaxy']
      expected = text

      assert_equal(expected, no_translate(text, words))
    end

    def test_no_translate_changes_content_inside_html_tags
      text = '<a>galaxy</a>'
      words = ['galaxy']
      expected = '<a><span class="notranslate">galaxy</span></a>'

      assert_equal(expected, no_translate(text, words))
    end
  end
else
  # This needs to be inside the if/else otherwise running this file with `ruby`
  # throws an error as `Liquid` is undefined.
  Liquid::Template.register_filter(Jekyll::NoTranslate)
end
