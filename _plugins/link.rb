# frozen_string_literal: true

module Jekyll
  # Replaces the built in link tag temporarily
  class CustomLinkTag < Liquid::Tag
    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def render(_context)
      # This is a workaround for https://github.com/jekyll/jekyll/issues/9179
      # We should remove it when 9179 is solved.
      #
      # Note that this does NOT support news posts with a date in the URL.
      "/training-material/#{@text.gsub(/\.md/, '.html')}"
    end
  end
end

Liquid::Template.register_tag('link', Jekyll::CustomLinkTag)
