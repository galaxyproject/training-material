# frozen_string_literal: true

require 'jekyll'

module Jekyll
  ##
  # This class modifies the page contents to replace {abbr} with the associated
  # abbreviation in a proper <abbr> tag. It's done as a generator because it's
  # easier to operate once per page and be able to easily see if we've generated
  # the abbreviation before.
  #
  # This could in theory operate as a simple tag, but I'm not sure how to keep
  # track of "# of times seen per page" there.
  class Abbreviate < Jekyll::Generator
    safe true

    def initialize(config) # :nodoc:
      super
      @config = config['abbreviate'] ||= {}
    end

    def generate(site) # :nodoc:
      site.pages
          .reject { |page| skip_layout?(page.data['layout']) }
          .each { |page| abbreviate page }
      site.posts.docs
          .reject { |post| skip_layout?(post.data['layout']) }
          .each { |post| abbreviate post }
    end

    private

    def abbreviate(page) # :nodoc:
      return unless page.data.key?('abbreviations')

      seen = {}
      page.data['abbreviations'].each do |abbr, definition|
        page.content = page.content.gsub(/\{(#{abbr})\}/) do
          if seen.key?(abbr)
            firstdef = false
          else
            firstdef = true
            seen[abbr] = true
          end

          if firstdef
            "#{definition} (#{abbr})"
          else
            "<abbr title=\"#{definition}\">#{abbr}</abbr>"
          end
        end
      end
    end

    def skip_layout?(layout)
      to_skip = @config['skip_layouts'] || []

      true if to_skip.empty?

      to_skip.include?(layout)
    end
  end
end
