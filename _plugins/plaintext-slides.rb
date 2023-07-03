# frozen_string_literal: true

module Jekyll
  # Convert our slides to plaintext
  # It's not a great convesion, the CSS classes are retained which are ugly
  # But there's no good way to parse those out since they use a wildly nonstandard syntax
  class PlaintextSlidesGenerator < Generator
    SLIDE_LAYOUTS = %w[
      tutorial_slides
      base_slides
      introduction_slides
      tutorial_slides_ai4life
    ].freeze

    ##
    # Generate a plaintext version of the slides
    # Params:
    # +site+:: The +Jekyll::Site+ object
    def generate(site)
      # layout: tutorial_slides
      # layout: base_slides

      site.pages.select { |page| SLIDE_LAYOUTS.include? page.data['layout'] }.each do |page|
        dir = File.dirname(File.join('.', page.url))
        page2 = Jekyll::Page.new(site, site.source, dir, page.name)
        page2.data['layout'] = 'slides-plain'
        page2.basename = if page2.data.key?('lang')
                           "slides-plain_#{page2.data['lang'].upcase}"
                         else
                           'slides-plain'
                         end
        page2.content = page2.content.gsub(/^name:\s*([^ ]+)\s*$/) do
          anchor = ::Regexp.last_match(1)

          "<span id=\"#{anchor.strip}\"><i class=\"fas fa-link\" aria-hidden=\"true\"></i> #{anchor}</span>"
        end
        if page2.data.key?('redirect_from')
          page2.data['redirect_from'].map { |x| x.gsub!(%r{/slides}, '/slides-plain') }
        end

        site.pages << page2
      end
    end
  end
end
