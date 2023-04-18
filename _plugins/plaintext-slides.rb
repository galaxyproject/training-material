module Jekyll
  class PlaintextSlidesGenerator < Generator
    SLIDE_LAYOUTS = [
      'tutorial_slides',
      'base_slides',
      'introduction_slides',
      'tutorial_slides_ai4life'
    ]

    ##
    # Generate a plaintext version of the slides
    # Params:
    # +site+:: The +Jekyll::Site+ object
    def generate(site)
      # layout: tutorial_slides
      # layout: base_slides

      site.pages.select{|page| SLIDE_LAYOUTS.include? page.data['layout'] }.each do |page|
        dir = File.dirname(File.join('.', page.url))
        page2 = Jekyll::Page.new(site, site.source, dir, page.name)
        page2.data['layout'] = 'slides-plain'
        if page2.data.has_key?('lang') then
          page2.basename = "slides-plain_#{page2.data['lang'].upcase}"
        else
          page2.basename = 'slides-plain'
        end
        page2.content = page2.content.gsub(/^name:\s*([^ ]+)\s*$/) {
          anchor = $1

          "<span id=\"#{anchor.strip}\"><i class=\"fas fa-link\" aria-hidden=\"true\"></i> #{anchor}</span>"
        }
        if page2.data.has_key?('redirect_from')
          page2.data['redirect_from'].map{|x| x.gsub!(/\/slides/, '/slides-plain') }
        end

        site.pages << page2
      end
    end
  end
end
