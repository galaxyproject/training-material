module Jekyll
  class PlaintextSlidesGenerator < Generator
    safe true

    def generate(site)
      # layout: tutorial_slides
      # layout: base_slides

      site.pages.select{|page| page.data['layout'] == 'tutorial_slides' or page.data['layout'] == 'base_slides' or page.data['layout'] == 'introduction_slides'}.each do |page|
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
        site.pages << page2
      end
    end
  end
end
