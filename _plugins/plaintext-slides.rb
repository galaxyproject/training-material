module Jekyll
  class PlaintextSlidesGenerator < Generator
    safe true

    def generate(site)
      # layout: tutorial_slides
      # layout: base_slides

      site.pages.select{|page| page.data['layout'] == 'tutorial_slides' or page.data['layout'] == 'base_slides'}.each do |page|
        dir = File.dirname(File.join('.', page.url))
        page2 = Jekyll::Page.new(site, site.source, dir, page.name)
        page2.data['layout'] = 'slides-plain'
        page2.basename = 'slides-plain'
        site.pages << page2
      end
    end
  end
end
