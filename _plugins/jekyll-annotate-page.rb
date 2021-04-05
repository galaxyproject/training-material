module Jekyll
  module AnnotateFilter
    def get_topic(page)
      # Arrays that will store all introduction slides and tutorials we discover.
      page['path'].split('/')[1]
    end

    def filter_slides(pages)
      # Arrays that will store all introduction slides and tutorials we discover.
      out = []
      for page in pages do
        if ['base_slides', 'introduction_slides', 'tutorial_slides'].include?(page['layout']) then
          out.push(page)
        end
      end
      return out
    end

  end
end

Liquid::Template.register_filter(Jekyll::AnnotateFilter)
