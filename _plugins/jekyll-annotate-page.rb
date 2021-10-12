module Jekyll
  module AnnotateFilter
    def get_topic(page)
      # Arrays that will store all introduction slides and tutorials we discover.
      page['path'].split('/')[1]
    end
  end
end

Liquid::Template.register_filter(Jekyll::AnnotateFilter)
