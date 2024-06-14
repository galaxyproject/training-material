# frozen_string_literal: true

require './_plugins/gtn'

module Jekyll
  ##
  # This class generates the GTN's author pags
  class RecordingPageGenerator < Generator
    safe true

    ##
    # This generates the recording pages, where needed.
    # Params
    # +site+:: The site object
    def generate(site)
      Jekyll.logger.info "[GTN/Videos] Generating recording pages"
      TopicFilter.list_all_materials(site).select{|m| m.has_key? 'recordings'}.each do |material|
        page2 = PageWithoutAFile.new(site, '', material['dir'], 'recordings/index.html')
        page2.content = nil
        page2.data['layout'] = 'recordings'
        page2.data['topic_name'] = material['topic_name']
        page2.data['tutorial_name'] = material['tutorial_name']
        page2.data['material'] = material
        site.pages << page2
      end
    end
  end
end
