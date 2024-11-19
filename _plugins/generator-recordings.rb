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
      materials = TopicFilter
        .list_all_materials(site)

      with_video = materials
        .select{|m| m.has_key? 'recordings' or m.has_key? 'slide_recordings'}

      Jekyll.logger.info "[GTN/Videos] #{with_video.length} materials with recordings found."
      materials.each do |material|
        page2 = PageWithoutAFile.new(site, '', material['dir'], 'recordings/index.html')
        page2.content = nil
        page2.data['layout'] = 'recordings'
        page2.data['topic_name'] = material['topic_name']
        page2.data['tutorial_name'] = material['tutorial_name']
        page2.data['material'] = material
        page2.data['title'] = 'Recordings for ' + material['title']
        site.pages << page2
      end
    end
  end
end
