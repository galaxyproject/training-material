# frozen_string_literal: true

require './_plugins/gtn'

module Jekyll
  ##
  # {Jekyll Generators}[https://jekyllrb.com/docs/plugins/generators/] are a way to let you generate files at runtime, without needing them to exist on disk.
  #
  # We use generators for lots of purposes, e.g.
  #
  # Real Generators, the way Jekyll intended:
  #
  # - Jekyll::Generators::PlaintextSlidesGenerator - turns slides.html into plain text non-JS versions.
  # - Jekyll::Generators::RecordingPageGenerator - emits a webpage for every tutorial that has recordings, in the GTN
  # - Jekyll::Generators::WorkflowPageGenerator - emits a webpage for every workflow in the GTN
  # - Jekyll::Generators::AuthorPageGenerator - emits a hall-of-fame entry for every contributor, organisation, and grant listed in our site metadata.
  # - Jekyll::Generators::RmarkdownGenerator - outputs the RMarkdown notebooks for tutorials that want them.
  # - Jekyll::Generators::SitemapGenerator2 - alternative for the jekyll-sitemap plugin that's a bit faster.
  # - Jekyll::Generators::SyntheticTopicGenerator - our synthetic tag based topics
  # - Jekyll::Generators::TagPageGenerator - topic page for every tag
  # - Jekyll::Generators::WorkflowPageGenerator
  #
  # Muck with page contents generators (probably should be hooks):
  #
  # - Jekyll::Generators::Abbreviate - turns +{something}+ into an abbreviation tag
  # - Jekyll::Generators::Figurify - handles our modified markdown for images
  # - Jekyll::Generators::Boxify - turns +<box-title>blah</box-title> ... {: .box}+ into GTN boxes.
  #
  # Other generators (also probably should be hooks):
  #
  # - Jekyll::Generators::EnvironmentVariablesGenerator - adds git revision, tags, and other environment variables to the site object
  module Generators

    ##
    # This class generates the GTN's workflow pages.
    class WorkflowPageGenerator < Generator
      safe true

      ##
      # Params
      # +site+:: The site object
      def generate(site)
        Jekyll.logger.info "[GTN/Workflows] Generating workflow pages"
        materials = Gtn::TopicFilter
          .list_all_materials(site)

        # [{"workflow"=>"Calling_variants_in_non-diploid_systems.ga",
        #   "tests"=>true,
        #   "url"=>"https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/non-dip/workflows/Calling_variants_in_non-diploid_systems.ga",
        #   "path"=>"topics/variant-analysis/tutorials/non-dip/workflows/Calling_variants_in_non-diploid_systems.ga",
        #   "wfid"=>"variant-analysis-non-dip",
        #   "wfname"=>"calling-variants-in-non-diploid-systems",
        #   "trs_endpoint"=>"https://training.galaxyproject.org/training-material/api/ga4gh/trs/v2/tools/variant-analysis-non-dip/versions/calling-variants-in-non-diploid-systems",
        #   "license"=>nil,
        #   "parent_id"=>"variant-analysis/non-dip",
        #   "topic_id"=>"variant-analysis",
        #   "tutorial_id"=>"non-dip",
        #   "creators"=>[],
        #   "name"=>"Calling variants in non-diploid systems",
        #   "title"=>"Calling variants in non-diploid systems",
        #   "test_results"=>nil,
        #   "modified"=>2024-03-18 12:38:46.293831071 +0100,
        #   "mermaid"=>
        #   "graph_dot"=>
        # }]
        # /api/workflows/#{topic_id}/#{tutorial_id}/#{wfid}/rocrate.zip
        shortlinks = site.data['shortlinks']['id'].invert

        materials.each do |material|
          (material['workflows'] || []).each do |workflow|
            page2 = PageWithoutAFile.new(site, '', '', "#{workflow['path'].gsub(/.ga$/, '.html')}")
            path = File.join('/', workflow['path'].gsub(/.ga$/, '.html'))
            page2.content = nil
            page2.data['title'] = workflow['title']
            page2.data['layout'] = 'workflow'
            page2.data['material'] = material
            page2.data['workflow'] = workflow
            page2.data['js_requirements'] = {'mathjax' => false, 'mermaid' => true}
            page2.data['short_id'] = shortlinks[path]
            page2.data['redirect_from'] = ["/short/#{shortlinks[path]}"]
            site.pages << page2
          end
        end
      end
    end
  end
end

