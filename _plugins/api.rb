# frozen_string_literal: true

require 'json'
require './_plugins/jekyll-topic-filter'
require './_plugins/gtn/metrics'
require './_plugins/gtn/scholar'
require './_plugins/gtn/git'
require './_plugins/gtn'

module Jekyll
  ##
  # This class generates the GTN's "api" by writing out a folder full of JSON files.
  class APIGenerator < Generator
    ##
    # Returns contributors, regardless of whether they are 'contributor' or 'contributions' style
    # Params:
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of contributor IDs
    def get_contributors(data)
      if data.key?('contributors')
        data['contributors']
      elsif data.key?('contributions')
        data['contributions'].keys.map { |k| data['contributions'][k] }.flatten
      end
    end

    ##
    # Use Jekyll's Markdown converter to convert text to HTML
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +text+:: +String+ of text to convert
    # Returns:
    # +String+ of markdown text
    def markdownify(site, text)
      site.find_converter_instance(
        Jekyll::Converters::Markdown
      ).convert(text.to_s)
    end

    ##
    # Recursively visit a hash and markdownify all strings inside
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +f+:: +Hash+ to visit
    # Returns:
    # +Hash+ with all strings markdownified
    def visitAndMarkdownify(site, f)
      case f
      when Array
        f.map! { |x| visitAndMarkdownify(site, x) }
      when Hash
        f = f.transform_values do |v|
          visitAndMarkdownify(site, v)
        end
      when String
        f = markdownify(site, f).strip.gsub(/<p>/, '').gsub(%r{</p>}, '')
      end
      f
    end

    ##
    # Map a contributor ID to a JSON object which includes links to their profile page and API endpoint
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +c+:: +String+ of contributor ID
    # Returns:
    # +Hash+ of contributor information
    def mapContributor(site, c)
      x = site.data['contributors']
              .fetch(c, {})
              .merge({
                       'id' => c,
                       'url' => site.config['url'] + site.config['baseurl'] + "/api/contributors/#{c}.json",
                       'page' => site.config['url'] + site.config['baseurl'] + "/hall-of-fame/#{c}/",
                     })
      visitAndMarkdownify(site, x)
    end

    ##
    # Generates /api/configuration.json
    # Params:
    # +site+:: +Jekyll::Site+ object
    # Returns:
    # nil
    def generateConfiguration(site)
      page2 = PageWithoutAFile.new(site, '', 'api/', 'configuration.json')
      site.config.update(Gtn::Git.discover)
      # Remove every key that starts with "cached_"
      conf = site.config.reject { |k, _v| k.to_s.start_with?('cached_') }
      page2.content = JSON.pretty_generate(conf)
      page2.data['layout'] = nil
      site.pages << page2
    end

    ##
    # Generates /api/version.json
    # Params:
    # +site+:: +Jekyll::Site+ object
    # Returns:
    # nil
    def generateVersion(site)
      page2 = PageWithoutAFile.new(site, '', 'api/', 'version.json')
      page2.content = JSON.pretty_generate(Gtn::Git.discover)
      page2.data['layout'] = nil
      site.pages << page2
    end

    ##
    # Generates /api/data-library.yaml
    # Params:
    # +site+:: +Jekyll::Site+ object
    # Returns:
    # nil
    def generateLibrary(site)
      puts '[GTN/API] Data Library'
      page2 = PageWithoutAFile.new(site, '', 'api/', 'data-library.yaml')
      data_libraries = Dir.glob('topics/**/data-library.yaml')
      data_libraries.map! { |x| YAML.load_file(x) }
      pp data_libraries
      page2.content = JSON.pretty_generate(Gtn::Git.discover)
      page2.data['layout'] = nil
      site.pages << page2
    end

    ##
    # Runs the generation process
    # Params:
    # +site+:: +Jekyll::Site+ object
    def generate(site)
      generateConfiguration(site)
      # For some reason the templating isn't working right here.
      # generateVersion(site)
      # TODO:
      # generateLibrary(site)

      # Full Bibliography
      Gtn::Scholar.load_bib(site)
      puts '[GTN/API] Bibliography'
      page3 = PageWithoutAFile.new(site, '', 'api/', 'gtn.bib')
      page3.content = site.config['cached_global_bib'].to_s
      page3.data['layout'] = nil
      site.pages << page3

      # Metrics endpoint, /metrics
      page2 = PageWithoutAFile.new(site, '', '', 'metrics')
      page2.content = "{% raw %}\n#{Gtn::Metrics.generate_metrics(site)}{% endraw %}"
      page2.data['layout'] = nil
      site.pages << page2

      # Contributors
      puts '[GTN/API] Contributors'
      page2 = PageWithoutAFile.new(site, '', 'api/', 'contributors.json')
      page2.content = JSON.pretty_generate(site.data['contributors'].map { |c, _| mapContributor(site, c) })
      page2.data['layout'] = nil
      site.pages << page2
      site.data['contributors'].each do |c, _|
        page4 = PageWithoutAFile.new(site, '', 'api/', "contributors/#{c}.json")
        page4.content = JSON.pretty_generate(mapContributor(site, c))
        page4.data['layout'] = nil
        site.pages << page4
      end

      page2 = PageWithoutAFile.new(site, '', 'api/', 'contributors.geojson')
      page2.content = JSON.pretty_generate(
        {
          'type' => 'FeatureCollection',
          'features' => site.data['contributors']
            .select { |_k, v| v.key? 'location' }
            .map do |k, v|
              {
                'type' => 'Feature',
                'geometry' => { 'type' => 'Point', 'coordinates' => [v['location']['lon'], v['location']['lat']] },
                'properties' => {
                  'name' => v.fetch('name', k),
                  'url' => "https://training.galaxyproject.org/training-material/hall-of-fame/#{k}/",
                  'joined' => v['joined'],
                  'orcid' => v['orcid'],
                  'id' => k,
                  'contact_for_training' => v.fetch('contact_for_training', false),
                }
              }
            end
        }
      )
      page2.data['layout'] = nil
      site.pages << page2

      # Trigger the topic cache to generate if it hasn't already
      puts '[GTN/API] Tutorials'
      TopicFilter.topic_filter(site, 'does-not-matter')
      TopicFilter.list_topics(site).map do |topic|
        out = site.data[topic].dup
        out['materials'] = TopicFilter.topic_filter(site, topic).map do |x|
          q = x.dup
          q['contributors'] = get_contributors(q).dup.map { |c| mapContributor(site, c) }

          q['urls'] = {}

          if !q['hands_on'].nil?
            q['urls']['hands_on'] = site.config['url'] + site.config['baseurl'] + "/api/topics/#{q['url'][8..-6]}.json"
          end

          if !q['slides'].nil?
            q['urls']['slides'] = site.config['url'] + site.config['baseurl'] + "/api/topics/#{q['url'][8..-6]}.json"
          end

          # Write out the individual page
          page6 = PageWithoutAFile.new(site, '', 'api/topics/', "#{q['url'][7..-6]}.json")
          # Delete the ref to avoid including it by accident
          q.delete('ref')
          page6.content = JSON.pretty_generate(q)
          page6.data['layout'] = nil
          site.pages << page6

          q
        end
        out['maintainers'] = out['maintainers'].map { |c| mapContributor(site, c) }

        page2 = PageWithoutAFile.new(site, '', 'api/topics/', "#{topic}.json")
        page2.content = JSON.pretty_generate(out)
        page2.data['layout'] = nil
        site.pages << page2
      end

      topics = {}
      puts '[GTN/API] Topics'
      # Individual Topic Indexes
      site.data.each_pair do |k, v|
        if v.is_a?(Hash) && v.key?('type') && v.key?('maintainers')

          topics[k] = {
            'name' => v['name'],
            'title' => v['title'],
            'summary' => v['summary'],
            'url' => site.config['url'] + site.config['baseurl'] + "/api/topics/#{k}.json",
            'maintainers' => v['maintainers'].map { |c| mapContributor(site, c) }
          }
        end
      end

      # Videos.json
      # {
      # "id": "transcriptomics/tutorials/mirna-target-finder/slides",
      # "topic": "Transcriptomics",
      # "title": "Whole transcriptome analysis of Arabidopsis thaliana"
      # },

      page2 = PageWithoutAFile.new(site, '', 'api/', 'videos.json')
      page2.content = JSON.pretty_generate(TopicFilter.list_videos(site).map do |m|
        {
          id: "#{m['topic_name']}/tutorials/#{m['tutorial_name']}/slides",
          topic: m['topic_name_human'],
          title: m['title']
        }
      end)
      page2.data['layout'] = nil
      site.pages << page2

      # Overall topic index
      page2 = PageWithoutAFile.new(site, '', 'api/', 'topics.json')
      page2.content = JSON.pretty_generate(topics)
      page2.data['layout'] = nil
      site.pages << page2

      puts '[GTN/API] Tutorial and Slide pages'

      TopicFilter.list_all_materials(site).each do |material|
        directory = material['dir']

        if material['slides']
          page5 = PageWithoutAFile.new(site, '', 'api/', "#{directory}/slides.json")
          p = material.dup
          p.delete('ref')
          p['contributors'] = get_contributors(p).dup.map { |c| mapContributor(site, c) }

          # Here we un-do the tutorial metadata priority, and overwrite with
          # slides metadata when available.
          slides_data = site.pages.select { |p2| p2.url == "/#{directory}/slides.html" }[0]
          p.update(slides_data.data) if slides_data&.data

          page5.content = JSON.pretty_generate(p)
          page5.data['layout'] = nil
          site.pages << page5
        end

        if material['hands_on']
          page5 = PageWithoutAFile.new(site, '', 'api/', "#{directory}/tutorial.json")
          p = material.dup
          p.delete('ref')
          p['contributors'] = get_contributors(p).dup.map { |c| mapContributor(site, c) }
          page5.content = JSON.pretty_generate(p)
          page5.data['layout'] = nil
          site.pages << page5
        end
      end

      # Deploy the feedback file as well
      page2 = PageWithoutAFile.new(site, '', 'api/', 'feedback.json')
      page2.content = JSON.pretty_generate(site.data['feedback'])
      page2.data['layout'] = nil
      site.pages << page2

      # Top Tools
      puts '[GTN/API] Top Tools'
      page2 = PageWithoutAFile.new(site, '', 'api/', 'top-tools.json')
      page2.content = JSON.pretty_generate(TopicFilter.list_materials_by_tool(site))
      page2.data['layout'] = nil
      site.pages << page2

      # Not really an API
      TopicFilter.list_materials_by_tool(site).each do |tool, tutorials|
        page2 = PageWithoutAFile.new(site, '', 'by-tool/', "#{tool.gsub('%20', ' ')}.html")
        page2.content = nil
        page2.data['layout'] = 'by_tool'
        page2.data['short_tool'] = tool
        page2.data['observed_tool_ids'] = tutorials['tool_id']
        page2.data['tutorial_list'] = tutorials['tutorials']
        site.pages << page2
      end

      # GA4GH TRS Endpoint
      # Please note that this is all a fun hack
      TopicFilter.list_all_materials(site).select { |m| m['workflows'] }.each do |material|
        material['workflows'].each do |workflow|
          wfid = workflow['wfid']
          wfname = workflow['wfname']

          page2 = PageWithoutAFile.new(site, '', "api/ga4gh/trs/v2/tools/#{wfid}/versions/", "#{wfname}.json")
          page2.content = JSON.pretty_generate(
            {
              'id' => wfname,
              'url' => site.config['url'] + site.config['baseurl'] + material['url'],
              'name' => 'v1',
              'author' => [],
              'descriptor_type' => ['GALAXY'],
            }
          )
          page2.data['layout'] = nil
          site.pages << page2

          page2 = PageWithoutAFile.new(site, '', "api/ga4gh/trs/v2/tools/#{wfid}/versions/#{wfname}/GALAXY",
                                       'descriptor.json')
          page2.content = JSON.pretty_generate(
            {
              'content' => File.read("#{material['dir']}/workflows/#{workflow['workflow']}"),
              'checksum' => [],
              'url' => nil,
            }
          )
          page2.data['layout'] = nil
          site.pages << page2
        end
      end
    end
  end
end
