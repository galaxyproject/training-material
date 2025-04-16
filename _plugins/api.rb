# frozen_string_literal: true

require 'json'

require './_plugins/jekyll-topic-filter'
require './_plugins/gtn/metrics'
require './_plugins/gtn/scholar'
require './_plugins/gtn/git'
require './_plugins/gtn/hooks'
require './_plugins/gtn/ro-crate'
require './_plugins/gtn'
require './_plugins/util'

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
  contrib_type, contrib = Gtn::Contributors.fetch(site, c)
  x = contrib
      .merge({
               'id' => c,
               'url' => site.config['url'] + site.config['baseurl'] + "/api/#{contrib_type}s/#{c}.json",
               'page' => site.config['url'] + site.config['baseurl'] + "/hall-of-fame/#{c}/",
             })
  visitAndMarkdownify(site, x)
end

module Jekyll
  module Generators
    ##
    # This class generates the GTN's "api" by writing out a folder full of JSON files.
    class APIGenerator

      def copy(site, source, dest)
        # It isn't unusual that some of these might not exist in dev envs.
        if File.exist?(site.in_source_dir(source))
          if ! Dir.exist?(File.dirname(site.in_dest_dir(dest)))
            FileUtils.mkdir_p(File.dirname(site.in_dest_dir(dest)))
          end

          FileUtils.cp(site.in_source_dir(source), site.in_dest_dir(dest))
        end
      end

      def write(site, dest, data, json: true, pretty: true)
        # Since we're doing this ourselves, need to be responsible for ensuring
        # the directory exists.
        if ! Dir.exist?(File.dirname(site.in_dest_dir(dest)))
          FileUtils.mkdir_p(File.dirname(site.in_dest_dir(dest)))
        end

        if json
          if pretty
            File.write(site.in_dest_dir(dest), JSON.pretty_generate(data))
          else
            File.write(site.in_dest_dir(dest), JSON.generate(data))
          end
        else
          # Pretty isn't relevant.
          File.write(site.in_dest_dir(dest), data)
        end
      end

      ##
      # Runs the generation process
      # Params:
      # +site+:: +Jekyll::Site+ object
      def generate(site)
        Jekyll.logger.info '[GTN/API] Generating API'

        write(site, 'api/configuration.json', site.config.reject { |k, _v| k.to_s.start_with?('cached_') })
        write(site, 'api/swagger.json', site.data['swagger'])
        write(site, 'api/version.json', Gtn::Git.discover)

        # Full Bibliography
        Jekyll.logger.debug '[GTN/API] Bibliography'
        write(site, 'api/gtn.bib', site.config['cached_global_bib'].to_s, json: false)

        # Metrics endpoint, /metrics
        write(site, 'api/metrics', Gtn::Metrics.generate_metrics(site), json: false)

        # Public tool listing
        write(site, 'api/psl.json', site.data['public-server-tools'], pretty: false)

        # Tool Categories
        write(site, 'api/toolcats.json', site.data['toolcats'], pretty: false)

        # Tool Categories
        write(site, 'api/toolshed-revisions.json', site.data['toolshed-revisions'], pretty: false)

        # Feedback Data
        write(site, 'api/feedback2.json', site.data['feedback2'], pretty: false)
        copy(site, 'metadata/feedback.csv', 'api/feedback.csv')
        copy(site, 'metadata/feedback2.yaml', 'api/feedback2.yaml')

        # Contributors
        Jekyll.logger.debug '[GTN/API] Contributors, Funders, Organisations'
        %w[contributors grants organisations].each do |type|
          pfo = site.data[type].map { |c, _| mapContributor(site, c) }
          write(site, "api/#{type}.json", pfo, pretty: false)

          site.data['contributors'].each do |c, _|
            write(site, "api/#{type}/#{c}.json", mapContributor(site, c))
          end
        end

        geojson = {
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
        write(site, "api/contributors.geojson", geojson)


        # Trigger the topic cache to generate if it hasn't already
        Jekyll.logger.debug '[GTN/API] Tutorials'
        Gtn::TopicFilter.topic_filter(site, 'does-not-matter')
        Gtn::TopicFilter.list_topics(site).map do |topic|
          out = site.data[topic].dup
          out['materials'] = Gtn::TopicFilter.topic_filter(site, topic).map do |x|
            q = x.dup
            q['contributors'] = Gtn::Contributors.get_contributors(q).dup.map do |c|
              mapContributor(site, c)
            end

            q['urls'] = {}

            if !q['hands_on'].nil?
              q['urls']['hands_on'] = site.config['url'] + site.config['baseurl'] + "/api/topics/#{q['url'][8..-6]}.json"
            end

            if !q['slides'].nil?
              q['urls']['slides'] = site.config['url'] + site.config['baseurl'] + "/api/topics/#{q['url'][8..-6]}.json"
            end

            # Write out the individual page
            # Delete the ref to avoid including it by accident
            q.delete('ref')
            q.delete('ref_tutorials')
            q.delete('ref_slides')
            write(site, "api/topics/#{q['url'][7..-6]}.json", q)

            q
          end
          out['editorial_board'] = out['editorial_board'].map do |c|
            mapContributor(site, c)
          end

          write(site, "api/topics/#{topic}.json", out)
        end

        topics = {}
        Jekyll.logger.debug '[GTN/API] Topics'
        # Individual Topic Indexes
        site.data.each_pair do |k, v|
          if v.is_a?(Hash) && v.key?('type') && v.key?('editorial_board')

            topics[k] = {
              'name' => v['name'],
              'title' => v['title'],
              'summary' => v['summary'],
              'url' => site.config['url'] + site.config['baseurl'] + "/api/topics/#{k}.json",
              'editorial_board' => v['editorial_board'].map { |c| mapContributor(site, c) }
            }
          end
        end

        # Videos.json
        # {
        # "id": "transcriptomics/tutorials/mirna-target-finder/slides",
        # "topic": "Transcriptomics",
        # "title": "Whole transcriptome analysis of Arabidopsis thaliana"
        # },

        videos = Gtn::TopicFilter.list_videos(site).map do |m|
          {
            id: "#{m['topic_name']}/tutorials/#{m['tutorial_name']}/slides",
            topic: m['topic_name_human'],
            title: m['title']
          }
        end
        write(site, "api/videos.json", videos)


        # Overall topic index
        write(site, "api/topics.json", topics)

        Jekyll.logger.debug '[GTN/API] Tutorial and Slide pages'

        # Deploy the feedback file as well
        write(site, "api/feedback.json", site.data['feedback'])

        # Top Tools
        Jekyll.logger.debug '[GTN/API] Top Tools'
        write(site, "api/top-tools.json", Gtn::TopicFilter.list_materials_by_tool(site))

        # GA4GH TRS Endpoint
        # Please note that this is all a fun hack
        Gtn::TopicFilter.list_all_materials(site).select { |m| m['workflows'] }.each do |material|
          material['workflows'].each do |workflow|
            wfid = workflow['wfid']
            wfname = workflow['wfname']

            ga4gh_blob =  {
              'id' => wfname,
              'url' => site.config['url'] + site.config['baseurl'] + material['url'],
              'name' => 'v1',
              'author' => [],
              'descriptor_type' => ['GALAXY'],
            }
            write(site, "api/ga4gh/trs/v2/tools/#{wfid}/versions/#{wfname}.json", ga4gh_blob)



            descriptor = {
              'content' => File.read("#{material['dir']}/workflows/#{workflow['workflow']}"),
              'checksum' => [],
              'url' => nil,
            }
            write(site, "api/ga4gh/trs/v2/tools/#{wfid}/versions/#{wfname}/GALAXY/descriptor.json", descriptor)
          end
        end
      end
    end
  end
end


Jekyll::Hooks.register :site, :post_read do |site|
  if Jekyll.env == 'production'
    Gtn::Hooks.by_tool(site)
  end
end

# Basically like `PageWithoutAFile`, we just write out the ones we'd created earlier.
Jekyll::Hooks.register :site, :post_write do |site|
  # No need to run this except in prod.
  if Jekyll.env == 'production'
    # Build our API
    api = Jekyll::Generators::APIGenerator.new
    api.generate(site)

    # Public tool listing: reorganised
    if site.data['public-server-tools'] && site.data['public-server-tools']['tools']
      site.data['public-server-tools']['tools'].each do |tool, version_data|
        path = File.join(site.dest, 'api', 'psl', "#{tool}.json")
        dir = File.dirname(path)
        FileUtils.mkdir_p(dir) unless File.directory?(dir)

        d = version_data.dup
        d.each_key do |k|
          # Replace the indexes with the server URLs from site['public-server-tools']['servers']
          d[k] = d[k].map { |v| site.data['public-server-tools']['servers'][v] }
        end

        File.write(path, JSON.generate(d))
      end
      Jekyll.logger.debug '[GTN/API/PSL] PSL written'
    else
      Jekyll.logger.debug '[GTN/API/PSL] PSL Dataset not available, are you in a CI environment?'
    end

    Gtn::TopicFilter.list_all_materials(site).each do |material|
      directory = material['dir']

      if material['slides']
        path = File.join(site.dest, 'api', directory, 'slides.json')
        p = material.dup
        p.delete('ref')
        p.delete('ref_tutorials')
        p.delete('ref_slides')
        p['contributors'] = Gtn::Contributors.get_contributors(p).dup.map { |c| mapContributor(site, c) }

        # Here we un-do the tutorial metadata priority, and overwrite with
        # slides metadata when available.
        slides_data = site.pages.select { |p2| p2.url == "/#{directory}/slides.html" }[0]
        p.update(slides_data.data) if slides_data&.data

        if !Dir.exist?(File.dirname(path))
          FileUtils.mkdir_p(File.dirname(path))
        end
        File.write(path, JSON.generate(p))
      end

      if material['hands_on']
        path = File.join(site.dest, 'api', directory, 'tutorial.json')
        p = material.dup
        p.delete('ref')
        p.delete('ref_tutorials')
        p.delete('ref_slides')
        p['contributors'] = Gtn::Contributors.get_contributors(p).dup.map { |c| mapContributor(site, c) }
        if !Dir.exist?(File.dirname(path))
          FileUtils.mkdir_p(File.dirname(path))
        end
        File.write(path, JSON.generate(p))
      end
    end

    # Import on-demand
    require 'securerandom'


    dir = File.join(site.dest, 'api', 'workflows')
    # ro-crate-metadata.json
    crate_start = Time.now
    count = 0
    Gtn::TopicFilter.list_all_materials(site).select { |m| m['workflows'] }.each do |material|
      material['workflows'].each do |workflow|
        Gtn::RoCrate.write(site, dir, material, workflow, site.config['url'], site.config['baseurl'])
        count += 1
      end
    end

    Jekyll.logger.debug "[GTN/API/WFRun] RO-Crate Metadata written in #{Time.now - crate_start} seconds for #{count} workflows"
  end
end
