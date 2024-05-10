# frozen_string_literal: true

require 'English'
require './_plugins/gtn/contributors'
require './_plugins/gtn/boxify'
require './_plugins/gtn/mod'
require './_plugins/gtn/images'
require './_plugins/gtn/synthetic'
require './_plugins/gtn/metrics'
require './_plugins/gtn/scholar'
require './_plugins/gtn/supported'
require './_plugins/gtn/toolshed'
require './_plugins/gtn/usegalaxy'
require './_plugins/jekyll-topic-filter'
require 'time'

puts "[GTN] You are running #{RUBY_VERSION} released on #{RUBY_RELEASE_DATE} for #{RUBY_PLATFORM}"
version_parts = RUBY_VERSION.split('.')
puts '[GTN] WARNING: This Ruby is pretty old, you might want to update.' if version_parts[0].to_i < 3

##
# This module contains functions that are used in the GTN, our internal functions that is.

module Jekyll
  # The main GTN function library
  module GtnFunctions
    def self.cache
      @@cache ||= Jekyll::Cache.new('GtnFunctions')
    end

    ##
    # List of elixir node country IDs (ISO 3166-1 alpha-2) and their names
    ELIXIR_NODES = {
      'au' => 'Australia',
      'be' => 'Belgium',
      'ch' => 'Switzerland',
      'cz' => 'Czechia',
      'de' => 'Germany',
      'dk' => 'Denmark',
      'ee' => 'Estonia',
      'es' => 'Spain',
      'fi' => 'Finland',
      'fr' => 'France',
      'gr' => 'Greece',
      'hu' => 'Hungary',
      'ie' => 'Ireland',
      'il' => 'Israel',
      'it' => 'Italy',
      'lu' => 'Luxembourg',
      'nl' => 'the Netherlands',
      'no' => 'Norway',
      'pt' => 'Portugal',
      'se' => 'Sweden',
      'si' => 'Slovenia',
      'uk' => 'United Kingdom',
    }.freeze

    ##
    # Returns the name of an elixir node, given its country ID
    # Params:
    # +name+:: The country ID of the node (ISO 3166-1 alpha-2)
    # Returns:
    # +String+:: The name of the node
    def elixirnode2name(name)
      ELIXIR_NODES[name]
    end

    ##
    # Obtain the most cited paper in the GTN
    # Params:
    # +citations+:: The citations to search through
    #
    # Returns:
    # +Hash+:: The papers including their text citation and citation count
    def top_citations(citations)
      if citations.nil?
        {}
      else
        citations.sort_by { |_k, v| v }.reverse.to_h.first(20).to_h do |k, v|
          [k, { 'count' => v, 'text' => Gtn::Scholar.render_citation(k) }]
        end
      end
    end

    ##
    # A slightly more unsafe slugify function
    # Params:
    # +text+:: The text to slugify
    # Returns:
    # +String+:: The slugified text
    #
    # Example:
    #  slugify_unsafe("Hello, World!") # => "Hello-World"
    def slugify_unsafe(text)
      # Gets rid of *most* things without making it completely unusable?
      text.gsub(%r{["'\\/;:,.!@#$%^&*()]}, '').gsub(/\s/, '-').gsub(/-+/, '-')
    end

    ##
    # Return human text for ruby types
    # Params:
    # +type+:: The type to humanize
    # Returns:
    # +String+:: The humanized type
    #
    # Example:
    #  humanize_types("seq") # => "List of Items"
    def humanize_types(type)
      data = {
        'seq' => 'List of Items',
        'str' => 'Free Text',
        'map' => 'A dictionary/map',
        'float' => 'Decimal Number',
        'int' => 'Integer Number',
        'bool' => 'Boolean'
      }
      data[type]
    end

    ##
    # Replaces newlines with newline + two spaces
    def replace_newline_doublespace(text)
      text.gsub(/\n/, "\n  ")
    end

    ##
    # Returns the publication date of a page, when it was merged into `main`
    # Params:
    # +page+:: The page to get the publication date of
    # Returns:
    # +String+:: The publication date of the page
    #
    def gtn_pub_date(path)
      # if it's not a string then log a warning
      path = path['path'] if !path.is_a?(String)
      # Automatically strips any leading slashes.
      Gtn::PublicationTimes.obtain_time(path.gsub(%r{^/}, ''))
    end

    ##
    # Returns the last modified date of a page
    # Params:
    # +page+:: The page to get the last modified date of
    # Returns:
    # +String+:: The last modified date of the page
    #
    # TODO: These two could be unified tbh
    def last_modified_at(page)
      Gtn::ModificationTimes.obtain_time(page['path'])
    end

    ##
    # Returns the last modified date of a page
    # Params:
    # +page+:: The page to get the last modified date of
    # Returns:
    # +String+:: The last modified date of the page
    #
    def gtn_mod_date(path)
      # Automatically strips any leading slashes.
      Gtn::ModificationTimes.obtain_time(path.gsub(%r{^/}, ''))
    end

    ##
    # How many times has a topic been mentioned in feedback?
    # Params:
    # +feedback+:: The feedback to search through
    # +name+:: The name of the topic to search for
    # Returns:
    # +Integer+:: The number of times the topic has been mentioned
    def how_many_topic_feedbacks(feedback, name)
      if feedback.nil?
        return 0
      end

      feedback.select { |x| x['topic'] == name }.length
    end

    ##
    # How many times has a tutorial been mentioned in feedback?
    # Params:
    # +feedback+:: The feedback to search through
    # +name+:: The name of the tutorial to search for
    # Returns:
    # +Integer+:: The number of times the tutorial has been mentioned
    def how_many_tutorial_feedbacks(feedback, name)
      if feedback.nil?
        return 0
      end

      feedback.select { |x| x['tutorial'] == name }.length
    end

    ##
    # Fix the titles of boxes in a page
    # Params:
    # +content+:: The content to fix
    # +lang+:: The language of the content
    # +key+:: The key of the content
    # Returns:
    # +String+:: The fixed content
    def fix_box_titles(content, lang, key)
      Gtn::Boxify.replace_elements(content, lang, key)
    end

    ##
    # Params:
    # +data+:: The page data
    # Returns:
    # +Array+:: The "authors" of the material (list of strings)
    #
    # Example:
    #  {% assign authors = page | filter_authors -%}
    def filter_authors(data)
      Gtn::Contributors.get_authors(data)
    end

    ##
    # Params:
    # +data+:: The site data
    # +string+:: The contributor id
    # Returns:
    # +Hash+:: The contributing entity
    #
    # Example:
    #  {% assign contrib = site | fetch_contributor: page.contributor -%}
    def fetch_contributor(site, id)
      Gtn::Contributors.fetch_contributor(site, id)
    end

    ##
    # Params:
    # +data+:: The contributor's data
    # Returns:
    # +String+:: The funding URL
    #
    # Example:
    #  {{ entity | fetch_funding_url }}
    def fetch_funding_url(entity)
      Gtn::Contributors.fetch_funding_url(entity)
    end

    ##
    # Params:
    # +data+:: The contributor's data
    # Returns:
    # +String+:: The avatar's URL
    #
    # Example:
    #  {{ entity | fetch_entity_avatar: 'alice', 120 }}
    def fetch_entity_avatar_url(entity, id, width)
      return 'ERROR_NO_ENTITY' if entity.nil?

      width.nil? ? '' : "width=\"#{width}\""
      if !entity['avatar'].nil?
        entity['avatar']
      elsif entity['github'] != false
        qp = width.nil? ? '' : "?s=#{width}"
        "https://avatars.githubusercontent.com/#{id}#{qp}"
      else
        '/training-material/assets/images/avatar.png'
      end
    end

    ##
    # Params:
    # +data+:: The contributor's data
    # Returns:
    # +String+:: The funding URL
    #
    # Example:
    #  {{ entity | fetch_entity_avatar: 'alice', 120 }}
    def fetch_entity_avatar(entity, id, width)
      if entity.nil?
        return '<img src="/training-material/assets/images/avatar.png" alt="ERROR_NO_ENTITY avatar" class="avatar"/>'
      end

      w = width.nil? ? '' : "width=\"#{width}\""
      url = fetch_entity_avatar_url(entity, id, width)
      %(<img src="#{url}" alt="#{entity['name']} avatar" #{w} class="avatar" />)
    end

    ##
    # Convert a fedi address to a link
    # Params:
    # +fedi_address+:: The fedi address to convert
    # Returns:
    # +String+:: The URL at which their profile is accessible
    #
    # Example:
    #  {{ contributors[page.contributor].fediverse | fedi2link }}
    #
    #  fedi2link("@hexylena@galaxians.garden") => "https://galaxians.garden/@hexylena"
    def fedi2link(fedi_address)
      fedi_address.gsub(/^@?(?<user>.*)@(?<host>.*)$/) { |_m| "https://#{$LAST_MATCH_INFO[:host]}/@#{$LAST_MATCH_INFO[:user]}" }
    end

    ##
    # Load an SVG file directly into the page
    # Params:
    # +path+:: The path of the SVG file (relative to GTN workspace root)
    # Returns:
    # +String+:: The SVG file contents
    #
    # Example:
    #  {{ "assets/images/mastodon.svg" | load_svg }}
    def load_svg(path)
      File.read(path).gsub(/\R+/, '')
    end

    def regex_replace(str, regex_search, value_replace)
      regex = /#{regex_search}/m
      str.gsub(regex, value_replace)
    end

    ##
    # This method does a single regex replacement
    #
    # = Example
    #
    #   {{ content | regex_replace: '<hr>', '' }}
    def regex_replace_once(str, regex_search, value_replace)
      regex = /#{regex_search}/m
      str.sub(regex, value_replace)
    end

    def convert_to_material_list(site, materials)
      # [{"name"=>"introduction", "topic"=>"admin"}]
      return [] if materials.nil?

      materials.map do |m|
        if m.key?('name') && m.key?('topic')
          found = TopicFilter.fetch_tutorial_material(site, m['topic'], m['name'])
          Jekyll.logger.warn "Could not find material #{m['topic']}/#{m['name']} in the site data" if found.nil?
          found
        elsif m.key?('external') && m['external']
          {
            'layout' => 'tutorial_hands_on',
            'name' => m['name'],
            'title' => m['name'],
            'hands_on' => 'external',
            'hands_on_url' => m['link'],
          }
        else
          Jekyll.logger.warn "[GTN] Unsure how to render #{m}"
        end
      end
    end

    ##
    # Convert a workflow path to a TRS path
    # Params:
    # +str+:: The workflow path
    # Returns:
    # +String+:: The TRS path
    #
    # Example:
    #  {{ "topics/metagenomics/tutorials/mothur-miseq-sop-short/workflows/workflow1_quality_control.ga" |
    #     convert_workflow_path_to_trs }}
    #  => "/api/ga4gh/trs/v2/tools/metagenomics-mothur-miseq-sop-short/versions/workflow1_quality_control"
    def convert_workflow_path_to_trs(str)
      return 'GTN_TRS_ERROR_NIL' if str.nil?

      m = str.match(%r{topics/(?<topic>.*)/tutorials/(?<tutorial>.*)/workflows/(?<workflow>.*)\.ga})
      return "/api/ga4gh/trs/v2/tools/#{m[:topic]}-#{m[:tutorial]}/versions/#{m[:workflow].downcase}" if m

      'GTN_TRS_ERROR'
    end

    def layout_to_human(layout)
      case layout
      when /slides/
        'Slides'
      when /tutorial_hands_on/
        'Hands-on'
      when 'faq'
        'FAQs'
      when 'news'
        'News'
      end
    end

    def get_version_number(page)
      Gtn::ModificationTimes.obtain_modification_count(page['path'])
    end

    def get_rating_histogram(site, material_id, recent: false)
      return {} if material_id.nil?

      feedbacks = recent ? get_recent_feedbacks_time(site, material_id) : get_feedbacks(site, material_id)

      return {} if feedbacks.nil? || feedbacks.empty?

      ratings = feedbacks.map { |f| f['rating'] }
      ratings.each_with_object(Hash.new(0)) { |w, counts| counts[w] += 1 }
    end

    def get_rating_histogram_chart(site, material_id)
      histogram = get_rating_histogram(site, material_id)
      return {} if histogram.empty?

      highest = histogram.map { |_k, v| v }.max
      histogram
        .map { |k, v| [k, [v, v / highest.to_f]] }
        .sort_by { |k, _v| -k }
        .to_h
    end

    def get_rating(site, material_id, recent: false)
      f = get_rating_histogram(site, material_id, recent: recent)
      rating = f.map { |k, v| k * v }.sum / f.map { |_k, v| v }.sum.to_f
      rating.round(1)
    end

    def get_rating_recent(site, material_id)
      r = get_rating(site, material_id, recent: true)
      r.nan? ? get_rating(site, material_id, recent: false) : r
    end

    # Only accepts an integer rating
    def to_stars(rating)
      if rating.nil? || (rating.to_i < 1) || (rating == '0') || rating.zero?
        %(<span class="sr-only">0 stars</span>) +
          '<i class="far fa-star" aria-hidden="true"></i>'
      elsif rating.to_i < 1
        '<span class="sr-only">0 stars</span><i class="far fa-star" aria-hidden="true"></i>'
      else
        %(<span class="sr-only">#{rating} stars</span>) +
          ('<i class="fa fa-star" aria-hidden="true"></i>' * rating.to_i)
      end
    end

    def get_feedbacks(site, material_id)
      return [] if material_id.nil?

      begin
        topic, tutorial = material_id.split('/')

        if tutorial.include?(':')
          language = tutorial.split(':')[1]
          tutorial = tutorial.split(':')[0]
          # If a language is supplied, then
          feedbacks = site.data['feedback2'][topic][tutorial]
                          .select { |f| (f['lang'] || '').downcase == language.downcase }
        else
          # English is the default
          feedbacks = site.data['feedback2'][topic][tutorial]
                          .select { |f| f['lang'].nil? }
        end
      rescue StandardError
        return []
      end

      return [] if feedbacks.nil? || feedbacks.empty?

      feedbacks
        .sort_by { |f| f['date'] }
        .reverse
        .map do |f|
        f['stars'] = to_stars(f['rating'])
        f
      end
    end

    def get_feedback_count(site, material_id)
      get_feedbacks(site, material_id).length
    end

    def get_feedback_count_recent(site, material_id)
      get_recent_feedbacks_time(site, material_id).length
    end

    def get_recent_feedbacks_time(site, material_id)
      feedbacks = get_feedbacks(site, material_id)
                  .select do |f|
                    f['pro']&.length&.positive? ||
                      f['con']&.length&.positive?
                  end
                  .map do |f|
        f['f_date'] = Date.parse(f['date']).strftime('%B %Y')
        f
      end

      feedbacks.select { |f| Date.parse(f['date']) > Date.today - 365 }
    end

    def get_recent_feedbacks(site, material_id)
      feedbacks = get_feedbacks(site, material_id)
                  .select do |f|
                    f['pro']&.length&.positive? ||
                      f['con']&.length&.positive?
                  end
                  .map do |f|
        f['f_date'] = Date.parse(f['date']).strftime('%B %Y')
        f
      end

      last_year = feedbacks.select { |f| Date.parse(f['date']) > Date.today - 365 }
      # If we have fewer than 20 in the last year, then extend further.
      if last_year.length < 20
        feedbacks
          .first(20)
          .group_by { |f| f['f_date'] }
      else
        # Otherwise just everything last year.
        last_year
          .group_by { |f| f['f_date'] }
      end
    end

    def tutorials_over_time_bar_chart(site)
      graph = Hash.new(0)
      TopicFilter.list_all_materials(site).each do |material|
        yymm = material['pub_date'].strftime('%Y-%m')
        graph[yymm] += 1
      end

      # Cumulative over time
      # https://stackoverflow.com/questions/71745593/how-to-do-a-single-line-cumulative-count-for-hash-values-in-ruby
      graph
        # Turns it into an array
        .sort_by { |k, _v| k }
        # Cumulative sum
        .each_with_object([]) { |(k, v), a| a << [k, v + a.last&.last.to_i] }.to_h
        .map { |k, v| { 'x' => k, 'y' => v } }
        .to_json
    end

    def list_usegalaxy_servers(_site)
      Gtn::Usegalaxy.servers.map { |x| x.transform_keys(&:to_s) }
    end

    def list_usegalaxy_servers_shuffle(_site)
      Gtn::Usegalaxy.servers.map { |x| x.transform_keys(&:to_s) }.shuffle
    end

    def topic_name_from_page(page, site)
      if page.key? 'topic_name'
        site.data[page['topic_name']]['title']
      else
        site.data.fetch(page['url'].split('/')[2], { 'title' => '' })['title']
      end
    end

    ##
    # Get the topic of a page's path
    # Params:
    # +page+:: The page to get the topic of, it will inspect page['path']
    # Returns:
    # +String+:: The topic of the page
    #
    # Example:
    #  {{ page | get_topic }}
    def get_topic(page)
      # Arrays that will store all introduction slides and tutorials we discover.
      page['path'].split('/')[1]
    end

    def shuffle(array)
      array.shuffle
    end

    def get_og_desc(site, page); end

    def get_og_title(site, page, reverse)
      og_title = []
      topic_id = page['path'].gsub(%r{^\./}, '').split('/')[1]

      if site.data.key?(topic_id)
        if site.data[topic_id].is_a?(Hash) && site.data[topic_id].key?('title')
          og_title = [site.data[topic_id]['title']]
        else
          Jekyll.logger.warn "Missing title for #{topic_id}"
        end
      end

      if page['layout'] == 'topic'
        og_title.push 'Tutorial List'
        return og_title.join(' / ')
      end

      material_id = page['path'].gsub(%r{^\./}, '').split('/')[3]
      material = nil
      material = fetch_tutorial_material(site, topic_id, material_id) if site.data.key? topic_id

      og_title.push material['title'] if !material.nil?

      case page['layout']
      when 'workflow-list'
        og_title.push 'Workflows'
      when 'faq-page', 'faqs'
        if page['path'] =~ %r{faqs/gtn}
          og_title.push 'GTN FAQs'
        elsif page['path'] =~ %r{faqs/galaxy}
          og_title.push 'Galaxy FAQs'
        else
          og_title.push 'FAQs'
        end
      when 'faq'
        og_title.push "FAQ: #{page['title']}"
      when 'learning-pathway'
        og_title.push "Learning Pathway: #{page['title']}"
      when 'tutorial_hands_on'
        og_title[-1]&.prepend 'Hands-on: '
      when /slides/
        og_title[-1]&.prepend 'Slide Deck: '
      else
        og_title.push page['title']
      end

      if reverse.to_s == 'true'
        og_title.compact.reverse.join(' / ').gsub(/Hands-on: Hands-on:/, 'Hands-on:')
      else
        og_title.compact.join(' / ').gsub(/Hands-on: Hands-on:/, 'Hands-on:')
      end
    end

    ##
    # Gets the 'default' link for a material, hands on if it exists, otherwise slides.
    # Params:
    # +material+:: The material to get the link for
    # Returns:
    # +String+:: The URL of the default link
    def get_default_link(material)
      return 'NO LINK' if material.nil?
      return 'NO LINK' if material == true

      url = nil

      url = "topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/slides.html" if material['slides']

      if material['hands_on'] && (material['hands_on'] != 'external' && material['hands_on'] != '')
        url = "topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/tutorial.html"
      end

      url
    end

    def group_icons(icons)
      icons.group_by { |_k, v| v }.transform_values { |v| v.map { |z| z[0] } }.invert
    end
  end
end

Liquid::Template.register_filter(Jekyll::GtnFunctions)

##
# This does post-modification to every page
# Mapping the authors to their human names, and copying the cover (when present) to 'image'
#
# This exists because the jekyll-feed plugin expects those fields to look like that.
Jekyll::Hooks.register :posts, :pre_render do |post, _out|
  post.data['author'] = Gtn::Contributors.get_authors(post.data).map do |c|
    Gtn::Contributors.fetch_name(post.site, c)
  end.join(', ')
  post.data['image'] = post.data['cover']
end

# We're going to do some find and replace, to replace `@gtn:contributorName` with a link to their profile.
Jekyll::Hooks.register :site, :pre_render do |site|
  site.posts.docs.each do |post|
    if post.content
      post.content = post.content.gsub(/@gtn:([a-zA-Z0-9_-]+)/) do |match|
        # Get first capture
        name = match.gsub('@gtn:', '')
        if site.data['contributors'].key?(name)
          "{% include _includes/contributor-badge-inline.html id=\"#{name}\" %}"
        else
          match
        end
      end
    end
  end
  site.pages.each do |page|
    if page.content
      page.content = page.content.gsub(/@gtn:([a-zA-Z0-9_-]+)/) do |match|
        name = match.gsub('@gtn:', '')
        if site.data['contributors'].key?(name)
          "{% include _includes/contributor-badge-inline.html id=\"#{name}\" %}"
        else
          match
        end
      end

      # This would also need to modify the box types themselves, not sure how is best to do that.
      page.content = page.content.gsub(/> \[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]/) do |match|
        if match =~ /(CAUTION|WARNING)/
          "> <warning-title></warning-title>"
        elsif match =~ /TIP/
          "> <tip-title></tip-title>"
        else
          "> <comment-title></comment-title>"
        end
      end
    end
  end
end

# Create back-refs for affiliations
Jekyll::Hooks.register :site, :post_read do |site|
  # Users list affiliations on their profile in site.data['contributors']
  # And we want to create a back-ref to the user from the affiliation
  site.data['contributors'].each do |name, contributor|
    if contributor.key?('affiliations')
      contributor['affiliations'].each do |affiliation|
        if site.data['organisations'].key?(affiliation)
          if !site.data['organisations'][affiliation].key?('members')
            site.data['organisations'][affiliation]['members'] = []
          end

          site.data['organisations'][affiliation]['members'] << name
        elsif site.data['funders'].key?(affiliation)
          site.data['funders'][affiliation]['members'] = [] if !site.data['funders'][affiliation].key?('members')

          site.data['funders'][affiliation]['members'] << name
        end
      end
    end

    if contributor.key?('former_affiliations')
      contributor['former_affiliations'].each do |affiliation|
        if site.data['organisations'].key?(affiliation)
          if !site.data['organisations'][affiliation].key?('former_members')
            site.data['organisations'][affiliation]['former_members'] = []
          end

          site.data['organisations'][affiliation]['former_members'] << name
        elsif site.data['funders'].key?(affiliation)
          if !site.data['funders'][affiliation].key?('former_members')
            site.data['funders'][affiliation]['former_members'] = []
          end

          site.data['funders'][affiliation]['former_members'] << name
        end
      end
    end
  end
end

if $PROGRAM_NAME == __FILE__
  result = Gtn::ModificationTimes.obtain_time(ARGV[0].gsub(%r{^/}, ''))
  puts "Modification time of #{ARGV[0].gsub(%r{^/}, '')} is #{result}"
end
