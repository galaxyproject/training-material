# frozen_string_literal: true

require 'json'
require 'yaml'
require './_plugins/gtn'

# The main GTN module to parse tutorials and topics into useful lists of things that can bes shown on topic pages
module TopicFilter
  ##
  # This function returns a list of all the topics that are available.
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # Returns:
  # +Array+:: The list of topics
  def self.list_topics(site)
    list_topics_h(site).keys
  end

  def self.list_topics_h(site)
    site.data.select { |_k, v| v.is_a?(Hash) && v.key?('editorial_board') }
  end

  ##
  # This function returns a list of all the topics that are available.
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # Returns:
  # +Array+:: The topic objects themselves
  def self.enumerate_topics(site)
    list_topics_h(site).values
  end

  ##
  # Fill the cache with all the topics
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # Returns:
  # +nil+
  def self.fill_cache(site)
    return if site.data.key?('cache_topic_filter')

    puts '[GTN/TopicFilter] Begin Cache Prefill'
    site.data['cache_topic_filter'] = {}

    # For each topic
    list_topics(site).each do |topic|
      site.data['cache_topic_filter'][topic] = filter_by_topic(site, topic)
    end
    puts '[GTN/TopicFilter] End Cache Prefill'
  end

  ##
  # This function returns a list of all the materials that are available for a specific topic.
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # +topic_name+:: The name of the topic
  # Returns:
  # +Array+:: The list of materials
  def self.topic_filter(site, topic_name)
    fill_cache(site)
    site.data['cache_topic_filter'][topic_name]
  end

  ##
  # This function returns a list of all the materials that are available for a
  # specific topic, but this time in a structured manner
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # +topic_name+:: The name of the topic
  # Returns:
  # +Hash+:: The subtopics and their materials
  #
  # Example:
  #  {
  #   "intro" => {
  #     "subtopic" => {"title" => "Introduction", "description" => "Introduction to the topic", "id" => "intro"},
  #     "materials" => [
  #       ...
  #     ]
  #   },
  #   "__OTHER__" => {
  #     "subtopic" => {"title" => "Other", "description" => "Other materials", "id" => "__OTHER__"},
  #     "materials" => [.. ]
  #   }
  #  ]
  def self.list_materials_structured(site, topic_name)
    # This method is built with the idea to replace the "topic_filter" command,
    # and instead of returning semi-structured data, we will immediately return
    # fully structured data for a specific "topic_name" query, like, "admin"
    #
    # Instead of returning a flat list of tutorials, instead we'll structure
    # them properly in subtopics (if they exist) or return the flat list
    # otherwise.
    #
    # This will let us generate new "views" into the tutorial lists, having
    # them arranged in new and exciting ways.

    fill_cache(site)

    # Here we want to either return data structured around subtopics

    if site.data[topic_name]['tag_based'].nil? && site.data[topic_name].key?('subtopics')
      # We'll construct a new hash of subtopic => tutorials
      out = {}
      seen_ids = []
      site.data[topic_name]['subtopics'].each do |subtopic, _v|
        specific_resources = filter_by_topic_subtopic(site, topic_name, subtopic['id'])
        out[subtopic['id']] = {
          'subtopic' => subtopic,
          'materials' => specific_resources
        }
        seen_ids += specific_resources.map { |x| x['id'] }
      end

      # And we'll have this __OTHER__ subtopic for any tutorials that weren't
      # in a subtopic.
      all_topics_for_tutorial = filter_by_topic(site, topic_name)
      out['__OTHER__'] = {
        'subtopic' => { 'title' => 'Other', 'description' => 'Assorted Tutorials', 'id' => 'other' },
        'materials' => all_topics_for_tutorial.reject { |x| seen_ids.include?(x['id']) }
      }
    elsif site.data[topic_name]['tag_based'] && site.data[topic_name]['custom_ordering']
      # TODO
      Jekyll.logger.error 'UNIMPLEMENTED'
      out = {}
    elsif site.data[topic_name]['tag_based'] # Tag based Topic
      # We'll construct a new hash of subtopic(parent topic) => tutorials
      out = {}
      seen_ids = []
      tn = topic_name.gsub('by_tag_', '')

      materials = filter_by_tag(site, tn)

      # Which topics are represented in those materials?
      seen_topics = materials.map { |x| x['topic_name'] }.sort

      # Treat them like subtopics, but fake subtopics.
      seen_topics.each do |parent_topic, _v|
        specific_resources = materials.select { |x| x['topic_name'] == parent_topic }
        out[parent_topic] = {
          'subtopic' => { 'id' => parent_topic, 'title' => site.data[parent_topic]['title'], 'description' => nil },
          'materials' => specific_resources
        }
        seen_ids += specific_resources.map { |x| x['id'] }
      end

      # And we'll have this __OTHER__ subtopic for any tutorials that weren't
      # in a subtopic.
      all_topics_for_tutorial = filter_by_tag(site, tn)
      out['__OTHER__'] = {
        'subtopic' => { 'title' => 'Other', 'description' => 'Assorted Tutorials', 'id' => 'other' },
        'materials' => all_topics_for_tutorial.reject { |x| seen_ids.include?(x['id']) }
      }
    else
      # Or just the list (Jury is still out on this one, should it really be a
      # flat list? Or in this identical structure.)
      out = {
        '__FLAT__' => {
          'subtopic' => nil,
          'materials' => filter_by_topic(site, topic_name)
        }
      }
    end

    # Cleanup empty sections
    out.delete('__OTHER__') if out.key?('__OTHER__') && out['__OTHER__']['materials'].empty?

    out.each do |_k, v|
      v['materials'].sort_by! { |m| [m.fetch('priority', 1), m['title']] }
    end

    out
  end

  ##
  # Fetch a specific tutorial material by topic and tutorial name
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # +topic_name+:: The name of the topic
  # +tutorial_name+:: The name of the tutorial
  # Returns:
  # +Hash+:: The tutorial material
  def self.fetch_tutorial_material(site, topic_name, tutorial_name)
    fill_cache(site)
    if site.data['cache_topic_filter'][topic_name].nil?
      Jekyll.logger.warn "Topic cache not filled, cannot fetch tutorial material for #{topic_name}"
    else
      site.data['cache_topic_filter'][topic_name].select { |p| p['tutorial_name'] == tutorial_name }[0]
    end
  end

  ##
  # Extract the list of tools used in a workflow
  # Params:
  # +data+:: The workflow data
  # Returns:
  # +Array+:: The list of tool IDs
  def self.extract_workflow_tool_list(data)
    out = data['steps'].select { |_k, v| v['type'] == 'tool' }.map { |_k, v| v['tool_id'] }.compact
    out += data['steps'].select do |_k, v|
             v['type'] == 'subworkflow'
           end.map { |_k, v| extract_workflow_tool_list(v['subworkflow']) }
    out
  end

  ##
  # Annotation of a path with topic and tutorial information
  # Params:
  # +path+:: The path to annotate
  # Returns:
  # +Hash+:: The annotation
  #
  # This is a bit of a hack, but it works for now.
  #
  # Example:
  #  /topics/assembly/tutorials/velvet-assembly/tutorial.md
  #  => {
  #    "topic" => "assembly",
  #    "topic_name" => "assembly",
  #    "material" => "assembly/velvet-assembly",
  #    "tutorial_name" => "velvet-assembly",
  #    "dir" => "topics/assembly/tutorials/velvet-assembly"
  #    "type" => "tutorial"
  #  }
  def self.annotate_path(path, layout)
    parts = path.split('/')
    parts.shift if parts[0] == '.'

    return nil if parts[0] != 'topics'

    return nil if parts[2] != 'tutorials'

    return nil if parts.length < 4

    material = {
      'topic' => parts[1], # Duplicate
      'topic_name' => parts[1],
      'material' => "#{parts[1]}/#{parts[3]}",
      'tutorial_name' => parts[3],
      'dir' => parts[0..3].join('/'),
    }

    return nil if path =~ %r{/faqs/}

    return nil if parts[-1] =~ /data[_-]library.yaml/ || parts[-1] =~ /data[_-]manager.yaml/

    if parts[4] =~ /tutorial.*\.md/ || layout == 'tutorial_hands_on'
      material['type'] = 'tutorial'
    elsif parts[4] =~ /slides.*\.html/ || %w[tutorial_slides base_slides introduction_slides].include?(layout)
      material['type'] = 'slides'
    elsif parts[4] =~ /ipynb$/
      material['type'] = 'ipynb'
    elsif parts[4] =~ /Rmd$/
      material['type'] = 'rmd'
    elsif parts[4] == 'workflows'
      material['type'] = 'workflow'
    elsif parts[4] == 'tours'
      material['type'] = 'tour'
    elsif parts[-1] == 'index.md'
      return nil
    else
      return nil
      # material['type'] = 'unknown'
    end

    material
  end

  ##
  # Get the list of posts from the site
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # Returns:
  # +Array+:: The list of posts
  #
  # This is a transition period function that can later be removed. It is added
  # because with the jekyll version we're using, site.posts is an iterable in
  # prod+dev (_config-dev.yml) modes, however! If we access site.posts.docs in
  # prod it's fine, while in dev mode, site.posts claims to be an Array (rather
  # than I guess a 'posts' object with a docs method). So we check if it has
  # docs and use that, otherwise just site.posts should be iterable.
  def self.get_posts(site)
    # Handle the transition period
    if site.posts.respond_to?(:docs)
      site.posts.docs
    else
      site.posts
    end
  end

  ##
  # Collate the materials into a large hash
  # Params:
  # +site+:: The +Jekyll::Site+ object
  # +pages+:: The list of pages to collate
  # Returns:
  # +Hash+:: The collated materials
  #
  # Example:
  # collate_materials(site, pages)
  # => {
  # "assembly/velvet-assembly" => {
  #  "topic" => "assembly",
  #  "topic_name" => "assembly",
  #  "material" => "assembly/velvet-assembly",
  #  "tutorial_name" => "velvet-assembly",
  #  "dir" => "topics/assembly/tutorials/velvet-assembly",
  #  "resources" => [
  #    {
  #    "type" => "slides",
  #    "url" => "/topics/assembly/tutorials/velvet-assembly/slides.html",
  #    "title" => "Slides",
  #    "priority" => 1
  #    },
  #    {
  #    "type" => "tutorial",
  #    "url" => "/topics/assembly/tutorials/velvet-assembly/tutorial.html",
  #    "title" => "Tutorial",
  #    "priority" => 2
  #    }
  #   ]
  #  }
  def self.collate_materials(site, pages)
    # In order to speed up queries later, we'll store a set of "interesting"
    # pages (i.e. things that are under `topic_name`)
    shortlinks = site.data['shortlinks']
    shortlinks_reversed = shortlinks['id'].invert

    get_posts(site).each do |post|
      post.data['short_id'] = shortlinks_reversed[post.url]
    end

    interesting = {}
    pages.each do |page|
      page.data['short_id'] = shortlinks_reversed[page.url]

      # Skip anything outside of topics.
      next if !page.url.include?('/topics/')

      # Extract the material metadata based on the path
      page.data['url'] = page.url
      material_meta = annotate_path(page.path, page.data['layout'])

      # If unannotated then we want to skip this material.
      next if material_meta.nil?

      mk = material_meta['material']

      if !interesting.key? mk
        interesting[mk] = material_meta.dup
        interesting[mk].delete('type') # Remove the type since it's specific, not generic
        interesting[mk]['resources'] = []
      end

      page.data['topic_name'] = material_meta['topic_name']
      page.data['tutorial_name'] = material_meta['tutorial_name']
      page.data['dir'] = material_meta['dir']
      page.data['short_id'] = shortlinks_reversed[page.data['url']]

      interesting[mk]['resources'].push([material_meta['type'], page])
    end

    interesting
  end

  def self.mermaid(wf)
    # We're converting it to Mermaid.js
    # flowchart TD
    #     A[Start] --> B{Is it?}
    #     B -- Yes --> C[OK]
    #     C --> D[Rethink]
    #     D --> B
    #     B -- No ----> E[End]

    output = "flowchart TD\n"
    wf['steps'].each_key do |id|
      step = wf['steps'][id]
      output += "  #{id}[\"#{step['name']}\"];\n"
      step = wf['steps'][id]
      step['input_connections'].each do |_, v|
        # if v is a list
        if v.is_a?(Array)
          v.each do |v2|
            output += "  #{v2['id']} -->|#{v2['output_name']}| #{id};\n"
          end
        else
          output += "  #{v['id']} -->|#{v['output_name']}| #{id};\n"
        end
      end
    end

    output
  end

  def self.resolve_material(site, material)
    # We've already
    # looked in every /topic/*/tutorials/* folder, and turn these disparate
    # resources into a page_obj as well. Most variables are copied directly,
    # either from a tutorial, or a slides (if no tutorial is available.) This
    # means we do not (cannot) support external_slides AND external_handson.
    # This is probably a sub-optimal situation we'll end up fixing someday.
    #
    tutorials = material['resources'].select { |a| a[0] == 'tutorial' }
    slides    = material['resources'].select { |a| a[0] == 'slides' }
    tours     = material['resources'].select { |a| a[0] == 'tours' }

    # Our final "page" object (a "material")
    page = nil

    slide_has_video = false
    slide_translations = []
    page_ref = nil

    if slides.length.positive?
      page = slides.min { |a, b| a[1].path <=> b[1].path }[1]
      slide_has_video = page.data.fetch('video', false)
      slide_translations = page.data.fetch('translations', [])
      page_ref = page
    end

    # No matter if there were slides, we override with tutorials if present.
    tutorial_translations = []
    if tutorials.length.positive?
      page = tutorials.min { |a, b| a[1].path <=> b[1].path }[1]
      tutorial_translations = page.data.fetch('translations', [])
      page_ref = page
    end

    if page.nil?
      Jekyll.logger.error '[GTN/TopicFilter] Could not process material'
      return {}
    end

    # Otherwise clone the metadata from it which works well enough.
    page_obj = page.data.dup
    page_obj['id'] = "#{page['topic_name']}/#{page['tutorial_name']}"
    page_obj['ref'] = page_ref

    id = page_obj['id']
    page_obj['video_library'] = {}

    if site.data.key?('video-library')
      page_obj['video_library']['tutorial'] = site.data['video-library']["#{id}/tutorial"]
      page_obj['video_library']['slides'] = site.data['video-library']["#{id}/slides"]
      page_obj['video_library']['demo'] = site.data['video-library']["#{id}/demo"]
      page_obj['video_library']['both'] = site.data['video-library'][id]
    end

    page_obj['video_library']['session'] = site.data['session-library'][id] if site.data.key?('session-library')

    # Sometimes `hands_on` is set to something like `external`, in which
    # case it is important to not override it. So we only do that if the
    # key isn't already set. Then we choose to set it to a test for the
    # tutorial being present. We probably don't need to test both, but it
    # is hard to follow which keys are which and safer to test for both in
    # case someone edits the code later. If either of these exist, we can
    # automatically set `hands_on: true`
    page_obj['hands_on'] = tutorials.length.positive? if !page_obj.key?('hands_on')

    # Same for slides, if there's a resource by that name, we can
    # automatically set `slides: true`
    page_obj['slides'] = slides.length.positive? if !page_obj.key?('slides')

    folder = material['dir']

    ymls = Dir.glob("#{folder}/quiz/*.yml") + Dir.glob("#{folder}/quiz/*.yaml")
    if ymls.length.positive?
      quizzes = ymls.map { |a| a.split('/')[-1] }
      page_obj['quiz'] = quizzes.map do |q|
        quiz_data = YAML.load_file("#{folder}/quiz/#{q}")
        {
          'id' => q,
          'path' => "#{folder}/quiz/#{q}",
          'title' => quiz_data['title'],
          'contributors' => quiz_data['contributors'],
        }
      end
    end

    # In dev configuration, this breaks for me. Not sure why config isn't available.
    domain = if !site.config.nil? && site.config.key?('url')
               "#{site.config['url']}#{site.config['baseurl']}"
             else
               'http://localhost:4000//training-material/'
             end
    # Similar as above.
    workflows = Dir.glob("#{folder}/workflows/*.ga") # TODO: support gxformat2
    if workflows.length.positive?
      workflow_names = workflows.map { |a| a.split('/')[-1] }
      page_obj['workflows'] = workflow_names.map do |wf|
        wfid = "#{page['topic_name']}-#{page['tutorial_name']}"
        wfname = wf.gsub(/.ga/, '').downcase
        trs = "api/ga4gh/trs/v2/tools/#{wfid}/versions/#{wfname}"
        wf_path = "#{folder}/workflows/#{wf}"
        wf_json = JSON.parse(File.read(wf_path))
        license = wf_json['license']
        creators = wf_json['creator'] || []
        wftitle = wf_json['name']

        # /galaxy-intro-101-workflow.eu.json
        workflow_test_results = Dir.glob(wf_path.gsub(/.ga$/, '.*.json'))
        workflow_test_outputs = {}
        workflow_test_results.each do |test_result|
          server = workflow_test_results[0].match(/\.(..)\.json$/)[1]
          workflow_test_outputs[server] = JSON.parse(File.read(test_result))
        end
        workflow_test_outputs = nil if workflow_test_outputs.empty?

        {
          'workflow' => wf,
          'tests' => Dir.glob("#{folder}/workflows/" + wf.gsub(/.ga/, '-test*')).length.positive?,
          'url' => "#{domain}/#{folder}/workflows/#{wf}",
          'path' => wf_path,
          'wfid' => wfid,
          'wfname' => wfname,
          'trs_endpoint' => "#{domain}/#{trs}",
          'license' => license,
          'creators' => creators,
          'name' => wf_json['name'],
          'title' => wftitle,
          'test_results' => workflow_test_outputs,
          'modified' => File.mtime(wf_path),
          'mermaid' => mermaid(wf_json),
        }
      end
    end

    # Really only used for tool list install for ephemeris, not general.
    page_obj['api'] = "#{domain}/api/topics/#{page['topic_name']}/tutorials/#{page['tutorial_name']}/tutorial.json"

    # Tool List
    #
    # This is exposed in the GTN API to help admins/devs easily get the tool
    # list for installation.
    page_obj['tools'] = []
    page_obj['tools'] += page.content.scan(/{% tool \[[^\]]*\]\(([^)]*)\)\s*%}/) if page_obj['hands_on']

    page_obj['workflows']&.each do |wf|
      wf_path = "#{folder}/workflows/#{wf['workflow']}"

      wf_data = JSON.parse(File.read(wf_path))
      page_obj['tools'] += extract_workflow_tool_list(wf_data)
    end
    page_obj['tools'] = page_obj['tools'].flatten.sort.uniq

    topic = site.data[page_obj['topic_name']]
    page_obj['supported_servers'] = if topic['type'] == 'use' || topic['type'] == 'basics'
                                      Gtn::Supported.calculate(site.data['public-server-tools'], page_obj['tools'])
                                    else
                                      []
                                    end

    topic_name_human = site.data[page_obj['topic_name']]['title']
    page_obj['topic_name_human'] = topic_name_human # TODO: rename 'topic_name' and 'topic_name' to 'topic_id'
    admin_install = Gtn::Toolshed.format_admin_install(site.data['toolshed-revisions'], page_obj['tools'],
                                                       topic_name_human, site.data['toolcats'])
    page_obj['admin_install'] = admin_install
    page_obj['admin_install_yaml'] = admin_install.to_yaml

    page_obj['tours'] = tours.length.positive?
    page_obj['video'] = slide_has_video
    page_obj['translations'] = {}
    page_obj['translations']['tutorial'] = tutorial_translations
    page_obj['translations']['slides'] = slide_translations
    page_obj['translations']['video'] = slide_has_video # Just demand it?
    # I feel less certain about this override, but it works well enough in
    # practice, and I did not find any examples of `type: <anything other
    # than tutorial>` in topics/*/tutorials/*/tutorial.md but that doesn't
    # make it future proof.
    page_obj['type'] = 'tutorial'

    if page_obj.key?('draft') && page_obj['draft']
      page_obj['tags'] = [] if !page_obj.key? 'tags'
      page_obj['tags'].push('work-in-progress')
    end

    page_obj
  end

  def self.process_pages(site, pages)
    # eww.
    return site.data['cache_processed_pages'] if site.data.key?('cache_processed_pages')

    materials = collate_materials(site, pages).map { |_k, v| resolve_material(site, v) }
    puts '[GTN/TopicFilter] Filling Materials Cache'
    site.data['cache_processed_pages'] = materials

    # Prepare short URLs
    shortlinks = site.data['shortlinks']
    mappings = Hash.new { |h, k| h[k] = [] }

    shortlinks.each_key do |kp|
      shortlinks[kp].each do |k, v|
        mappings[v].push("/short/#{k}")
      end
    end
    # Update the materials with their short IDs + redirects
    pages.select { |p| mappings.keys.include? p.url }.each do |p|
      # Set the short id on the material
      if p['ref']
        # Initialise redirects if it wasn't set
        p['ref'].data['redirect_from'] = [] if !p['ref'].data.key?('redirect_from')
        p['ref'].data['redirect_from'].push(*mappings[p.url])
        p['ref'].data['redirect_from'].uniq!
      else
        p.data['redirect_from'] = [] if !p.data.key?('redirect_from')

        p.data['redirect_from'].push(*mappings[p.url])
        p.data['redirect_from'].uniq!
      end
    end
    # Same for news
    get_posts(site).select { |p| mappings.keys.include? p.url }.each do |p|
      # Set the short id on the material
      p.data['redirect_from'] = [] if !p.data.key?('redirect_from')
      p.data['redirect_from'].push(*mappings[p.url])
      p.data['redirect_from'].uniq!
    end

    materials
  end

  ##
  # This is a helper function to get all the materials in a site.
  def self.list_all_materials(site)
    process_pages(site, site.pages)
  end

  ##
  # This is a helper function to get all the materials in a site.
  def self.list_videos(site)
    materials = process_pages(site, site.pages)
    materials.select { |x| x['video'] == true }
  end

  ##
  # List every tag used across all materials.
  # This is used to generate the tag cloud.
  #
  # Parameters:
  # +site+:: The +Jekyll::Site+ object, used to get the list of pages.
  # Returns:
  # +Array+:: An array of strings, each string is a tag. (sorted and unique)
  #
  def self.list_all_tags(site)
    materials = process_pages(site, site.pages)
    (materials.map { |x| x['tags'] || [] }.flatten + list_topics(site)).sort.uniq
  end

  def self.filter_by_topic(site, topic_name)
    # Here we make a (cached) call to load materials into memory and sort them
    # properly.
    materials = process_pages(site, site.pages)

    # Select out the materials by topic:
    resource_pages = materials.select { |x| x['topic_name'] == topic_name }

    # If there is nothing with that topic name, try generating it by tags.
    resource_pages = materials.select { |x| (x['tags'] || []).include?(topic_name) } if resource_pages.empty?

    # The complete resources we'll return is the introduction slides first
    # (EDIT: not anymore, we rely on prioritisation!)
    # and then the rest of the pages.
    resource_pages = resource_pages.sort_by { |k| k.fetch('priority', 1) }

    Jekyll.logger.error "Error? Could not find any relevant pages for #{topic_name}" if resource_pages.empty?

    resource_pages
  end

  def self.filter_by_tag(site, topic_name)
    # Here we make a (cached) call to load materials into memory and sort them
    # properly.
    materials = process_pages(site, site.pages)

    # Select those with that topic ID or that tag
    resource_pages = materials.select { |x| x['topic_name'] == topic_name }
    resource_pages += materials.select { |x| (x['tags'] || []).include?(topic_name) }

    # The complete resources we'll return is the introduction slides first
    # (EDIT: not anymore, we rely on prioritisation!)
    # and then the rest of the pages.
    resource_pages = resource_pages.sort_by { |k| k.fetch('priority', 1) }

    Jekyll.logger.error "Error? Could not find any relevant tagged pages for #{topic_name}" if resource_pages.empty?

    resource_pages
  end

  ##
  # Filter a list of materials by topic and subtopic.
  def self.filter_by_topic_subtopic(site, topic_name, subtopic_id)
    resource_pages = filter_by_topic(site, topic_name)

    # Select out materials with the correct subtopic
    resource_pages = resource_pages.select { |x| x['subtopic'] == subtopic_id }

    if resource_pages.empty?
      Jekyll.logger.error "Error? Could not find any relevant pages for #{topic_name} / #{subtopic_id}"
    end

    resource_pages
  end

  ##
  # Get a list of contributors for a list of materials
  # Parameters:
  # +materials+:: An array of materials
  # Returns:
  # +Array+:: An array of contributors as strings.
  def self.identify_contributors(materials, site)
    materials
      .map { |_k, v| v['materials'] }.flatten
      # Not 100% sure why this flatten is needed? Probably due to the map over hash
      .map { |mat| Gtn::Contributors.get_contributors(mat) }.flatten.uniq.shuffle
      .reject { |c| Gtn::Contributors.funder?(site, c) }
  end

  ##
  # Get a list of funders for a list of materials
  # Parameters:
  # +materials+:: An array of materials
  # Returns:
  # +Array+:: An array of funders as strings.
  def self.identify_funders(materials, site)
    materials
      .map { |_k, v| v['materials'] }.flatten
      # Not 100% sure why this flatten is needed? Probably due to the map over hash
      .map { |mat| Gtn::Contributors.get_contributors(mat) }.flatten.uniq.shuffle
      .select { |c| Gtn::Contributors.funder?(site, c) }
  end

  ##
  # Get the version of a tool.
  # Parameters:
  # +tool+:: A tool string
  # Returns:
  # +String+:: The version of the tool.
  #
  # Examples:
  # get_version("toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.0") => "1.0.0"
  def self.get_version(tool)
    if tool.count('/') > 4
      tool.split('/')[-1]
    else
      tool
    end
  end

  ##
  # Get a short version of a tool.
  # Parameters:
  # +tool+:: A tool string
  # Returns:
  # +String+:: The short version of the tool.
  #
  # Examples:
  # short_tool("toolshed.g2.bx.psu.edu/repos/galaxyp/regex_find_replace/regex1/1.0.0") => "galaxyp/regex1"
  def self.short_tool(tool)
    if tool.count('/') > 4
      "#{tool.split('/')[2]}/#{tool.split('/')[4]}"
    else
      tool
    end
  end

  ##
  # List materials by tool
  # Parameters:
  # +site+:: The +Jekyll::Site+ object, used to get the list of pages.
  # Returns:
  # +Hash+:: A hash of tool_id => {
  #   "tool_id" => [tool_id, version],
  #   "tutorials" => [tutorial_id, tutorial_title, topic_title, tutorial_url]
  # }
  def self.list_materials_by_tool(site)
    tool_map = {}

    list_all_materials(site).each do |m|
      m.fetch('tools', []).each do |tool|
        sid = short_tool(tool)
        tool_map[sid] = { 'tool_id' => [], 'tutorials' => [] } if !tool_map.key?(sid)

        tool_map[sid]['tool_id'].push([tool, get_version(tool)])
        tool_map[sid]['tutorials'].push([
                                          m['id'], m['title'], site.data[m['topic_name']]['title'], m['url']
                                        ])
      end
    end

    # Uniqueify/sort
    t = tool_map.to_h do |k, v|
      v['tool_id'].uniq!
      v['tool_id'].sort_by! { |k2| k2[1] }
      v['tool_id'].reverse!

      v['tutorials'].uniq!
      v['tutorials'].sort!
      [k, v]
    end

    # Order by most popular tool
    t.sort_by { |_k, v| v['tutorials'].length }.reverse.to_h
  end
end

module Jekyll
  # The "implementation" of the topic filter as liquid accessible filters
  module ImplTopicFilter
    ##
    # List the most recent contributors to the GTN.
    # Parameters:
    # +contributors+:: A hash of contributors
    # +count+:: The number of contributors to return
    # Returns:
    # +Hash+:: A hash of contributors
    #
    # Example:
    # most_recent_contributors(contributors, 5)
    # => {
    #  "hexylena" => {
    #  "name" => "Hexylena",
    #  "avatar" => "https://avatars.githubusercontent.com/u/458683?v=3",
    #  ...
    #  }
    # }
    def most_recent_contributors(contributors, count)
      # Remove non-hof
      hof = contributors.reject { |_k, v| v.fetch('halloffame', 'yes') == 'no' }
      # Get keys + sort by joined date
      hof_k = hof.keys.sort do |x, y|
        hof[y].fetch('joined', '2016-01') <=> hof[x].fetch('joined', '2016-01')
      end

      # Transform back into hash
      hof_k.slice(0, count).to_h { |k| [k, hof[k]] }
    end

    ##
    # Find the most recently modified tutorials
    # Parameters:
    # +site+:: The +Jekyll::Site+ object, used to get the list of pages.
    # +exclude_recently_published+:: Do not include ones that were recently
    #                                published in the slice, to make it look a bit nicer.
    # Returns:
    # +Array+:: An array of the 10 most recently modified pages
    # Example:
    #  {% assign latest_tutorials = site | recently_modified_tutorials %}
    def recently_modified_tutorials(site, exclude_recently_published: true)
      tutorials = site.pages.select { |page| page.data['layout'] == 'tutorial_hands_on' }

      latest = tutorials.sort do |x, y|
        Gtn::ModificationTimes.obtain_time(y.path) <=> Gtn::ModificationTimes.obtain_time(x.path)
      end

      latest_published = recently_published_tutorials(site)
      latest = latest.reject { |x| latest_published.include?(x) } if exclude_recently_published

      latest.slice(0, 10)
    end

    ##
    # Find the most recently published tutorials
    # Parameters:
    # +site+:: The +Jekyll::Site+ object, used to get the list of pages.
    # Returns:
    # +Array+:: An array of the 10 most recently published modified pages
    # Example:
    #  {% assign latest_tutorials = site | recently_modified_tutorials %}
    def recently_published_tutorials(site)
      tutorials = site.pages.select { |page| page.data['layout'] == 'tutorial_hands_on' }

      latest = tutorials.sort do |x, y|
        Gtn::PublicationTimes.obtain_time(y.path) <=> Gtn::PublicationTimes.obtain_time(x.path)
      end

      latest.each { |x| puts [x.path, Gtn::PublicationTimes.obtain_time(x.path)] }
      latest.slice(0, 10)
    end

    def topic_count(resources)
      # Count lines in the table except introduction slides
      resources.length
    end

    ##
    # Fetch a tutorial material's metadata
    # Parameters:
    # +site+:: The +Jekyll::Site+ object, used to get the list of pages.
    # +topic_name+:: The name of the topic
    # +page_name+:: The name of the page
    # Returns:
    # +Hash+:: The metadata for the tutorial material
    #
    # Example:
    #  {% assign material = site | fetch_tutorial_material:page.topic_name,page.tutorial_name%}
    def fetch_tutorial_material(site, topic_name, page_name)
      TopicFilter.fetch_tutorial_material(site, topic_name, page_name)
    end

    def list_topics_ids(site)
      ['introduction'] + TopicFilter.list_topics(site).filter { |k| k != 'introduction' }
    end

    def list_topics_h(site)
      TopicFilter.list_topics(site)
    end

    def list_topics_by_category(site, category)
      q = TopicFilter.list_topics(site).map do |k|
        [k, site.data[k]]
      end

      # Alllow filtering by a category, or return "all" otherwise.
      if category == 'non-tag'
        q = q.select { |_k, v| v['tag_based'].nil? }
      elsif category == 'science'
        q = q.select { |_k, v| %w[use basics].include? v['type'] }
      elsif category == 'technical'
        q = q.select { |_k, v| %w[admin-dev data-science instructors].include? v['type'] }
      elsif category == 'science-technical'
        q = q.select { |_k, v| %w[use basics admin-dev data-science instructors].include? v['type'] }
      elsif category != 'all'
        q = q.select { |_k, v| v['type'] == category }
      end

      # Sort alphabetically by titles
      q.sort { |a, b| a[1]['title'] <=> b[1]['title'] }
    end

    def to_keys(arr)
      arr.map { |k| k[0] }
    end

    def to_vals(arr)
      arr.map { |k| k[1] }
    end

    def list_materials_by_tool(site)
      TopicFilter.list_materials_by_tool(site)
    end

    def list_materials_structured(site, topic_name)
      TopicFilter.list_materials_structured(site, topic_name)
    end

    def list_all_tags(site)
      TopicFilter.list_all_tags(site)
    end

    def topic_filter(site, topic_name)
      TopicFilter.topic_filter(site, topic_name)
    end

    def topic_filter_tutorial_count(site, topic_name)
      TopicFilter.topic_filter(site, topic_name).length
    end

    def identify_contributors(materials, site)
      TopicFilter.identify_contributors(materials, site)
    end

    def identify_funders(materials, site)
      TopicFilter.identify_funders(materials, site)
    end
  end
end

Liquid::Template.register_filter(Jekyll::ImplTopicFilter)
