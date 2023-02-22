require 'json'
require 'yaml'
require './_plugins/gtn.rb'


module TopicFilter
  def self.list_topics(site)
    site.data.select{|k, v| v.is_a?(Hash) && v.has_key?('maintainers')}.map{|k, v| k}
  end

  def self.fill_cache(site)
    if not site.data.has_key?('cache_topic_filter')
      puts "[GTN/TopicFilter] Begin Cache Prefill"
      site.data['cache_topic_filter'] = Hash.new

      # For each topic
      self.list_topics(site).each{|topic|
        site.data['cache_topic_filter'][topic] = self.filter_by_topic(site, topic)
      }
      puts "[GTN/TopicFilter] End Cache Prefill"
    end
  end

  def self.topic_filter(site, topic_name)
    self.fill_cache(site)
    site.data['cache_topic_filter'][topic_name]
  end

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

    self.fill_cache(site)

    # Here we want to either return data structured around subtopics
    if site.data[topic_name].has_key?("subtopics")
      # We'll construct a new hash of subtopic => tutorials
      out = Hash.new
      seen_ids = []
      site.data[topic_name]['subtopics'].each{|subtopic, v|
        specific_resources = self.filter_by_topic_subtopic(site, topic_name, subtopic['id'])
        out[subtopic['id']] = {
          "subtopic" => subtopic,
          "materials" => specific_resources
        }
        seen_ids += specific_resources.map{|x| x['id'] }
      }

      # And we'll have this __OTHER__ subtopic for any tutorials that weren't
      # in a subtopic.
      all_topics_for_tutorial = self.filter_by_topic(site, topic_name)
      out["__OTHER__"] = {
        "subtopic" => {"title" => "Other", "description" => "Assorted Tutorials", "id" => "other"},
        "materials" => all_topics_for_tutorial.select{|x| ! seen_ids.include?(x['id']) }
      }
    elsif site.data[topic_name]['tag_based'] and site.data[topic_name]['custom_ordering']
      # TODO
      puts "UNIMPLEMENTED"
      out = Hash.new
    elsif site.data[topic_name]['tag_based'] # Tag based Topic
      # We'll construct a new hash of subtopic(parent topic) => tutorials
      out = Hash.new
      seen_ids = []

      materials = self.filter_by_topic(site, topic_name)

      # Which topics are represented in those materials?
      seen_topics = materials.map{|x| x['topic_name']}.sort

      # Treat them like subtopics, but fake subtopics.
      seen_topics.each{|parent_topic, v|
        specific_resources = materials.select{|x| x['topic_name'] == parent_topic}
        out[parent_topic] = {
          "subtopic" => {"id" => parent_topic, "title" => site.data[parent_topic]['title'], "description" => nil},
          "materials" => specific_resources
        }
        seen_ids += specific_resources.map{|x| x['id'] }
      }

      # And we'll have this __OTHER__ subtopic for any tutorials that weren't
      # in a subtopic.
      all_topics_for_tutorial = self.filter_by_topic(site, topic_name)
      out["__OTHER__"] = {
        "subtopic" => {"title" => "Other", "description" => "Assorted Tutorials", "id" => "other"},
        "materials" => all_topics_for_tutorial.select{|x| ! seen_ids.include?(x['id']) }
      }
    else
      # Or just the list (Jury is still out on this one, should it really be a
      # flat list? Or in this identical structure.)
      out = {
        "__FLAT__" => {
          "subtopic" => nil,
          "materials" => self.filter_by_topic(site, topic_name)
        }
      }
    end

    # Cleanup empty sections
    if out.has_key?("__OTHER__") and out["__OTHER__"]["materials"].length == 0
      out.delete("__OTHER__")
    end

    out.each{|k, v|
      v['materials'].sort_by!{|m| [m.fetch('priority', 1), m['title']] }
    }

    out
  end

  def self.fetch_tutorial_material(site, topic_name, tutorial_name)
    self.fill_cache(site)
    site.data['cache_topic_filter'][topic_name].select{|p| p['tutorial_name'] == tutorial_name}[0]
  end

  def self.extract_workflow_tool_list(data)
    out = data['steps'].select{|k, v| v['type'] == "tool"}.map{|k, v| v['tool_id']}.select{|x| ! x.nil?}
    out += data['steps'].select{|k, v| v['type'] == "subworkflow"}.map{|k, v| self.extract_workflow_tool_list(v['subworkflow'])}
    out
  end

  def self.annotate_path(path)
    parts = path.split('/')
    if parts[0] == '.'
      parts.shift
    end

    if parts[0] != 'topics'
      return nil
    end

    if parts[2] != 'tutorials'
      return nil
    end

    if parts.length < 4
      return nil
    end

    material = {
      "topic" => parts[1], # Duplicate
      "topic_name" => parts[1],
      "material" => parts[1] + '/' + parts[3],
      "tutorial_name" => parts[3],
      "dir" => parts[0..3].join('/'),
    }

    if path =~ /\/faqs\//
      return nil
    end

    if parts[-1] =~ /data[_-]library.yaml/ || parts[-1] =~ /data[_-]manager.yaml/
      return nil
    end

    if parts[-1] == 'tools.yaml'
      return nil
    end

    if parts[4] =~ /tutorial.*\.md/
      material['type'] = 'tutorial'
    elsif parts[4] =~ /slides.*\.html/
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

    return material

  end

  def self.collate_materials(pages)
    # In order to speed up queries later, we'll store a set of "interesting"
    # pages (i.e. things that are under `topic_name`)
    interesting = {}
    for page in pages do
      page_parts = page.url.split('/')
      # Skip anything outside of topics.
      if not page.url.include?('/topics/') then next end

      # Extract the material metadata based on the path
      page.data['url'] = page.url
      material_meta = self.annotate_path(page.path)

      # If unannotated then we want to skip this material.
      if material_meta.nil? then next end

      mk = material_meta['material']

      if not interesting.has_key? mk
        interesting[mk] = material_meta.dup
        interesting[mk].delete("type") # Remove the type since it's specific, not generic
        interesting[mk]["resources"] = []
      end

      page.data['topic_name'] = material_meta['topic_name']
      page.data['tutorial_name'] = material_meta['tutorial_name']
      page.data['dir'] = material_meta['dir']

      interesting[mk]["resources"].push([material_meta["type"], page])
    end

    interesting
  end

  def self.resolve_material(site, material)
    # We've already
    # looked in every /topic/*/tutorials/* folder, and turn these disparate
    # resources into a page_obj as well. Most variables are copied directly,
    # either from a tutorial, or a slides (if no tutorial is available.) This
    # means we do not (cannot) support external_slides AND external_handson.
    # This is probably a sub-optimal situation we'll end up fixing someday.
    #
    tutorials = material["resources"].select{|a| a[0] == 'tutorial'}
    slides    = material["resources"].select{|a| a[0] == 'slides'}
    tours     = material["resources"].select{|a| a[0] == 'tours'}

    # Our final "page" object (a "material")
    page = nil

    slide_has_video = false
    slide_translations = []
    if slides.length > 0 then
      page = slides.sort{|a, b| a[1].path <=> b[1].path}[0][1]
      slide_has_video = page.data.fetch('video', false)
      slide_translations = page.data.fetch('translations', [])
    end

    # No matter if there were slides, we override with tutorials if present.
    tutorial_translations = []
    if tutorials.length > 0 then
      page = tutorials.sort{|a, b| a[1].path <=> b[1].path}[0][1]
      tutorial_translations = page.data.fetch('translations', [])
    end

    if page.nil? then
      puts "[GTN/TopicFilter] Could not process material"
      return {}
    end

    # Otherwise clone the metadata from it which works well enough.
    page_obj = page.data.dup
    page_obj['id'] = page['topic_name'] + '/' + page['tutorial_name']

    id = page_obj['id']
    page_obj['video_library'] = Hash.new

    if site.data.has_key?("video-library") then
      page_obj['video_library']["tutorial"] = site.data['video-library'][id + "/tutorial"]
      page_obj['video_library']["slides"] = site.data['video-library'][id + "/slides"]
      page_obj['video_library']["demo"] = site.data['video-library'][id + "/demo"]
      page_obj['video_library']["both"] = site.data['video-library'][id]
    end

    if site.data.has_key?("session-library") then
      page_obj['video_library']["session"] = site.data['session-library'][id]
    end

    # Sometimes `hands_on` is set to something like `external`, in which
    # case it is important to not override it. So we only do that if the
    # key isn't already set. Then we choose to set it to a test for the
    # tutorial being present. We probably don't need to test both, but it
    # is hard to follow which keys are which and safer to test for both in
    # case someone edits the code later. If either of these exist, we can
    # automatically set `hands_on: true`
    if not page_obj.has_key?("hands_on") then
      page_obj['hands_on'] = tutorials.length > 0
    end

    # Same for slides, if there's a resource by that name, we can
    # automatically set `slides: true`
    if not page_obj.has_key?("slides") then
      page_obj['slides'] = slides.length > 0
    end

    folder = material["dir"]

    ymls = Dir.glob("#{folder}/quiz/*.yml") + Dir.glob("#{folder}/quiz/*.yaml")
    if ymls.length > 0
      quizzes = ymls.map{ |a| a.split('/')[-1] }
      page_obj['quiz'] = quizzes.map{|q|
        quiz_data = YAML.load_file("#{folder}/quiz/#{q}")
        {
          "id" => q,
          "path" => "#{folder}/quiz/#{q}",
          "title" => quiz_data['title'],
          "contributors" => quiz_data['contributors'],
        }
      }
    end

    # Similar as above.
    workflows = Dir.glob("#{folder}/workflows/*.ga") # TODO: support gxformat2
    if workflows.length > 0
      workflow_names = workflows.map{ |a| a.split('/')[-1] }
      page_obj['workflows'] = workflow_names.map{|wf|
        x = {
          "workflow" => wf,
          "tests" => Dir.glob("#{folder}/workflows/" + wf.gsub(/.ga/, '-test*')).length > 0,
          "url" => "#{site.config['url']}#{site.config['baseurl']}/#{folder}/workflows/#{wf}",
        }
        x
      }
    end

    # Tool List
    #
    # This is exposed in the GTN API to help admins/devs easily get the tool
    # list for installation.
    page_obj['tools'] = []
    if page_obj['hands_on']
      page_obj['tools'] += page.content.scan(/{% tool \[[^\]]*\]\(([^)]*)\)\s*%}/)
    end

    if page_obj['workflows']
      page_obj['workflows'].each{|wf|
        wf_path = "#{folder}/workflows/#{wf['workflow']}"

        wf_data = JSON.parse(File.open(wf_path).read)
        page_obj['tools'] += self.extract_workflow_tool_list(wf_data)
      }
    end
    page_obj['tools'] = page_obj['tools'].flatten.sort.uniq

    page_obj['tours'] = tours.length > 0
    page_obj['video'] = slide_has_video
    page_obj['translations'] = Hash.new
    page_obj['translations']["tutorial"] = tutorial_translations
    page_obj['translations']["slides"] = slide_translations
    page_obj['translations']['video'] = slide_has_video # Just demand it?
    # I feel less certain about this override, but it works well enough in
    # practice, and I did not find any examples of `type: <anything other
    # than tutorial>` in topics/*/tutorials/*/tutorial.md but that doesn't
    # make it future proof.
    page_obj['type'] = 'tutorial'

    if page_obj.has_key?("draft") and page_obj['draft'] then
      if ! page_obj.has_key? 'tags'
        page_obj['tags'] = Array.new
      end
      page_obj['tags'].push('work-in-progress')
    end

    page_obj
  end

  def self.process_pages(site, pages)
    # eww.
    if site.data.has_key?('cache_processed_pages') then 
      return site.data['cache_processed_pages']
    end

    materials = self.collate_materials(pages).map{|k, v| self.resolve_material(site, v) }
    puts "[GTN/TopicFilter] Filling Materials Cache"
    site.data['cache_processed_pages'] = materials

    materials
  end

  def self.list_all_materials(site)
    self.process_pages(site, site.pages)
  end

  def self.list_all_tags(site)
    materials = self.process_pages(site, site.pages)
    materials.map{|x| x.fetch('tags', [])}.flatten.sort.uniq
  end

  def self.filter_by_topic(site, topic_name)
    # Here we make a (cached) call to load materials into memory and sort them
    # properly.
    materials = self.process_pages(site, site.pages)

    # Select out the materials by topic:
    resource_pages = materials.select{|x| x['topic_name'] == topic_name}

    # If there is nothing with that topic name, try generating it by tags.
    if resource_pages.length == 0
      resource_pages = materials.select{|x| x.fetch('tags', []).include?(topic_name)}
    end

    # The complete resources we'll return is the introduction slides first
    # (EDIT: not anymore, we rely on prioritisation!)
    # and then the rest of the pages.
    resource_pages = resource_pages.sort_by{ |k| k.fetch('priority', 1)}

    if resource_pages.length == 0 then
      puts "Error? Could not find any relevant pages for #{topic_name}"
    end

    resource_pages
  end

  def self.filter_by_topic_subtopic(site, topic_name, subtopic_id)
    resource_pages = self.filter_by_topic(site, topic_name)

    # Select out materials with the correct subtopic
    resource_pages = resource_pages.select{|x| x['subtopic'] == subtopic_id }

    if resource_pages.length == 0 then
      puts "Error? Could not find any relevant pages for #{topic_name} / #{subtopic_id}"
    end

    resource_pages
  end

  def self.get_contributors(material)
    if material.has_key?("contributors")
      material['contributors']
    else
      material['contributions'].map{|k, v| v}.flatten
    end
  end

  def self.identify_contributors(materials)
    materials
      .map{|k, v| v['materials']}.flatten # Not 100% sure why this flatten is needed? Probably due to the map over hash
      .map{|mat| get_contributors(mat) }.flatten.uniq.shuffle
  end

  def self.get_version(tool)
    if tool.count('/') > 4 then
      tool.split('/')[-1]
    else
      tool
    end
  end

  def self.short_tool(tool)
    if tool.count('/') > 4 then
      short_tool = tool.split('/')[2] + '/' + tool.split('/')[4]
    else
      short_tool = tool
    end
    short_tool
  end

  def self.list_materials_by_tool(site)
    tool_map = Hash.new

    self.list_all_materials(site).each do |m|
      m.fetch('tools', []).each do |tool|
        sid = self.short_tool(tool)
        if ! tool_map.has_key?(sid)
          tool_map[sid] = {"tool_id" => [], "tutorials" => []}
        end

        tool_map[sid]["tool_id"].push([tool, self.get_version(tool)])
        tool_map[sid]["tutorials"].push([
          m['id'], m['title'], site.data[m['topic_name']]['title'], m['url']
        ])
      end
    end

    # Uniqueify/sort
    t = tool_map.map{|k, v| 
      v["tool_id"].uniq!
      v["tool_id"].sort_by!{|k| k[1]}
      v["tool_id"].reverse!

      v["tutorials"].uniq!
      v["tutorials"].sort!
      [k, v]
    }.to_h

    # Order by most popular tool
    t.sort_by{|k, v| v["tutorials"].length}.reverse.to_h
  end
end


module Jekyll
  module ImplTopicFilter

    def most_recent_contributors(contributors, count)
      # Remove non-hof
      hof = contributors.select{ |k, v| v.fetch("halloffame", "yes") != "no" }
      # Get keys + sort by joined date
      hof_k = hof.keys.sort{ |x, y|
        hof[y].fetch('joined', '2016-01') <=> hof[x].fetch('joined', '2016-01')
      }

      # Transform back into hash
      Hash[hof_k.slice(0, count).collect{|k| [k, hof[k]]}]
    end

    def recently_modified_tutorials(site)
      tutorials = site.pages.select{|page| page.data['layout'] == 'tutorial_hands_on' }

      latest = tutorials.sort{ |x, y|
        Gtn::ModificationTimes.obtain_time(y.path) <=> Gtn::ModificationTimes.obtain_time(x.path)
      }
      latest.slice(0, 10)
    end

    def topic_count(resources)
      # Count lines in the table except introduction slides
      resources.length
    end

    def fetch_tutorial_material(site, topic_name, page_name)
      TopicFilter.fetch_tutorial_material(site, topic_name, page_name)
    end

    def list_topics_by_category(site, category)
      q = TopicFilter.list_topics(site).map{|k|
        [k, site.data[k]]
      }

      # Alllow filtering by a category, or return "all" otherwise.
      if category == "non-tag" then
        q = q.select{|k, v| v['tag_based'].nil? }
      elsif category != "all" then
        q = q.select{|k, v| v['type'] == category }
      end

      # Sort alphabetically by titles
      q = q.sort{|a, b| a[1]['title'] <=> b[1]['title'] }

      q
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

    def identify_contributors(materials)
      TopicFilter.identify_contributors(materials)
    end

  end
end

Liquid::Template.register_filter(Jekyll::ImplTopicFilter)
