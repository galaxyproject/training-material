require 'json'
require 'yaml'
require './_plugins/gtn.rb'


module TopicFilter
  def self.list_topics(site)
    site.data.select{|k, v| v.is_a?(Hash) && v.has_key?('maintainers')}.map{|k, v| k}
  end

  def self.topic_filter(site, topic_name)
    if not site.data.has_key?('cache_topic_filter')
      site.data['cache_topic_filter'] = Hash.new

      # For each topic
      self.list_topics(site).each{|topic|
        site.data['cache_topic_filter'][topic] = self.run_topic_filter(site.pages, topic)
      }
    end

    site.data['cache_topic_filter'][topic_name]
  end

  def self.fetch_tutorial_material(site, topic_name, tutorial_name)
    if not site.data.has_key?('cache_topic_filter')
      site.data['cache_topic_filter'] = Hash.new

      # For each topic
      self.list_topics(site).each{|topic|
        site.data['cache_topic_filter'][topic] = self.run_topic_filter(site.pages, topic)
      }
    end

    site.data['cache_topic_filter'][topic_name].select{|p| p['tutorial_name'] == tutorial_name}[0]
  end

  def self.extract_workflow_tool_list(data)
    out = data['steps'].select{|k, v| v['type'] == "tool"}.map{|k, v| v['tool_id']}.select{|x| ! x.nil?}
    out += data['steps'].select{|k, v| v['type'] == "subworkflow"}.map{|k, v| self.extract_workflow_tool_list(v['subworkflow'])}
    out
  end

  def self.run_topic_filter(pages, topic_name)
    # Arrays that will store all introduction slides and tutorials we discover.
    resource_intro = []
    resource_pages = []

    # In order to speed up queries later, we'll store a set of "interesting"
    # pages (i.e. things that are under `topic_name`)
    interesting = {}
    for page in pages do
      page_parts = page.url.split('/')
      # Skip anything outside of topics.
      if not page.url.include?('/topics/') then next end
      page.data['url'] = page.url
      # If the path is long enough and it is under the topic name
      if page_parts.length > 3 and page_parts[2] == topic_name
        # Automate the tutorial-name thing. This writes back to the shared
        # data structure.
        page.data['topic_name'] = page_parts[2]

        # Slides are one directory level shorter and re-use name for their identity.
        if page_parts[3] == 'slides' then
          page.data['tutorial_name'] = page_parts[4].sub(/\.html$/, '')
        else
          page.data['tutorial_name'] = page_parts[4]
        end

        # And then store in our interesting stuff
        key = page.url.sub(/^\//, '') # strip leading slash since later queries don't have it.
        interesting[key] = page
      end
    end

    # Theory is as follows:
    #
    # for each folder we need to extract
    #   if any slides/workflows/tutorials/etc
    #
    # and we need to update this object with the tutorial's title in order to
    # make the object complete. It is called a page_obj or a resource within
    # the context of the script but basically it's the same as the old
    # entries from the metadata.yaml file. A thing that has slides:yes/no,
    # hands_on:yes/no, and all of the other associated data so we can
    # construct the single line in the topic page with the tutorial's
    # resources.

    # All of the files in the intro slides and tutorials directories
    intro_slides = Dir.glob("topics/#{topic_name}/slides/*").sort
    tutorial_folders = Dir.glob("topics/#{topic_name}/tutorials/*").sort

    # Turn these into an object. Override the slides to true (since they're
    # slides) and the type to introduction, to prevent people from having to
    # set those variables. All other variables will be copied directly with
    # the `page.data.dup` so people can use external hands ons and similar.
    for intro_slide in intro_slides do
      page = interesting[intro_slide]
      page_obj = page.data.dup
      page_obj['slides'] = true
      page_obj['type'] = 'introduction'
      resource_intro.push(page_obj)
    end

    # Look in every /topic/*/tutorials/* folder, and turn these disparate
    # resources into a page_obj as well. Most variables are copied directly,
    # either from a tutorial, or a slides (if no tutorial is available.) This
    # means we do not (cannot) support external_slides AND external_handson.
    # This is probably a sub-optimal situation we'll end up fixing someday.
    for folder in tutorial_folders do
      # Discover resources, defined as any file in a given tutorial folder
      resources = Dir.glob("#{folder}/*").map{ |a| a.split('/')[-1] }

      # We will compare this against the 'interesting' hash we built earlier,
      # let's find everything from interesting that is in this specific topic
      # AND tutorial (interesting only contained things in this topic).
      # We add the trailing slash because tutorial directories sharing a
      # common prefix would cause issues.
      known_pages = interesting.select{|a| a.include?("#{folder}/")}

      # We need to extract metadata from some of the pages to determine the
      # proper 'title' attribute for a tutorial. We can only pull this if we
      # have the page object in jekyll to work with. So we'll look for the
      # specific keys in `known_pages` (keyed on filesystem paths) and select
      # for the things we're interested in.
      tutorial_page_keys = known_pages.keys.select{|a| a.include?('tutorial.html')}
      slides_page_keys   = known_pages.keys.select{|a| a.include?('slides.html')}

      # We'll handle slides first and have hands-on override.
      page = false
      slide_has_video = false
      slide_translations = []
      if slides_page_keys.length == 1 then
        page = interesting[slides_page_keys[0]]
        slide_has_video = page.data.fetch('video', false)
        slide_translations = page.data.fetch('translations', [])
      end

      tutorial_translations = []
      if tutorial_page_keys.length == 1 then
        page = interesting[tutorial_page_keys[0]]
        tutorial_translations = page.data.fetch('translations', [])
      end

      # If no tutorial OR slides are found, then we have an issue
      if page == false then
        puts "Error? No tutorial OR slides found in #{folder}. We saw #{known_pages.keys}"
        next
      end

      # Otherwise clone the metadata from it which works well enough.
      page_obj = page.data.dup

      # Sometimes `hands_on` is set to something like `external`, in which
      # case it is important to not override it. So we only do that if the
      # key isn't already set. Then we choose to set it to a test for the
      # tutorial being present. We probably don't need to test both, but it
      # is hard to follow which keys are which and safer to test for both in
      # case someone edits the code later. If either of these exist, we can
      # automatically set `hands_on: true`
      if not page_obj.has_key?("hands_on") then
        page_obj['hands_on'] = resources.include?('tutorial.md') or resources.include?('tutorial.html')
      end

      # Same for slides, if there's a resource by that name, we can
      # automatically set `slides: true`
      if not page_obj.has_key?("slides") then
        page_obj['slides'] = resources.include?('slides.html')
      end

      if resources.include?("quiz") then
        ymls = Dir.glob("#{folder}/quiz/*.yml") + Dir.glob("#{folder}/quiz/*.yaml")
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
      if resources.include?('workflows')
        workflow_files = Dir.glob("#{folder}/workflows/*.ga").map{ |a| a.split('/')[-1] }
        page_obj['workflows'] = workflow_files.map{|wf|
          x = {
            "workflow" => wf,
            "tests" => Dir.glob("#{folder}/workflows/" + wf.gsub(/.ga/, '-test*')).length > 0,
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


      page_obj['tours'] = resources.include?('tours')
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

      if page_obj.has_key?("enable") and !page_obj['enable'] then
        if ! page_obj.has_key? 'tags'
          page_obj['tags'] = Array.new
        end
        page_obj['tags'].push('work-in-progress')
      end

      # Push onto our stack.
      resource_pages.push(page_obj)
    end

    # The complete resources we'll return is the introduction slides first
    # (regardless of alphabetisation), and then the rest of the pages.
    resource_pages = resource_intro + resource_pages.sort_by{ |k| k.fetch('title', '').downcase}

    if resource_pages.length == 0 then
      puts "Error? Could not find any relevant pages for #{topic_name}"
    end

    # Apparently return is optional?
    resource_pages
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

    def last_modified_at(page)
      Gtn::ModificationTimes.obtain_time(page['path'])
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
      resources.select{ |a| a['type'] != 'introduction' }.length
    end

    def fetch_tutorial_material(site, topic_name, page_name)
      TopicFilter.fetch_tutorial_material(site, topic_name, page_name)
    end

    def list_topics(site, category)
      q = TopicFilter.list_topics(site).map{|k|
        [k, site.data[k]]
      }

      # Alllow filtering by a category, or return "all" otherwise.
      if category != "all"
        q = q.select{|k, v| v['type'] == category }
      end

      # Sort alphabetically by titles
      q = q.sort{|a, b| a[1]['title'] <=> b[1]['title'] }

      # But move introduction to the start
      q = q.select{|k, v| k == "introduction"} + q.select{|k, v| k != "introduction"}

      q
    end

    def topic_filter(site, topic_name)
      TopicFilter.topic_filter(site, topic_name)
    end

    ELIXIR_NODES = {
      "au" => "Australia",
      "be" => "Belgium",
      "ch" => "Switzerland",
      "cz" => "Czechia",
      "de" => "Germany",
      "dk" => "Denmark",
      "ee" => "Estonia",
      "es" => "Spain",
      "fi" => "Finland",
      "fr" => "France",
      "gr" => "Greece",
      "hu" => "Hungary",
      "ie" => "Ireland",
      "il" => "Israel",
      "it" => "Italy",
      "lu" => "Luxembourg",
      "nl" => "the Netherlands",
      "no" => "Norway",
      "pt" => "Portugal",
      "se" => "Sweden",
      "si" => "Slovenia",
      "uk" => "United Kingdom",
    }

    def elixirnode2name(name)
      ELIXIR_NODES[name]
    end

    def slugify_unsafe(text)
      # Gets rid of *most* things without making it completely unusable?
      text.gsub(/["'\\\/-;:,.!@#$%^&*()-]/, '').gsub(/\s/, '-')
    end

    def humanize_types(type)
      data = {
        "seq" => "List of Items",
        "str" => "Free Text",
        "map" => "A dictionary/map",
        "float" => "Decimal Number",
        "int" => "Integer Number",
        "bool" => "Boolean"
      }
      data[type]
    end

    def replace_newline_doublespace(text)
      text.gsub(/\n/, "\n  ")
    end
  end
end

Liquid::Template.register_filter(Jekyll::ImplTopicFilter)
